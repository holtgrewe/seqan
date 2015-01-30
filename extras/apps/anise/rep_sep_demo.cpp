#include <iostream>

#include "rep_sep/clique_enum.h"
#include "rep_sep/clique_set.h"
#include "rep_sep/cluster_linking.h"
#include "rep_sep/feature_map.h"
#include "rep_sep/local_variation_store.h"
#include "rep_sep/merge_read_set.h"
#include "rep_sep/pair_based_sep.h"
#include "rep_sep/read.h"
#include "rep_sep/read_set.h"
#include "rep_sep/read_cover.h"
#include "rep_sep/separating_columns.h"
#include "rep_sep/string_packed_pop_count.h"

using namespace rep_sep;

FeatureMap buildFeatureMap(LocalVariationStore const & varStore)
{
    FeatureMap featureMap;

    for (unsigned columnID = 0; columnID < length(varStore.columns); ++columnID)
        featureMap.insert(FeatureDescription(/*id=*/columnID,
                                             /*kind=*/FeatureDescription::COLUMN,
                                             /*contigID=*/varStore.positions[columnID].first,
                                             /*pos=*/varStore.positions[columnID].second,
                                             /*coverage=*/length(varStore.columns[columnID])));

    featureMap.refresh();
    return featureMap;
}

class ClusterLinkingDriver
{
    struct ReadFeature
    {
        unsigned readID { (unsigned)-1 };
        int value { -1 };

        bool operator<(ReadFeature const & other) const
        { return std::make_pair(readID, value) < std::make_pair(other.readID, other.value); }

        ReadFeature() = default;
        ReadFeature(unsigned readID, int value) : readID(readID), value(value) {}
    };

public:
    ClusterLinkingDriver(FeatureReadSet & readSet,
                         FeatureReadSet const & atomicReadSet,
                         FeatureMap const & featureMap,
                         TFragmentStore /*const*/ & inStore,
                         ReadSeparatorOptions const & options) :
            readSet(readSet), atomicReadSet(atomicReadSet), featureMap(featureMap), /*inStore(inStore),*/
            options(options),
            contigID(atomicReadSet.numActualReads, (unsigned)-1),
            globalizer(/*minOverlap=*/3, (options.verbosity >= 2))
    {
        globalizer.init(inStore);  // for mate ids
    }

    void run();

private:

    // Build map from feature id to read.
    void buildMappings();
    // Perform cluster linking step.
    void performLinking();
    // Build resulting FeatureReadSet.
    void buildResult();

    // Output / Input

    FeatureReadSet & readSet;
    FeatureReadSet const & atomicReadSet;
    FeatureMap const & featureMap;
    // Input fragment store.
    // TFragmentStore /*const*/ & inStore;
    // Configuration for repeat separation.
    ReadSeparatorOptions const & options;

    // State

    std::map<unsigned, std::set<ReadFeature>> readsForFeature;
    std::vector<unsigned> contigID;  // read to contig
    ClusterGlobalizer globalizer;
};

void ClusterLinkingDriver::buildMappings()
{
    // Build map from feature id to read.
    unsigned readSetEntryID = 0;
    for (auto const & read : atomicReadSet)
    {
        for (unsigned readID : read.subReads)
            contigID.at(readID) = read.contigID;

        for (auto feature : read.features)
        {
            auto it = featureMap.find(feature.id);
            if (it == featureMap.end() || it->kind != FeatureDescription::COLUMN)
                continue;
            for (unsigned readID : read.subReads)
                readsForFeature[feature.id].insert(ReadFeature(readID, feature.value));
        }

        ++readSetEntryID;
    }
}

void ClusterLinkingDriver::performLinking()
{
    // Perform cluster linking steps.
    //
    // The mapping values is for mapping the actual values to 0..(k-1).
    for (auto const & featureReads : readsForFeature)
    {
        auto const & featureValues = featureReads.second;
        std::map<unsigned, unsigned> mapping;  // read -> value (first actual, later 0..(k-1)).
        for (auto const & pair : featureValues)
            mapping[pair.readID] = pair.value;
        std::map<int, int> values;  // actual -> 0..(k-1)
        for (auto const & pair : mapping)
            if (!values.count(pair.second))
            {
                unsigned k = values.size();
                values[pair.second] = k;
            }
        if (options.verbosity >= 3)
        {
            std::cerr << "VALUES\n";
            for (auto const & pair : values)
                std::cerr << pair.first << " -> " << pair.second << "\n";
        }
        for (auto & pair : mapping)
            pair.second = values[pair.second];
        globalizer.run(values.size(), mapping);
    }
}

void ClusterLinkingDriver::buildResult()
{
    // Finalize and build resulting FeatureReadSet.
    readSet.reads.clear();
    globalizer.stop();

    for (auto const & cluster : globalizer.partition())
    {
        if (cluster.empty())
            continue;
        Read read;
        bool first = true;
        for (auto readID : cluster)
            if (first)
            {
                read = atomicReadSet.reads.at(readID);
                SEQAN_ASSERT_EQ(read.id, readID);
                read.id = readSet.reads.size();
                first = false;
            }
            else
            {
                read.mergeWithThis(atomicReadSet.reads.at(readID), -1, true, true);
            }
        readSet.reads.emplace_back(std::move(read));
    }
}

void ClusterLinkingDriver::run()
{
    buildMappings();
    performLinking();
    buildResult();

    if (options.verbosity >= 2)
        std::cerr << "Done Cluster Linking.\nComputing read cover.\n";
}


int main()
{
    using seqan::Dna5;
    using seqan::Standard;
    using rep_sep::LocalVariationStore;

    LocalVariationStore lvs;
    lvs.consensus = "CGAA";

    resize(lvs.positions, 4);
    lvs.positions[0] = std::make_pair(0, 1);
    lvs.positions[1] = std::make_pair(0, 2);
    lvs.positions[2] = std::make_pair(0, 3);
    lvs.positions[3] = std::make_pair(0, 4);

    resize(lvs.profile, 4);
    lvs.profile[0].count[ordValue(Dna5('C'))] = 3;
    lvs.profile[0].count[ordValue(Dna5('T'))] = 1;
    lvs.profile[1].count[ordValue(Dna5('G'))] = 2;
    lvs.profile[1].count[ordValue(Dna5('A'))] = 2;
    lvs.profile[2].count[ordValue(Dna5('A'))] = 3;
    lvs.profile[2].count[ordValue(Dna5('G'))] = 1;
    lvs.profile[2].count[ordValue(Dna5('C'))] = 1;
    lvs.profile[3].count[ordValue(Dna5('C'))] = 3;
    lvs.profile[3].count[ordValue(Dna5('T'))] = 2;

    resize(lvs.columns, 4);
    resize(lvs.columns[0], 4);
    resize(lvs.columns[1], 4);
    resize(lvs.columns[2], 5);
    resize(lvs.columns[3], 5);

    // TODO(holtgrew): fill columns
    typedef LocalVariationStore::TColumnEntry TColumnEntry;
    typedef LocalVariationStore::TConsensusAlphabet C;
    lvs.columns[0][0] = TColumnEntry(0u, C('C'), 'I');
    lvs.columns[0][1] = TColumnEntry(2u, C('C'), 'I');
    lvs.columns[0][2] = TColumnEntry(4u, C('C'), 'I');
    lvs.columns[0][3] = TColumnEntry(6u, C('T'), 'I');

    lvs.columns[1][0] = TColumnEntry(2u, C('G'), 'I');
    lvs.columns[1][1] = TColumnEntry(4u, C('G'), 'I');
    lvs.columns[1][2] = TColumnEntry(6u, C('A'), 'I');
    lvs.columns[1][3] = TColumnEntry(8u, C('A'), 'I');

    lvs.columns[2][0] = TColumnEntry(1u, C('A'), 'I');
    lvs.columns[2][1] = TColumnEntry(3u, C('A'), 'I');
    lvs.columns[2][2] = TColumnEntry(5u, C('A'), 'I');
    lvs.columns[2][3] = TColumnEntry(7u, C('G'), 'I');
    lvs.columns[2][4] = TColumnEntry(9u, C('C'), 'I');

    lvs.columns[3][0] = TColumnEntry(1u, C('C'), 'I');
    lvs.columns[3][1] = TColumnEntry(4u, C('C'), 'I');
    lvs.columns[3][2] = TColumnEntry(5u, C('C'), 'I');
    lvs.columns[3][3] = TColumnEntry(7u, C('T'), 'I');
    lvs.columns[3][4] = TColumnEntry(9u, C('T'), 'I');

    // fill coveringReads
    resize(lvs.coveringReads, 4);
    for (int i = 0; i < 4; ++i) {
        for (unsigned j = 0; j < length(lvs.columns[i]); ++j)
            appendValue(lvs.coveringReads[i], lvs.columns[i][j].i1);
        std::sort(begin(lvs.coveringReads[i], Standard()), end(lvs.coveringReads[i], Standard()));
    }

    // fill compression map
    for (int i = 0; i < 10; ++i)
        appendValue(lvs.compressionMap[i], i);

    std::cerr << "CONSENSUS\t" << lvs.consensus << "\n";
    std::cerr << "PROFILE\t";
    for (int pos = 0; pos < 4; ++pos)
        for (int i = 0; i <= 5; ++i)
            std::cerr << "lvs.profile[" << pos << "].count[" << i << "] == " << lvs.profile[pos].count[i] << "\n";
    
    TFragmentStore store;
    for (unsigned i = 0; i < 5; ++i)
        appendMatePair(store, "", "", "", "");
    for (unsigned i = 0; i < 10; ++i)
        appendAlignedRead(store, i, 0, 0, 0, i / 2);

    FeatureMap featureMap = buildFeatureMap(lvs);
    std::cerr << "\nfeatureMap\n";
    featureMap.print(std::cerr);

    FeatureReadSet readSet, atomicReadSet;
    buildAtomicReadSet(atomicReadSet, store, lvs);
    buildFeatureReadSet(readSet, store, lvs);

    FeatureReadSet initialReadSet(readSet);
    ReadSeparatorOptions options;
    options.verbosity = 11;
    ClusterLinkingDriver helper(readSet, atomicReadSet, featureMap, store, options);
    helper.run();

    std::cerr << "readSet\n";
    readSet.print(std::cerr);

    return 0;
}
