// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNER_H_
#define EXTRAS_INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNER_H_

#include <map>
#include <vector>

#include <seqan/index.h>
#include <seqan/index/find_pigeonhole.h>
#include <seqan/store.h>

#include "consensus_alignment_options.h"

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class OverlapInfo_
// ----------------------------------------------------------------------------

// Store an overlap candidate between sequences seq0 and seq1 (seq0 left of seq1), starting at pos1 in seq0:
//
// Examples (independent of relative read lengths):
//
// (1)   XXXXXXXXXXXXXXXXXX      seq0
//         XXXXXXXXXXXXXXXXXXXX  seq1, pos1 = 2
//
// (2)   XXXXXXXXXXXXXXXXXX      seq0
//       XXXXXXXXXXXXXXXXXXXX  seq1, pos1 = 0
//
// (3)   XXXXXXXXXXXXXXXXXX      seq0
//           XXXXXXXXXXXXXXXXXXXX  seq1, pos1 = 4

struct OverlapInfo_
{
    unsigned seq0;
    unsigned seq1;
    int len0, pos1;
    int numErrors;  // number of alignment errors

    OverlapInfo_() : seq0(-1), seq1(-1), len0(0), pos1(0), numErrors(0) {}
    OverlapInfo_(unsigned seq0, unsigned seq1, int len0, int pos1, int numErrors = -1) :
            seq0(seq0), seq1(seq1), len0(len0), pos1(pos1), numErrors(numErrors) {}

    bool operator<(OverlapInfo_ const & other) const
    {
        return (std::make_pair(std::make_pair(seq0, seq1), pos1) <
                std::make_pair(std::make_pair(other.seq0, other.seq1), other.pos1));
    }
};

// ----------------------------------------------------------------------------
// Class ConsensusAligner_
// ----------------------------------------------------------------------------

// Implementation helper of the consensusAlignment() function.

template <typename TFragmentStore>
class ConsensusAligner_
{
public:
    ConsensusAligner_(TFragmentStore & store, ConsensusAlignmentOptions const & options) :
            store(store), options(options)
    {}

    void run();

private:

    // Compute overlap info using banded alignment.
    OverlapInfo_ computeOverlapInfo(unsigned lhs, unsigned rhs, int lDiag, int uDiag) const;

    // Build position-based all-to-all overlap infos.
    void buildPositionBasedOverlapInfos(std::vector<OverlapInfo_> & infos) const;
    // Build contig-wise all-to-all overlap infos.
    void buildContigAllToAllOverlapInfos(std::vector<OverlapInfo_> & infos) const;
    // Build all-to-all overlap infos for all reads.
    void buildAllToAllOverlapInfos(std::vector<OverlapInfo_> & infos) const;
    // Build all-to-all overlap infos for the given read ids.
    void buildAllToAllOverlapInfos(std::vector<OverlapInfo_> & infos,
                                        std::vector<unsigned> const & readIDs) const;

    // The FragmentStore to use for consensus computation.
    TFragmentStore & store;
    // The configuration of the consensus alignment.
    ConsensusAlignmentOptions const & options;
};

template <typename TFragmentStore>
inline OverlapInfo_ ConsensusAligner_<TFragmentStore>::computeOverlapInfo(
        unsigned lhs, unsigned rhs, int lDiag, int uDiag) const
{
    typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
    typedef typename Value<TReadSeqStore>::Type TReadSeq;

    Align<TReadSeq> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), store.readSeqStore[lhs]);  // TODO(holtgrew): setSource cannot use infix!
    assignSource(row(align, 1), store.readSeqStore[rhs]);

    Score<int, Simple> scoringScheme(1, -1, -1);

    AlignConfig<true, true, true, true> alignConfig;

    // TODO(holtgrew): Fix bands, add bands.
    globalAlignment(align, scoringScheme, alignConfig);

    if (options.verbosity >= 2)
        std::cerr << "before swapping\n" << align << "\n";

    // Create overlap candidate.
    if (isGap(row(align, 0), 0))  // lhs is right and rhs is left, swap roles
    {
        using std::swap;
        swap(lhs, rhs);
        swap(row(align, 0), row(align, 1));
    }

    if (options.verbosity >= 2)
        std::cerr << "after swapping\n" << align << "\n";

    // Count errors in overlap of alignment.
    int beginPos = countGaps(begin(row(align, 1), Standard()));
    int endPosH = length(row(align, 0));
    for (; endPosH > beginPos && isGap(row(align, 0), endPosH - 1); --endPosH)
        continue;
    int endPosV = length(row(align, 1));
    for (; endPosV > beginPos && isGap(row(align, 1), endPosV - 1); --endPosV)
        continue;
    int endPos = std::min(endPosH, endPosV);

    typedef typename Row<Align<TReadSeq> >::Type TRow;
    typedef typename Iterator<TRow, Standard>::Type TIterator;
    TIterator itH = iter(row(align, 0), beginPos, Standard());
    TIterator itHEnd = iter(row(align, 0), endPos, Standard());
    TIterator itV = iter(row(align, 1), beginPos, Standard());
    TIterator itVEnd = iter(row(align, 1), endPos, Standard());
    int numErrors = 0;
    for (; itH != itHEnd; ++itH, ++itV)
        numErrors += (isGap(itH) || isGap(itV) || (*itH != *itV));
    SEQAN_ASSERT(itV == itVEnd);

    // Build result.
    return OverlapInfo_(lhs, rhs, length(store.readSeqStore[lhs]), beginPos, numErrors);
}

template <typename TFragmentStore>
inline void ConsensusAligner_<TFragmentStore>::buildPositionBasedOverlapInfos(
        std::vector<OverlapInfo_> & infos) const
{
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadStoreIter;

    // Obtain sorted copy of store's aligned read store.
    TAlignedReadStore sortedAlignedReads = store.alignedReadStore;
    sortAlignedReads(sortedAlignedReads, SortBeginPos());
    sortAlignedReads(sortedAlignedReads, SortContigId());

    // Iterate over aligned read store, compute overlaps indicated in this multi-read alignment, and collect them.
    std::vector<OverlapInfo_> overlaps;  // will be expanded to alignments later
    TAlignedReadStoreIter it = begin(store.alignedReadStore, Standard());
    TAlignedReadStoreIter itEnd = end(store.alignedReadStore, Standard());
    for (; it != itEnd; ++it)
        // Consider all overlaps with it->readId right of it but overlapping.
        for (TAlignedReadStoreIter it2 = it; it2 != itEnd; ++it2)
        {
            if (it == it2)
                continue;  // no overlap with self
            if (it2->contigId != it->contigId || it2->beginPos > it->endPos + options.posDelta)
                break;  // traversed contig or went too far to the right on same contig
            int pos = it2->beginPos - it->beginPos;
            infos.push_back(computeOverlapInfo(it->readId, it2->readId,
                                               pos - options.posDelta, pos + options.posDelta));
            if (options.verbosity >= 2)
                std::cerr << "OVERLAP\t" << store.readSeqStore[infos.back().seq0] << "\tseq0=" << infos.back().seq0
                          << "\t" << store.readSeqStore[infos.back().seq0] << "\tseq1=" << infos.back().seq1
                          << "\tlen0=" << infos.back().len0 << "\tpos1=" << infos.back().pos1
                          << "\tnumErrors=" << infos.back().numErrors << "\n";
        }
}

template <typename TFragmentStore>
inline void ConsensusAligner_<TFragmentStore>::buildContigAllToAllOverlapInfos(
        std::vector<OverlapInfo_> & infos) const
{
    std::map<unsigned, std::vector<unsigned> > readIDs;  // factorized by contig ID
    for (unsigned i = 0; i < length(store.alignedReadStore); ++i)
        readIDs[store.alignedReadStore[i].contigId].push_back(store.alignedReadStore[i].readId);

    for (std::map<unsigned, std::vector<unsigned> >::const_iterator it = readIDs.begin(); it != readIDs.end(); ++it)
        buildAllToAllOverlapInfos(infos, it->second);
}

template <typename TFragmentStore>
inline void ConsensusAligner_<TFragmentStore>::buildAllToAllOverlapInfos(
        std::vector<OverlapInfo_> & infos) const
{
    std::vector<unsigned> readIDs;
    for (unsigned readID = 0; readID < length(store.readSeqStore); ++readID)
        readIDs.push_back(readID);
    buildAllToAllOverlapInfos(infos, readIDs);
}

template <typename TFragmentStore>
inline void ConsensusAligner_<TFragmentStore>::buildAllToAllOverlapInfos(
        std::vector<OverlapInfo_> & infos,
        std::vector<unsigned> const & readIDs) const
{
    typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
    typedef typename Value<TReadSeqStore>::Type TReadSeq;
    typedef typename Value<TReadSeq>::Type TAlphabet;
    typedef StringSet<TReadSeq, Dependent<> > TStringSet;
    typedef typename Iterator<TStringSet, Standard>::Type TStringSetIter;

    // Get copy of fragment store's read seq store.
    TStringSet superSet;
    for (unsigned i = 0; i < length(store.readSeqStore); ++i)
        appendValue(superSet, store.readSeqStore[i]);

    // Obtain subset of fragment store's read sequences.
    TStringSet subSet;
    for (std::vector<unsigned>::const_iterator it = readIDs.begin(); it != readIDs.end(); ++it)
        assignValueById(subSet, superSet, *it);

    // Build q-gram index over the read subset.
    typedef Shape<TAlphabet, OneGappedShape>        TShape;
    typedef IndexQGram<TShape, OpenAddressing>      TIndexSpec;
    typedef Index<TStringSet const, TIndexSpec>     TIndex;
    typedef Pattern<TIndex, Pigeonhole<> >          TFilterPattern;
    typedef Finder<Dna5String const, Pigeonhole<> > TFilterFinder;

    double maxErrorRate = 1.0 / options.kMerSize;

    TIndex index(subSet);
    TFilterPattern filterPattern(index);
    _patternInit(filterPattern, maxErrorRate);

    // Perform the pigeonhole-based search.
    for (TStringSetIter it = begin(subSet, seqan::Rooted()); !atEnd(it); ++it)
    {
        unsigned seq0 = position(it);
        TFilterFinder filterFinder(*it);
        while (find(filterFinder, filterPattern, maxErrorRate))
        {
            if (length(countOccurrencesMultiple(index, filterPattern.shape)) > (unsigned)options.kMerMaxOcc)
                continue;  // Ignore, too many matching sequences.

            unsigned seq1 = filterFinder.curHit->ndlSeqNo;
            // TODO(holtgrew): Care about dupes, i.e. use >=.
            if (seq1 == seq0)
                continue;  // Skip hits with self.
            __int64 lDiag = filterFinder.curHit->hstkPos;
            __int64 uDiag = filterFinder.curHit->hstkPos + filterFinder.curHit->bucketWidth - length(subSet[seq1]);
            SEQAN_ASSERT_GEQ(uDiag, lDiag);

            infos.push_back(computeOverlapInfo(positionToId(subSet, seq0), positionToId(subSet, seq1),
                                               lDiag, uDiag));

            if (options.verbosity >= 2)
                std::cerr << "OVERLAP\t" << store.readSeqStore[infos.back().seq0] << "\tseq0=" << infos.back().seq0
                          << "\t" << store.readSeqStore[infos.back().seq0] << "\tseq1=" << infos.back().seq1
                          << "\tlen0=" << infos.back().len0 << "\tpos1=" << infos.back().pos1
                          << "\tnumErrors=" << infos.back().numErrors << "\n";
        }
    }
}

template <typename TFragmentStore>
inline void ConsensusAligner_<TFragmentStore>::run()
{
    // Build overlap infos, based on position and contig ID if configured to do so.
    std::vector<OverlapInfo_> overlapInfos;
    if (!options.useContigID)
        buildAllToAllOverlapInfos(overlapInfos);
    else if (!options.usePositions)
        buildContigAllToAllOverlapInfos(overlapInfos);
    else
        buildPositionBasedOverlapInfos(overlapInfos);
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNER_H_
