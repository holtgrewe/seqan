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
// Class OverlapCandidate_
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

struct OverlapCandidate_
{
    unsigned seq0;
    unsigned seq1;
    int pos1;

    OverlapCandidate_() : seq0(-1), seq1(-1), pos1(0) {}
    OverlapCandidate_(unsigned seq0, unsigned seq1, int pos1) : seq0(seq0), seq1(seq1), pos1(pos1) {}
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

    // Computes an overlap candidate by using global alignment.  Returns default constructed OverlapCandidate_ in case
    // of insufficient similarity.
    OverlapCandidate_ buildCandidate(unsigned lhs, unsigend rhs) const;

    // Build contig-wise all-to-all overlap candidates.
    void buildContigAllToAllOverlapCandidates(std::vector<OverlapCandidate_> & candidates) const;
    // Build all-to-all overlap candidates for the given read ids.
    void buildAllToAllOverlapCandidates(std::vector<OverlapCandidate_> & candidates,
                                        std::vector<unsigned> const & readIDs) const;
    // Build all-to-all overlap candidates for all reads.
    void buildAllToAllOverlapCandidates(std::vector<OverlapCandidate_> & candidates) const;

    // The FragmentStore to use for consensus computation.
    TFragmentStore & store;
    // The configuration of the consensus alignment.
    ConsensusAlignmentOptions const & options;
};

// TODO(holtgrew): Use q-gram index to speed up this step.
inline OverlapCandidate_ ConsensusAligner_<TFragmentStore>::buildCandidate(unsigned lhs, unsigend rhs) const
{
    OverlapCandidate_ result;

    return result;
}

template <typename TFragmentStore>
inline void ConsensusAligner_<TFragmentStore>::buildContigAllToAllOverlapCandidates(
        std::vector<OverlapCandidate_> & candidates) const
{
    std::map<unsigned, std::vector<unsigned> > readIDs;  // factorized by contig ID
    for (unsigned i = 0; i < length(store.alignedReadStore); ++i)
        readIDs[store.alignedReadStore[i].contigId].push_back(store.alignedReadStore[i].readId);

    for (std::map<unsigned, std::vector<unsigned> >::const_iterator it = readIDs.begin(); it != readIDs.end(); ++it)
        buildAllToAllOverlapCandidates(candidates, it->second);
}

template <typename TFragmentStore>
inline void ConsensusAligner_<TFragmentStore>::buildAllToAllOverlapCandidates(
        std::vector<OverlapCandidate_> & candidates,
        std::vector<unsigned> const & readIDs) const
{
    typedef std::vector<unsigned>::const_iterator TIter;
    for (TIter it = readIDs.begin(); it != readIDs.end(); ++it)
        for (TIter it2 = it; it2 != readIDs.end(); ++it2)
            if (it != it2)
            {
                OverlapCandidate_ cand = buildCandidate(*it, *it2);
                if (cand.seq0 == (unsigned)-1)  // invalid
                    continue;
                candidates.push_back(cand);
            }
}

template <typename TFragmentStore>
inline void ConsensusAligner_<TFragmentStore>::buildAllToAllOverlapCandidates(
        std::vector<OverlapCandidate_> & candidates) const
{
    std::vector<unsigned> readIDs;
    for (unsigned readID = 0; readID < length(store.readSeqStore); ++readID)
        readIDs.push_back(readID);
    buildAllToAllOverlapCandidates(candidates, readIDs);
}

template <typename TFragmentStore>
inline void ConsensusAligner_<TFragmentStore>::run()
{
    // Build overlap candidates, based on position and contig ID if configured to do so.
    std::vector<OverlapCandidate_> overlapCandidates;
    if (!options.useContigID)
        buildAllToAllOverlapCandidates(overlapCandidates);
    else if (!options.usePositions)
        buildContigAllToAllOverlapCandidates(overlapCandidates);
    else
        buildPositionBasedOverlapCandidates(overlapCandidates);
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNER_H_
