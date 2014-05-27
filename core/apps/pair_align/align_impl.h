// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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

#ifndef CORE_APPS_PAIR_ALIGN_ALIGN_IMPL_H_
#define CORE_APPS_PAIR_ALIGN_ALIGN_IMPL_H_

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ---------------------------------------------------------------------------
// Function _loadSequences()
// ---------------------------------------------------------------------------

template <typename TSeqSet, typename TNameSet>
bool _loadSequences(TSeqSet & sequences, 
                    TNameSet & fastaIDs,
                    char const * fileName)
{
    MultiFasta multiFasta;
    if (!open(multiFasta.concat, fileName, OPEN_RDONLY))
        return false;
    AutoSeqFormat format;
    guessFormat(multiFasta.concat, format); 
    split(multiFasta, format);
    unsigned seqCount = length(multiFasta);
    resize(sequences, seqCount, Exact());
    resize(fastaIDs, seqCount, Exact());
    for (unsigned i = 0; i < seqCount; ++i) 
    {
        assignSeqId(fastaIDs[i], multiFasta[i], format);
        assignSeq(sequences[i], multiFasta[i], format);
    }
    return (seqCount > 0);
}

// ---------------------------------------------------------------------------
// Function globalAlignment()
// ---------------------------------------------------------------------------

// TODO(holtgrew): Make publically available.
template<typename TStringSet, typename TCargo, typename TSpec>
inline int
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
                Lcs)
{
    return globalAlignment(g, stringSet(g), Lcs());
}

// ---------------------------------------------------------------------------
// Function performPairwiseAlignment()
// ---------------------------------------------------------------------------

template <typename TAlphabet, typename TAlignConfig, typename TScore, typename TSeqFile,
          typename TMethod, typename TDiag, typename TOutputFormat, typename TOutfile>
inline void
performPairwiseAlignment(TScore const & sc,
                         TSeqFile & seqfile,
                         TMethod method,
                         TDiag low,
                         TDiag high,
                         bool banded,
                         TOutputFormat outputFormat,
                         TOutfile & outfile) 
{
    // Load the 2 sequences
    typedef String<TAlphabet> TSequence;
    StringSet<TSequence, Owner<> > sequenceSet;
    StringSet<String<char> > sequenceNames;
    _loadSequences(sequenceSet, sequenceNames, seqfile.c_str());

    // Fix low and high diagonal.
    low = _max(low, -1 * (int) length(sequenceSet[1]));
    high = _min(high, (int) length(sequenceSet[0]));

    // Align the sequences
    Graph<Alignment<StringSet<TSequence, Dependent<> >, void, WithoutEdgeId> > gAlign(sequenceSet);
    
    int aliScore = 0;

    // Banded alignment?
    if (!banded)
    {
        if (method == 0)
            aliScore = globalAlignment(gAlign, sc, TAlignConfig(), NeedlemanWunsch());
        else if (method == 1)
            aliScore = globalAlignment(gAlign, sc, TAlignConfig(), Gotoh());
        else if (method == 2)
            aliScore = localAlignment(gAlign, sc);
        else if (method == 3)
            aliScore = globalAlignment(gAlign, Lcs());
    }
    else
    {
        if (method == 0)
            aliScore = globalAlignment(gAlign, sc, TAlignConfig(), low, high, NeedlemanWunsch());
        else if (method == 1)
            aliScore = globalAlignment(gAlign, sc, TAlignConfig(), low, high, Gotoh());
    }
    
    // Alignment output
    std::cout << "Alignment score: " << aliScore << std::endl;
    if (outputFormat == 0)
    {
        FILE* strmWrite = fopen(outfile.c_str(), "w");
        write(strmWrite, gAlign, sequenceNames, FastaFormat());
        fclose(strmWrite);
    }
    else if (outputFormat == 1)
    {
        FILE* strmWrite = fopen(outfile.c_str(), "w");
        write(strmWrite, gAlign, sequenceNames, MsfFormat());
        fclose(strmWrite);
    }
}

#endif  // #ifndef CORE_APPS_PAIR_ALIGN_ALIGN_IMPL_H_
