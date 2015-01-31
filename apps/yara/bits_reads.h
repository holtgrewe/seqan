// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains functions to map read sequences to read/pair ids.
// ==========================================================================

#ifndef APP_YARA_BITS_READS_H_
#define APP_YARA_BITS_READS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

template <typename TReadSeqs>
inline typename Size<TReadSeqs>::Type
getReadSeqsCount(TReadSeqs const & readSeqs)
{
    return length(readSeqs);
}

template <typename TReadSeqs>
inline typename Size<TReadSeqs>::Type
getReadsCount(TReadSeqs const & readSeqs)
{
    return length(readSeqs) / 2;
}

template <typename TReadSeqs>
inline typename Size<TReadSeqs>::Type
getPairsCount(TReadSeqs const & readSeqs)
{
    return length(readSeqs) / 4;
}

template <typename TReadSeqs, typename TReadSeqID>
inline bool isFwdReadSeq(TReadSeqs const & readSeqs, TReadSeqID readSeqID)
{
    SEQAN_ASSERT_LT(readSeqID, getReadSeqsCount(readSeqs));
    return readSeqID < getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqID>
inline bool isRevReadSeq(TReadSeqs const & readSeqs, TReadSeqID readSeqID)
{
    SEQAN_ASSERT_LT(readSeqID, getReadSeqsCount(readSeqs));
    return !isFwdReadSeq(readSeqs, readSeqID);
}

template <typename TReadSeqs, typename TReadSeqID>
inline bool isFirstMate(TReadSeqs const & readSeqs, TReadSeqID readSeqID)
{
    SEQAN_ASSERT_LT(readSeqID, getReadSeqsCount(readSeqs));
    return isFwdReadSeq(readSeqs, readSeqID) ? readSeqID < getPairsCount(readSeqs) : readSeqID < getPairsCount(readSeqs) + getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqID>
inline bool isSecondMate(TReadSeqs const & readSeqs, TReadSeqID readSeqID)
{
    SEQAN_ASSERT_LT(readSeqID, getReadSeqsCount(readSeqs));
    return !isFirstMate(readSeqs, readSeqID);
}

template <typename TReadSeqs, typename TPairID>
inline typename Size<TReadSeqs>::Type
getFirstMateFwdSeqID(TReadSeqs const & /* readSeqs */, TPairID pairID)
{
    return pairID;
}

template <typename TReadSeqs, typename TPairID>
inline typename Size<TReadSeqs>::Type
getSecondMateFwdSeqID(TReadSeqs const & readSeqs, TPairID pairID)
{
    SEQAN_ASSERT_LT(pairID, getPairsCount(readSeqs));
    return pairID + getPairsCount(readSeqs);
}

template <typename TReadSeqs, typename TPairID>
inline typename Size<TReadSeqs>::Type
getFirstMateRevSeqID(TReadSeqs const & readSeqs, TPairID pairID)
{
//    SEQAN_ASSERT_LT(pairID, getPairsCount(readSeqs));
    return getFirstMateFwdSeqID(readSeqs, pairID) + getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TPairID>
inline typename Size<TReadSeqs>::Type
getSecondMateRevSeqID(TReadSeqs const & readSeqs, TPairID pairID)
{
    SEQAN_ASSERT_LT(pairID, getPairsCount(readSeqs));
    return getSecondMateFwdSeqID(readSeqs, pairID) + getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqID>
inline typename Size<TReadSeqs>::Type
getReadID(TReadSeqs const & readSeqs, TReadSeqID readSeqID)
{
    return isFwdReadSeq(readSeqs, readSeqID) ? readSeqID : readSeqID - getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqID>
inline typename Size<TReadSeqs>::Type
getPairID(TReadSeqs const & readSeqs, TReadSeqID readSeqID)
{
    SEQAN_ASSERT_LT(readSeqID, getReadSeqsCount(readSeqs));

    typename Size<TReadSeqs>::Type pairID = readSeqID;

    if (isRevReadSeq(readSeqs, readSeqID))
        pairID -= getReadsCount(readSeqs);

    if (isSecondMate(readSeqs, readSeqID))
        pairID -= getPairsCount(readSeqs);

    SEQAN_ASSERT_LT(pairID, getPairsCount(readSeqs));
    
    return pairID;
}

template <typename TReadSeqs, typename TReadID>
inline typename Size<TReadSeqs>::Type
getMateID(TReadSeqs const & readSeqs, TReadID readID)
{
    typename Size<TReadSeqs>::Type pairID = getPairID(readSeqs, readID);

    return isFirstMate(readSeqs, readID) ? getSecondMateFwdSeqID(readSeqs, pairID) : getFirstMateFwdSeqID(readSeqs, pairID);
}

template <typename TReadSeqs, typename TReadSeqID>
inline typename Size<TReadSeqs>::Type
getMateSeqID(TReadSeqs const & readSeqs, TReadSeqID readSeqID)
{
    SEQAN_ASSERT_LT(readSeqID, getReadSeqsCount(readSeqs));

    typename Size<TReadSeqs>::Type pairID = getPairID(readSeqs, readSeqID);

    if (isFirstMate(readSeqs, readSeqID))
    {
        if (isFwdReadSeq(readSeqs, readSeqID))
            return getSecondMateRevSeqID(readSeqs, pairID);
        else
            return getSecondMateFwdSeqID(readSeqs, pairID);
    }
    else
    {
        if (isFwdReadSeq(readSeqs, readSeqID))
            return getFirstMateRevSeqID(readSeqs, pairID);
        else
            return getFirstMateFwdSeqID(readSeqs, pairID);
    }
}

#endif  // #ifndef APP_YARA_BITS_READS_H_
