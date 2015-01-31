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

#ifndef APP_YARA_BITS_SEEDS_H_
#define APP_YARA_BITS_SEEDS_H_

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function getHostPos()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TPos>
inline typename StringSetPosition<THost>::Type
getHostPos(StringSet<THost, Segment<TSpec> > const & me, TPos pos)
{
    return me.positions[pos];
}

// ----------------------------------------------------------------------------
// Function getSeedIDs()
// ----------------------------------------------------------------------------
// getSeqNosByHostSeqNo(seqNo)

template <typename THost, typename TSpec, typename TReadID>
inline Pair<typename ID<StringSet<THost, Segment<TSpec> > const>::Type>
getSeedIDs(StringSet<THost, Segment<TSpec> > const & seeds, TReadID readID)
{
    typedef StringSet<THost, Segment<TSpec> > const         TStringSet;
    typedef typename ID<TStringSet>::Type                   TID;
    typedef typename StringSetPositions<THost>::Type const  TPositions;
    typedef typename Iterator<TPositions, Standard>::Type   TPositionsIter;
    typedef typename Value<TPositions>::Type                TPos;

    TPositionsIter posBegin = begin(seeds.positions, Standard());
    TPositionsIter posEnd = end(seeds.positions, Standard());
    TPos pos(readID, 0);

    TPositionsIter seedsBegin = std::lower_bound(posBegin, posEnd, TPos(readID, 0));
    TPositionsIter seedsEnd = std::lower_bound(posBegin, posEnd, TPos(readID + 1, 0));

    return Pair<TID>(position(seedsBegin, seeds.positions), position(seedsEnd, seeds.positions));
}

// ----------------------------------------------------------------------------
// Function getReadSeqID()
// ----------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TSeedID>
inline typename ID<THost>::Type
getReadSeqID(StringSet<THost, Segment<TSpec> > const & seeds, TSeedID seedID)
{
    return getSeqNo(getHostPos(seeds, seedID));
}

// ----------------------------------------------------------------------------
// Function getPosInRead()
// ----------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TSeedID>
inline Pair<typename Position<typename Value<THost>::Type>::Type>
getPosInRead(StringSet<THost, Segment<TSpec> > const & seeds, TSeedID seedID)
{
    typedef typename Position<typename Value<THost>::Type>::Type TPos;

    TPos seedBegin = getSeqOffset(getHostPos(seeds, seedID));
    TPos seedEnd = seedBegin + length(value(seeds, seedID));

    return Pair<TPos>(seedBegin, seedEnd);
}

#endif  // #ifndef APP_YARA_BITS_SEEDS_H_
