// ==========================================================================
//                         Mason - A Read Simulator
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

#ifndef APPS_MASON2_GAP_PROJECTION_H_
#define APPS_MASON2_GAP_PROJECTION_H_

#include <seqan/align.h>
#include <seqan/bam_io.h>

namespace seqan {

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

// move it0 and it1 to point to the first non-gap character in it0 and return number of increments.
template <typename TGapsIter>
unsigned _moveToNextNonGap(TGapsIter & it0, TGapsIter & it1, TGapsIter const & it0End) {
    unsigned result = 0;

    do {
        ++result;
        ++it0;
        ++it1;
    } while (it0 != it0End && isGap(it0));

    return result;
}

// append CIGAR operation, no canceling out of I/D
template <typename TCigarString>
inline void _appendCigarOperation(TCigarString & cigar, CigarElement<> el) {
    if (!el.count)
        return;

    if (empty(cigar)) {
        appendValue(cigar, el);
    } else {
        if (back(cigar).operation == el.operation)
            back(cigar).count += el.count;
        else
            appendValue(cigar, el);
    }
}

/*!
 * @fn projectGapsTransitively
 * @brief Given two pairwise alignments of A to B and B to C, infer pairwise alignment of A to C.
 *
 * The sequence of <tt>outRefGaps</tt> and <tt>outReadGaps</tt> has to be set, clipping is added
 * from <tt>refGaps1</tt> and <tt>readGaps2</tt>.
 */
template <typename TGapsOut, typename TGapsIn1, typename TGapsIn2>
void projectGapsTransitively(
        TGapsOut & outRefGaps,
        TGapsOut & outReadGaps,
        TGapsIn1 const & refGaps1,
        TGapsIn1 const & haploGaps1,
        TGapsIn2 const & haploGaps2,
        TGapsIn2 const & readGaps2)
{
    String<CigarElement<> > cigar;

    typedef typename Iterator<TGapsOut, Standard>::Type TGapsIterOut;
    typedef typename Iterator<TGapsIn1 const, Standard>::Type TGapsIter1;
    typedef typename Iterator<TGapsIn2 const, Standard>::Type TGapsIter2;

    // Handle trailing and leading gaps.
    int leadingRefGaps1Count = countLeadingGaps(refGaps1);
    int trailingRefGaps1Count = countTrailingGaps(refGaps1);
    int leadingHaploGaps1Count = countLeadingGaps(haploGaps1);
    int trailingHaploGaps1Count = countTrailingGaps(haploGaps1);
    int leadingHaploGaps2Count = countLeadingGaps(haploGaps2);
    int trailingHaploGaps2Count = countTrailingGaps(haploGaps2);
    int leadingReadGaps2Count = countLeadingGaps(readGaps2);
    int trailingReadGaps2Count = countTrailingGaps(readGaps2);

    // std::cerr << "GAPS\t" << leadingRefGaps1Count << "\t"
    //           << trailingRefGaps1Count << "\t"
    //           << leadingHaploGaps1Count << "\t"
    //           << trailingHaploGaps1Count << "\t"
    //           << leadingHaploGaps2Count << "\t"
    //           << trailingHaploGaps2Count << "\t"
    //           << leadingReadGaps2Count << "\t"
    //           << trailingReadGaps2Count << "\n";

    // Below, we move the iterators right of the leading gaps, adjust cigar for this here.
    if (leadingHaploGaps1Count > 0)
        _appendCigarOperation(cigar, CigarElement<>('D', leadingHaploGaps1Count));
    if (leadingHaploGaps2Count > 0)
        _appendCigarOperation(cigar, CigarElement<>('I', leadingHaploGaps2Count));
    if (countGaps(refGaps1, 0) > countGaps(readGaps2, 0))
        _appendCigarOperation(cigar, CigarElement<>('I', countLeadingGaps(refGaps1) - countLeadingGaps(readGaps2)));
    if (countLeadingGaps(readGaps2) > countLeadingGaps(refGaps1))
        _appendCigarOperation(cigar, CigarElement<>('D', countLeadingGaps(readGaps2) - countLeadingGaps(refGaps1)));

    // Create end iterators of the four gaps.
    TGapsIter1 refGaps1End = end(refGaps1, Standard()) - std::max(trailingRefGaps1Count, trailingHaploGaps1Count);
    TGapsIter1 haploGaps1End = end(haploGaps1, Standard()) - std::max(trailingRefGaps1Count, trailingHaploGaps1Count);
    TGapsIter2 haploGaps2End = end(haploGaps2, Standard()) - std::max(trailingHaploGaps2Count, trailingReadGaps2Count);
    TGapsIter2 readGaps2End = end(readGaps2, Standard()) - std::max(trailingHaploGaps2Count, trailingReadGaps2Count);

    // Create iterators of the four gaps.
    int leadingShift = std::max(leadingRefGaps1Count, leadingReadGaps2Count);
    TGapsIter1 refGaps1Iter = begin(refGaps1, Standard()) + leadingShift + leadingHaploGaps1Count;
    TGapsIter1 haploGaps1Iter = begin(haploGaps1, Standard()) + leadingShift + leadingHaploGaps1Count;
    TGapsIter2 haploGaps2Iter = begin(haploGaps2, Standard()) + leadingShift + leadingHaploGaps2Count;
    TGapsIter2 readGaps2Iter = begin(readGaps2, Standard()) + leadingShift + leadingHaploGaps2Count;

    // These flags indicate whether the iterator aligns against a gap before the next round.
    bool aliGap1 = false, aliGap2 = false;

    while (refGaps1Iter != refGaps1End) {
        SEQAN_ASSERT(haploGaps1Iter != haploGaps1End);
        SEQAN_ASSERT(haploGaps2Iter != haploGaps2End);
        SEQAN_ASSERT(readGaps2Iter != readGaps2End);

        unsigned placed1 = _moveToNextNonGap(haploGaps1Iter, refGaps1Iter, haploGaps1End);
        SEQAN_ASSERT_GEQ(placed1, 1u);
        unsigned placed2 = _moveToNextNonGap(haploGaps2Iter, readGaps2Iter, haploGaps2End);
        SEQAN_ASSERT_GEQ(placed2, 1u);

        // std::cerr << "placed1=" << placed1 << "\tplaced2=" << placed2;
        // std::cerr << "\taliGap1=" << aliGap1 << "\taliGap2=" << aliGap2 << "\n";
        // std::cerr << "SEQS\t";
        // if (isGap(refGaps1Iter))
        //     std::cerr << "-";
        // else
        //     std::cerr << *refGaps1Iter;
        // std::cerr << " ";
        // if (isGap(haploGaps1Iter))
        //     std::cerr << "-";
        // else
        //     std::cerr << *haploGaps1Iter;
        // std::cerr << "\t";
        // if (isGap(haploGaps2Iter))
        //     std::cerr << "-";
        // else
        //     std::cerr << *haploGaps2Iter;
        // std::cerr << " ";
        // if (isGap(readGaps2Iter))
        //     std::cerr << "-";
        // else
        //     std::cerr << *readGaps2Iter;
        // std::cerr << "\n";

        // Place character for alignment.
        if (!aliGap1 && !aliGap2)
            _appendCigarOperation(cigar, CigarElement<>('M', 1));
        else if (!aliGap1 && aliGap2)
            _appendCigarOperation(cigar, CigarElement<>('D', 1));
        else if (aliGap1 && !aliGap2)
            _appendCigarOperation(cigar, CigarElement<>('I', 1));

        // std::cerr << "CIGAR\t";
        // for (unsigned i = 0; i < length(cigar); ++i)
        //     std::cerr << cigar[i].count << cigar[i].operation << " ";
        // std::cerr << "\n";

        // Add insertion/deletion characters.
        if (placed1 > 0u)
            _appendCigarOperation(cigar, CigarElement<>('D', placed1 - 1));
        if (placed2 > 0u)
            _appendCigarOperation(cigar, CigarElement<>('I', placed2 - 1));

        aliGap1 = !atEnd(refGaps1Iter) && isGap(refGaps1Iter);
        aliGap2 = !atEnd(readGaps2Iter) && isGap(readGaps2Iter);
    }

    // Adjust cigar for trailing gaps.
    if (trailingHaploGaps1Count > 0)
        _appendCigarOperation(cigar, CigarElement<>('D', trailingHaploGaps1Count));
    if (trailingHaploGaps2Count > 0)
        _appendCigarOperation(cigar, CigarElement<>('I', trailingHaploGaps2Count));
    if (countTrailingGaps(refGaps1) > countTrailingGaps(readGaps2))
        _appendCigarOperation(cigar, CigarElement<>('I', countTrailingGaps(refGaps1) - countTrailingGaps(readGaps2)));
    if (countTrailingGaps(readGaps2) > countTrailingGaps(refGaps1))
        _appendCigarOperation(cigar, CigarElement<>('D', countTrailingGaps(readGaps2) - countTrailingGaps(refGaps1)));

    // for (unsigned i = 0; i < length(cigar); ++i)
    //     std::cerr << cigar[i].count << cigar[i].operation << " ";
    // std::cerr << "\n";

    // Build resulting gaps.
    clearGaps(outRefGaps);
    clearClipping(outRefGaps);
    setClippedEndPosition(outRefGaps, endPosition(refGaps1));
    setClippedBeginPosition(outRefGaps, beginPosition(refGaps1));

    clearGaps(outReadGaps);
    clearClipping(outReadGaps);
    setClippedEndPosition(outReadGaps, endPosition(readGaps2));
    setClippedBeginPosition(outReadGaps, beginPosition(readGaps2));

    // std::cerr << "REF \t" << outRefGaps << "\n";
    // std::cerr << "READ\t" << outReadGaps << "\n";

    TGapsIterOut itRefOut = begin(outRefGaps, Standard());
    TGapsIterOut itReadOut = begin(outReadGaps, Standard());
    for (unsigned i = 0; i < length(cigar); ++i) {
        CigarElement<> const el = cigar[i];
        //std::cerr << el.count << el.operation;
        switch (el.operation) {
            case 'I':
                insertGaps(itRefOut, el.count);
                break;
            case 'D':
                insertGaps(itReadOut, el.count);
                break;
        }
        itRefOut += el.count;
        itReadOut += el.count;
    }
    //std::cerr << "\n";
}

}  // namespace seqan

#endif  // #ifndef APPS_MASON2_GAP_PROJECTION_H_
