// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STORE_ALIGN_H
#define SEQAN_HEADER_STORE_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Aligned Read Store
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class AlignedReadStoreElement
 * @headerfile <seqan/store.h>
 * @brief Represents an alignment of a read to a contig.
 *
 * @signature template <typename TPos, typename TGapAnchors[, typename TSpec]>
 *            struct AlignedReadStoreElement;
 *
 * @tparam TPos        The position type to use.
 * @tparam TGapAnchors Type of a read @link GapAnchor @endlink.
 * @tparam TSpec       The specializing type.  Default: <tt>void</tt>.
 *
 * Value type of the @link FragmentStore::alignedReadStore @endlink.  In contrast to all other @link FragmentStore
 * @endlink stores, the @link AlignedReadStoreElement::id @endlink of an aligned read is explicitely stored as a member
 * to allow for reordering the @link FragmentStore::alignedReadStore @endlink and still recover the id.
 */

/*!
 * @fn AlignedReadStoreElement::AlignedReadStoreElement
 * @brief Constructor.
 *
 * @signature AlignedReadStoreElement::AlignedReadStoreElement();
 * @signature AlignedReadStoreElement::AlignedReadStoreElement(id, readID, contigID, beginPos, endPos[, gaps]);
 *
 * @param[in] id       The alignment id refers to associated alignment information in @link
 *                     FragmentStore::alignQualityStore @endlink or @link FragmentStore::alignedReadTagStore @endlink.
 * @param[in] readID   Refers to the aligned read in the @link FragmentStore::readStore @endlink.
 * @param[in] contigID Refers to the contig in the @link FragmentStore::contigStore @endlink the read is aligned to.
 * @param[in] beginPos Begin position of the alignment in gap-space.
 * @param[in] endPos   End position of the alignment in gap-space.
 * @param[in] gaps     A @link String @endlink of @link GapAnchor @endlink objects.
 *
 * The default constructors sets all ids to <tt>AlignedReadStoreElement::INVALID_ID</tt> and @link
 * AlignedReadStoreElement::beginPos @endlink and @link AlignedReadStoreElement::endPos @endlink to 0.
 */

 /*!
 * @typedef AlignedReadStoreElement::TID
 * @signature typedef (..) TAlignedReadStoreElement::TID;
 * @brief Type of all store ids.
 *
 * <tt>TID</tt> is the result of <tt>ID&lg;AlignedReadStoreElement&lgt;&gt; &gt;::Type</tt>, see @link ID @endlink.
 *
 * @typedef AlignedReadStoreElement::TPos
 * @signature typedef (..) TAlignedReadStoreElement::TPos;
 * @brief Type of the @link AlignedReadStoreElement::beginPos @endlink and @link AlignedReadStoreElement::endPos
 *        @endlink.
 *
 * @typedef AlignedReadStoreElement::TGapAnchors
 * @signature typedef (..) TAlignedReadStoreElement::TGapAnchors;
 * @brief Type of the @link AlignedReadStoreElement::gaps @endlink member.
 *
 * @typedef AlignedReadStoreElement::TSpec
 * @signature typedef (..) TAlignedReadStoreElement::TSpec;
 * @brief The specializing type.
 */

/*!
 * @var TID AlignedReadStoreElement::INVALID_ID;
 * @brief Constant expression that is the value of an invalid id.
 *
 * @var TID AlignedReadStoreElement::id;
 * @brief The alignment id refers to associated alignment information in @link FragmentStore::alignQualityStore
 *        @endlink or @link FragmentStore::alignedReadTagStore @endlink.
 *
 * @var TID AlignedReadStoreElement::readID;
 * @brief Refers to the aligned read in @link FragmentStore::readStore @endlink.
 *
 * @var TID AlignedReadStoreElement::contigID;
 * @brief Refers to the contig in the @link FragmentStore::contigStore @endlink that the read is aligned with.
 *
 * @var TID AlignedReadStoreElement::pairMatchID;
 * @brief Two read alignments having the same pairMatchID form a valid pair match.  If it equals <tt>INVALID_ID</tt>
 *        then the read is either not paired or could not be aligned as part of a pair.
 *
 * @var TPos AlignedReadStoreElement::beginPos;
 * @brief Begin position in gap-space.
 *
 * If beginPos &lt; endPos then the read is aligned to the reverse strand where beginPos and endPos are the
 * corresponding positions on the forward strand.
 *
 * @var TPos AlignedReadStoreElement::endPos;
 * @brief End position in gap-space.
 *
 * If beginPos &lt; endPos then the read is aligned to the reverse strand where beginPos and endPos are the
 * corresponding positions on the forward strand.
 *
 * @var TGaps AlignedReadStoreElement::gaps;
 * @brief String of read @link GapAnchor @endlink objects.  Can be used to create a @link AnchorGaps @endlink
 *        alignment row.
 */

template <typename TPos_, typename TGapAnchor_, typename TSpec_ = void>
struct AlignedReadStoreElement
{
	typedef typename ID<AlignedReadStoreElement>::Type	TID;
	typedef TPos_										TPos;
	typedef TGapAnchor_									TGapAnchor;
	typedef TSpec_										TSpec;
	typedef String<TGapAnchor>							TGapAnchors;

	static const TID INVALID_ID;
	
	TID			id;
	TID			readID;
	TID			contigID;
	TID			pairMatchID;	// unique id. for multiple mate-pair matches (not matePairID)
	TPos		beginPos;		// begin position of the gapped sequence in gapped contig sequence
	TPos		endPos;			// end position of ..., for reverse aligned reads holds end < begin
	TGapAnchors	gaps;

	AlignedReadStoreElement() : id(INVALID_ID), readID(INVALID_ID), contigID(INVALID_ID), pairMatchID(INVALID_ID), beginPos(0), endPos(0) {}

	AlignedReadStoreElement(TID _id, TID _readID, TID _contigID, TPos _beginPos, TPos _endPos) :
		id(_id), 
		readID(_readID), 
		contigID(_contigID), 
		pairMatchID(INVALID_ID), 
		beginPos(_beginPos), 
		endPos(_endPos) {}

	AlignedReadStoreElement(TID _id, TID _readID, TID _contigID, TPos _beginPos, TPos _endPos, TGapAnchors const &_gaps) :
		id(_id), 
		readID(_readID), 
		contigID(_contigID), 
		pairMatchID(INVALID_ID), 
		beginPos(_beginPos), 
		endPos(_endPos),
		gaps(_gaps) {}

    inline bool operator==(AlignedReadStoreElement const & other) const
    {
        return id == other.id &&
                readID == other.readID &&
                contigID == other.contigID &&
                pairMatchID == other.pairMatchID &&
                beginPos == other.beginPos &&
                endPos == other.endPos &&
                gaps == other.gaps;
    }
};


// TODO(holtgrew): I find this useful for debugging purposes. Keep it?
template <typename TStream, typename TPos, typename TGapAnchor, typename TSpec>
TStream & operator<<(TStream & stream, AlignedReadStoreElement<TPos, TGapAnchor, TSpec> const & e) {
    return stream << "AlignedReadStore(id=" << e.id << ", readID=" << e.readID << ", contigID=" << e.contigID << ", pairMatchID=" << e.pairMatchID << ", beginPos=" << e.beginPos << ", endPos=" << e.endPos << ", {gaps})";
}


//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec> 
const typename ID<AlignedReadStoreElement<TPos, TGapAnchor, TSpec> >::Type 
AlignedReadStoreElement<TPos, TGapAnchor, TSpec>::INVALID_ID = MaxValue<typename ID<AlignedReadStoreElement<TPos, TGapAnchor, TSpec> >::Type>::VALUE; 

//////////////////////////////////////////////////////////////////////////////

/*!
 * @class AlignQualityStoreElement
 * @headerfile <seqan/store.h>
 * @brief Stores the alignment qualities.
 *
 * @signature template <typename TScore[, typename TSpec]>
 *            struct AlignQualityStoreElement;
 *
 * @tparam TScore Type to store align and pair score values.
 * @tparam TSpec  The specializing type.  Default: <tt>void</tt>.
 *
 * Value type of @link FragmentStore::alignQualityStore @endlink string.
 *
 *
 * @fn AlignQualityStoreElement::AlignQualityStoreElement
 * @brief Constructor, sets all members to 0.
 *
 * @signature AlignQualityStoreElement::AlignQualityStoreElement();
 * 
 *
 * @var TScore AlignQualityStoreElement::pairScore;
 * @brief Combined score of both alignmetns and pair match.
 *
 * @var TScore AlignQualityStoreElement::score;
 * @brief Score of the alignment.
 *
 * @var TCount AlignQualityStoreElement::errors;
 * @brief Absolute number of errors in the alignment (<tt>unsigned char</tt>).
 */

template <typename TScore, typename TSpec = void>
struct AlignQualityStoreElement
{
	TScore				pairScore;		// score of the mate-pair alignment (this read is part of)
	TScore				score;			// score of the single read alignment
	unsigned char		errors;			// absolute number of errors (Hamming or edit distance)
	
	AlignQualityStoreElement():
		pairScore(0),
		score(0),
		errors(MaxValue<unsigned char>::VALUE) {}

	AlignQualityStoreElement(TScore _pairScore, TScore _score, unsigned char _errors):
		pairScore(_pairScore),
		score(_score),
		errors(_errors) {}

    inline bool operator==(AlignQualityStoreElement const & other)
    {
        return pairScore == other.pairScore &&
                score == other.score &&
                errors == other.errors;
    }
};


//////////////////////////////////////////////////////////////////////////////
// Sorting tags
//////////////////////////////////////////////////////////////////////////////

/*!
 * @defgroup SortAlignedReadTags Tags for sortAlignedReads
 * @brief Tags to select a specific field to sort the @link FragmentStore::alignedReadStore @endlink by.
 *
 * @see sortAlignedReads
 * @see lowerBoundAlignedReads
 * @see upperBoundAlignedReads
 *
 * @tag SortAlignedReadTags#SortContigID
 * @headerfile <seqan/store.h>
 * @brief Sort aligned reads by @link AlignedReadStoreElement::contigID @endlink.
 *
 * @signature typedef Tag<SortContigID_> const SortContigID;
 *
 * @tag SortAlignedReadTags#SortID
 * @headerfile <seqan/store.h>
 * @brief Sort aligned reads by @link AlignedReadStoreElement::id @endlink.
 *
 * @signature typedef Tag<SortID_> const SortID;
 *
 * @tag SortAlignedReadTags#SortBeginPos
 * @headerfile <seqan/store.h>
 * @brief Sort aligned reads by @link AlignedReadStoreElement::beginPos @endlink.
 *
 * @signature typedef Tag<SortBeginPos_> const SortBeginPos;
 *
 * @tag SortAlignedReadTags#SortEndPos
 * @headerfile <seqan/store.h>
 * @brief Sort aligned reads by @link AlignedReadStoreElement::endPos @endlink.
 *
 * @signature typedef Tag<SortEndPos_> const SortEndPos;
 *
 * @tag SortAlignedReadTags#SortPairMatchID
 * @headerfile <seqan/store.h>
 * @brief Sort aligned reads by @link AlignedReadStoreElement::pairMatchID @endlink.
 *
 * @signature typedef Tag<SortPairMatchID_> const SortPairMatchID;
 *
 * @tag SortAlignedReadTags#SortReadID
 * @headerfile <seqan/store.h>
 * @brief Sort aligned reads by @link AlignedReadStoreElement::readID @endlink.
 *
 * @signature typedef Tag<SortReadID_> const SortReadID;
 */

struct SortContigID_;
typedef Tag<SortContigID_> const SortContigID;

struct SortID_;
typedef Tag<SortID_> const SortID;

struct SortBeginPos_;
typedef Tag<SortBeginPos_> const SortBeginPos;

struct SortEndPos_;
typedef Tag<SortEndPos_> const SortEndPos;

struct SortPairMatchID_;
typedef Tag<SortPairMatchID_> const SortPairMatchID;

struct SortReadID_;
typedef Tag<SortReadID_> const SortReadID;


//////////////////////////////////////////////////////////////////////////////
// Sorting functors
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec>
inline void swap(
    AlignedReadStoreElement<TPos, TGapAnchor, TSpec> &a,
    AlignedReadStoreElement<TPos, TGapAnchor, TSpec> &b)
{
    std::swap(a.id, b.id);
    std::swap(a.readID, b.readID);
    std::swap(a.contigID, b.contigID);
    std::swap(a.pairMatchID, b.pairMatchID);
    std::swap(a.beginPos, b.beginPos);
    std::swap(a.endPos, b.endPos);
    swap(a.gaps, b.gaps);
}

template <typename TAlignedRead, typename TTag>
struct _LessAlignedRead;

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortID> :
	public std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return (a1.id) < (a2.id);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortContigID> :
	public std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.contigID < a2.contigID;
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortBeginPos> :
	public std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return _min(a1.beginPos, a1.endPos) < _min(a2.beginPos, a2.endPos);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortEndPos> :
	public std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return _max(a1.beginPos, a1.endPos) < _max(a2.beginPos, a2.endPos);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortPairMatchID> :
	public std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.pairMatchID < a2.pairMatchID;
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortReadID> :
	public std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.readID < a2.readID;
	}
};

//////////////////////////////////////////////////////////////////////////////
// Sorting function
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): Maybe add equalRangeAlignedReads?

/*!
 * @fn sortAlignedReads
 * @headerfile <seqan/store.h>
 * @brief Stably read alignments, e.g. in @link FragmentStore::alignedReadStore @endlink.
 *
 * @signature void sortAlignedReads(alignStore, sortTag);
 * @signature void sortAlignedReads(alignStore, lessFunctor);
 *
 * @param[in,out] alignStore  The @link ContainerConcept sequence @endlink of @link AlignedReadStoreElement @endlink
 *                            to be sorted, e.g. @link FragmentStore::alignedReadStore @endlink.
 * @param[in]     sortTag     Tag for selecting the member to sort by.  See @link SortAlignedReadTags @endlink.
 * @param[in]     lessFunctor A functor to pass to <tt>std::stable_sort</tt> for sorting the sequence.
 *
 * @see SortAlignedReadTags
 * @see lowerBoundAlignedReads
 * @see upperBoundAlignedReads
 */

/*!
 * @fn lowerBoundAlignedReads
 * @headerfile <seqan/store.h>
 * @brief Performs a binary lower bound search on read alignments.
 *
 * @signature TIter1 lowerBoundAlignedReads(alignStore, value, sortTag);
 * @signature TIter2 lowerBoundAlignedReads(itBegin, itEnd, value, sortTag);
 *
 * @param[in,out] alignStore  The @link ContainerConcept sequence @endlink of @link AlignedReadStoreElement @endlink
 *                            to be searched, e.g. @link FragmentStore::alignedReadStore @endlink.
 * @param[in]     itBegin     Iterator to the begin of the sequence to search.
 * @param[in]     itEnd       Iterator to the end of the sequence to search.
 * @param[in]     value       The value to search for.
 * @param[in]     sortTag     Tag for selecting the member to compare by.  See @link SortAlignedReadTags @endlink.
 *
 * @return TIter1 Iterator to the lower bound item.  If <tt>TAlignStore</tt> is the type of <tt>alignStore</tt> then
 *                TIter1 is the result of <tt>Iterator&lt;TAlignStore, Standard&gt;::Type</tt>.
 * @return TIter2 Iterator to the lower bound item.  Has the same type as <tt>itBegin</tt> and <tt>itEnd</tt>.
 *
 * @see SortAlignedReadTags
 * @see sortAlignedReads
 * @see upperBoundAlignedReads
 */

/*!
 * @fn upperBoundAlignedReads
 * @headerfile <seqan/store.h>
 * @brief Performs a binary upper bound search on read alignments.
 *
 * @signature TIter1 upperBoundAlignedReads(alignStore, value, sortTag);
 * @signature TIter2 upperBoundAlignedReads(itBegin, itEnd, value, sortTag);
 *
 * @param[in,out] alignStore  The @link ContainerConcept sequence @endlink of @link AlignedReadStoreElement @endlink
 *                            to be searched, e.g. @link FragmentStore::alignedReadStore @endlink.
 * @param[in]     itBegin     Iterator to the begin of the sequence to search.
 * @param[in]     itEnd       Iterator to the end of the sequence to search.
 * @param[in]     value       The value to search for.
 * @param[in]     sortTag     Tag for selecting the member to compare by.  See @link SortAlignedReadTags @endlink.
 *
 * @return TIter1 Iterator to the upper bound item.  If <tt>TAlignStore</tt> is the type of <tt>alignStore</tt> then
 *                TIter1 is the result of <tt>Iterator&lt;TAlignStore, Standard&gt;::Type</tt>.
 * @return TIter2 Iterator to the upper bound item.  Has the same type as <tt>itBegin</tt> and <tt>itEnd</tt>.
 *
 * @see SortAlignedReadTags
 * @see sortAlignedReads
 * @see lowerBoundAlignedReads
 */

template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign& alignStore, Tag<TSortSpec> const &) 
{
	std::stable_sort(
		begin(alignStore, Standard() ), 
		end(alignStore, Standard() ), 
		_LessAlignedRead<typename Value<TAlign>::Type, Tag<TSortSpec> const>() );
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign const & alignStore, Tag<TSortSpec> const &) 
{
	std::stable_sort(
		begin(const_cast<TAlign&>(alignStore), Standard() ), 
		end(const_cast<TAlign&>(alignStore), Standard() ), 
		_LessAlignedRead<typename Value<TAlign>::Type, Tag<TSortSpec> const>() );
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign & alignStore, TFunctorLess const &less) 
{
	std::stable_sort(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		less);
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign const & alignStore, TFunctorLess const &less) 
{
	std::stable_sort(
		begin(const_cast<TAlign&>(alignStore), Standard()), 
		end(const_cast<TAlign&>(alignStore), Standard()), 
		less);
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
lowerBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortID const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
upperBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortID const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
lowerBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortContigID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigID = val;
	return std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigID const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortContigID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigID = val;
	return std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
upperBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortContigID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigID = val;
	return std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigID const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortContigID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigID = val;
	return std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
lowerBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
upperBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
lowerBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
upperBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
lowerBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortPairMatchID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchID = val;
	return std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchID const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortPairMatchID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchID = val;
	return std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
upperBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortPairMatchID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchID = val;
	return std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchID const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortPairMatchID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchID = val;
	return std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
lowerBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortReadID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readID = val;
	return std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadID const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortReadID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readID = val;
	return std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
upperBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortReadID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readID = val;
	return std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadID const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortReadID) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readID = val;
	return std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadID const>() );
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortID) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortID const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortID) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortContigID) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.contigID = val;
	return std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortContigID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortContigID) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.contigID = val;
	return std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortContigID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortBeginPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortBeginPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortEndPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortEndPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortPairMatchID) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.pairMatchID = val;
	return std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortPairMatchID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortPairMatchID) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.pairMatchID = val;
	return std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortPairMatchID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortReadID) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.readID = val;
	return std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortReadID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortReadID) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.readID = val;
	return std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortReadID const>() );
}

//////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortID) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortID) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortContigID) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.contigID = val;
	return std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortContigID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortContigID) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.contigID = val;
	return std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortContigID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortBeginPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortBeginPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortEndPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortEndPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortPairMatchID) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchID = val;
	return std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortPairMatchID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortPairMatchID) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchID = val;
	return std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortPairMatchID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortReadID) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.readID = val;
	return std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortReadID const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortReadID) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.readID = val;
	return std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortReadID const>() );
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
