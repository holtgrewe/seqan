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

#ifndef SEQAN_HEADER_STORE_CONTIG_H
#define SEQAN_HEADER_STORE_CONTIG_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Contig Store
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class ContigStoreElement
 * @headerfile <seqan/store.h>
 * @brief Represents a single contig.
 *
 * @signature template <typename TContigSeq, typename TGapAnchor[, typename TSpec]>
 *            class ContigStoreElement;
 *
 * @tparam TContigSeq Type to store the contig sequence.
 * @tparam TGapAnchor Type of the @link GapAnchor @endlink.
 * @tparam TSpec      Specializing type.  Default: <tt>void</tt>.
 *
 * Value type of the @link FragmentStore::contigStore @endlink string.
 *
 *
 * @fn ContigStoreElement::ContigStoreElement
 * @brief Constructor.
 *
 * @signature ContigStoreElement::ContigStoreElement();
 *
 * Set fileID to INVALID_ID and usage, fileBeginPos, fileEndPos to 0.
 */

/*!
 * @typedef ContigStoreElement::TContigSeq
 * @brief Type fo the seq member.
 *
 * @signature typedef (...) ContigStoreElement::TContigSeq;
 *
 * @typedef ContigStoreElement::TGapAnchors
 * @brief Type of the gaps member.
 *
 * @signature typedef (...) ContigStoreElement::TGapAnchors;
 *
 *
 * @typedef ContigStoreElement::TPos
 * @brief Type of the fileBeginPos and fileEndPos emmbers.
 *
 * @signature typedef (...) ContigStoreElement::TPos;
 *
 *
 * @typedef ContigStoreElement::TSpec
 * @brief The specializing type.
 *
 * @signature typedef (...) ContigStoreElement::TSpec;
 */

/*!
 * @var TContigSeq ContigStoreElement::seq;
 * @brief Contig sequence.
 *
 * @var TGaps ContigStoreElement::gaps;
 * @brief String of contig @link GapAnchor @endlink objects.  CAn be used to create a @link AnchorGaps @endlink
 *        alignment row.
 *
 * @var unsigned ContigStoreElement::usage;
 * @brief Count the number of locks.
 *
 * @see FragmentStore#lockContig
 *
 * @var TID ContigStoreElement::fileID;
 * @brief Refers to a file in the @link FragmentStore::contigFileStore @endlink or is INVALID_ID if the contig has no
 *        file association.
 *
 * @var TPos ContigStoreElement::fileBeginPos
 * @brief Begin position of the contig sequence fragment in the file.
 *
 * @var TPos ContigStoreElement::fileEndPos
 * @brief End position of the contig sequence fragment in the file.
 *
 * @var TID ContigStoreElement::INVALID_ID;
 * @brief Constant to represent an invalid.
 */

template <typename TContigSeq_, typename TGapAnchor_, typename TSpec_ = void>
struct ContigStoreElement
{
	typedef typename ID<ContigStoreElement>::Type	TID;
	
	typedef TContigSeq_			TContigSeq;
	typedef TGapAnchor_			TGapAnchor;
	typedef TSpec_				TSpec;
	typedef __int64				TPos;
	typedef String<TGapAnchor>	TGapAnchors;

	static const TID INVALID_ID;

	TContigSeq	seq;
	TGapAnchors	gaps;
	
// dynamic loading and disposing of contigs
	unsigned	usage;			// number of threads,... using this contig
	TID			fileID;
	TPos		fileBeginPos;
	TPos		fileEndPos;

	ContigStoreElement() : usage(0), fileID(INVALID_ID), fileBeginPos(0), fileEndPos(0) {}

    inline bool operator==(ContigStoreElement const & other) const
    {
        return usage == other.usage &&
                fileID == other.fileID &&
                fileBeginPos == other.fileBeginPos &&
                fileEndPos == other.fileEndPos &&
                seq == other.seq &&
                gaps == other.gaps;
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename TContigSeq_, typename TGapAnchor_, typename TSpec_> 
const typename ID<ContigStoreElement<TContigSeq_, TGapAnchor_, TSpec_> >::Type 
ContigStoreElement<TContigSeq_, TGapAnchor_, TSpec_>::INVALID_ID = MaxValue<typename ID<ContigStoreElement<TContigSeq_, TGapAnchor_, TSpec_> >::Type>::VALUE; 

//////////////////////////////////////////////////////////////////////////////

/*!
 * @class ContigFile
 * @headerfile <seqan/store.h>
 * @brief Represents a file containing contigs.
 *
 * @signature template <[typename TSpec]>
 *            class ContigFile;
 *
 * @tparam TSpec The specializing type.  Default: <tt>void</tt>.
 *
 * Value type of @link FragmentStore::contigFileStore @endlink string.
 */

/*!
 * @var CharString ContigFile::fileName;
 * @brief Contig file name.
 *
 * @var AutoSeqFormat ContigFile::format;
 * @brief Stores the contig file format, auto-detect in @link FragmentStore#loadContigs @endlink.
 *
 * @var TID ContigFile::firstContigID;
 * @brief The contigID of the first sequence in the file.  Subsequent contig sequences have an increasing contigID.
 */

template <typename TSpec_ = void>
struct ContigFile
{
//IOREV instead of storing filename and format store TFile ?
	typedef typename ID<ContigFile>::Type	TID;

	static const TID INVALID_ID;

	CharString		fileName;
	AutoSeqFormat	format;
	TID				firstContigID;	// first sequence of the file corresponds to this contigID

    inline bool operator==(ContigFile const & other) const
    {
        return fileName == other.fileName &&
                format == other.format &&
                firstContigID == other.firstContigID;
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec_> 
const typename ID<ContigFile<TSpec_> >::Type 
ContigFile<TSpec_>::INVALID_ID = MaxValue<typename ID<ContigFile<TSpec_> >::Type>::VALUE; 

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

//  LocalWords:  ContigFile
