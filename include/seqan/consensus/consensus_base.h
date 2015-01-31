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
// Author: Tobias Rausch <rausch@embl.de>
// ==========================================================================

#ifndef SEQAN_HEADER_CONSENSUS_BASE_H
#define SEQAN_HEADER_CONSENSUS_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Segment Match Generation tag
//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): This apparently belongs into graph_msa?

struct OverlapLibrary_;
typedef Tag<OverlapLibrary_> const OverlapLibrary;



//////////////////////////////////////////////////////////////////////////////
// Consensus tag
//////////////////////////////////////////////////////////////////////////////

/*
 * @defgroup ConsensusCallingTags Consensus Calling Tags
 * @brief Tags for consensus calling.
 *
 * @tag ConsensusCallingTags#MajorityVote
 * @headerfile <seqan/consensus.h>
 * @brief Consensus based on the most common character.
 *
 * @tag ConsensusCallingTags#Bayesian
 * @headerfile <seqan/consensus.h>
 * @brief Consensus based on a bayesian probability.
 */


struct MajorityVote_;
typedef Tag<MajorityVote_> const MajorityVote;

struct Bayesian_;
typedef Tag<Bayesian_> const Bayesian;



//////////////////////////////////////////////////////////////////////////////
// Read alignment and Consensus Generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////

struct ConsensusOptions {
public:
	// Method
	// 0: graph-based multiple sequence alignment
	// 1: realign
	int method;

	// ReAlign Method
	// 0: Needleman-Wunsch
	// 1: Gotoh
	int rmethod;

	// Bandwidth of overlap alignment
	int bandwidth;

	// Number of computed overlaps per read (at the beginning and end of a read)
	int overlaps;

	// Minimum match length of a computed overlap
	int matchlength;

	// Minimum quality (in percent identity) of a computed overlap
	int quality;

	// Window size, only relevant for insert sequencing
	// If window == 0, no insert sequencing is assumed
	int window;
	
	// Output
	// 0: seqan style
	// 1: afg output format
	// 2: frg output format
	// 3: cgb output format
    // 4: Sam output format
	int output;

	// Multi-read alignment
	bool noalign;

	// Offset all reads, so the first read starts at position 0
	bool moveToFront;

	// Include reference genome
	bool include;

	// Scoring object for overlap alignments
	Score<int> sc;

	// Various input and output files
	std::string readsfile;				// File of reads in FASTA format
	std::string afgfile;				// AMOS afg file input
	std::string samfile;				// Sam file input
	std::string contigsfile;			// FASTA reference file for Sam input
	std::string outfile;				// Output file name
	
	// Initialization
	ConsensusOptions() :
	    method(0), rmethod(0), bandwidth(0), overlaps(0), matchlength(0), quality(0),
	    window(0), output(0), noalign(false), moveToFront(0), include(0), sc(2, -6, -4, -9)
	{}
};


//////////////////////////////////////////////////////////////////////////////////

// Copy reads for the whole contig out of fragStore and into strSet.  The start and end
// positions of the alignments go into startEndPos.

template<typename TValue, typename TStrSpec, typename TPosPair, typename TStringSpec, typename TSpec, typename TConfig, typename TID>
inline void 
_loadContigReads(StringSet<TValue, Owner<TStrSpec> > & strSet,
				 String<TPosPair, TStringSpec> & startEndPos,
				 FragmentStore<TSpec, TConfig> const & fragStore,
				 TID const contigID)
{
	typedef FragmentStore<TSpec, TConfig> const TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadPos TReadPos;


	// Sort aligned reads according to contig id
	sortAlignedReads(fragStore.alignedReadStore, SortContigID());
	resize(strSet, length(fragStore.alignedReadStore));

	// Retrieve all reads, limit them to the clear range and if required reverse complement them
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore const>::Type TAlignIter;
	TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigID, SortContigID());
	TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigID, SortContigID());
	TSize numRead = 0;
	TReadPos begClr = 0;
	TReadPos endClr = 0;
	TSize lenRead = 0;
	TSize offset = 0;
	for(;alignIt != alignItEnd; ++alignIt) {
		offset = _min(alignIt->beginPos, alignIt->endPos);
		getClrRange(fragStore, *alignIt, begClr, endClr);
		strSet[numRead] = infix(fragStore.readSeqStore[alignIt->readID], begClr, endClr);
		lenRead = endClr - begClr;
		if (alignIt->beginPos < alignIt->endPos) appendValue(startEndPos, TPosPair(offset, offset + lenRead), Generous());
		else {
			reverseComplement(strSet[numRead]);
			appendValue(startEndPos, TPosPair(offset + lenRead, offset), Generous());
		}
		++numRead;
	}
	resize(strSet, numRead, Exact());
}




//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TMatrix, typename TSize2, typename TSize, typename TReadSlot> 
inline bool
convertAlignment(FragmentStore<TSpec, TConfig>& fragStore,
				 TMatrix& mat,
				 TSize2 contigID,
				 TSize& coverage,
				 TReadSlot& slot)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename TFragmentStore::TContigPos TContigPos;

	// Gap char
	TValue gapChar = gapValue<TValue>();

	// Sort according to contigID
	sortAlignedReads(fragStore.alignedReadStore, SortContigID());
	
	// Find range of the given contig
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore, Standard>::Type TAlignIter;
	TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigID, SortContigID());
	TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigID, SortContigID());

	// Sort the reads according to the begin position
	sortAlignedReads(infix(fragStore.alignedReadStore, alignIt - begin(fragStore.alignedReadStore, Standard()), alignItEnd - begin(fragStore.alignedReadStore, Standard())), SortBeginPos());
	TAlignIter alignItBegin = alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigID, SortContigID());
	alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigID, SortContigID());

	// Get the maximum coverage and the slot for each read
	typedef String<TSize> TFirstFreePos;
	typedef typename Iterator<TFirstFreePos, Standard>::Type TPosIter;
	TFirstFreePos freePos;
	TSize pos = 0;
	TSize maxTmp = 0;
	TSize numCol = 0;
	reserve(slot, alignItEnd - alignIt);
	for(;alignIt != alignItEnd; ++alignIt) {
		TPosIter itPos = begin(freePos, Standard());
		TPosIter itPosEnd = end(freePos, Standard());
		pos = 0;
		for(;itPos != itPosEnd; ++itPos, ++pos) 
			if ((TContigPos)*itPos < _min(alignIt->beginPos, alignIt->endPos)) break;
		if (pos + 1 > length(freePos)) resize(freePos, pos+1, Generous());
		maxTmp = _max(alignIt->beginPos, alignIt->endPos);
		freePos[pos] = maxTmp;
		if (maxTmp > numCol) numCol = maxTmp;
		appendValue(slot, pos);
	}
	coverage = length(freePos);
	clear(freePos);

	// Fill the matrix
	typedef typename Iterator<TMatrix, Standard>::Type TMatIter;
	resize(mat, coverage * numCol, '.');
	alignIt = alignItBegin;
	TSize readPos = 0;
	TMatIter matIt = begin(mat, Standard());
	typename TFragmentStore::TReadSeq myRead;
	for(;alignIt != alignItEnd; ++alignIt, ++readPos) {
		typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor>, Standard>::Type TReadGapsIter;
		TReadGapsIter itGaps = begin(alignIt->gaps, Standard());
		TReadGapsIter itGapsEnd = end(alignIt->gaps, Standard());
		
		// Place each read inside the matrix
		myRead = fragStore.readSeqStore[alignIt->readID];
		TSize lenRead = length(myRead);
		TSize offset = alignIt->beginPos;
		if (alignIt->beginPos > alignIt->endPos) {
			reverseComplement(myRead);
			offset = alignIt->endPos;
		}
		matIt = begin(mat, Standard());
		matIt += (slot[readPos] * numCol + offset);

		typedef typename Iterator<typename TFragmentStore::TReadSeq, Standard>::Type TReadIter;
		TReadIter seqReadIt = begin(myRead, Standard());

		// First clear range
		TSize mySeqPos = 0;
		int diff = 0;
		if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) {
			mySeqPos = itGaps->seqPos;
			diff = -1 * mySeqPos;
			seqReadIt += mySeqPos;
		}
		TSize clr2 = lenRead;
		TSize stop = 0;
		for(;itGaps != itGapsEnd; ++itGaps) {
			// Any clipped sequence at the end
			stop =  itGaps->seqPos;
			if (diff - ((int) itGaps->gapPos - (int) itGaps->seqPos) > 0) 
				clr2 = stop = lenRead - (diff - ((int) itGaps->gapPos - (int) itGaps->seqPos));
			
			for(;mySeqPos < stop; ++matIt, ++seqReadIt, ++mySeqPos) 
				*matIt = *seqReadIt;

			for(int i = 0; i < ((int) itGaps->gapPos - (int) itGaps->seqPos) - diff; ++i, ++matIt) 
				*matIt = gapChar;
	
			diff = (itGaps->gapPos - itGaps->seqPos);
		}
		for(;mySeqPos < clr2; ++mySeqPos, ++seqReadIt, ++matIt) 
			*matIt = *seqReadIt;
	}
	//for(TSize row = 0; row < coverage; ++row) {
	//	for(TSize col = 0; col<numCol; ++col) {
	//		std::cout << mat[row * numCol + col];
	//	}
	//	std::cout << std::endl;
	//}
	return true;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TMatrix, typename TSize2, typename TSize> 
inline bool
convertAlignment(FragmentStore<TSpec, TConfig>& fragStore,
				 TMatrix& mat,
				 TSize2 contigID,
				 TSize& coverage)
{
	String<TSize> slot;
	return convertAlignment(fragStore, mat, contigID, coverage, slot);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TMatrix> 
inline bool
convertAlignment(FragmentStore<TSpec, TConfig>& fragStore,
				 TMatrix& mat)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	TSize coverage;
	return convertAlignment(fragStore, mat, 0, coverage);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TGappedConsensus, typename TSize> 
inline void
getGappedConsensus(FragmentStore<TSpec, TConfig>& fragStore,
				   TGappedConsensus& gappedConsensus,
				   TSize contigID)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<TGappedConsensus>::Type TValue;
	typedef typename TFragmentStore::TContigPos TContigPos;
	
	TValue gapChar = gapValue<TValue>();
	typedef typename Iterator<typename TFragmentStore::TContigSeq, Standard>::Type TContigIter;
	TContigIter seqContigIt = begin(fragStore.contigStore[contigID].seq, Standard());
	TContigIter seqContigItEnd = end(fragStore.contigStore[contigID].seq, Standard());
	typedef typename Iterator<String<typename TFragmentStore::TContigGapAnchor>, Standard>::Type TGapsIter;
	TGapsIter itGaps = begin(fragStore.contigStore[contigID].gaps, Standard());
	TGapsIter itGapsEnd = end(fragStore.contigStore[contigID].gaps, Standard());
	int diff = 0;
	TContigPos mySeqPos = 0;
	for(;itGaps != itGapsEnd; goNext(itGaps)) {
		for(;mySeqPos < itGaps->seqPos; ++seqContigIt, ++mySeqPos) 
			appendValue(gappedConsensus, *seqContigIt, Generous());
			
		for(int i = 0; i < ((int) itGaps->gapPos - (int) itGaps->seqPos) - diff; ++i) 
			appendValue(gappedConsensus, gapChar, Generous());
			diff = (itGaps->gapPos - itGaps->seqPos);
	}
	for(;seqContigIt != seqContigItEnd; ++seqContigIt) 
		appendValue(gappedConsensus, *seqContigIt, Generous());
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TGappedConsensus, typename TSize> 
inline void
assignGappedConsensus(FragmentStore<TSpec, TConfig>& fragStore,
					  TGappedConsensus& gappedCons,
					  TSize contigID)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<TGappedConsensus>::Type TValue;
	TValue gapChar = gapValue<TValue>();

	// Update the contig
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigStoreElement;
	TContigStoreElement& contigEl = fragStore.contigStore[contigID];
	clear(contigEl.gaps);
	clear(contigEl.seq);

	// Create the sequence and the gap anchors
	typedef typename Iterator<TGappedConsensus, Standard>::Type TStringIter;
	TStringIter seqIt = begin(gappedCons, Standard());
	TStringIter seqItEnd = end(gappedCons, Standard());
	typedef typename TFragmentStore::TReadPos TReadPos;
	typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
	TReadPos ungappedPos = 0;
	TReadPos gappedPos = 0;
	bool gapOpen = false;
	for(;seqIt != seqItEnd; ++seqIt, ++gappedPos) {
		if (*seqIt == gapChar) gapOpen = true;				
		else {
			if (gapOpen) {
				appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());
				gapOpen = false;
			}
			Dna5Q letter = *seqIt;
			assignQualityValue(letter, 'D');
			appendValue(contigEl.seq, letter);
			++ungappedPos;
		}
	}
	if (gapOpen) 
		appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());
}



//////////////////////////////////////////////////////////////////////////////////

/*
 * @fn consensusAlignment
 * @headerfile <seqan/consensus.h>
 * @brief Compute consensus alignment.
 * 
 * @signature void consensusAlignment(alignmentGraph, beginEndPos[, options])
 * 
 * @param[out] alignmentGraph  Alignment graph to build.
 * @param[in]  options         Optional settings for the consenus alignment, type: <tt>ConsensusOptions</tt>.
 * @param[in]  beginEndPos     Interval start and end position for the read's alignment; <tt>String&lt;Pair&lt;TPos, TPos&gt; &gt;</tt>.
 * 
 * @section Example
 * 
 * @code{.cpp}
 * #include <seqan/sequence.h>
 * #include <seqan/graph_align.h>
 * #include <seqan/consensus.h>
 *  
 * int main()
 * {
 *     using namespace seqan;
 *  
 *     typedef StringSet<Dna5String> TStringSet;
 *     typedef Graph<Alignment<TStringSet, void, WithoutEdgeID> > TAlignGraph;
 *  
 *     TStringSet readSet;
 *     String<Pair<TSize> > begEndPos;
 *  
 *     appendValue(readSet, "CCCAGTGA");
 *     appendValue(begEndPos, Pair<TSize>(0, 5));
 *     appendValue(readSet, "AGGGACTGT");
 *     appendValue(begEndPos, Pair<TSize>(3, 9));
 *  
 *     TAlignGraph alignmentGraph(readSet);
 *     consensusAlignment(alignmentGraph, begEndPos);
 *  
 *     return 0;
 * }
 * @endcode
 */

template<typename TStringSet, typename TCargo, typename TSpec, typename TSize, typename TConfigOptions>
inline void
consensusAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gOut,
				   String<Pair<TSize, TSize> >& begEndPos,
				   TConfigOptions const& consOpt) 
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TOutGraph;
	typedef typename ID<TOutGraph>::Type TID;

	// Initialization
	TStringSet& seqSet = stringSet(gOut);

	// Select all overlapping reads and record the diagonals of the band
	String<Pair<TID, TID> > pList;
	String<Pair<int, int> > diagList;
	if (consOpt.window == 0)
        selectPairsAssembly(seqSet, begEndPos, consOpt.bandwidth, pList, diagList);
	else
        selectPairsAllAgainstAll(seqSet, begEndPos, consOpt.window, pList, diagList);

	// Set-up a sparse distance matrix
	Graph<Undirected<double> > pairGraph;
	
	// Containers for segment matches and corresponding scores 
	typedef String<Fragment<> > TFragmentString;
	TFragmentString matches;
	typedef String<int> TScoreValues;
	TScoreValues scores;

	// Compute segment matches from global pairwise alignments
	appendSegmentMatches(seqSet, pList, diagList, begEndPos, consOpt.sc, consOpt.matchlength, consOpt.quality, consOpt.overlaps, matches, scores, pairGraph, OverlapLibrary() );
	clear(pList);
	clear(diagList);

	// If there are no alignment matches, return
	if (!length(matches)) return;

	// Use these segment matches for the initial alignment graph
	typedef Graph<Alignment<TStringSet, TSize> > TGraph;
	TGraph g(seqSet);
	buildAlignmentGraph(matches, scores, g, consOpt.sc, ReScore() );
	clear(matches);
	clear(scores);

	// Guide Tree
	Graph<Tree<double> > guideTree;
	upgmaTree(pairGraph, guideTree);
	clear(pairGraph);

	// Triplet library extension
	graphBasedTripletLibraryExtension(g);

	// Perform a progressive alignment
	progressiveAlignment(g, guideTree, gOut);
	clear(g);
	clear(guideTree);
}

//////////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TSize>
inline void
consensusAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gOut,
				   String<Pair<TSize, TSize> >& begEndPos) 
{
	ConsensusOptions consOpt;
	consensusAlignment(gOut, begEndPos, consOpt);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TFragSpec, typename TConfig, typename TStringSet, typename TCargo, typename TSpec, typename TContigID>
inline void
updateContig(FragmentStore<TFragSpec, TConfig>& fragStore,
			 Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
			 TContigID contigID)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef std::map<TSize, TSize> TComponentLength;
	typedef char TValue;

	// Initialization
	TStringSet& strSet = stringSet(g);
	TSize nseq = length(strSet);
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';
	TSize maxCoverage = 0;
	TSize len = 0;
	String<TValue> mat;

	// Store for each read the begin position, the end position and the row in the alignment matrix
	String<TSize> readBegEndRowPos;
	resize(readBegEndRowPos, 3*nseq);

	// Strongly Connected Components, topological sort, and length of each component
	String<TSize> component;
	String<TSize> order;
	TComponentLength compLength;
	if (convertAlignment(g, component, order, compLength)) {
		TSize numOfComponents = length(order);
		
		// Assign to each sequence the start and end (in terms of component ranks)
		typedef String<std::pair<TSize, TSize> > TComponentToRank;
		TComponentToRank compToRank;
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) 
			appendValue(compToRank, std::make_pair(order[compIndex], compIndex), Generous());
		std::sort(begin(compToRank, Standard()), end(compToRank, Standard()));

		typedef Pair<TSize, TSize> TRankPair;
		typedef String<TRankPair> TSequenceToRanks;
		TSequenceToRanks seqToRank;
		resize(seqToRank, nseq);
		typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
		TVertexIterator itVertex(g);
		for(;!atEnd(itVertex);++itVertex) {
			TVertexDescriptor vert = value(itVertex);
			TSize seq = idToPosition(strSet, sequenceID(g, vert));
			if (fragmentBegin(g, vert) == 0) 
				seqToRank[seq].i1 = std::lower_bound(begin(compToRank, Standard()), end(compToRank, Standard()), std::make_pair((TSize) component[vert], (TSize) 0))->second;
			if (fragmentBegin(g, vert) + fragmentLength(g, vert) == length(strSet[seq]))
				seqToRank[seq].i2 = std::lower_bound(begin(compToRank, Standard()), end(compToRank, Standard()), std::make_pair((TSize) component[vert], (TSize) 0))->second;
		}
		clear(compToRank);

		// Assign the sequences to rows
		String<TSize> seqToRow;
		resize(seqToRow, nseq);
		maxCoverage = 0;
		typedef String<bool> TLeftOver;
		typedef typename Iterator<TLeftOver, Standard>::Type TLeftOverIter;
		TLeftOver leftOver;
		resize(leftOver, nseq, true);
		typedef String<std::pair<TSize, TSize> > TSeqToBegin;
		typedef typename Iterator<TSeqToBegin, Standard>::Type TSeqToBeginIter;
		TSeqToBegin seqToBegin;
		TSize finishedSeq = 0;
		while(finishedSeq < nseq) {
			TLeftOverIter itL = begin(leftOver, Standard());
			TLeftOverIter itLEnd = end(leftOver, Standard());
			for(TSize pos = 0; itL != itLEnd; ++itL, ++pos) 
				if (*itL) appendValue(seqToBegin, std::make_pair((seqToRank[pos]).i1, pos), Generous());
			std::sort(begin(seqToBegin, Standard()), end(seqToBegin, Standard()));
			
			TSize endPos = 0;
			TSeqToBeginIter itSB = begin(seqToBegin, Standard());
			TSeqToBeginIter itSBEnd = end(seqToBegin, Standard());
			for(;itSB != itSBEnd;++itSB) {
				if (endPos <= (*itSB).first) {
					TSize currentSeq = (*itSB).second;
					seqToRow[currentSeq] = maxCoverage;
					endPos = (seqToRank[currentSeq]).i2 + 2;
					leftOver[currentSeq] = false;
					++finishedSeq;
				}	
			}
			clear(seqToBegin);
			++maxCoverage;
		}
		clear(leftOver);

		// Create the matrix
		len = 0;
		String<TSize> compOffset;
		resize(compOffset, numOfComponents);
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			compOffset[order[compIndex]] = len;
			len+=compLength[order[compIndex]];
		}
		resize(mat, len * maxCoverage, gapChar);

		// Fill in the segments
		typedef typename Infix<TString>::Type TInfix;
		typedef typename Iterator<TInfix, Standard>::Type TInfixIter;
		typedef typename TGraph::TPosToVertexMap_ TPosToVertexMap;
		for(typename TPosToVertexMap::const_iterator it = g.data_pvMap.begin();it != g.data_pvMap.end(); ++it) {
			TInfix str = label(g,it->second);
			TSize c = property(component, it->second);
			TSize row = seqToRow[idToPosition(strSet, it->first.first)];
			//if (row == 0) {
			//	std::cout << sequenceID(g, it->second) << ':' << str << ',' << strSet[sequenceID(g, it->second)] << std::endl;
			//	std::cout << getProperty(component, it->second) << ',' << order[compIndex] << std::endl;
			//	std::cout << (seqToRank[sequenceID(g, it->second)]).i1 << ',' << (seqToRank[sequenceID(g, it->second)]).i2 << std::endl;
			//}
			TInfixIter sIt = begin(str, Standard());
			TInfixIter sItEnd = end(str, Standard());
			TSize i = compOffset[c];
			for(TSize pCol = i;sIt!=sItEnd;++sIt, ++pCol, ++i) 
				mat[row * len + pCol] = *sIt;
		}
		String<bool> active;
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			TSize offset = compOffset[order[compIndex]];
			TSize currentCompLength = compLength[order[compIndex]];

			clear(active);
			resize(active, maxCoverage, false);

			// Find the empty rows
			for(TSize i=0;i<nseq; ++i) {
				if (((seqToRank[i]).i1 <= compIndex) && ((seqToRank[i]).i2 >= compIndex)) 
					active[(seqToRow[i])] = true;
			}
			
			// Substitute false gaps with special gap character
			for(TSize i = 0; i < maxCoverage; ++i) {
				if (!(active[i])) {
					for(TSize pCol = offset;pCol < offset + currentCompLength;++pCol) 
						mat[i * len + pCol] = specialGap;
				}
			}
		}

		// Get the new begin and end positions
		for(TSize i=0;i<nseq; ++i) {
			TVertexDescriptor lastVertex = findVertex(const_cast<TGraph&>(g), positionToID(strSet, i), length(strSet[i]) - 1);
			TSize readBegin = compOffset[getProperty(component, findVertex(const_cast<TGraph&>(g), positionToID(strSet, i), 0))];
			TSize readEnd = compOffset[getProperty(component, lastVertex)] + fragmentLength(const_cast<TGraph&>(g), lastVertex);
			readBegEndRowPos[3*i] = readBegin;
			readBegEndRowPos[3*i+1] = readEnd;
			readBegEndRowPos[3*i+2] = seqToRow[i];
		}

	}
	clear(component);
	clear(order);
	compLength.clear();

	
	//// Debug code
	//for(TSize row = 0; row<maxCoverage; ++row) {
	//	for(TSize col = 0; col<len; ++col) {
	//		std::cout << mat[row * len + col];			
	//	}
	//	std::cout << std::endl;
	//}

	// Create the new consensus
	typedef typename Value<TString>::Type TAlphabet;
	String<TValue> gappedCons;
	consensusCalling(mat, gappedCons, maxCoverage, TAlphabet(), MajorityVote());

	// Assign new consensus
	assignGappedConsensus(fragStore, gappedCons, contigID);

	// Update all aligned reads
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename TFragmentStore::TReadPos TReadPos;
	typedef typename TFragmentStore::TContigPos TContigPos;
	typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
	//typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	sortAlignedReads(fragStore.alignedReadStore, SortContigID());
	TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigID, SortContigID());
	TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigID, SortContigID());
	TReadPos ungappedPos = 0;
	TReadPos gappedPos = 0;
	bool gapOpen;
	for(TSize i = 0;alignIt != alignItEnd; ++alignIt, ++i) {
		TSize lenRead = length(fragStore.readSeqStore[alignIt->readID]);
		TReadPos begClr = 0;
		TReadPos endClr = 0;
		getClrRange(fragStore, *alignIt, begClr, endClr);
		clear(alignIt->gaps);
		ungappedPos = begClr;
		if (alignIt->beginPos > alignIt->endPos) ungappedPos = lenRead - endClr;
		if (ungappedPos != 0) appendValue(alignIt->gaps, TContigGapAnchor(ungappedPos, 0));
		gappedPos = 0;
		gapOpen = false;
		for(TSize column = readBegEndRowPos[3*i]; column<readBegEndRowPos[3*i + 1]; ++column, ++gappedPos) {
			if (mat[readBegEndRowPos[3*i + 2] * len + column] == gapChar) gapOpen = true;				
			else {
				if (gapOpen) {
					appendValue(alignIt->gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());
					gapOpen = false;
				}
				++ungappedPos;
			}
		}
		if (gapOpen) appendValue(alignIt->gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());
		if (alignIt->beginPos < alignIt->endPos) {
			if (endClr != (TContigPos)lenRead) 
				appendValue(alignIt->gaps, TContigGapAnchor(lenRead, lenRead + (gappedPos - ungappedPos) - (lenRead - endClr)), Generous());
		} else {
			if (begClr != 0) 
				appendValue(alignIt->gaps, TContigGapAnchor(lenRead, lenRead + (gappedPos - ungappedPos) - begClr), Generous());
		}

		// Set new begin and end position
		if (alignIt->beginPos < alignIt->endPos) {
			alignIt->beginPos = readBegEndRowPos[3*i];
			alignIt->endPos = readBegEndRowPos[3*i+1];
		} else {
			alignIt->beginPos = readBegEndRowPos[3*i+1];
			alignIt->endPos = readBegEndRowPos[3*i];
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TCounters, typename TSize, typename TAlphabet>
inline void
_countLetters(String<TValue, TSpec> const & mat,
              TCounters & counterValues,
              TSize alignDepth,
              TAlphabet)
{
	typedef String<TValue, TSpec> const TMatrix;
	typedef typename Iterator<TMatrix, Standard>::Type TMatIter;

	// Initialization
	TSize len = length(mat) / alignDepth;
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';
	TSize alphabetSize = ValueSize<TAlphabet>::VALUE;

	// Set-up counter values
	typedef typename Value<TCounters>::Type TCounter;
	typedef typename Iterator<TCounters, Standard>::Type TCounterIt;
	resize(counterValues, len);
	for(TSize i=0;i<len; ++i) {
		TCounter counter;
		resize(counter, alphabetSize + 1, 0);
		counterValues[i] = counter;
	}

	// Count all 
	TMatIter matIt = begin(mat, Standard());
	TMatIter matItEnd = end(mat, Standard());
	TCounterIt countIt = begin(counterValues, Standard());
	TSize pos = 0;
	for(; matIt != matItEnd; ++matIt, ++countIt, ++pos) {
		if (pos % len == 0) countIt = begin(counterValues, Standard());
		if (*matIt != specialGap) {
			if (*matIt == gapChar) ++value(*countIt, alphabetSize);
			else ++value(*countIt, ordValue((TAlphabet) *matIt));
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TGappedCons, typename TAlignDepth, typename TAlphabet>
inline void
consensusCalling(String<TValue, TSpec> const& mat,
				 TGappedCons& gappedConsensus,
				 TAlignDepth maxCoverage,
				 TAlphabet,
				 Bayesian)
{
	typedef double TProbability;
	typedef String<TProbability> TProbabilityDistribution;
	typedef String<TProbabilityDistribution> TPositionalPrDist;
	typedef typename Iterator<TPositionalPrDist, Standard>::Type TPosPrDistIter;

	typedef typename Size<String<TValue, TSpec> >::Type TSize;
	TSize alphabetSize = ValueSize<TAlphabet>::VALUE;
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';

	// Set-up the counters
	typedef String<TSize> TCounter;
	typedef String<TCounter> TCounters;
	TCounters counterValues;
	_countLetters(mat, counterValues, maxCoverage, TAlphabet() );


	// Initialization
	TSize len = length(mat) / maxCoverage;
	TProbabilityDistribution backroundDist;
	resize(backroundDist, alphabetSize + 1, ((TProbability) 1 / (TProbability) (alphabetSize + 1)));
	
	// Get an initial consensus
	typedef typename Iterator<TCounters, Standard>::Type TCounterIt;
	TCounterIt countIt = begin(counterValues, Standard());
	TCounterIt countItEnd = end(counterValues, Standard());
	TPositionalPrDist posPrDist;
	TValue c = TAlphabet();
	for(;countIt != countItEnd; ++countIt) {
		TSize max = 0;
		typedef typename Iterator<TCounter, Standard>::Type TCIt;
		TCIt cIt = begin(*countIt, Standard());
		TCIt cItEnd = end(*countIt, Standard());
		TSize pos = 0;
		for(;cIt != cItEnd; ++cIt, ++pos) {
			if (*cIt > max) {
				max = *cIt;
				c = (pos == alphabetSize) ? gapChar : (TValue) TAlphabet(pos);
			}
		}
		TProbabilityDistribution prDist;
		resize(prDist, alphabetSize + 1, 0);
		if (c == gapChar) prDist[alphabetSize] = 1;
		else prDist[ordValue((TAlphabet) c)] = 1;
		appendValue(posPrDist, prDist, Generous());
	}

	TSize run = 1;
	TProbabilityDistribution pI;
	TProbabilityDistribution pIJ;
	TProbabilityDistribution pIOld;
	TProbabilityDistribution pIJOld;
	while (run) {
		// Store the values from the last iteration
		pIOld = pI;
		pIJOld = pIJ;

		// Count all letters in the consensus
		TProbabilityDistribution nI;
		resize(nI, alphabetSize + 1, 0);
		TPosPrDistIter itPosPrDist = begin(posPrDist, Standard());
		TPosPrDistIter itPosPrDistEnd = end(posPrDist, Standard());
		for(;itPosPrDist!=itPosPrDistEnd; ++itPosPrDist) 
			for(TSize i = 0; i<(alphabetSize + 1); ++i) 
				nI[i] += (*itPosPrDist)[i];
	
		// Composition probabilities
		clear(pI);
		resize(pI, alphabetSize + 1);
		TProbability lenPosPrDist = (TProbability) length(posPrDist);
		for(TSize i = 0; i<length(pI); ++i) 
			pI[i] = nI[i] / lenPosPrDist;
		

		// Count all letters that agree / disagree with the consensus
		TProbabilityDistribution nIJ;
		resize(nIJ, (alphabetSize + 1) * (alphabetSize + 1), 0);
		typedef String<TValue, TSpec> TMatrix;
		typedef typename Iterator<TMatrix, Standard>::Type TMatIter;
		TMatIter matIt = begin(mat, Standard());
		TMatIter matItEnd = end(mat, Standard());
		itPosPrDist = begin(posPrDist, Standard());
		TSize pos = 0;
		for(; matIt != matItEnd; ++matIt, ++itPosPrDist, ++pos) {
			if (pos % len == 0) itPosPrDist = begin(posPrDist, Standard());
			TValue c = *matIt;
			if (c != specialGap) {
				TSize fragJ = (c != gapChar) ? ordValue(TAlphabet(c)) : alphabetSize;
				for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) 
					nIJ[consI * (alphabetSize + 1) + fragJ] += (*itPosPrDist)[consI];	
			}
		}

		// Sequencing error probabilities
		clear(pIJ);
		resize(pIJ, (alphabetSize + 1) * (alphabetSize + 1));
		TProbability sumIJ = 0;
		for(TSize diag = 0; diag<(alphabetSize + 1); ++diag) sumIJ += nIJ[diag * (alphabetSize + 1) + diag];
		for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) 
			for(TSize fragJ = 0; fragJ<(alphabetSize + 1); ++fragJ)
				pIJ[consI * (alphabetSize + 1) + fragJ] = nIJ[consI * (alphabetSize + 1) + fragJ] / sumIJ;
	
		// Recompute positional probability distribution
		itPosPrDist = begin(posPrDist, Standard());
		TSize col = 0;
		for(;itPosPrDist!=itPosPrDistEnd; ++itPosPrDist, ++col) {
			TProbabilityDistribution prDist;
			resize(prDist, alphabetSize + 1);
			for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) {
				TProbability numerator = pI[consI];
				TProbability denominator = 0;
				for(TSize allI = 0; allI<(alphabetSize + 1); ++allI) {
					TProbability denominatorSub = value(pI, allI);
					for(TSize row = 0; row < maxCoverage; ++row) {
						TValue c = mat[row * len + col];
						if (c != specialGap) {
							TSize fragJ = (c != gapChar) ? ordValue(TAlphabet(c)) : alphabetSize;
							if (allI == consI) 
								numerator *= pIJ[allI * (alphabetSize + 1) + fragJ]; 
							denominatorSub *= pIJ[allI * (alphabetSize + 1) + fragJ]; 
						}
					}
					denominator += denominatorSub;
				}
				prDist[consI] = numerator / denominator;
			}
			*itPosPrDist = prDist;
		}	

		// Check termination criterion
		TProbability eps = 0.00001;
		typedef typename Iterator<TProbabilityDistribution, Standard>::Type TProbIter;
		TProbIter pIter = begin(pIOld, Standard());
		TProbIter pIterCompare = begin(pI, Standard());
		TProbIter pIterEnd = end(pIOld, Standard());
		TSize runOld = run;
		for(;pIter != pIterEnd; ++pIter, ++pIterCompare) {
			if (*pIter > *pIterCompare) {
				if (*pIter - *pIterCompare > eps) {
					++run;
					break;
				}
			} else {
				if (*pIterCompare - *pIter > eps) {
					++run;
					break;
				}
			}
		}
		if (runOld == run) {
			pIter = begin(pIJOld, Standard());
			pIterCompare = begin(pIJ, Standard());
			pIterEnd = end(pIJOld, Standard());
			for(;pIter != pIterEnd; ++pIter, ++pIterCompare) {
				if (*pIter > *pIterCompare) {
					if (*pIter - *pIterCompare > eps) {
						++run;
						break;
					}
				} else {
					if (*pIterCompare - *pIter > eps) {
						++run;
						break;
					}
				}
			}
		}

		if (runOld == run) {
			std::cout << "Iterations: " << run << std::endl;
			run = 0;
		}
	}
	
	// Compute the most likely consensus
	TPosPrDistIter itPosPrDist = begin(posPrDist, Standard());
	TPosPrDistIter itPosPrDistEnd = end(posPrDist, Standard());
	clear(gappedConsensus);
	for(;itPosPrDist!=itPosPrDistEnd; ++itPosPrDist) {
		TProbability max = 0;
		TSize ind = 0;
		for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) {
			if ((*itPosPrDist)[consI] > max) {
				max = (*itPosPrDist)[consI];
				ind = consI;
			}
		}
		if (ind == alphabetSize) appendValue(gappedConsensus, gapChar);
		else appendValue(gappedConsensus, TAlphabet(ind));
	}

}

//////////////////////////////////////////////////////////////////////////////


template <typename TFragSpec, typename TConfig, typename TContigID>
inline void
consensusCalling(FragmentStore<TFragSpec, TConfig>& fragStore,
				 TContigID contigID,
				 Bayesian)
{
	SEQAN_CHECKPOINT

	typedef FragmentStore<TFragSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadSeq TReadSeq;
	typedef typename Value<TReadSeq>::Type TAlphabet;
	typedef char TValue;

	// Convert the contig to an alignment matrix
	typedef String<TValue> TAlignMat;
	TAlignMat mat;
	TSize maxCoverage;
	convertAlignment(fragStore, mat, contigID, maxCoverage);

	// Call the consensus
	String<TValue> gappedConsensus;
	consensusCalling(mat, gappedConsensus, maxCoverage, TAlphabet(), Bayesian());

	// Assign the new consensus
	assignGappedConsensus(fragStore, gappedConsensus, contigID);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TGappedCons, typename TAlignDepth, typename TAlphabet>
inline void
consensusCalling(String<TValue, TSpec> const& mat,
				 TGappedCons& gappedConsensus,
				 TAlignDepth maxCoverage,
				 TAlphabet,
				 MajorityVote)
{
	typedef typename Size<String<TValue, TSpec> >::Type TSize;
	TSize alphabetSize = ValueSize<TAlphabet>::VALUE;
	TValue gapChar = gapValue<TValue>();
	
	// Set-up the counters
	typedef String<TSize> TCounter;
	typedef String<TCounter> TCounters;
	TCounters counterValues;
	_countLetters(mat, counterValues, maxCoverage, TAlphabet() );
	
	// Get the consensus
	typedef typename Iterator<TCounters, Standard>::Type TCounterIt;
	TCounterIt countIt = begin(counterValues, Standard());
	TCounterIt countItEnd = end(counterValues, Standard());
	clear(gappedConsensus);
	TSize max = 0;
	TValue c = TValue();
	TSize pos = 0;
	for(;countIt != countItEnd; ++countIt) {
		max = 0;	
		typedef typename Iterator<TCounter, Standard>::Type TCIt;
		TCIt cIt = begin(*countIt, Standard());
		TCIt cItEnd = end(*countIt, Standard());
		pos = 0;
		for(;cIt != cItEnd; ++cIt, ++pos) {
			if (*cIt > max) {
				max = *cIt;
				c = (pos == alphabetSize) ? gapChar : (TValue) TAlphabet(pos);
			}
		}
		appendValue(gappedConsensus, c);
	}
}


//////////////////////////////////////////////////////////////////////////////


template <typename TFragSpec, typename TConfig, typename TContigID>
inline void
consensusCalling(FragmentStore<TFragSpec, TConfig>& fragStore,
				 TContigID contigID,
				 MajorityVote)
{
	typedef FragmentStore<TFragSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadSeq TReadSeq;
	typedef typename Value<TReadSeq>::Type TAlphabet;
	typedef char TValue;

	// Convert the contig to an alignment matrix
	typedef String<TValue> TAlignMat;
	TAlignMat mat;
	TSize maxCoverage;
	convertAlignment(fragStore, mat, contigID, maxCoverage);

	// Call the consensus
	String<TValue> gappedConsensus;
	consensusCalling(mat, gappedConsensus, maxCoverage, TAlphabet(), MajorityVote());

	// Assign the new consensus
	assignGappedConsensus(fragStore, gappedConsensus, contigID);
}


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Old proprietary FastaReadFormat
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

struct FastaReadFormat_;
typedef Tag<FastaReadFormat_> const FastaReadFormat;

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec, typename TConfig>
inline void
write(
    TFile & file,
    FragmentStore<TSpec, TConfig>& fragStore,
    FastaReadFormat)
{
//IOREV _nodoc_ what has this got to do with consensus? why is it here?
	// Basic types
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TContigPos TContigPos;
	//typedef typename Value<TFile>::Type TValue;
	typedef char TMultiReadChar;
	TMultiReadChar gapChar = gapValue<TMultiReadChar>();

    typename DirectionIterator<TFile, Output>::Type target = directionIterator(file, Output());

	typedef typename Iterator<typename TFragmentStore::TContigStore, Standard>::Type TContigIter;
	TContigIter contigIt = begin(fragStore.contigStore, Standard() );
	TContigIter contigItEnd = end(fragStore.contigStore, Standard() );
	for(TSize idCount = 0;contigIt != contigItEnd; ++contigIt, ++idCount) {
		// Alignment matrix
		typedef String<TMultiReadChar> TAlignMat;
		TAlignMat mat;
		TSize maxCoverage;
		String<TSize> readSlot;
		convertAlignment(fragStore, mat, idCount, maxCoverage, readSlot);
		TSize len = length(mat) / maxCoverage;
		
		// Gapped consensus sequence
		typedef String<TMultiReadChar> TGappedConsensus;
		TGappedConsensus gappedConsensus;
		getGappedConsensus(fragStore, gappedConsensus, idCount);

		// Print the alignment matrix
		String<TSize> coverage;
		resize(coverage, len, 0);
		typedef typename Iterator<TGappedConsensus, Standard>::Type TConsIter;
		TConsIter itCons = begin(gappedConsensus, Standard());
		TSize winSize = 60;
		int offset = 2;
		TSize column = 0;
		while (column<len) {
			TSize window_end = column + winSize;
			if (window_end >= len) window_end = len;
			// Position
			for(int i = 0; i<offset - 2; ++i) writeValue(target, ' ');
			write(target, "Pos: ");
			appendNumber(target, column);
			writeValue(target, '\n');
			// Ruler
			for(int i = 0; i<offset + 3; ++i) writeValue(target, ' ');
			for(TSize local_col = 1; local_col<window_end - column + 1; ++local_col) {
				if ((local_col % 10)==0)
                    writeValue(target, ':');
				else if ((local_col % 5)==0)
                    writeValue(target, '.');
				else
                    writeValue(target, ' ');
			}
			writeValue(target, '\n');
			// Matrix
			for(TSize row = 0; row<maxCoverage; ++row) {
				TSize tmp = row;
				int off = 0;
				while (tmp / 10 != 0) {
					tmp /= 10;
					++off;
				}
				for(int i = 0; i<offset - off; ++i) writeValue(target, ' ');
				appendNumber(target, row);
				writeValue(target, ':');
				writeValue(target, ' ');
				for(TSize local_col = column; local_col<window_end; ++local_col) {
					writeValue(target, mat[row * len + local_col]);
					if (mat[row * len + local_col] != '.') ++coverage[local_col];
				}
				writeValue(target, '\n');
			}
			writeValue(target, '\n');
	
			// Consensus
			for(int i = 0; i<offset; ++i)
                writeValue(target, ' ');
			write(target, "C: ");
			for(unsigned int local_col = column; local_col<window_end; ++local_col, ++itCons) 
				writeValue(target, *itCons);
			writeValue(target, '\n');
			for(int i = 0; i<offset-1; ++i)
                writeValue(target, ' ');
			write(target, ">2: ");
			for(unsigned int local_col = column; local_col<window_end; ++local_col) {
				if (coverage[local_col] > 2)
                    writeValue(target, gappedConsensus[local_col]);
				else
                    writeValue(target, gapChar);
			}
			writeValue(target, '\n');
			writeValue(target, '\n');
			column+=winSize;
		}
		writeValue(target, '\n');
		writeValue(target, '\n');

		// Print all aligned reads belonging to this contig

		// Sort according to contigID
		sortAlignedReads(fragStore.alignedReadStore, SortContigID());
	
		// Find range of the given contig
		typedef typename Iterator<typename TFragmentStore::TAlignedReadStore, Standard>::Type TAlignIter;
		TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigID());
		TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigID());

		// Sort the reads according to the begin position
		sortAlignedReads(infix(fragStore.alignedReadStore, alignIt - begin(fragStore.alignedReadStore, Standard()), alignItEnd - begin(fragStore.alignedReadStore, Standard())), SortBeginPos());
		TAlignIter alignItTmp = lowerBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigID());
		TAlignIter alignItTmpEnd = upperBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigID());
		String<std::pair<TSize, TSize> > idToPos;
		reserve(idToPos, alignItTmpEnd - alignItTmp);
		for(TSize iCount = 0; alignItTmp!=alignItTmpEnd; ++iCount, ++alignItTmp) 
			appendValue(idToPos, std::make_pair(alignItTmp->id, readSlot[iCount]));
		std::sort(begin(idToPos, Standard()), end(idToPos, Standard()));

		// Sort the reads according to the id
		sortAlignedReads(infix(fragStore.alignedReadStore, alignIt - begin(fragStore.alignedReadStore, Standard()), alignItEnd - begin(fragStore.alignedReadStore, Standard())), SortID());
		alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigID());
		alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigID());

		bool noNamesPresent = (length(fragStore.readNameStore) == 0);
		for(TSize iCount = 0;alignIt != alignItEnd; ++alignIt, ++iCount) {

			// Print all reads
			write(target, "typ:");
			if (!noNamesPresent)
            {
				writeValue(target, 'R');
				appendNumber(target, iCount);
			}
            else
            {
                write(target, fragStore.readNameStore[alignIt->readID]);
            }
			write(target, "\nseq:");
			write(target, fragStore.readSeqStore[alignIt->readID]);
			write(target, "\nPos:");
			appendNumber(target, alignIt->beginPos);
			writeValue(target, ',');
			appendNumber(target, alignIt->endPos);
			writeValue(target, '\n');
#ifndef CELERA_OFFSET
			TSize begClr = 0;
			TSize endClr = 0;
			getClrRange(fragStore, *alignIt, begClr, endClr);
			write(target, "clr:");
			appendNumber(target, begClr);
			writeValue(target, ',');
			appendNumber(target, endClr);
			writeValue(target, '\n');
#endif
			std::stringstream gapCoords;
			TSize letterCount = 0;
			TSize gapCount = 0;
			for(TContigPos column = _min(alignIt->beginPos, alignIt->endPos); column < _max(alignIt->beginPos, alignIt->endPos); ++column) {
				if (mat[idToPos[iCount].second * len + column] == gapChar) {
					++gapCount;
					gapCoords << letterCount << ' ';
				} else ++letterCount;
			}
			write(target, "dln:");
			appendNumber(target, gapCount);
			writeValue(target, '\n');
			write(target, "del:");
			write(target, gapCoords.str());
			writeValue(target, '\n');
			writeValue(target, '\n');
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// Read simulator format: Simple fasta read file with positions
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TSpec, typename TConfig, typename TFilePath>
inline void
_convertSimpleReadFile(TFile& /*file*/,
					   FragmentStore<TSpec, TConfig>& /*fragStore*/,
					   TFilePath& /*filePath*/,
					   bool /*moveToFront*/)
{
return;
////IOREV _nodoc_ huge undocumented function, uses custom IO based on iostream and FILE* :S
//	// Basic types
//	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
//	typedef typename ID<TFragmentStore>::Type TID;
//	typedef typename Size<TFragmentStore>::Type TSize;
//	//typedef typename Value<TFile>::Type TValue;
//	typedef typename TFragmentStore::TContigPos TPos;
//	typedef typename TFragmentStore::TReadSeq TReadSeq;
//
//	// All fragment store element types
//	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
//	typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
//	typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
//	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
//	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
//	
//
//	// All maps to mirror file ids to our internal ids
//	typedef std::map<TID, TID> TIDMap;
//	TIDMap libIDMap;
//	TIDMap frgIDMap;
//	TIDMap readIDMap;
//
//    // Create RecordReader object.
//    typename DirectionIterator<TFile, Input>::Type reader = directionIterator(file, Input());
//
//	// Parse the file and convert the internal ids
//	TPos maxPos = 0;
//	TPos minPos = MaxValue<TPos>::VALUE;
//	TID count = 0;
//    if (atEnd(reader))
//        return false;
//    CharString buffer;
//	while (!atEnd(reader))
//    {
//		// New read?
//		if (value(reader) == '>')
//        {
//			TAlignedElement alignEl;
//			TID id = count;
//			TID fragID = TReadStoreElement::INVALID_ID;
//			TID repeatID = 0;
//
//            goNext(reader);
//            if (skipWhitespaces(reader) != 0)
//                return 1;
//
//			// Get the layout positions
//            clear(buffer);
//            if (readDigits(buffer, reader) != 0)
//                return 1;
//            if (!lexicalCast2(alignEl.beginPos, buffer))
//                return 1;
//            goNext(reader);
//            if (skipWhitespaces(reader) != 0)
//                return 1;
//            clear(buffer);
//            if (readDigits(buffer, reader) != 0)
//                return 1;
//            if (!lexicalCast2(alignEl.endPos, buffer))
//                return 1;
//			
//			// Any attributes?
//			String<char> eid;
//			String<char> qlt;
//			TReadSeq seq;
//			if (value(reader) == '[')
//            {
//				String<char> fdIdentifier;
//				while (value(reader) != ']')
//                {
//                    goNext(reader);
//                    if (skipWhitespaces(reader) != 0)
//                        return 1;
//					clear(fdIdentifier);
//                    if (readAlphaNums(fdIdentifier, reader) != 0)
//                        return 1;
//                    goNext(reader);  // Skip "="
//					if (fdIdentifier == "id")
//                    {
//                        clear(buffer);
//                        if (readDigits(buffer, reader) != 0)
//                            return 1;
//                        if (!lexicalCast2(id, buffer))
//                            return 1;
//					} else if (fdIdentifier == "fragID") {
//                        clear(buffer);
//                        if (readDigits(buffer, reader) != 0)
//                            return 1;
//                        if (!lexicalCast2(fragID, buffer))
//                            return 1;
//					} else if (fdIdentifier == "repeatID") {
//                        clear(buffer);
//                        if (readDigits(buffer, reader) != 0)
//                            return 1;
//                        if (!lexicalCast2(repeatID, buffer))
//                            return 1;
//					} else if (fdIdentifier == "eid") {
//                        if (readUntilOneOf(eid, reader, ',', ']') != 0)
//                            return 1;
//					} else if (fdIdentifier == "qlt") {
//                        if (readUntilOneOf(qlt, reader, ',', ']') != 0)
//                            return 1;
//					} else {
//						// Jump to next attribute
//                        // TODO(holtgrew): Add skipUntilOneOf()?
//                        if (readUntilOneOf(buffer, reader, ',', ']') != 0)
//                            return 1;
//					}
//				}
//			}
//            if (skipLine(reader) != 0)
//                return 1;
//            if (skipWhitespaces(reader) != 0)
//                return 1;
//            while (!atEnd(reader) && value(reader) != '>')
//            {
//                if (readLetters(seq, reader) != 0)
//                    return 1;
//                int res = skipWhitespaces(reader);
//                if (res != 0 && res != EOF_BEFORE_SUCCESS)
//                    return 1;
//            }
//			
//			// Set quality
//			typedef typename Iterator<TReadSeq, Standard>::Type TReadIter;
//			typedef typename Iterator<String<char> >::Type TQualIter;
//			TReadIter begIt = begin(seq, Standard() );
//			TReadIter begItEnd = begin(seq, Standard() );
//			if (length(qlt)) {
//				TQualIter qualIt = begin(qlt);
//				TQualIter qualItEnd = end(qlt);
//				for(;qualIt != qualItEnd; goNext(qualIt), goNext(begIt)) assignQualityValue(value(begIt), value(qualIt));
//			} else {
//				for(;begIt != begItEnd; goNext(begIt)) assignQualityValue(value(begIt), 'D');
//			}
//
//			// Set eid if not given
//			if (empty(eid)) {
//				std::stringstream input;
//				input << "R" << id;
//				input << "-" << repeatID;
//				eid = input.str().c_str();
//			}
//
//			// Insert the read
//			readIDMap.insert(std::make_pair(id, static_cast<TID>(length(fragStore.readStore))));
//			appendRead(fragStore, seq, fragID);
//			appendValue(fragStore.readNameStore, eid, Generous());
//
//			// Insert an aligned read
//			TSize readLen = length(seq);
//			if (alignEl.beginPos < alignEl.endPos) {
//				if ((TPos)readLen != alignEl.endPos - alignEl.beginPos) {
//					alignEl.endPos = alignEl.beginPos + readLen;
//				}
//				if (alignEl.beginPos < minPos) minPos = alignEl.beginPos;
//				if (alignEl.endPos > maxPos) maxPos = alignEl.endPos;
//			} else {
//				if ((TPos)readLen != alignEl.beginPos - alignEl.endPos) {
//					alignEl.beginPos = alignEl.endPos + readLen;
//				}
//				if (alignEl.endPos < minPos) minPos = alignEl.endPos;
//				if (alignEl.beginPos > maxPos) maxPos = alignEl.beginPos;
//			}
//			alignEl.readID = id;
//			alignEl.pairMatchID =  fragID;
//			alignEl.contigID = 0;
//			alignEl.id = length(fragStore.alignedReadStore);
//			appendValue(fragStore.alignedReadStore, alignEl, Generous());
//			++count;
//		} else {
//            if (skipLine(reader) != 0)
//                return 1;
//		}
//	}
//
//	// Read contig or reference sequence
//	TContigElement contigEl;
//	std::string fileName = filePath + 'S';
//	SeqFileIn strmRef;
//	String<char> contigEid = "C0";
//    if (open(strmRef, fileName.c_str()))
//    {
//        typename DirectionIterator<TFile, Input>::Type readerRef = directionIterator(strmRef, Input());
//        clear(contigEid);
//        readRecord(contigEid, contigEl.seq, strmRef);
//		close(strmRef);
//	}
//	if (empty(contigEl.seq)) {
//		typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
//		if (moveToFront) appendValue(contigEl.gaps, TContigGapAnchor(0, maxPos - minPos), Generous());
//		else appendValue(contigEl.gaps, TContigGapAnchor(0, maxPos), Generous());
//	}
//	appendValue(fragStore.contigStore, contigEl, Generous());
//	appendValue(fragStore.contigNameStore, contigEid, Generous());
//
//
//	// Read fragments
//	fileName = filePath + 'F';
//	FILE* strmFrag = fopen(fileName.c_str(), "rb");
//    if (!strmFrag)
//        return 1;
//    RecordReader<FILE *, SinglePass<> > readerFrag(strmFrag);
//    while (!atEnd(readerFrag))
//    {
//        if (value(readerFrag) == '>')
//        {
//            TMatePairElement matePairEl;
//            goNext(readerFrag);
//            if (skipWhitespaces(readerFrag) != 0)
//                return 1;
//
//            // Get the fragment id
//            clear(buffer);
//            if (readDigits(buffer, readerFrag) != 0)
//                return 1;
//            TID id = 0;
//            if (!lexicalCast2(id, buffer))
//                return 1;
//			
//            // Any attributes?
//            std::stringstream input;
//            input << "F" << id;
//            String<char> eid(input.str().c_str());
//            if (value(readerFrag) == '[')
//            {
//                String<char> fdIdentifier;
//                while (value(readerFrag) != ']')
//                {
//                    goNext(readerFrag);
//                    if (skipWhitespaces(readerFrag) != 0)
//                        return 1;
//                    clear(fdIdentifier);
//                    if (readAlphaNums(fdIdentifier, readerFrag) != 0)
//                        return 1;
//                    goNext(readerFrag);
//                    if (fdIdentifier == "libID")
//                    {
//                        clear(buffer);
//                        if (readDigits(buffer, readerFrag) != 0)
//                            return 1;
//                        if (!lexicalCast2(matePairEl.libID, buffer))
//                            return 1;
//                    } else if (fdIdentifier == "eid") {
//                        clear(eid);
//                        if (readUntilOneOf(eid, readerFrag, ',', ']') != 0)
//                            return 1;
//                    } else {
//                        // Jump to next attribute
//                        if (readUntilOneOf(buffer, readerFrag, ',', ']') != 0)
//                            return 1;
//                    }
//                }
//            }
//            if (skipLine(readerFrag) != 0)
//                return 1;
//
//            // Read the two reads belonging to this mate pair
//            clear(buffer);
//            for (int i = 0; i < 2; ++i)
//            {
//                if (skipWhitespaces(readerFrag) != 0)
//                    return 1;
//                clear(buffer);
//                if (readDigits(buffer, readerFrag) != 0)
//                    return 1;
//                if (!lexicalCast2(matePairEl.readID[i], buffer))
//                    return 1;
//                if (!i)  // Skip ','.
//                    goNext(readerFrag);
//            }
//            int res = skipLine(readerFrag);
//            if (res != 0 && res != EOF_BEFORE_SUCCESS)
//                return 1;
//
//            // Insert mate pair
//            if (matePairEl.readID[0] != matePairEl.readID[1]) {
//                frgIDMap.insert(std::make_pair(id, static_cast<TID>(length(fragStore.matePairStore))));
//                appendValue(fragStore.matePairStore, matePairEl, Generous());
//                appendValue(fragStore.matePairNameStore, eid, Generous());
//            }
//        } else {
//            int res = skipLine(readerFrag);
//            if (res != 0 && res != EOF_BEFORE_SUCCESS)
//                return 1;
//        }
//    }
//    fclose(strmFrag);
//
//	// Read libraries
//	fileName = filePath + 'L';
//	FILE* strmLib = fopen(fileName.c_str(), "rb");
//    if (!strmLib)
//        return 1;
//    RecordReader<FILE *, SinglePass<> > readerLib(strmLib);
//    while (!atEnd(readerLib))
//    {
//        if (value(readerLib) == '>')
//        {
//            TLibraryStoreElement libEl;
//            goNext(readerLib);
//            if (skipWhitespaces(readerLib) != 0)
//                return 1;
//
//            // Get the fragment id
//            clear(buffer);
//            if (readDigits(buffer, readerLib) != 0)
//                return 1;
//            TID id = 0;
//            if (!lexicalCast2(id, buffer))
//                return 1;
//			
//            // Any attributes?
//            std::stringstream input;
//            input << "L" << id;
//            String<char> eid(input.str().c_str());
//            if (value(readerLib) == '[')
//            {
//                String<char> fdIdentifier;
//                while (value(readerLib) != ']')
//                {
//                    goNext(readerLib);
//                    if (skipWhitespaces(readerLib) != 0)
//                        return 1;
//                    clear(fdIdentifier);
//                    if (readAlphaNums(fdIdentifier, readerLib) != 0)
//                        return 1;
//                    if (fdIdentifier == "eid")
//                    {
//                        clear(eid);
//                        if (readUntilOneOf(eid, readerLib, ',', ']') != 0)
//                            return 1;
//                    } else {
//                        // Jump to next attribute
//                        // TODO(holtgrew): skipUntilOneOf()?
//                        if (readUntilOneOf(buffer, readerLib, ',', ']') != 0)
//                            return 1;
//                    }
//                }
//            }
//            if (skipLine(readerLib) != 0)
//                return 1;
//            if (skipWhitespaces(readerLib) != 0)
//                return 1;
//
//            // Read the mean and standard deviation
//            clear(buffer);
//            if (readDigits(buffer, readerLib) != 0)
//                return 1;
//            if (!lexicalCast2(libEl.mean, buffer))
//                return 1;
//            if (skipWhitespaces(readerLib) != 0)
//                return 1;
//            goNext(readerLib);
//            clear(buffer);
//            if (readDigits(buffer, readerLib) != 0)
//                return 1;
//            if (!lexicalCast2(libEl.std, buffer))
//                return 1;
//            int res = skipLine(readerLib);
//            if (res != 0 && res != EOF_BEFORE_SUCCESS)
//                return 1;
//
//            // Insert mate pair
//            libIDMap.insert(std::make_pair(id, static_cast<TID>(length(fragStore.libraryStore))));
//            appendValue(fragStore.libraryStore, libEl, Generous());
//            appendValue(fragStore.libraryNameStore, eid, Generous());
//        } else {
//            int res = skipLine(readerLib);
//            if (res != 0 && res != EOF_BEFORE_SUCCESS)
//                return 1;
//        }
//    }
//    fclose(strmLib);
//	
//	// Renumber all ids
//	typedef typename TIDMap::const_iterator TIDMapIter;
//	typedef typename Iterator<typename TFragmentStore::TMatePairStore>::Type TMateIter;
//	TMateIter mateIt = begin(fragStore.matePairStore);
//	TMateIter mateItEnd = end(fragStore.matePairStore);
//	for(;mateIt != mateItEnd; goNext(mateIt)) {
//		if (mateIt->libID != TMatePairElement::INVALID_ID) {
//			TIDMapIter libIDPos = libIDMap.find(mateIt->libID);
//			if (libIDPos != libIDMap.end()) mateIt->libID = libIDPos->second;
//			else mateIt->libID = TMatePairElement::INVALID_ID;
//		}
//		if (mateIt->readID[0] != TMatePairElement::INVALID_ID) {
//			TIDMapIter readIDPos = readIDMap.find(mateIt->readID[0]);
//			if (readIDPos != readIDMap.end()) mateIt->readID[0] = readIDPos->second;
//			else mateIt->readID[0] = TMatePairElement::INVALID_ID;
//		}
//		if (mateIt->readID[1]!= TMatePairElement::INVALID_ID) {
//			TIDMapIter readIDPos = readIDMap.find(mateIt->readID[1]);
//			if (readIDPos != readIDMap.end()) mateIt->readID[1] = readIDPos->second;
//			else mateIt->readID[0] = TMatePairElement::INVALID_ID;
//		}
//	}
//	typedef typename Iterator<typename TFragmentStore::TReadStore>::Type TReadIter;
//	TReadIter readIt = begin(fragStore.readStore);
//	TReadIter readItEnd = end(fragStore.readStore);
//	for(;readIt != readItEnd; goNext(readIt)) {
//		if (readIt->matePairID != TReadStoreElement::INVALID_ID) {
//			TIDMapIter mateIDPos = frgIDMap.find(readIt->matePairID);
//			if (mateIDPos != frgIDMap.end()) readIt->matePairID = mateIDPos->second;
//			else readIt->matePairID = TReadStoreElement::INVALID_ID;
//		}
//	}
//	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
//	TAlignIter alignIt = begin(fragStore.alignedReadStore);
//	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
//	for(;alignIt != alignItEnd; goNext(alignIt)) {
//		if (alignIt->readID != TAlignedElement::INVALID_ID) {
//			TIDMapIter readIDPos = readIDMap.find(alignIt->readID);
//			if (readIDPos != readIDMap.end()) alignIt->readID = readIDPos->second;
//			else alignIt->readID = TAlignedElement::INVALID_ID;
//		}
//		if (moveToFront) {
//			alignIt->beginPos -= minPos;
//			alignIt->endPos -= minPos;
//		}
//	}
//	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Rudimentary write functions for CeleraFrg and Celera Cgb
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec, typename TConfig>
inline void 
_writeCeleraFrg(TFile& file,
				FragmentStore<TSpec, TConfig>& fragStore) 
{
//IOREV _nodoc_ 
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadPos TReadPos;

    typename DirectionIterator<TFile, Output>::Type target = directionIterator(file, Output());

	// Iterate over all aligned reads to get the clear ranges
	typedef Pair<TReadPos, TReadPos> TClrRange;
	String<TClrRange> clearStr;
	resize(clearStr, length(fragStore.readStore));
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	for(;alignIt != alignItEnd; goNext(alignIt))
    {
		TReadPos begClr = 0;
		TReadPos endClr = 0;
		getClrRange(fragStore, value(alignIt), begClr, endClr);
		value(clearStr, alignIt->readID) = TClrRange(begClr, endClr);
	}

	// Write Reads
	typedef typename Iterator<typename TFragmentStore::TReadStore>::Type TReadIter;
	TReadIter readIt = begin(fragStore.readStore);
	TReadIter readItEnd = end(fragStore.readStore);
	bool noNamesPresent = (length(fragStore.readNameStore) == 0);
	for (TSize idCount = 0; readIt != readItEnd; goNext(readIt), ++idCount)
    {
		write(target, "{FRG\nact:A\nacc:");
		appendNumber(target, idCount + 1);
		write(target, "\ntyp:R\n");
		if (!noNamesPresent)
        {
			write(target, "src:\n");
			write(target, fragStore.readNameStore[idCount]);
			write(target, "\n.\n");
		}
		write(target, "etm:0\nseq:\n");
        writeWrappedString(target, fragStore.readSeqStore[idCount], 70);
		write(target, ".\nqlt:\n");
        typedef typename Value<typename TFragmentStore::TReadSeqStore>::Type TReadSeq;
        typedef QualityExtractor<typename Value<TReadSeq>::Type> TQualityExtractor;
        ModifiedString<TReadSeq const, ModView<TQualityExtractor> > quals(fragStore.readSeqStore[idCount]);
        writeWrappedString(target, quals, 70);
		// Note: Clear range does not have to be ordered, e.g. no indication for reverse complemented reads, this is happening in cgb records
		write(target, ".\nclr:");
		appendNumber(target, clearStr[idCount].i1);
		writeValue(target, ',');
		appendNumber(target, clearStr[idCount].i2);
		write(target, "\n}\n");
	}
}


//////////////////////////////////////////////////////////////////////////////

// Writes out the first contig only.

// TODO(holtgrew): This is only used in seqcons, should it go into app?
template<typename TFile, typename TSpec, typename TConfig>
inline void
_writeCeleraCgb(TFile& file,
				FragmentStore<TSpec, TConfig>& fragStore)
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename ID<TFragmentStore>::Type TID;
	//typedef typename TFragmentStore::TReadPos TReadPos;

    typename DirectionIterator<TFile, Output>::Type target = directionIterator(file, Output());

	// Write the first contig
	TID contigID = 0;

	// Sort the reads according to position
	sortAlignedReads(fragStore.alignedReadStore, SortBeginPos());

	// Write Header
	write(target, "{IUM\nacc:0\nsrc:\ngen> @@ [0,0]\n.\ncov:0.000\nsta:X\nfur:X\nabp:0\nbbp:0\nlen:");
	appendNumber(target, length((value(fragStore.contigStore, contigID)).seq));
	write(target, "\ncns:\n.\nqlt:\n.\nfor:0\nnfr:");
	appendNumber(target, length(fragStore.readStore));
	writeValue(target, '\n');

	// Write reads
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	TSize offsetLeft = _min(alignIt->beginPos, alignIt->endPos);
	for (;alignIt != alignItEnd; goNext(alignIt))
    {
		if (contigID != alignIt->contigID)
            continue;
		write(target, "{IMP\ntyp:R\nmid:");
		appendNumber(target, alignIt->readID + 1);
		write(target, "\ncon:0\npos:");
		appendNumber(target, alignIt->beginPos - offsetLeft);
		writeValue(target, ',');
		appendNumber(target, alignIt->endPos - offsetLeft);
		write(target, "\ndln:0\ndel:\n}\n");
	}
	write(target, "}\n");
}

}  // namespace seqan

#endif //#ifndef SEQAN_HEADER_...
