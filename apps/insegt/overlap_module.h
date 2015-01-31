/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $ID$
 ==========================================================================*/

#ifndef SEQAN_HEADER_OVERLAP_MODULE_H
#define SEQAN_HEADER_OVERLAP_MODULE_H
// #define DEBUG_OVERLAP_MODULE

#include <seqan/bam_io.h>
#include <seqan/store.h>

namespace SEQAN_NAMESPACE_MAIN
{

// Assign ContigID of children to parents    
template<typename TSpec, typename TConfig>
inline void
adjustParentEntries(FragmentStore<TSpec, TConfig> & fragStore)
{
    // Iterate over all parents and assign contigID of first child
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	                TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type 				                        TAnnotationStoreElement;
    typedef typename Iterator<FragmentStore<TSpec, TConfig>, AnnotationTree<> >::Type   TIter;

    TIter it = begin(fragStore, AnnotationTree<>());
    // Go down from root to the first parent
    goDown(it);
    do
    { // Iterate over all parents
        if (getAnnotation(it).contigID == TAnnotationStoreElement::INVALID_ID)
        {
            getAnnotation(it).contigID = getAnnotation(nodeDown(it)).contigID;  // If no child exists: nothing happens
            getAnnotation(it).beginPos = TAnnotationStoreElement::INVALID_POS;
            getAnnotation(it).endPos = TAnnotationStoreElement::INVALID_POS;
        }
    } while (goRight(it));
}
    
//////////////////////////////////////////////////////////////////////////////
// getIDsFroRead
//////////////////////////////////////////////////////////////////////////////

template<typename TAnnoIDs, typename TSpec, typename TConfig, typename TIntervalTree, typename TIntervals>
inline void
getIDsForRead(TAnnoIDs & ids, FragmentStore<TSpec, TConfig> & fragStore, TIntervalTree & intervalTree, TIntervals & alignIntervals, unsigned offsetInterval)
{
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos 		TContigPos;
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TID 				TID;
	//typedef typename Value<TIntervals>::Type 				TInterval;
	typedef 	 String<TID> 						TResult;
	typedef typename Iterator<TIntervals >::Type 				TIntervalIter;
	typedef typename Iterator<StringSet<TResult > >::Type			TResultIter;
	//typedef typename Iterator<TResult >::Type				TIDIter;
	
	static const TID INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	resize(ids, length(alignIntervals));
	
	TIntervalIter itI = begin(alignIntervals);
	TIntervalIter itIEnd = end(alignIntervals);
	TResultIter itR = begin(ids);
	// TResultIter itREnd = end(ids);

	// search mapped annotations for each interval of the aligned read and store them in the corresponding list 'ids'
	for ( ; itI != itIEnd; goNext(itI), goNext(itR))
	{
		findIntervalsForInterval(value(itR), intervalTree, getValue(itI), offsetInterval);	
	}
	
	// check for each aligment-interval, if the inner interval-borders fit to the borders of the annotation id:
	itI = begin(alignIntervals);
	itR = begin(ids);
	TID currentID;
	TContigPos beginPos;
	TContigPos endPos;
	
	for ( ; itI != itIEnd; goNext(itI), goNext(itR))
	{
		for (unsigned i = 0; i < length(getValue(itR)); ++i)
		{
			currentID = getValue(getValue(itR), i);
			beginPos = getValue(fragStore.annotationStore, currentID).beginPos;
			endPos = getValue(fragStore.annotationStore, currentID).endPos;
			
			if (beginPos > endPos)
			{
				beginPos = endPos;
				endPos = beginPos;
			}
			// begin of read
			if (itR == begin(ids) && length(ids) > 1)
			{
				if (static_cast<TContigPos>(getValue(itI).i2 + offsetInterval) < endPos) // if the borders don't fit: delete annotation-id
				{
					erase(value(itR), i);
					--i;
				}
			}
			// end of read
			else if (position(itR, ids) == endPosition(ids) - 1 && length(ids) > 1u)
			{
				if (static_cast<TContigPos>(getValue(itI).i1 - offsetInterval) > beginPos)
				{
					erase(value(itR), i);
					--i;
				}
			}
			// in the middle of the read
			else if (length(ids) > 2)
			{
				if (static_cast<TContigPos>(getValue(itI).i2 + offsetInterval) < endPos)
				{
					erase(value(itR), i);
					--i;
				}
				else if (static_cast<TContigPos>(getValue(itI).i1 - offsetInterval) > beginPos)
				{
					erase(value(itR), i);
					--i;
				}
			}
		}
		if (empty(getValue(itR)) )  // if aligment-interval doesn't fit to any annotation, append INVALID_ID to mark this
			appendValue(value(itR), INVALID_ID, Generous());
	}
}


//////////////////////////////////////////////////////////////////////////////
////// ReadAnnoStoreELement
//////////////////////////////////////////////////////////////////////////////
template <typename TID>
struct ReadAnnoStoreElement
{
	typedef StringSet<String<TID> > TAnnoIDs;
	
	TAnnoIDs	annoIDs;
	String<TID> 	parentIDs;     // not only for exon-annotations -> more than one parentID possible, only if whole read mapped in parent
	TID		contigID;
};


//////////////////////////////////////////////////////////////////////////////
////// assign IDs to ReadAnnoStore
//////////////////////////////////////////////////////////////////////////////
template<typename TReadAnnoStore, typename TSpec, typename TConfig, typename TID, typename TAnnoIDs>
inline void
assignToReadAnnoStore(TReadAnnoStore &readAnnoStore, FragmentStore<TSpec, TConfig> & fragStore, TID readID, TAnnoIDs &annoIDs)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename Iterator<TAnnoIDs>::Type				TAnnoIDsIter;	
	typedef typename Value<TAnnoIDs>::Type					TIDs;
	typedef typename Iterator<TIDs>::Type					TIDsIter;
	
	static const TID INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	// assign annotationIDs:
	value(readAnnoStore, readID).annoIDs = annoIDs;	
		
	// assign coresponding parentIDs:
	clear (value(readAnnoStore, readID).parentIDs);
	if(!empty(annoIDs))
	{
		TIDsIter itID = begin(front(annoIDs));
		TIDsIter itIDEnd = end(front(annoIDs));
		for ( ; itID != itIDEnd; goNext(itID))		// read maps in gene, if all intervals map in gene: at least one exon of the gene has to occur in the id-list of the first interval
			if (getValue(itID) != INVALID_ID && !isElement_unsorted(getValue(fragStore.annotationStore, getValue(itID)).parentID, getValue(readAnnoStore, readID).parentIDs))	
				appendValue(value(readAnnoStore, readID).parentIDs, getValue(fragStore.annotationStore, getValue(itID)).parentID, Generous() );
	
		if (!empty(getValue(readAnnoStore, readID).parentIDs))
		{
			TAnnoIDsIter itA = begin(annoIDs);
			TAnnoIDsIter itAEnd = end(annoIDs);
			goNext(itA);
			for ( ; itA != itAEnd; goNext(itA))	// not only for exon-annotations -> more than one parentID possible
			{
				itID = begin(*itA);		// for each interval of read:
				itIDEnd = end(*itA);
				for (unsigned i = 0; i < length(getValue(readAnnoStore, readID).parentIDs); ++i) // check if at least one child of the parentID occurs 
				{						
					for ( ; itID != itIDEnd; goNext(itID))
					{
						if (getValue(itID) != INVALID_ID && getValue(fragStore.annotationStore, getValue(itID)).parentID == getValue(getValue(readAnnoStore, readID).parentIDs, i) )
							break;
					}
					if (itID == itIDEnd)			 // if not, delete parentID
					{
						erase(value(readAnnoStore, readID).parentIDs, i);
						--i;
					}
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
////// buildTupleCountStore
//////////////////////////////////////////////////////////////////////////////
template <typename TID>
struct TupleCountStoreElement
{
	typedef String<TID>			TTuple;				
	typedef String<TTuple > 		TTupleList;
	typedef String<unsigned>		TTupleCounts;
	typedef String<double>			TTupleNorm;

	TTupleList		readConnections;
	TTupleCounts		readConnectionCounts;
	TTupleNorm		readConnectionNorm;
	TTupleList		matePairConnections; 
	TTupleCounts		matePairConnectionCounts;
	TTupleNorm		matePairConnectionNorm;
};


//////////////////////////////////////////////////////////////////////////////
template<typename TTupleCountStore, typename TSpec, typename TConfig, typename TReadAnnoStore>
inline void
buildTupleCountStore(TTupleCountStore & tupleCountStore, 
		     FragmentStore<TSpec, TConfig> &  fragStore, 
		     TReadAnnoStore & readAnnoStore, 
		     unsigned n, 
		     bool exact_nTuple)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos		TPos;
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore 		TReadStore;
	typedef typename Value<TReadStore>::Type 				TReadStoreElement;
	typedef typename TReadStoreElement::TID					TReadID;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TID 				TID; 
	typedef typename Iterator<TReadAnnoStore>::Type 			TReadIter;
	typedef typename Value<TReadAnnoStore>::Type				TReadAnnoStoreElement;
	typedef typename TReadAnnoStoreElement::TAnnoIDs			TAnnoIDs;
	typedef typename Iterator<TAnnoIDs>::Type 				TAnnoIDsIter;
	typedef typename Value<TAnnoIDs>::Type					TIDs;
	typedef typename Iterator<TIDs>::Type					TIDsIter;
	
	static const TReadID INVALID_READ_ID = TReadStoreElement::INVALID_ID;
	static const TID INVALID_ANNO_ID = TAnnotationStoreElement::INVALID_ID;
	
	resize(tupleCountStore, length(fragStore.annotationStore)); 
	
	bool validMate;
	TReadIter itRead = begin(readAnnoStore);
	TReadIter itReadEnd = end(readAnnoStore);
	TIDsIter itP;
	TIDsIter itPEnd;
	TAnnoIDs annoIDs;
	TAnnoIDs tupleSet;
	TReadID readID;
	TReadID matePairID;
	TReadID secReadID;
	TAnnoIDs secTupleSet;
	TID firstAnnoID1;
	TID firstAnnoID2;
	TAnnoIDsIter itTuple;
	TAnnoIDsIter itTupleEnd;
	TAnnoIDsIter itSecTuple;
	TAnnoIDsIter itSecTupleEnd;
	TAnnoIDsIter itAnnoIDs;
	TAnnoIDsIter itAnnoIDsEnd;
	TIDs matePairTuple;
	TPos beginPos1;	
	TPos endPos1;
	TPos beginPos2;
	TPos endPos2;
	unsigned pos;
	
	for ( ; itRead != itReadEnd; goNext(itRead))
	{
		if (!empty(getValue(itRead).parentIDs) )
		{ 
			itP = begin(getValue(itRead).parentIDs);
			itPEnd = end(getValue(itRead).parentIDs);
			for ( ; itP != itPEnd; goNext(itP) )
			{
				validMate = false;
				// create list of all possible tuples for current read:
				annoIDs = getValue(itRead).annoIDs;
				clear(tupleSet);
				// create all Tuple of length n:
				if (exact_nTuple && n <= length(annoIDs)) create_nTuple(tupleSet, fragStore, annoIDs, getValue(itP), n);
				// create all max-Tuple (whole read) for current parentID:
				else if (!exact_nTuple && n == 0) create_nTuple(tupleSet, fragStore, annoIDs, getValue(itP), length(annoIDs));	
				// create all tuple >= n for current parentID:
				else if (!exact_nTuple) create_Tuple(tupleSet, fragStore, annoIDs, getValue(itP), n);	
				if (!empty(tupleSet))
				{
					// create if necessary list of all possible tuples for second matepair-read:
					readID = position(itRead, readAnnoStore);
					matePairID = getValue(fragStore.readStore, readID).matePairID;
					clear(secTupleSet);
					if (matePairID != INVALID_READ_ID)
					{
						if(getValue(getValue(fragStore.matePairStore, matePairID).readID, 0) == readID)
							secReadID = getValue(getValue(fragStore.matePairStore, matePairID).readID, 1);
						else
							secReadID = getValue(getValue(fragStore.matePairStore, matePairID).readID, 0);
				
						if ( secReadID != INVALID_READ_ID )	
						{
							//if (!empty(getValue(readAnnoStore, secReadID).annoIDs)) 
							if ( isElement_unsorted(getValue(itP), getValue(readAnnoStore, secReadID).parentIDs) )	// p in parents of matepair? -> annoIDs is not empty
							{
								validMate = true;
								annoIDs = getValue(readAnnoStore, secReadID).annoIDs;
								firstAnnoID1 = front(front(tupleSet));	// ids necessary to check positions in aligment  
								firstAnnoID2 = front(front(annoIDs));	// can't be INVALID_ID, because parents was checked
							
								// check if current read-position is smaller than the position of the second read -> tuple are ordered by position
								if ( (getValue(fragStore.annotationStore, firstAnnoID1).beginPos <= 
									getValue(fragStore.annotationStore,firstAnnoID1).endPos && 
								      getValue(fragStore.annotationStore, firstAnnoID1).beginPos < 
								      	getValue(fragStore.annotationStore, firstAnnoID2).endPos) ||
								     (getValue(fragStore.annotationStore, firstAnnoID1).beginPos > 
								     	getValue(fragStore.annotationStore, firstAnnoID1).endPos && 
								      getValue(fragStore.annotationStore, firstAnnoID1).endPos < 
								      	getValue(fragStore.annotationStore, firstAnnoID2).beginPos)  ) 
								{	
									if (exact_nTuple && n <= length(annoIDs)) create_nTuple(secTupleSet, fragStore, annoIDs, getValue(itP), n);
									else if (!exact_nTuple && n == 0) create_nTuple(secTupleSet, fragStore, annoIDs, getValue(itP), length(annoIDs));		
									else if (!exact_nTuple) create_Tuple(secTupleSet, fragStore, annoIDs, getValue(itP), n);
								}
							}
						}
					}
					else validMate = true;
			
					// access to tupleCountStore for all tuple of current read:
					if (validMate)
					{
						itTuple = begin(tupleSet);
						itTupleEnd = end(tupleSet);
						for ( ; itTuple != itTupleEnd; goNext(itTuple))	
						{
							firstAnnoID1 = front(getValue(itTuple));
							erase(value(itTuple), 0);			// first id is not stored; is know by position in tupleCountStore

							// readConnections:
							if (!empty(getValue(itTuple)))
							{
								if (searchValue(pos, getValue(itTuple), getValue(tupleCountStore, firstAnnoID1).readConnections)) 
									++value(value(tupleCountStore, firstAnnoID1).readConnectionCounts, pos);
								else 
								{
									if (pos != endPosition(getValue(tupleCountStore, firstAnnoID1).readConnections) )
									{
										resizeSpace(value(tupleCountStore, firstAnnoID1).readConnections, 1, pos, pos, Generous());
										assignValue(value(tupleCountStore, firstAnnoID1).readConnections, pos, getValue(itTuple));
										insertValue(value(tupleCountStore, firstAnnoID1).readConnectionCounts, pos, 1, Generous());
									}
									else
									{
										appendValue(value(tupleCountStore, firstAnnoID1).readConnections, getValue(itTuple), Generous());
										appendValue(value(tupleCountStore, firstAnnoID1).readConnectionCounts, 1, Generous());
									}
								}
							}
							// matePairConnections: 
							if (!empty(secTupleSet))
							{
								itSecTuple = begin(secTupleSet);
								itSecTupleEnd = end(secTupleSet);
								for ( ; itSecTuple != itSecTupleEnd; goNext(itSecTuple) )
								{
									matePairTuple = getValue(itTuple);
									// INVALID_ID: sign for connection by matepair (apart from that, there are no INVALID_IDs in the list)
									appendValue(matePairTuple, INVALID_ANNO_ID, Generous());				
									if (!empty(getValue(itTuple)) && back(getValue(itTuple)) == front(getValue(itSecTuple)) )		// no id 2x allowed
									{	
										if (exact_nTuple == 0 && n == 0) erase(value(itSecTuple), 0);
										else continue;							// tupel would be created double or tupel wouldn't have the length n anymore
									}
									append(matePairTuple, getValue(itSecTuple), Generous());
					
									if (empty(getValue(itTuple))) 							
									{
										beginPos1 = getValue(fragStore.annotationStore, firstAnnoID1).beginPos;
										endPos1 = getValue(fragStore.annotationStore, firstAnnoID1).endPos;
									}
									else
									{
										beginPos1 = getValue(fragStore.annotationStore, back(getValue(itTuple))).beginPos;
										endPos1 = getValue(fragStore.annotationStore, back(getValue(itTuple))).endPos;
									}
									// begin position of first annotation in tuple of second read
									beginPos2 = getValue(fragStore.annotationStore, front(getValue(itSecTuple))).beginPos; 
									endPos2 = getValue(fragStore.annotationStore, front(getValue(itSecTuple))).endPos;
									if ( (beginPos1 <= endPos1 && endPos1 < beginPos2) ||			// no overlapping annotations allowed
									     (endPos1 < beginPos1 && beginPos1 < endPos2) )
									{
										if (searchValue(pos, matePairTuple, getValue(tupleCountStore, firstAnnoID1).matePairConnections))
											++value(value(tupleCountStore, firstAnnoID1).matePairConnectionCounts, pos);
										else 
										{
											if (pos != endPosition(getValue(tupleCountStore, firstAnnoID1).matePairConnections) )
											{
												resizeSpace(value(tupleCountStore, firstAnnoID1).matePairConnections, 1, pos, pos, Generous());
												assignValue(value(tupleCountStore, firstAnnoID1).matePairConnections, pos, matePairTuple);
												insertValue(value(tupleCountStore, firstAnnoID1).matePairConnectionCounts, pos, 1, Generous());
											}
											else
											{
												appendValue(value(tupleCountStore, firstAnnoID1).matePairConnections, matePairTuple, Generous());
												appendValue(value(tupleCountStore, firstAnnoID1).matePairConnectionCounts, 1, Generous());
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
////// buildAnnoCountStoreg++  -I../seqan/projects/library/ -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -pedantic -lrt  main.cpp   -o main
//////////////////////////////////////////////////////////////////////////////
template<typename TAnnoCountStore, typename TSpec, typename TConfig, typename TReadAnnoStore>
inline void
buildAnnoCountStore(TAnnoCountStore & annoCountStore, FragmentStore<TSpec, TConfig> & fragStore, TReadAnnoStore & readAnnoStore)
{
	typedef typename Iterator<TReadAnnoStore>::Type 			TReadIter;
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore 		TReadStore;
	typedef typename Value<TReadStore>::Type 				TReadStoreElement;
	typedef typename Value<TReadAnnoStore>::Type				TReadAnnoStoreElement;
	typedef typename TReadAnnoStoreElement::TAnnoIDs			TAnnoIDs;
	typedef typename Iterator<TAnnoIDs>::Type				TAnnoIDsIter;
	typedef typename Value<TAnnoIDs>::Type					TIDs;
	typedef typename Iterator<TIDs>::Type 					TIDsIter;
	typedef typename Value<TIDs>::Type					TID;
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	
	static const TID INVALID_READ_ID = TReadStoreElement::INVALID_ID;
	static const TID INVALID_ANNO_ID = TAnnotationStoreElement::INVALID_ID;
	
	resize(annoCountStore, length(fragStore.annotationStore), 0);
	
	TReadIter itRead = begin(readAnnoStore);
	TReadIter itReadEnd = end(readAnnoStore);
	TID readID;
	TID matePairID;
	TID secReadID = 0;
	TIDs interSecIDs;
	TIDsIter itP;
	TIDsIter itPEnd;
	TAnnoIDsIter itAnnoIDs;
	TAnnoIDsIter itAnnoIDsEnd;
	TIDsIter itID;
	TIDsIter itIDEnd;
	
	// increment for each read respective to its mapped ids the count in the annoCountStore 
	for ( ; itRead != itReadEnd; goNext(itRead))
	{
		if (!empty(getValue(itRead).annoIDs) )
		{
			readID = position(itRead, readAnnoStore);
			matePairID = getValue(fragStore.readStore, readID).matePairID;
			if (matePairID != INVALID_READ_ID)
			{
				if (getValue(getValue(fragStore.matePairStore, matePairID).readID, 0) == readID)
					secReadID = getValue(getValue(fragStore.matePairStore, matePairID).readID, 1);
				else
					secReadID = getValue(getValue(fragStore.matePairStore, matePairID).readID, 0);
			}
			// for each parentID: we just want to count annotations, in which the read mapped
			if (!empty(getValue(itRead).parentIDs))
			{
				itP = begin(getValue(itRead).parentIDs);
				itPEnd = end(getValue(itRead).parentIDs);
				for (; itP != itPEnd; goNext(itP) )
				{
					// check mate-read to prevent double counts
					if (matePairID != INVALID_READ_ID)
					{
						if (!isElement_unsorted(getValue(itP), getValue(readAnnoStore, secReadID).parentIDs) )
							continue; // if matepair read doesn't map in same parentID: no count (go to next parentID)
					
						// count annotations, which occur in both reads, shouldn't be increment for the read with the bigger readID 
						if (secReadID < readID )	 	
						{
							clear(interSecIDs);								
							// just check the periphery annotations
							// if the end of the current read mapped in a same annotation as the start of the second read:  
							if ( interSec(interSecIDs, back(getValue(itRead).annoIDs), front(getValue(readAnnoStore, secReadID).annoIDs)) ) 
							{
								for (unsigned i = 0; i < length(interSecIDs); ++i)
								{
									if (getValue(fragStore.annotationStore, getValue(interSecIDs, i) ).parentID ==  getValue(itP))
										--value(annoCountStore, getValue(interSecIDs, i));
									// decrement the corresponding count
								}
							}
							// or if the start of the current read mapped in a same annotation as the end of the second read: 
							else if ( interSec(interSecIDs, front(getValue(itRead).annoIDs), 
								  back(getValue(readAnnoStore, secReadID).annoIDs)) )	
							{
								for (unsigned i = 0; i < length(interSecIDs); ++i)
								{
									if (getValue(fragStore.annotationStore, getValue(interSecIDs, i) ).parentID ==  getValue(itP))
										--value(annoCountStore, getValue(interSecIDs, i));
								}
							}	
							if (getValue(itP) != INVALID_ANNO_ID) --value(annoCountStore, getValue(itP));
						}
					}
			
					// count for all annoIDs
					itAnnoIDs = begin(getValue(itRead).annoIDs);
					itAnnoIDsEnd = end(getValue(itRead).annoIDs);
					for ( ; itAnnoIDs != itAnnoIDsEnd; goNext(itAnnoIDs))
					{
						itID = begin(*itAnnoIDs);
						itIDEnd = end(*itAnnoIDs);
						for ( ; itID != itIDEnd; goNext(itID))
							if (getValue(itID) != INVALID_ANNO_ID && getValue(fragStore.annotationStore, getValue(itID)).parentID == getValue(itP) )
								++value(annoCountStore, getValue(itID));
					}
		
					// count for parentIDs (already selected)
					if (getValue(itP) != INVALID_ANNO_ID) ++value(annoCountStore, getValue(itP));
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
////// Overlap Module
//////////////////////////////////////////////////////////////////////////////
template<typename TReadAnnoStore, typename TAnnoCountStore, typename TTupleCountStore, typename TSpec, typename TConfig>
inline void
getResults(TReadAnnoStore & readAnnoStore,
	   TAnnoCountStore & annoCountStore,
	   TTupleCountStore & tupleCountStore, 
	   FragmentStore<TSpec, TConfig> & fragStore, 
	   unsigned tupelSize,
	   bool exact_nTuple,
	   unsigned offsetInterval,
	   unsigned thresholdGaps,
	   bool unknownO)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TID 				TID;
	typedef typename FragmentStore<TSpec, TConfig>::TIntervalTreeStore 	TIntervalTreeStore;
	typedef typename Iterator<TIntervalTreeStore>::Type			TIntervalTree;
	typedef typename Value<TReadAnnoStore>::Type 				TReadAnnoStoreElement;
	typedef typename TReadAnnoStoreElement::TAnnoIDs 			TAnnoIDs;
	
	typedef typename FragmentStore<TSpec, TConfig>::TAlignedReadStore	TAlignedReadStore;
	typedef typename Position<TAlignedReadStore>::Type 			TAlignPos;
	typedef 	 String<AlignIntervalsStoreElement<> > 			TAlignIntervalsStore;
	typedef typename Iterator<TAlignIntervalsStore>::Type 			TAlignIntervalsStoreIter;
	
	resize(readAnnoStore, length(fragStore.readStore));
	
	TIntervalTree intervalTree;
	
	// extract intervals from alignedReadStore and store them in AlignIntervalsStore:
	TAlignIntervalsStore alignIntervalsStore;
	buildAlignIntervalsStore(alignIntervalsStore, fragStore, thresholdGaps);

	if (!empty(alignIntervalsStore))
	{
		TAlignPos alignPos;
		TID contigID;
		TID readID;
		TAnnoIDs  ids;
	
		TAlignIntervalsStoreIter it = begin(alignIntervalsStore);
		TAlignIntervalsStoreIter itEnd = end(alignIntervalsStore);
		// for each item in alignIntervalsStore:
		for ( ; it != itEnd; goNext(it))
		{
			// get ids from alignedReadStore (same position as in alignIntervalsStore):
			alignPos = position(it, alignIntervalsStore);
			contigID = getValue(fragStore.alignedReadStore, alignPos).contigID;
			readID = getValue(fragStore.alignedReadStore, alignPos).readID;
			// get respective intervalTree
			if (unknownO ||  getValue(fragStore.alignedReadStore, alignPos).beginPos <= getValue(fragStore.alignedReadStore, alignPos).endPos)
				intervalTree = begin(fragStore.intervalTreeStore_F, Standard()) + contigID; 	//getValue(fragStore.intervalTreeStore_F, contigID);
			else 
				intervalTree = begin(fragStore.intervalTreeStore_R, Standard()) + contigID;        //getValue(fragStore.intervalTreeStore_R, contigID);
			
			// get annotationStore-IDs for these intervals:
			clear(ids);
			if ((*intervalTree).interval_counter != 0)
				getIDsForRead(ids, fragStore, *intervalTree, getValue(it).intervals, offsetInterval);
			// assign IDs from mapped annotations to readAnnoStore:
			value(readAnnoStore, readID).contigID = contigID;
			assignToReadAnnoStore(readAnnoStore, fragStore, readID, ids);
		}
	}
	buildAnnoCountStore(annoCountStore, fragStore, readAnnoStore);
	buildTupleCountStore(tupleCountStore, fragStore, readAnnoStore, tupelSize, exact_nTuple);
}


//////////////////////////////////////////////////////////////////////////////
/// get normalized values for annotations And get Map for Gene orientations (necessary for annotation Output)
//////////////////////////////////////////////////////////////////////////////
template<typename TAnnoNormStore, typename TMapO, typename TAnnoCountStore, typename TSpec, typename TConfig>
inline void
normalizeAnnoCounts(TAnnoNormStore &annoNormStore, TMapO &mapO, TAnnoCountStore &annoCountStore, FragmentStore<TSpec, TConfig> &fragStore)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TID				TID;
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos		TPos;
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore		TReadStore;
	typedef typename Size<TReadStore>::Type					TReadStoreSize; 
	typedef typename Iterator<TAnnotationStore>::Type			TAnnoIter;
	typedef typename Iterator<TAnnoCountStore>::Type 			TCountIter;
	typedef typename Iterator<TAnnoNormStore>::Type				TNormIter;
	typedef typename Size<TPos>::Type					TSize;
	typedef 	 String<TSize>						TChildrenLengths;
	typedef typename Iterator<String<TSize> >::Type				TLengthIter;	
	typedef 	 Pair<TID, TChildrenLengths>				TPair;
	typedef typename Value<TMapO>::Type					TPairO;
	typedef 	 Map<TPair>						TMap;
	typedef typename Iterator<TMap>::Type					TMapIter;
	
	static const TID INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	static const TPos INVALID_POS = TAnnotationStoreElement::INVALID_POS;
	
	resize(annoNormStore, length(annoCountStore), 0);
	
	TReadStoreSize readNo = length(fragStore.readStore) - length(fragStore.matePairStore);
	
	TMap map;
	clear(map);
	clear(mapO);
	
	if(!empty(fragStore.annotationStore))
	{
		TAnnoIter itA = begin(fragStore.annotationStore);
		TAnnoIter itAEnd = end(fragStore.annotationStore);
		TCountIter itC = begin(annoCountStore);
		TNormIter itN = begin(annoNormStore);
		TChildrenLengths childrenLengths;
		TPair pair;
		TPairO pairO;
		TSize length;
		
		for ( ; itA != itAEnd; goNext(itA), goNext(itC), goNext(itN) )
		{
			if (getValue(itA).beginPos == INVALID_POS && getValue(itA).parentID == INVALID_ID) 	// make entry for each gene/parent in map:
			{
				clear(childrenLengths);
				pair.i1 = position(itA, fragStore.annotationStore);
				pair.i2 = childrenLengths;
				insert(map, pair);								// for lengths
		
				pairO.i1 = position(itA, fragStore.annotationStore);
				pairO.i2 = 0;
				insert(mapO, pairO);								// for orientation
			}
			else if (getValue(itA).beginPos != INVALID_POS)						// for each exon/child: 
			{
				if (getValue(itA).beginPos <= getValue(itA).endPos)
					length = getValue(itA).endPos - getValue(itA).beginPos;
				else
					length = getValue(itA).beginPos - getValue(itA).endPos;
		
				value(itN) = ((double)1000000000 * (double)getValue(itC))/((double)readNo * (double)length);		// calculate normalized expression-value
		
				if (getValue(itA).parentID != INVALID_ID)					// append length to gene/parent lengths
				{
					appendValue(mapValue(map, getValue(itA).parentID), length, Generous());
					if (getValue(itA).beginPos > getValue(itA).endPos)
						mapValue(mapO, getValue(itA).parentID) = 1;
				}
			}
		}
	}
	
	if (!empty(map))
	{
		TMapIter itM = begin(map);
		TMapIter itMEnd = end(map);
		TSize length;
		TLengthIter itL;
		TLengthIter itLEnd;
		for ( ; itM != itMEnd; goNext(itM))	// calculate normalized gene/parent expression-values
		{
			length = 0;
			itL = begin(value(itM).i2);
			itLEnd = end(value(itM).i2);
			for ( ; itL != itLEnd; goNext(itL))
				length += getValue(itL);
			value(annoNormStore, value(itM).i1) = ((double)1000000000 * (double)getValue(annoCountStore, value(itM).i1) )/((double)readNo * (double)length);
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
/// get normalized values for tuple
//////////////////////////////////////////////////////////////////////////////
template<typename TTupleCountStore, typename TSpec, typename TConfig>
inline void
normalizeTupleCounts(TTupleCountStore &tupleCountStore, FragmentStore<TSpec, TConfig> &fragStore)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TID				TID;
	typedef typename Value<TTupleCountStore>::Type				TTupleCountStoreElement;
	typedef typename TTupleCountStoreElement::TTupleList			TTupleList;
	typedef typename TTupleCountStoreElement::TTupleCounts			TTupleCounts;
	typedef typename TTupleCountStoreElement::TTupleNorm			TTupleNorm;
	typedef typename TTupleCountStoreElement::TTuple			TTuple;
	typedef typename Iterator<TTupleCountStore>::Type 			TStoreIter;
	typedef typename Iterator<TTupleList>::Type				TTupleListIter;
	typedef typename Iterator<TTupleCounts>::Type				TCountIter;
	typedef typename Iterator<TTupleNorm>::Type				TNormIter;
	typedef typename Iterator<TTuple>::Type					TTupleIter;
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos		TPos;
	typedef typename Size<TPos>::Type					TSize;
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore		TReadStore;
	typedef typename Size<TReadStore>::Type					TReadStoreSize; 
	
	static const TID INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	TReadStoreSize readNo = length(fragStore.readStore) - length(fragStore.matePairStore);
	
	if (!empty(tupleCountStore))
	{
		TStoreIter itS = begin(tupleCountStore);
		TStoreIter itSEnd = end(tupleCountStore);
		TTupleListIter itT;
		TTupleListIter itTEnd;
		TCountIter itC;
		TNormIter itN;
		TSize tupleLength;
		TTupleIter itID;
		TTupleIter itIDEnd;
		for ( ; itS != itSEnd; goNext(itS))
		{
			// readConnections:
			resize(value(itS).readConnectionNorm, length(getValue(itS).readConnections));
			if (!empty(getValue(itS).readConnections))
			{
				itT = begin(itS->readConnections);
				itTEnd = end(itS->readConnections);
				itC = begin(itS->readConnectionCounts);
				itN = begin(itS->readConnectionNorm);
				for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
				{
					tupleLength = 0;
					itID = begin(*itT);
					itIDEnd = end(*itT);
					for ( ; itID != itIDEnd; goNext(itID))
					{	
						if (getValue(fragStore.annotationStore, *itID).beginPos <= getValue(fragStore.annotationStore, *itID).endPos)
							tupleLength += getValue(fragStore.annotationStore, *itID).endPos - getValue(fragStore.annotationStore, *itID).beginPos;
						else
							tupleLength += getValue(fragStore.annotationStore, *itID).beginPos - getValue(fragStore.annotationStore, *itID).endPos;
					}
					value(itN) = ((double)1000000000 * (double)*itC) / ((double)readNo * (double)tupleLength);
				}
			}
			// matePairConnections:
			resize(value(itS).matePairConnectionNorm, length(itS->matePairConnections));
			if (!empty(itS->matePairConnections))
			{
				itT = begin(itS->matePairConnections);
				itTEnd = end(itS->matePairConnections);
				itC = begin(itS->matePairConnectionCounts);
				itN = begin(itS->matePairConnectionNorm);
				for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
				{
					tupleLength = 0;
					itID = begin(*itT);
					itIDEnd = end(*itT);
					for ( ; itID != itIDEnd; goNext(itID))
					{
						if (*itID != INVALID_ID)
						{
							if (getValue(fragStore.annotationStore, *itID).beginPos <= getValue(fragStore.annotationStore, *itID).endPos)
								tupleLength += getValue(fragStore.annotationStore, *itID).endPos - getValue(fragStore.annotationStore, *itID).beginPos;
							else
								tupleLength += getValue(fragStore.annotationStore, getValue(itID)).beginPos - getValue(fragStore.annotationStore, *itID).endPos;
						}
					}
					value(itN) = ((double)1000000000 * (double)*itC) / ((double)readNo * (double)tupleLength);
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
/// NGS Overlapper main function
//////////////////////////////////////////////////////////////////////////////
template<typename TOptions>
inline bool
ngsOverlapper(TOptions &options)
{	
	FragmentStore<> fragStore;
#ifdef DEBUG_OVERLAP_MODULE
	SEQAN_PROTIMESTART(find1_time);	
#endif 
	// build contigStore from FASTA file
#ifdef DEBUG_OVERLAP_MODULE
	std::cout << "load Sam..." << std::endl;
#endif 
	// read aligned reads in FragmentStore from Sam files
    BamFileIn inSam(toCString(options.nameSAM));
    readRecords(fragStore, inSam);
#ifdef DEBUG_OVERLAP_MODULE
	std::cout << "  loaded " << length(fragStore.readSeqStore) << " read records" << std::endl;
#endif 

#ifdef DEBUG_OVERLAP_MODULE
	SEQAN_PROTIMESTART(find2_time);
#endif 
    // read annotations from GFF or GTF
    GffFileIn inGff(toCString(options.nameGFF));
    readRecords(fragStore, inGff);

	adjustParentEntries(fragStore);

	// create IntervalTreeStore:
#ifdef DEBUG_OVERLAP_MODULE
	SEQAN_PROTIMESTART(find3_time);
#endif 
	createIntervalTreeStore(fragStore, options.unknownO);
#ifdef DEBUG_OVERLAP_MODULE
	std::cout << "create intervalTreeStores from annotationStore took: \t" << SEQAN_PROTIMEDIFF(find3_time) << " seconds" << std::endl;
#endif 
	
	// build stores for results:
	String<ReadAnnoStoreElement<unsigned> >			readAnnoStore; 	
	String<unsigned> 					annoCountStore;
	String<TupleCountStoreElement<unsigned> >		tupleCountStore;
	String<TupleCountStoreElement_Fusion<unsigned> >	tupleCountStore_Fusion;			// additional Store, if fusion genes should be checked

	// get results with additional check for transfusion genes (will be changed later additionally) 
	if (options.fusion == 1)
		getResults_Fusion(readAnnoStore, annoCountStore, tupleCountStore, tupleCountStore_Fusion, fragStore, options.nTuple, options.exact_nTuple, options.offsetInterval, options.thresholdGaps, options.unknownO);
	else // get normal results:
		getResults(readAnnoStore, annoCountStore, tupleCountStore, fragStore, options.nTuple, options.exact_nTuple, options.offsetInterval, options.thresholdGaps, options.unknownO);
	
	
	// normalize:
	String<double>			annoNormStore;
	Map<Pair<unsigned, bool> > 	mapO;
	normalizeAnnoCounts(annoNormStore, mapO, annoCountStore, fragStore);
	normalizeTupleCounts(tupleCountStore, fragStore);
	if (options.fusion == 1)
		normalizeTupleCounts_Fusion(tupleCountStore_Fusion, fragStore);


	// output:
	std::fstream readOutput;
	readOutput.open(toCString(options.readOutputFileName), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    if (!readOutput.good())
        std::cerr << "ERROR:  Could not open output file " << options.readOutputFileName << "\n";
	createReadCountGFF(readOutput, readAnnoStore, fragStore);
	readOutput.close();
	
	std::fstream annoOutput;
	annoOutput.open(toCString(options.annoOutputFileName), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    if (!annoOutput.good())
        std::cerr << "ERROR:  Could not open output file " << options.annoOutputFileName << "\n";
	createAnnoCountGFF(annoOutput, annoCountStore, annoNormStore, fragStore, mapO);
	annoOutput.close();
	
	std::fstream tupleOutput;
	tupleOutput.open(toCString(options.tupleOutputFileName), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    if (!tupleOutput.good())
        std::cerr << "ERROR:  Could not open output file " << options.tupleOutputFileName << "\n";
	createTupleCountGFF(tupleOutput, tupleCountStore, fragStore, options.thresholdCount, options.thresholdRPKM);
	tupleOutput.close();

	// additional output, if fusion genes were checked
	if (options.fusion == 1)
	{
		std::fstream tupleOutput_Fusion;
		tupleOutput_Fusion.open(toCString(options.tupleFusionOutputFileName), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
		createTupleCountGFF_Fusion(tupleOutput_Fusion, tupleCountStore_Fusion, fragStore, options.thresholdCount, options.thresholdRPKM);
		tupleOutput_Fusion.close();
	}

	
#ifdef DEBUG_OVERLAP_MODULE
	std::cout << "ngsOverlapper-function took: \t" << SEQAN_PROTIMEDIFF(find1_time) << " seconds" << std::endl;
	std::cout << "ngsOverlapper-function without reading Sam took: \t" << SEQAN_PROTIMEDIFF(find2_time) << " seconds" << std::endl;
	std::cout << "ngsOverlapper-function and create IntervalTreeStore without reading Sam took:\t" << SEQAN_PROTIMEDIFF(find2_time) - SEQAN_PROTIMEDIFF(find3_time) << " seconds" << std::endl;
#endif 
	return true;
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
