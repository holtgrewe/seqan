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
  fragStoreRCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $ID$
 ==========================================================================*/

#ifndef SEQAN_HEADER_FUSION_H
#define SEQAN_HEADER_FUSION_H
//#define DEBUG_OVERLAP_MODULE

namespace SEQAN_NAMESPACE_MAIN
{ 

//////////////////////////////////////////////////////////////////////////////
////// buildTupleCountStore_Fusion 
//////////////////////////////////////////////////////////////////////////////
template <typename TID>
struct TupleCountStoreElement_Fusion
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
template<typename TTupleCountStore, typename TTupleCountStore_Fusion, typename TSpec, typename TConfig, typename TReadAnnoStore>
inline void
buildTupleCountStore_Fusion(TTupleCountStore & tupleCountStore, 
		     TTupleCountStore_Fusion & tupleCountStore_Fusion, 
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
	resize(tupleCountStore_Fusion, length(fragStore.annotationStore)); 
	
	bool validMate;
	TReadIter itRead = begin(readAnnoStore);
	TReadIter itReadEnd = end(readAnnoStore);
	TIDsIter itP;
	TIDsIter itPEnd;
	TIDsIter itP2;
	TIDsIter itP2End;
	TAnnoIDs annoIDs;
	TAnnoIDs tupleSet;
	TReadID readID;
	TReadID matePairID;
	TReadID secReadID;
	TAnnoIDs secTupleSet;
	TAnnoIDs tempSecTupleSet;
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
							else if (!empty(getValue(readAnnoStore, secReadID).parentIDs)) 	// parent of matepair different -> possible transfusion
							{
								//validMate = true;
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
									itP2 = begin(value(readAnnoStore, secReadID).parentIDs);
									itP2End = end(value(readAnnoStore, secReadID).parentIDs);
									clear(tempSecTupleSet);
									for ( ; itP2 != itP2End; goNext(itP2))
									{
										if (exact_nTuple && n <= length(annoIDs)) create_nTuple(tempSecTupleSet, fragStore, annoIDs, getValue(itP2), n);
										else if (!exact_nTuple && n == 0) create_nTuple(tempSecTupleSet, fragStore, annoIDs, getValue(itP2), length(annoIDs));		
										else if (!exact_nTuple) create_Tuple(tempSecTupleSet, fragStore, annoIDs, getValue(itP2), n);

										itSecTuple = begin(tempSecTupleSet);
										itSecTupleEnd = end(tempSecTupleSet);
										for ( ; itSecTuple != itSecTupleEnd; goNext(itSecTuple))
										{										
											appendValue(secTupleSet, getValue(itSecTuple));
										}
									}
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
					else if (!empty(secTupleSet))	// Transfusion: matepairs with different parents
					{
						itTuple = begin(tupleSet);
						itTupleEnd = end(tupleSet);
						for ( ; itTuple != itTupleEnd; goNext(itTuple))	
						{
							firstAnnoID1 = front(getValue(itTuple));
							erase(value(itTuple), 0);			// first id is not stored; is know by position in tupleCountStore
			
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
									if (searchValue(pos, matePairTuple, getValue(tupleCountStore_Fusion, firstAnnoID1).matePairConnections))
									{
										++value(value(tupleCountStore_Fusion, firstAnnoID1).matePairConnectionCounts, pos);
									}
									else 
									{
										if (pos != endPosition(getValue(tupleCountStore_Fusion, firstAnnoID1).matePairConnections) )
										{
											resizeSpace(value(tupleCountStore_Fusion, firstAnnoID1).matePairConnections, 1, pos, pos, Generous());
											assignValue(value(tupleCountStore_Fusion, firstAnnoID1).matePairConnections, pos, matePairTuple);
											insertValue(value(tupleCountStore_Fusion, firstAnnoID1).matePairConnectionCounts, pos, 1, Generous());
										}
										else
										{
											appendValue(value(tupleCountStore_Fusion, firstAnnoID1).matePairConnections, matePairTuple, Generous());
											appendValue(value(tupleCountStore_Fusion, firstAnnoID1).matePairConnectionCounts, 1, Generous());
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
//create tupleCountGFF_Fusion
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TTupleCountStore_Fusion, typename TSpec, typename TConfig>
inline void
createTupleCountGFF_Fusion(TFile & tupleOutput_Fusion, TTupleCountStore_Fusion & tupleCountStore_Fusion, FragmentStore<TSpec, TConfig> & fragStore, unsigned thresholdCount, double thresholdRPKM)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TID 				TID;
	
	typedef typename Iterator<TTupleCountStore_Fusion>::Type 		TCountStoreIter;
	typedef typename Value<TTupleCountStore_Fusion>::Type			TTupleCountStoreElement_Fusion;
	typedef typename TTupleCountStoreElement_Fusion::TTupleList		TTupleList;
	typedef typename TTupleCountStoreElement_Fusion::TTupleCounts		TTupleCounts;
	typedef typename TTupleCountStoreElement_Fusion::TTupleNorm		TTupleNorm;
	typedef typename Value<TTupleList>::Type				TTupel;
	typedef typename Iterator<TTupleList>::Type				TTupleListIter;
	typedef typename Iterator<TTupleCounts>::Type				TCountIter;
	typedef typename Iterator<TTupleNorm>::Type				TNormIter;
	typedef typename Iterator<TTupel>::Type					TTupelIter;

	static const TID INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	if (!empty(tupleCountStore_Fusion))
	{
		TCountStoreIter itCountStore = begin(tupleCountStore_Fusion);
		TCountStoreIter itCountStoreEnd = end(tupleCountStore_Fusion);
		TAnnotationStoreElement currentElement;
		TTupleListIter itT;
		TTupleListIter itTEnd;
		TCountIter itC;
		TNormIter itN;
		TTupelIter itID;
		TTupelIter itIDEnd;
		for ( ; itCountStore != itCountStoreEnd; goNext(itCountStore))
		{
			currentElement = getValue(fragStore.annotationStore, position(itCountStore, tupleCountStore_Fusion));

			/*
			itT = begin(getValue(itCountStore).readConnections);
			itTEnd = end(getValue(itCountStore).readConnections);
			itC = begin(getValue(itCountStore).readConnectionCounts);
			itN = begin(getValue(itCountStore).readConnectionNorm);
			// read connections:
			for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
			{
				
				if (getValue(itC) >= thresholdCount && getValue(itN) >= thresholdRPKM)
				{
					// contig-name
					streamPut(tupleOutput_Fusion, getValue(fragStore.contigNameStore, currentElement.contigID));
					itID = begin(getValue(itT));
					itIDEnd = end(getValue(itT));
					for ( ; itID != itIDEnd; goNext(itID))
					{
						if (getValue(itID) == INVALID_ID)
						{
							goNext(itID);
							streamPut(tupleOutput_Fusion, "~");
							streamPut(tupleOutput_Fusion, getValue(fragStore.contigNameStore, getValue(fragStore.annotationStore, getValue(itID)).contigID));
						}
					}
					streamPut(tupleOutput_Fusion, '\t');					

					// parent-names
					streamPut(tupleOutput_Fusion, getValue(fragStore.annotationNameStore, currentElement.parentID));
					itID = begin(getValue(itT));
					for ( ; itID != itIDEnd; goNext(itID))
					{
						if (getValue(itID) == INVALID_ID)
						{
							goNext(itID);
							streamPut(tupleOutput_Fusion, "~");
							streamPut(tupleOutput_Fusion, getValue(fragStore.annotationNameStore, getValue(fragStore.annotationStore, getValue(itID)).parentID));
						}
					}
					streamPut(tupleOutput_Fusion, '\t');
					
					// orientation
					if ( currentElement.beginPos <= currentElement.endPos )
					{
						streamPut(tupleOutput_Fusion, "+");
					}
					else
					{
						streamPut(tupleOutput_Fusion, "-");
					}
					itID = begin(getValue(itT));
					for ( ; itID != itIDEnd; goNext(itID))
					{
						if (getValue(itID) == INVALID_ID)
						{
							goNext(itID);
							streamPut(tupleOutput_Fusion, "~");
							if ( getValue(fragStore.annotationStore, getValue(itID)).beginPos <= getValue(fragStore.annotationStore, getValue(itID)).endPos )
							{
								streamPut(tupleOutput_Fusion, "+");
							}
							else
							{
								streamPut(tupleOutput_Fusion, "-");
							}
						}
					}
					streamPut(tupleOutput_Fusion, "\t");

					// first annotationID of tuple (store implicit)
					streamPut(tupleOutput_Fusion, getValue(fragStore.annotationNameStore, position(itCountStore, tupleCountStore_Fusion)));
					// other annotationIDs
					itID = begin(getValue(itT));
					itIDEnd = end(getValue(itT));
					for ( ; itID != itIDEnd ; goNext(itID))
					{
						if (getValue(itID) != INVALID_ID)
						{
							streamPut(tupleOutput_Fusion, ":");
							streamPut(tupleOutput_Fusion, getValue(fragStore.annotationNameStore, getValue(itID)));
						}
						else
						{
							streamPut(tupleOutput_Fusion, "~");
							goNext(itID);
							streamPut(tupleOutput_Fusion, getValue(fragStore.annotationNameStore, getValue(itID)));
						}
					}
					streamPut(tupleOutput_Fusion, '\t');
					// tuple count
					streamPut(tupleOutput_Fusion, getValue(itC));
					streamPut(tupleOutput_Fusion, '\t');
					// normalized tuple count
					_streamPutDouble(tupleOutput_Fusion, getValue(itN));
					streamPut(tupleOutput_Fusion, '\n');
				}
			}
			*/
	
			// matePairConnections:
			itT = begin(getValue(itCountStore).matePairConnections);
			itTEnd = end(getValue(itCountStore).matePairConnections);
			itC = begin(getValue(itCountStore).matePairConnectionCounts);
			itN = begin(getValue(itCountStore).matePairConnectionNorm);
			for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
			{if (getValue(itC) >= thresholdCount && getValue(itN) >= thresholdRPKM)
				{
					// contig-name
					tupleOutput_Fusion << getValue(fragStore.contigNameStore, currentElement.contigID);
					itID = begin(getValue(itT));
					itIDEnd = end(getValue(itT));
					for ( ; itID != itIDEnd; goNext(itID))
					{
						if (getValue(itID) == INVALID_ID)
						{
							goNext(itID);
							tupleOutput_Fusion << '^'
                                               << getValue(fragStore.contigNameStore,
                                                           getValue(fragStore.annotationStore,
                                                                    getValue(itID)).contigID);
						}
					}
					tupleOutput_Fusion << '\t';
				
					// parent-name
					if (currentElement.parentID == INVALID_ID )
						tupleOutput_Fusion << "NO_PARENT";
					else
						tupleOutput_Fusion << getValue(fragStore.annotationNameStore, currentElement.parentID);
					itID = begin(getValue(itT));
					for ( ; itID != itIDEnd; goNext(itID))
					{
						if (getValue(itID) == INVALID_ID)
						{
							goNext(itID);
							if (getValue(fragStore.annotationStore, getValue(itID)).parentID == INVALID_ID)
								tupleOutput_Fusion << "^NO_PARENT";
							else
								tupleOutput_Fusion << "^"
                                                   << getValue(fragStore.annotationNameStore, getValue(fragStore.annotationStore, getValue(itID)).parentID);
						}
					}
					tupleOutput_Fusion << '\t';

					// orientation
					if ( currentElement.beginPos <= currentElement.endPos )
						tupleOutput_Fusion << "+";
					else
						tupleOutput_Fusion << "-";
					itID = begin(getValue(itT));
					for ( ; itID != itIDEnd; goNext(itID))
					{
						if (getValue(itID) == INVALID_ID)
						{
							goNext(itID);
							tupleOutput_Fusion << "^";
							if ( getValue(fragStore.annotationStore, getValue(itID)).beginPos <= getValue(fragStore.annotationStore, getValue(itID)).endPos )
								tupleOutput_Fusion << "+";
							else
								tupleOutput_Fusion << "-";
						}
					}
					tupleOutput_Fusion << "\t";
				
					// first annotationID of tuple
					tupleOutput_Fusion << getValue(fragStore.annotationNameStore, position(itCountStore, tupleCountStore_Fusion));
			
					// other annotationIDs of first read
					itID = begin(getValue(itT));
					itIDEnd = end(getValue(itT));
					for ( ; itID != itIDEnd && getValue(itID) != INVALID_ID; goNext(itID))
					{
						tupleOutput_Fusion << ":" << getValue(fragStore.annotationNameStore, getValue(itID));
					}
					goNext(itID);
					tupleOutput_Fusion << "^";
		
					// annotationIDs of second read
					tupleOutput_Fusion << getValue(fragStore.annotationNameStore, getValue(itID));
					goNext(itID);
					for ( ; itID != itIDEnd; goNext(itID))
						tupleOutput_Fusion << ":" << getValue(fragStore.annotationNameStore, getValue(itID));
					tupleOutput_Fusion << '\t';
					
					// tuple count
					tupleOutput_Fusion << *itC << '\t';
					// normalized tuple count
					tupleOutput_Fusion << formattedNumber("%f", *itN) << '\n';
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
/// get normalized values for tuple of transfusion genes
//////////////////////////////////////////////////////////////////////////////
template<typename TTupleCountStore_Fusion, typename TSpec, typename TConfig>
inline void
normalizeTupleCounts_Fusion(TTupleCountStore_Fusion &tupleCountStore_Fusion, FragmentStore<TSpec, TConfig> &fragStore)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TID			TID;
	typedef typename Value<TTupleCountStore_Fusion>::Type			TTupleCountStoreElement_Fusion;
	typedef typename TTupleCountStoreElement_Fusion::TTupleList		TTupleList;
	typedef typename TTupleCountStoreElement_Fusion::TTupleCounts		TTupleCounts;
	typedef typename TTupleCountStoreElement_Fusion::TTupleNorm		TTupleNorm;
	typedef typename TTupleCountStoreElement_Fusion::TTuple			TTuple;
	typedef typename Iterator<TTupleCountStore_Fusion>::Type 		TStoreIter;
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

	if (!empty(tupleCountStore_Fusion))
	{
		TStoreIter itS = begin(tupleCountStore_Fusion);
		TStoreIter itSEnd = end(tupleCountStore_Fusion);
		TTupleListIter itT;
		TTupleListIter itTEnd;
		TCountIter itC;
		TNormIter itN;
		TSize tupleLength;
		TTupleIter itID;
		TTupleIter itIDEnd;
		for ( ; itS != itSEnd; goNext(itS))
		{
			/*
			// readConnections:
			resize(value(itS).readConnectionNorm, length(getValue(itS).readConnections));
			if (!empty(getValue(itS).readConnections))
			{
				itT = begin(getValue(itS).readConnections);
				itTEnd = end(getValue(itS).readConnections);
				itC = begin(getValue(itS).readConnectionCounts);
				itN = begin(getValue(itS).readConnectionNorm);
				for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
				{
					tupleLength = 0;
					itID = begin(getValue(itT));
					itIDEnd = end(getValue(itT));
					for ( ; itID != itIDEnd; goNext(itID))
					{	
						if (getValue(itID) != INVALID_ID)
						{
							if (getValue(fragStore.annotationStore, getValue(itID)).beginPos <= getValue(fragStore.annotationStore, getValue(itID)).endPos)
								tupleLength += getValue(fragStore.annotationStore, getValue(itID)).endPos - getValue(fragStore.annotationStore, getValue(itID)).beginPos;
							else
								tupleLength += getValue(fragStore.annotationStore, getValue(itID)).beginPos - getValue(fragStore.annotationStore, getValue(itID)).endPos;
						}
					}
					value(itN) = ((double)1000000000 * (double)getValue(itC)) / ((double)readNo * (double)tupleLength);
				}
			}
			*/
			// matePairConnections:
			resize(value(itS).matePairConnectionNorm, length(getValue(itS).matePairConnections));
			if (!empty(getValue(itS).matePairConnections))
			{
				itT = begin(getValue(itS).matePairConnections);
				itTEnd = end(getValue(itS).matePairConnections);
				itC = begin(getValue(itS).matePairConnectionCounts);
				itN = begin(getValue(itS).matePairConnectionNorm);
				for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
				{
					tupleLength = 0;
					itID = begin(getValue(itT));
					itIDEnd = end(getValue(itT));
					for ( ; itID != itIDEnd; goNext(itID))
					{
						if (getValue(itID) != INVALID_ID)
						{
							if (getValue(fragStore.annotationStore, getValue(itID)).beginPos <= getValue(fragStore.annotationStore, getValue(itID)).endPos)
								tupleLength += getValue(fragStore.annotationStore, getValue(itID)).endPos - getValue(fragStore.annotationStore, getValue(itID)).beginPos;
							else
								tupleLength += getValue(fragStore.annotationStore, getValue(itID)).beginPos - getValue(fragStore.annotationStore, getValue(itID)).endPos;
						}
					}
					value(itN) = ((double)1000000000 * (double)getValue(itC)) / ((double)readNo * (double)tupleLength);
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
////// Overlap Module
//////////////////////////////////////////////////////////////////////////////
template<typename TReadAnnoStore, typename TAnnoCountStore, typename TTupleCountStore, typename TTupleCountStore_Fusion, typename TSpec, typename TConfig>
inline void
getResults_Fusion(TReadAnnoStore & readAnnoStore,
	   TAnnoCountStore & annoCountStore,
	   TTupleCountStore & tupleCountStore, 
	   TTupleCountStore_Fusion & tupleCountStore_Fusion, 
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
	buildTupleCountStore_Fusion(tupleCountStore, tupleCountStore_Fusion,  fragStore, readAnnoStore, tupelSize, exact_nTuple);
}



//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
