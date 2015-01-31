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

#ifndef SEQAN_HEADER_CREATE_GFF_H
#define SEQAN_HEADER_CREATE_GFF_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//create readCountGFF
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TReadAnnoStore, typename TSpec, typename TConfig>
inline void
createReadCountGFF(TFile & readOutput, TReadAnnoStore & readAnnoStore, FragmentStore<TSpec, TConfig> & fragStore)
{	
	typedef typename Iterator<TReadAnnoStore>::Type 			TCountIter;
	typedef typename Value<TReadAnnoStore>::Type				TReadAnnoStoreElement;
	typedef typename TReadAnnoStoreElement::TAnnoIDs			TAnnoIDs;
	typedef typename Iterator<TAnnoIDs>::Type				TAnnoIDsIter;
	typedef typename Value<TAnnoIDs>::Type					TIDs;
	typedef typename Iterator<TIDs>::Type					TIDsIter;	
					
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TID 				TID;
	
	static const TID INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	TCountIter itCountStore = begin(readAnnoStore);
	TCountIter itCountStoreEnd = end(readAnnoStore);
	TAnnoIDsIter itAnnoIDs;
	TAnnoIDsIter itAnnoIDsEnd;
	TID firstID;
	TIDs allParentIDs;
	TIDsIter itID;
	TIDsIter itIDEnd;
	TIDsIter itP;
	TIDsIter itPEnd;
	bool help;
	bool invalid;
	for ( ; itCountStore != itCountStoreEnd; goNext(itCountStore))
	{
		// read-name:
        readOutput << getValue(fragStore.readNameStore, position(itCountStore, readAnnoStore))
                   << '\t';
		
		if (empty(getValue(itCountStore).annoIDs) )
		{
            readOutput << ".\t.\t.\t.\t.";
		}
		else
		{
			// contig-name:
			readOutput << getValue(fragStore.contigNameStore, getValue(itCountStore).contigID)
                                   << '\t';
			
			itAnnoIDs = begin(getValue(itCountStore).annoIDs);
			itAnnoIDsEnd = end(getValue(itCountStore).annoIDs);
			while (itAnnoIDs != itAnnoIDsEnd && front(getValue(itAnnoIDs)) == INVALID_ID)
			{
				goNext(itAnnoIDs);
			}
	
			if (itAnnoIDs != itAnnoIDsEnd) // not only INVALID_IDS
			{
				firstID = front(getValue(itAnnoIDs));
				
				// orientation:
				if (getValue(fragStore.annotationStore, firstID).beginPos <= getValue(fragStore.annotationStore, firstID).endPos)
				{
                    readOutput << "+\t";
				}
				else
                {
                    readOutput << "-\t";
                }
				
				// Annotation-IDs:
				
				allParentIDs =  getValue(itCountStore).parentIDs;
				// output for first parentID if possible; get other parentIDs, in which the read doesn't map (entirely with all of his intervals)
				itAnnoIDs = begin(getValue(itCountStore).annoIDs);
				itAnnoIDsEnd = end(getValue(itCountStore).annoIDs);
				for ( ; itAnnoIDs != itAnnoIDsEnd; goNext(itAnnoIDs))
				{
					help = false;
					itID = begin(*itAnnoIDs);
					itIDEnd = end(*itAnnoIDs);
					for ( ; itID != itIDEnd; goNext(itID))
					{
						if (!empty(getValue(itCountStore).parentIDs) && getValue(itID) != INVALID_ID && 	// if current parentID == first id of parentIDs (of read entry in readAnnoStore)
						    getValue(fragStore.annotationStore, getValue(itID)).parentID == front(allParentIDs))
						{
							if (help) 	// not the first annotation for this read-interval -> ";" sign for overlapping annotations
                                readOutput << ';';
                            readOutput << getValue(fragStore.annotationNameStore, getValue(itID));
							help = true;
						}
						else if (getValue(itID) != INVALID_ID && getValue(fragStore.annotationStore, getValue(itID)).parentID != INVALID_ID &&	// get other parentIDs 
							!isElement_unsorted(getValue(fragStore.annotationStore, getValue(itID)).parentID, allParentIDs) ) //?
						{
							appendValue(allParentIDs, getValue(fragStore.annotationStore, getValue(itID)).parentID, Generous() );
						}
					}
					if ( !empty(getValue(itCountStore).parentIDs) && position(itAnnoIDs, getValue(itCountStore).annoIDs) != endPosition(getValue(itCountStore).annoIDs) - 1)
                        readOutput << ':';
				}
				if (!empty(getValue(itCountStore).parentIDs))
					readOutput << '\t' << getValue(fragStore.annotationNameStore, front(allParentIDs)) << '\t';
				// outputs for all other parentIDs
				itP = begin(allParentIDs);
				itPEnd = end(allParentIDs);
				if (!empty(getValue(itCountStore).parentIDs)) goNext(itP);
				for ( ; itP != itPEnd; goNext(itP))
				{
					itAnnoIDs = begin(getValue(itCountStore).annoIDs);
					itAnnoIDsEnd = end(getValue(itCountStore).annoIDs);
					for ( ; itAnnoIDs != itAnnoIDsEnd; goNext(itAnnoIDs))
					{
						invalid = true;				// if no annotation for the current parent in  interval -> UNKOWN_REGION
						itID = begin(*itAnnoIDs);	
                        itIDEnd = end(*itAnnoIDs);
						for ( ; itID != itIDEnd; goNext(itID))
						{
							if (getValue(itID) != INVALID_ID && getValue(fragStore.annotationStore, getValue(itID)).parentID == getValue(itP))
							{
								if (!invalid)	// not the first annotation for this interval -> ";" sign for overlapping annotations
                                    readOutput << ';';
								readOutput << getValue(fragStore.annotationNameStore, getValue(itID));
								invalid = false;
							}
						}
						if (invalid)
                            readOutput << "UNKNOWN_REGION";
						if (position(itAnnoIDs, getValue(itCountStore).annoIDs) != endPosition(getValue(itCountStore).annoIDs) - 1)
                            readOutput << ':';
					}
                    readOutput << '\t';
					if (getValue(itP) != INVALID_ID)
                        readOutput << getValue(fragStore.annotationNameStore, getValue(itP));
					else
                        readOutput << "NO_PARENT";
                    readOutput << '\t';
				}
			}
			else  // only INVALID_IDS
			{	
                readOutput << ".\t";
						
				// invalid_ids for each interval
				itAnnoIDs = begin(getValue(itCountStore).annoIDs);
				itAnnoIDsEnd = end(getValue(itCountStore).annoIDs);
				for ( ; itAnnoIDs != itAnnoIDsEnd; goNext(itAnnoIDs))
				{
                    readOutput << "UNKNOWN_REGION";
					if (position(itAnnoIDs, getValue(itCountStore).annoIDs) != endPosition(getValue(itCountStore).annoIDs) - 1)
                        readOutput << ':';
				}
                readOutput << "\tUNKNOWN_REGION";
			}
		}
        readOutput << "\n";
	}	
}



//////////////////////////////////////////////////////////////////////////////
//create AnnoCountGFF
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TAnnoCountStore, typename TAnnoNormStore, typename TSpec, typename TConfig, typename TMap>
inline void
createAnnoCountGFF(TFile & annoOutput, TAnnoCountStore & annoCountStore, TAnnoNormStore &annoNormStore, FragmentStore<TSpec, TConfig> & fragStore, TMap &mapO)
{
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos 		TContigPos;
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TID 				TID;
	typedef typename Iterator<TAnnoCountStore>::Type 			TCountIter;
	typedef typename Iterator<TAnnoNormStore>::Type				TNormIter;
	typedef typename Iterator<TAnnotationStore>::Type 			TAnnoIter;
	
	static const TID INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	static const TContigPos INVALID_POS = TAnnotationStoreElement::INVALID_POS;

	TCountIter itCount = begin(annoCountStore);
	TCountIter itCountEnd = end(annoCountStore);
	TAnnoIter itAnno = begin(fragStore.annotationStore);
	TNormIter itNorm = begin(annoNormStore);
	
	for ( ; itCount != itCountEnd; goNext(itCount), goNext(itAnno), goNext(itNorm))
	{
        if (getValue(itAnno).typeID != INVALID_ID)
	        if (fragStore.annotationTypeStore[getValue(itAnno).typeID] == "<root>") continue;  
		// contig-name
		if (getValue(itAnno).contigID == INVALID_ID )
			annoOutput << "INVALID_ID\t";
		else
			annoOutput << getValue(fragStore.contigNameStore, getValue(itAnno).contigID) << '\t';
		annoOutput << "Annotation_Count\tregion\t";
		// startposition endposition orientation . 
		if (getValue(itAnno).beginPos == INVALID_POS)
		{
			annoOutput << ".\t.\t" << getValue(itCount);
			if (getValue(itAnno).parentID == INVALID_ID || (fragStore.annotationStore[getValue(itAnno).parentID].typeID != INVALID_ID && fragStore.annotationTypeStore[fragStore.annotationStore[getValue(itAnno).parentID].typeID] == "<root>"))
			{
				if (mapValue(mapO, position(itAnno, fragStore.annotationStore)) == 0)
					annoOutput << "\t+\t.\t";
				else
					annoOutput << "\t-\t.\t";
			}
			else
            {
                annoOutput << "\t.\t.\t";
            }
		}
		else
		{
			if (getValue(itAnno).beginPos <= getValue(itAnno).endPos)
				annoOutput << getValue(itAnno).beginPos + 1
				           << '\t'
				           << getValue(itAnno).endPos
				           << '\t'
				           << getValue(itCount)
				           << "\t+\t.\t";
			else
				annoOutput << getValue(itAnno).endPos + 1
				           << '\t'
				           << getValue(itAnno).beginPos
				           << '\t'
				           << getValue(itCount)
				           << "\t-\t.\t";
		}
		// annotation-name (parent annotation-name)
		if (getValue(itAnno).parentID == INVALID_ID || (fragStore.annotationStore[getValue(itAnno).parentID].typeID != INVALID_ID && fragStore.annotationTypeStore[fragStore.annotationStore[getValue(itAnno).parentID].typeID] == "<root>"))
			annoOutput << "ID="
			           << getValue(fragStore.annotationNameStore, position(itAnno, fragStore.annotationStore))
			           << ';';
		else
			annoOutput << "ID="
			           << getValue(fragStore.annotationNameStore, position(itAnno, fragStore.annotationStore))
			           << ";ParentID="
			           << getValue(fragStore.annotationNameStore, getValue(itAnno).parentID)
			           << ';';
        annoOutput << formattedNumber("%f", *itNorm) << ";\n";
	}
}


//////////////////////////////////////////////////////////////////////////////
//create tupleCountGFF
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TTupleCountStore, typename TSpec, typename TConfig>
inline void
createTupleCountGFF(TFile & tupleOutput, TTupleCountStore & tupleCountStore, FragmentStore<TSpec, TConfig> & fragStore, unsigned thresholdCount, double thresholdRPKM)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TID 				TID;
	
	typedef typename Iterator<TTupleCountStore>::Type 			TCountStoreIter;
	typedef typename Value<TTupleCountStore>::Type				TTupleCountStoreElement;
	typedef typename TTupleCountStoreElement::TTupleList			TTupleList;
	typedef typename TTupleCountStoreElement::TTupleCounts			TTupleCounts;
	typedef typename TTupleCountStoreElement::TTupleNorm			TTupleNorm;
	typedef typename Value<TTupleList>::Type				TTupel;
	typedef typename Iterator<TTupleList>::Type				TTupleListIter;
	typedef typename Iterator<TTupleCounts>::Type				TCountIter;
	typedef typename Iterator<TTupleNorm>::Type				TNormIter;
	typedef typename Iterator<TTupel>::Type					TTupelIter;
	
	static const TID INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	if (!empty(tupleCountStore))
	{
		TCountStoreIter itCountStore = begin(tupleCountStore);
		TCountStoreIter itCountStoreEnd = end(tupleCountStore);
		TAnnotationStoreElement currentElement;
		TTupleListIter itT;
		TTupleListIter itTEnd;
		TCountIter itC;
		TNormIter itN;
		TTupelIter itID;
		TTupelIter itIDEnd;
		for ( ; itCountStore != itCountStoreEnd; goNext(itCountStore))
		{
			currentElement = getValue(fragStore.annotationStore, position(itCountStore, tupleCountStore));
	
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
					tupleOutput << getValue(fragStore.contigNameStore, currentElement.contigID) << '\t';
					// parent-name
					if (currentElement.parentID == INVALID_ID )
						tupleOutput << "NO_PARENT\t";
					else
						tupleOutput << getValue(fragStore.annotationNameStore, currentElement.parentID) << '\t';
					// orientation
					if ( currentElement.beginPos <= currentElement.endPos )
						tupleOutput << "+\t";
					else
						tupleOutput << "-\t";
					// first annotationID of tuple (store implicit)
					tupleOutput << getValue(fragStore.annotationNameStore, position(itCountStore, tupleCountStore));
					// other annotationIDs
					itID = begin(getValue(itT));
					itIDEnd = end(getValue(itT));
					for ( ; itID != itIDEnd; goNext(itID))
						tupleOutput << ":" << getValue(fragStore.annotationNameStore, getValue(itID));
					tupleOutput << '\t';
					// tuple count
					tupleOutput << getValue(itC) << '\t';
					// normalized tuple count
                    tupleOutput << formattedNumber("%f", *itN) << '\n';
				}
			}
			//matepair connections:
			itT = begin(getValue(itCountStore).matePairConnections);
			itTEnd = end(getValue(itCountStore).matePairConnections);
			itC = begin(getValue(itCountStore).matePairConnectionCounts);
			itN = begin(getValue(itCountStore).matePairConnectionNorm);
			for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
			{
				if (getValue(itC) >= thresholdCount && getValue(itN) >= thresholdRPKM)
				{
					// contig-name
					tupleOutput << getValue(fragStore.contigNameStore, currentElement.contigID) << '\t';
					// parent-name
					if (currentElement.parentID == INVALID_ID )
						tupleOutput << "NO_PARENT\t";
					else
						tupleOutput << getValue(fragStore.annotationNameStore, currentElement.parentID) << '\t';
					// orientation
					if ( currentElement.beginPos <= currentElement.endPos )
						tupleOutput << "+\t";
					else
						tupleOutput << "-\t";
					// first annotationID of tuple
					tupleOutput << getValue(fragStore.annotationNameStore, position(itCountStore, tupleCountStore));
					// other annotationIDs of first read
					itID = begin(getValue(itT));
					itIDEnd = end(getValue(itT));
					for ( ; itID != itIDEnd && getValue(itID) != INVALID_ID; goNext(itID))
						tupleOutput << ":" << getValue(fragStore.annotationNameStore, getValue(itID));
					goNext(itID);
                    tupleOutput << '^';
					// annotationIDs of second read
					tupleOutput << getValue(fragStore.annotationNameStore, getValue(itID));
					goNext(itID);
					for ( ; itID != itIDEnd; goNext(itID))
						tupleOutput << ':' << getValue(fragStore.annotationNameStore, getValue(itID));
                    tupleOutput << '\t';
					// tuple count
					tupleOutput << *itC << '\t';
					// normalized tuple count
                    tupleOutput << formattedNumber("%f", *itN) << '\n';
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
