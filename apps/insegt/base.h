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

#ifndef SEQAN_HEADER_BASE_H
#define SEQAN_HEADER_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// create possible tuples with length n:
//////////////////////////////////////////////////////////////////////////////
template<typename TStringSet, typename TSpec, typename TConfig, typename TID>
inline void
create_nTuple(TStringSet &tupleSet, FragmentStore<TSpec, TConfig> &fragStore, const TStringSet &annoIDs, const TID &parentID, const unsigned &n)	
{
	typedef typename Iterator<TStringSet>::Type 				TSetIter;
	typedef typename Value<TStringSet>::Type 				TString;
	typedef typename Iterator<TString const>::Type			TStringIter;
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	
	static const TID INVALID_ANNO_ID = TAnnotationStoreElement::INVALID_ID;
	
	clear(tupleSet);
	TStringIter itStr;
	TStringIter itStrEnd;
	TSetIter itSet;
	TSetIter itSetEnd;
	TStringSet tempAnnoIDs;
	resize(tempAnnoIDs, n);
	TStringSet tempTupleSet;
	
	for (unsigned i = 0; i <= length(annoIDs) - n; ++i)
	{
		for ( unsigned j = i; j < i + n; ++j)
		{
			clear(value(tempAnnoIDs, j - i));
			itStr = begin(annoIDs[j]);
			itStrEnd = end(annoIDs[j]);
			for ( ; itStr != itStrEnd; goNext(itStr))
			{
				if (getValue(itStr) != INVALID_ANNO_ID && getValue(fragStore.annotationStore, getValue(itStr)).parentID == parentID)
				{
					appendValue(value(tempAnnoIDs, j - i), getValue(itStr), Generous() );
				}
			}
		}
		
		clear(tempTupleSet);
		createCombinations(tempTupleSet, tempAnnoIDs);
		itSet = begin(tempTupleSet);
		itSetEnd = end(tempTupleSet);
		for ( ; itSet != itSetEnd; goNext(itSet))
		{
			appendValue(tupleSet, getValue(itSet), Generous());
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// create all possible tuples: for each n <= length
//////////////////////////////////////////////////////////////////////////////
template<typename TStringSet, typename TSpec, typename TConfig, typename TID>
inline void
create_Tuple(TStringSet &tupleSet, FragmentStore<TSpec, TConfig> &fragStore, const TStringSet &annoIDs, const TID &parentID, const unsigned &n)	
{
	typedef typename Iterator<TStringSet>::Type TIter;
	
	clear(tupleSet);
	TIter it;
	TIter itEnd;
	TStringSet tempTupleSet;
	for (unsigned i = 1; i <= n && i <= length(annoIDs); ++i)
	{
		clear(tempTupleSet);
		create_nTuple(tempTupleSet, fragStore, annoIDs, parentID, i);
		
		it = begin(tempTupleSet);
		itEnd = end(tempTupleSet);
		for ( ; it != itEnd; goNext(it))
		{
			appendValue(tupleSet, getValue(it), Generous());
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// create possible tuples: 
//////////////////////////////////////////////////////////////////////////////
template<typename TStringSet>
inline void
createCombinations(TStringSet &tupleSet, TStringSet &annoIDs)
{
	//typedef typename Value<TStringSet>::Type 	TString;
	//typedef typename Value<TString>::Type		TID;
	typedef typename Iterator<TStringSet>::Type 	TStringSetIter;
	//typedef typename Iterator<TString>::Type	TStringIter;
	
	if (!empty(annoIDs))
	{
		TStringSetIter itStrSet = begin(annoIDs);
		TStringSetIter itStrSetEnd = end(annoIDs);
	
		unsigned n = 1;
		for ( ; itStrSet != itStrSetEnd; goNext(itStrSet))
		{
			n = n * length(getValue(itStrSet));
		}
	
		clear(tupleSet);
		resize(tupleSet, n);
		/*
		TStringSetIter itT = begin(tupleSet);			
		TStringSetIter itTEnd = end(tupleSet);
	
		for ( ; itT != itTEnd; goNext(itT))
		{
			resize(value(itT), length(annoIDs));
		}
		*/
	
		// TID id;
		unsigned pos;
		unsigned help = n;
		for (unsigned i = 0; i < length(annoIDs); ++i)
		{
			unsigned m = length(getValue(annoIDs, i));
			help = help / m;
			for (unsigned j = 0; j< length(getValue(annoIDs, i)); ++j)
			{
				// id = getValue(getValue(annoIDs, i), j);
				for (unsigned k = 0; k < unsigned (n/ (help*m)) ; ++k)
				{
					for (unsigned l = 0; l < help; ++l)
					{
						pos = j*help + k*m*help + l;
						// assignValue(value(tupleSet, pos), i, getValue(getValue(annoIDs, i), j) );
						appendValue(value(tupleSet, pos), getValue(getValue(annoIDs, i), j), Generous() );
					}
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// getDifference
//////////////////////////////////////////////////////////////////////////////
template<typename TString>
inline void
getDifference(TString &difference, const TString &string1, const TString &string2)
{
	typedef typename Iterator<TString>::Type TIter;
	
	resize(difference, length(string1));
	TString str1 = string1;
	TString str2 = string2;
	std::stable_sort(begin(str1), end(str1));
	std::stable_sort(begin(str2), end(str2));
	
	TIter it = std::set_difference(begin(str1), end(str1), begin(str2), end(str2), begin(difference));
	erase(difference, position(it, difference), endPosition(difference));
}

//////////////////////////////////////////////////////////////////////////////
// searchValue
//////////////////////////////////////////////////////////////////////////////
template<typename TPos, typename TValue, typename TString>
inline bool
searchValue(TPos & pos, TValue const & value, TString const & string)
{
	typedef typename Iterator<TString const>::Type TIter;
	
	TIter it = std::lower_bound(begin(string), end(string), value);
	pos = position(it, string);
	
	if (it == end(string) || getValue(it) != value) return false;
	else return true;
}


//////////////////////////////////////////////////////////////////////////////
// intersection
//////////////////////////////////////////////////////////////////////////////

template<typename TList>
inline bool
interSec(TList &result, const TList &list1, const TList  &list2)
{
	//typedef typename Value<TList>::Type TValue;
	typedef typename Iterator<TList>::Type TIter;
	
	resize(result, length(list1));
	TList sort1 = list1;
	TList sort2 = list2;
	std::stable_sort(begin(sort1), end(sort1));
	std::stable_sort(begin(sort2), end(sort2));
	
	TIter itR = std::set_intersection(begin(sort1), end(sort1), begin(sort2), end(sort2), begin(result));
	
	erase(result, position(itR, result), endPosition(result));
	if (empty(result)) return false;
	else return true; 
}



//////////////////////////////////////////////////////////////////////////////
// isElement
//////////////////////////////////////////////////////////////////////////////

template<typename TValue>
inline bool
isElement(const TValue &value, const String<TValue> &string)
{
	return std::binary_search(begin(string), end(string), value);	
}


//////////////////////////////////////////////////////////////////////////////
// isElement
//////////////////////////////////////////////////////////////////////////////

template<typename TValue>
inline bool
isElement_unsorted(TValue const & value, String<TValue> const & string)
{
	typedef typename Iterator<String<TValue> const>::Type TIter;
	
	TIter it = std::find(begin(string), end(string), value);
	
	if (atEnd(it, string)) return false;
	else return true;
}


//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
