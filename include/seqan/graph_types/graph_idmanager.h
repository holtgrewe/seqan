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

// TODO(holtgrew): Move to misc?

#ifndef SEQAN_HEADER_GRAPH_IDMANAGER_H
#define SEQAN_HEADER_GRAPH_IDMANAGER_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// IDManager
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class IDManager
 * @headerfile <seqan/graph_types.h>
 * @brief ID manager that provides uniq ids for vertices and edges.
 *
 * @signature template <typename TIDType, typename TSpec>
 *            class IDManager;
 *
 * @tparam TIDType The type to use for storing IDs.  Use the Value metafunctio to get the id type managed by
 *                 a given IDManager.
 * @tparam TSpec   Specializing type.
 */

template <typename TIDType, typename TSpec>
class IDManager 
{
	public:
		String<TIDType> data_freeIDs;  
		String<bool> data_in_use;   //1 = in use, 0 = not in use

//____________________________________________________________________________	
	public:
		IDManager()
		{
			SEQAN_CHECKPOINT
			clear(data_in_use);
			clear(data_freeIDs);
		}

		~IDManager() 
		{
			SEQAN_CHECKPOINT
		}

		IDManager(IDManager const & _other)
		{
			SEQAN_CHECKPOINT
			data_freeIDs = _other.data_freeIDs;
			data_in_use = _other.data_in_use;
		}

		IDManager const& 
		operator = (IDManager const& _other) 
		{
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			data_freeIDs = _other.data_freeIDs;
			data_in_use = _other.data_in_use;
			return *this;
		}

//____________________________________________________________________________
};
	

//////////////////////////////////////////////////////////////////////////////
// IDManager - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn IDManager#Value
 * @brief Return the value type of an IDManager.
 *
 * @signature Value<TIDManager>::Type;
 */

template<typename TIDType, typename TSpec> 
struct Value<IDManager<TIDType, TSpec> > 
{
	typedef TIDType Type;
};

template<typename TIDType, typename TSpec> 
struct Value<IDManager<TIDType, TSpec> const> 
{
	typedef TIDType Type;
};

//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn IDManager#Spec
 * @brief Return the specializing type of an IDManager.
 *
 * @signature Spec<TIDManager>::Type;
 */

template<typename TIDType, typename TSpec> 
struct Spec<IDManager<TIDType, TSpec> > 
{
	typedef TSpec Type;
};

template<typename TIDType, typename TSpec> 
struct Spec<IDManager<TIDType, TSpec> const> 
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn IDManager#obatinID
 * @brief Obatains a new id from the id manager.
 *
 * @signature TID obtainID(manager);
 *
 * @param[in,out] manager IDManager to get ids from.
 *
 * @return TID A new unique id.
 */

template<typename TIDType, typename TSpec>
inline typename Value<IDManager<TIDType, TSpec> >::Type 
obtainID(IDManager<TIDType, TSpec>& idm) 
{
	SEQAN_CHECKPOINT

	TIDType _id;
	if (!empty(idm.data_freeIDs)) {
		_id = getValue(idm.data_freeIDs, length(idm.data_freeIDs) - 1);
		resize(idm.data_freeIDs, length(idm.data_freeIDs) - 1, Generous());
		assignValue(idm.data_in_use, _id, true);
	} else {
		if (empty(idm.data_in_use)) _id = 0;
		else _id = (TIDType) length(idm.data_in_use);
		resize(idm.data_in_use, _id + 1, Generous());
		assignValue(idm.data_in_use, _id, true);
	}
	return _id;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn IDManager#releaseID
 * @brief Releases a given id so it can be reused later on.
 *
 * @signature void releaseID(manager, id);
 *
 * @param[in,out] manager The IDManager to use.
 * @param[in]     id      The ID to release.
 */

template<typename TIDType, typename TSpec, typename TID>
inline void 
releaseID(IDManager<TIDType, TSpec>& idm, 
		  TID const id) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(idm,id));
	if (id == (TID) length(idm.data_in_use) - 1) {
		resize(idm.data_in_use, length(idm.data_in_use) - 1, Generous());
	} else {
		assignValue(idm.data_in_use, id, false);
		appendValue(idm.data_freeIDs, id, Generous());
	}
	if (idCount(idm)==0) {
		releaseAll(idm);
	}
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn IDManager#releaseAll
 * @brief Releases all ids handled by this IDManager.
 *
 * @signature void releaseAll(manager);
 *
 * @param[in,out] manager The IDManager to release the ids from.
 */


template<typename TIDType, typename TSpec>
inline void 
releaseAll(IDManager<TIDType, TSpec>& idm) 
{
	SEQAN_CHECKPOINT
	clear(idm.data_freeIDs);
	clear(idm.data_in_use);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn IDManager#getIDUpperBound
 * @brief Returns the largest distributed id plus 1.
 *
 * @signature TID getIDUpperBound(manager);
 *
 * @param[in] manager The IDManager to query.
 *
 * @return TID The return value is guaranteed to be an upper bound on all distributed ids.
 */

template<typename TIDType, typename TSpec>
inline typename Value<IDManager<TIDType, TSpec> >::Type 
getIDUpperBound(IDManager<TIDType, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	return (empty(idm.data_in_use)) ? 0 : (typename Value<IDManager<TIDType, TSpec> >::Type) length(idm.data_in_use);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn IDManager#getIDLowerBound
 * @brief Returns the smalles distributed id.
 *
 * @signature TID getIDLowerBound(manager);
 *
 * @param[in] manager The IDManager to query.
 *
 * @return TID The smallest obtained id so far.
 */

template<typename TIDType, typename TSpec>
inline typename Value<IDManager<TIDType, TSpec> >::Type 
getIDLowerBound(IDManager<TIDType, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	for(TIDType it = 0; it < length(idm.data_in_use); ++it) {
		if (getValue(idm.data_in_use, it)) return it;
	}
	return 0;
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn IDManager#idCount
 * @brief Return the number of obtained ids.
 *
 * @signature TID idCount(manager);
 *
 * @param[in] manager The IDManager to query.
 *
 * @return TID The number of ids obtained so far.
 */

template<typename TIDType, typename TSpec>
inline typename Value<IDManager<TIDType, TSpec> >::Type 
idCount(IDManager<TIDType, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	return (length(idm.data_in_use) - length(idm.data_freeIDs));
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn IDManager#idInUse
 * @brief Return whether an id is in use.
 *
 * @signature bool idInUse(manager, id);
 *
 * @param[in] manager The IDManager to query.
 * @param[in] id      The id to check for.
 *
 * @return bool true if the id is already in use and false otherwise.
 */

template<typename TIDType, typename TSpec, typename TID>
inline bool 
idInUse(IDManager<TIDType, TSpec> const& idm, 
		TID const id)
{
	SEQAN_CHECKPOINT
	return (id < static_cast<TID>(length(idm.data_in_use))) ? idm.data_in_use[id] : false;
}


//////////////////////////////////////////////////////////////////////////////
// Dummy IDManager
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class CountingIDManager
 * @extends IDManager
 * @headerfile <seqan/graph_types.h>
 * @brief IDManager that just counts the number of ids in use.
 *
 * @signature template <typename TSpec>
 *            class IDManager<void, TSpec>;
 *
 * @tparam TSpec Specializing type.
 */

template<typename TSpec>
class IDManager<void, TSpec> 
{
	public:
		typedef typename ID<IDManager>::Type TIDType;
		TIDType data_idCount;

//____________________________________________________________________________	
	public:
		IDManager() : data_idCount(0) 
		{
			SEQAN_CHECKPOINT
		}

		~IDManager() 
		{
			SEQAN_CHECKPOINT
		}

		IDManager(IDManager const & _other) : data_idCount(_other.data_idCount) 
		{
			SEQAN_CHECKPOINT
		}

		IDManager const& 
		operator = (IDManager const& _other) 
		{
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			data_idCount = _other.data_idCount;
			return *this;
		}

//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Dummy IDManager - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec> 
struct Value<IDManager<void, TSpec> > {
	typedef typename Size<IDManager<void, TSpec> >::Type Type;
};

template<typename TSpec> 
struct Value<IDManager<void, TSpec> const> {
	typedef typename Size<IDManager<void, TSpec> const>::Type Type;
};


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Value<IDManager<void, TSpec> >::Type 
obtainID(IDManager<void, TSpec>& idm) 
{
	SEQAN_CHECKPOINT
	++idm.data_idCount;
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TID>
inline void 
releaseID(IDManager<void, TSpec>& idm, 
		  TID const) 
{
	SEQAN_CHECKPOINT
	if (idm.data_idCount > 0) --idm.data_idCount;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline void 
releaseAll(IDManager<void, TSpec>& idm) 
{
	SEQAN_CHECKPOINT
	idm.data_idCount = 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Value<IDManager<void, TSpec> >::Type 
getIDUpperBound(IDManager<void, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	// Must be data_idCount in order to resize property maps!!!
	// Don't change to 0
	return idm.data_idCount;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Value<IDManager<void, TSpec> >::Type 
getIDLowerBound(IDManager<void, TSpec> const&)
{
	SEQAN_CHECKPOINT
	return 0;
}

//////////////////////////////////////////////////////////////////////////////


template <typename TSpec>
inline typename Value<IDManager<void, TSpec> >::Type 
idCount(IDManager<void, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	return idm.data_idCount;
}

//////////////////////////////////////////////////////////////////////////////


template <typename TSpec, typename TID>
inline bool 
idInUse(IDManager<void, TSpec> const&, 
		TID const) 
{
	SEQAN_CHECKPOINT
	return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
