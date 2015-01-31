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

#ifndef SEQAN_HEADER_STORE_ANNOTATION_H
#define SEQAN_HEADER_STORE_ANNOTATION_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Annotation Store
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class AnnotationStoreElement
 * @headerfile <seqan/store.h>
 * @brief Represents an annotation of a genomic feature.
 *
 * @signature template <typename TPos[, typename TSpec]>
 *            struct AnnotationStoreElement;
 *
 * @tparam TPos  The position type to use.
 * @tparam TSpec The specializing type.
 *
 * Value type of the @link FragmentStore::annotationStore @endlink stirng.  The id of an annotation (aka annotationID)
 * is not stored explicitly as it is implicitly given byt he position in the @link FragmentStore::annotationStore
 * @endlink.
 *
 *
 * @fn AnnotationStoreElement::AnnotationStoreElement()
 * @brief Constructor.
 *
 * @signature AnnotationStoreElement::AnnotationStoreElement();
 *
 * The default constructor sets all members to <tt>INVALID_ID</tt> and beginPos and endPos to <tt>INVALID_POS</tt>.
 */

/*!
 * @typedef AnnotationStoreElement::TID
 * @brief Type of annotationID and contigID.
 *
 * @signature typedef (...) TAnnotationStoreElement::TID;
 *
 * TID is the result of <tt>ID&lt;AnnotationStoreElement&lt;&gt; &gt;::Type</tt>, see @link ID @endlink.
 *
 *
 * @typedef AnnotationStoreElement::TPos
 * @brief Type of the @link AnnotationStoreElement::beginPos @endlink and @link AnnotationStoreElement::endPos
 *        @endlink.
 *
 * @signature typedef (...) TAnnotationStoreElement::TPos;
 *
 * @typedef AnnotationStoreElement::TValues
 * @brief @link StringSet @endlink type of the @link AnnotationStoreElement::values @endlink member.
 *
 * @signature typedef (...) TAnnotationStoreElement::TValues;
 */

/*!
 * @var TID AnnotationStoreElement::INVALID_ID;
 * @brief Constant to represent an invalid id.
 *
 * @var TID AnnotationStoreElement::INVALID_POS;
 * @brief Constant to represent an invalid position.
 *
 * @var TID AnnotationStoreElement::contigID;
 * @brief Refers to the contig in @link FragmentStore::contigStore @endlink that the annotation is part of.
 *
 * @var TID AnnotationStoreElement::typeID;
 * @brief Refers to an entry in the @link FragmentStore::annotationTypeStore @endlink.
 *
 * There are some type ids predefined for commonly used types, e.g. $ANNO_GENE$.  See @link
 * FragmentStore::PredefinedAnnotationTypes @endlink.
 *
 * @var TPos AnnotationStoreElement::beginPos;
 * @brief Begin position of the annotation in gap space.
 *
 * If endPos &lt; beginPos then the annotation feature is located on the reverse strand where beginPos and endPos are
 * the corresponding positions on the forward strand.
 *
 * @var TPos AnnotationStoreElement::endPos;
 * @brief End position of the annotation in gap space.
 *
 * If endPos &lt; beginPos then the annotation feature is located on the reverse strand where beginPos and endPos are
 * the corresponding positions on the forward strand.
 *
 * @var TValues AnnotationStoreElement::values;
 * @brief @link StringSet @endlink that stores additional values addressed by keyID.
 *
 * The GFF and FTF file formats allow to define user-specific key-value pairs.  The set of all keys addressed by keyID
 * are stored in @link FragmentStore::annotationKeyStore @endlink.
 *
 * @var TID AnnotationStoreElement::parentID;
 * @brief The id of the parent annotation.
 *
 * @var TID AnnotationStoreElement::nextSiblingID;
 * @brief The id of the right sibling annotation.
 *
 * @var TID AnnotationStoreElement::lastChildID;
 * @brief The id of the parent annotation.
 */

template <typename TPos_, typename TSpec = void>
struct AnnotationStoreElement
{
	typedef typename ID<AnnotationStoreElement>::Type		TID;
	typedef TPos_											TPos;
	typedef StringSet<CharString, Owner< ConcatDirect<> > >	TValues;

	static const TID  INVALID_ID;
	static const TPos INVALID_POS;
	
	TID					parentID;
	TID					contigID;
	TID					countID;
	TID					typeID;			// gene, intron, ...

	TPos				beginPos;		// begin position of the annotation in the gapped contig sequence (i.e. in gap-space)
	TPos				endPos;			// end position of ..., for annotations on the reverse strand holds end < begin
	
	TID					lastChildID;	// generated back links to child
	TID					nextSiblingID;	// and sibling
	
	TValues				values;			// stores values for each keyID of (key,value) pairs

	AnnotationStoreElement() : 
		parentID(INVALID_ID), contigID(INVALID_ID), countID(INVALID_ID), typeID(INVALID_ID), 
		beginPos(INVALID_POS), endPos(INVALID_POS),
		lastChildID(INVALID_ID), nextSiblingID(INVALID_ID) {}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TSpec>
const typename ID<AnnotationStoreElement<TPos, TSpec> >::Type
AnnotationStoreElement<TPos, TSpec>::INVALID_ID = MaxValue<typename ID<AnnotationStoreElement<TPos, TSpec> >::Type>::VALUE;

template <typename TPos, typename TSpec>
const TPos
AnnotationStoreElement<TPos, TSpec>::INVALID_POS = MaxValue<TPos>::VALUE;

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec = void>
struct AnnotationTree {};

/*!
 * @class AnnotationTreeIterator
 * @extends Iter
 * @headerfile <seqan/store.h>
 * @brief Iterator of the annotation tree represetned by a @link FragmentStore @endlink.
 *
 * This iterator can move @link AnnotationTreeIterator#goDown down @endlink, @link AnnotationTreeIterator#goRight right
 * @endlink, and @link AnnotationTreeIterator#goUp up @endlink in the tree and supports a preorder DFS traversal using
 * the functions @link RootedRandomAccessIteratorConcept#goBegin @endlink, @link InputIteratorConcept#goNext @endlink,
 * and @link RootedIteratorConcept#atEnd @endlink.
 *
 * Preorder means that the iterator visits the node befor it visits its children.
 *
 * To obtain the type of the AnnotationTreeIterator for a FragmentStore specializiation <tt>TFragmentStore</tt>, you can
 * use the metafunction @link ContainerConcept#Iterator Iterator @endlink as follows:
 * <tt>Iterator&lt;TFragmentStore&gt;::Type</tt>.
 *
 * @signature template <typename TFragmentStore>
 *            class Iter<TFragmentStore, AnnotationTree<> >;
 *
 * @tparam TFragmentStore The FragmentStore to iterate over.
 *
 * @section Remarks
 *
 * To access the annotation, the iterator points to, use @link AnnotationTreeIterator#getAnnotation @endlink.  The
 * annotation id is returned by @link IteratorAssociatedTypesConcept#value @endlink.
 *
 * @section Example
 *
 * <img src="AnnotationTree.png" title="Typical annotation tree hierarchy." />
 *
 * A new annotation iterator can be instantiated as follows:
 *
 * @code{.cpp}
 * FragmentStore<> store;
 * Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
 * it = begin(store, AnnotationTree<>());
 * @endcode
 *
 * Or shorter:
 *
 * @code{.cpp}
 * FragmentStore<> store;
 * Iterator<FragmentStore<>, AnnotationTree<> >::Type it(store);
 * @endcode
 *
 *
 * @fn AnnotationTreeIterator::AnnotationTreeIterator
 * @brief Constructor
 *
 * @signature Iter::Iter();
 * @signature Iter::Iter(store[, startInNode]);
 *
 * @param[in] store       The FragmentStore with the annotation tree to iterate.
 * @param[in] startInNode Annotation id of the ndoe the iterator should start at.  Default: 0, the root node id.
 *
 * The @link ContainerConcept#begin @endlink function can also be used to create a tree iterator that starts in the root
 * node:
 *
 * @code{.cpp}
 * FragmentStore<> store;
 * Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
 * it = begin(store, AnnotationTree<>());
 * @endcode
 */

template <typename TFragmentStore, typename TSpec>
class Iter<TFragmentStore, AnnotationTree<TSpec> >
{
public:
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TID					TID;

	TFragmentStore *store;
	TID _id;
	
	Iter():
		store(NULL),
		_id(TAnnotation::INVALID_ID) {}

	Iter(TFragmentStore &_store):
		store(&_store),
		_id(0) {}

	Iter(TFragmentStore &_store, TID startInNode):
		store(&_store),
		_id(startInNode) {}

	Iter(TFragmentStore &_store, MinimalCtor):
		store(&_store),
		_id(TAnnotation::INVALID_ID) {}

	inline Iter const &
	operator = (Iter const &_origin)
	{
		store = &container(_origin);
		_id = _origin._id;
		return *this;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
struct Iterator< TFragmentStore, AnnotationTree<TSpec> > {
	typedef Iter< TFragmentStore, AnnotationTree<TSpec> > Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
struct Value< Iter< TFragmentStore, AnnotationTree<TSpec> > >:
	VertexDescriptor<TFragmentStore> {};

template <typename TFragmentStore, typename TSpec>
struct Size< Iter< TFragmentStore, AnnotationTree<TSpec> > > :
	Size<TFragmentStore> {};

template <typename TFragmentStore, typename TSpec>
struct Position< Iter< TFragmentStore, AnnotationTree<TSpec> > > :
	Position<TFragmentStore> {};

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline typename VertexDescriptor<TFragmentStore>::Type & 
value(Iter< TFragmentStore, AnnotationTree<TSpec> > &it) { 
	return it._id;
}

template <typename TFragmentStore, typename TSpec>
inline typename VertexDescriptor<TFragmentStore>::Type const & 
value(Iter< TFragmentStore, AnnotationTree<TSpec> > const &it) { 
	return it._id;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline TFragmentStore & 
container(Iter< TFragmentStore, AnnotationTree<TSpec> > &it) { 
	return *it.store;
}

template <typename TFragmentStore, typename TSpec>
inline TFragmentStore & 
container(Iter< TFragmentStore, AnnotationTree<TSpec> > const &it) { 
	return *it.store;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AnnotationTreeIterator#getAnnotation
 * @brief Returns the current annotation.
 *
 * @signature TAnnotation getAnnotation(it);
 *
 * @param[in] it The AnnotationTreeIterator to query for its annotation.
 *
 * @return TAnnotation A reference to the AnnotationStoreElement the iterator points to.
 */

// TODO(holtgrew): Better use Reference in both cases?

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationStore>::Type
getAnnotation(Iter<TFragmentStore const, AnnotationTree<TSpec> > const & it)
{
	return getValue(it.store->annotationStore, it._id);
}

template <typename TFragmentStore, typename TSpec>
inline typename Reference<typename TFragmentStore::TAnnotationStore>::Type
getAnnotation(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return value(it.store->annotationStore, it._id);
}

/*!
 * @fn AnnotationTreeIterator#getName
 * @brief Returns the identifier of the current annotation.
 *
 * @signature TName getName(it);
 *
 * @param[in] it An AnnotationTreeIterator to query.
 *
 * @return TName The name of the current annotation.  This is a reference to the corresponding position in @link
 *               FragmentStore::annotationNameStore @endlink.
 *
 * @see AnnotationTreeIterator#setName
 */

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationNameStore>::Type
getName(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return getAnnoName(*it.store, it._id);
}

/*!
 * @fn AnnotationTreeIterator#setName
 * @brief Sets the identifier of the current annotation.
 *
 * @signature void setName(it, name);
 *
 * @param[in,out] it   Iterator to the annotation to set the name for.
 * @param[in]     name The new identifier of the annotation pointed to by it.
 *
 * @see AnnotationTreeIterator#getName
 */

template <typename TFragmentStore, typename TSpec, typename TName>
inline void
setName(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it, TName const & name)
{
	if (length(it.store->annotationNameStore) <= it._id)
		resize(it.store->annotationNameStore, it._id + 1);
	it.store->annotationNameStore[it._id] = name;
}

/*!
 * @fn AnnotationTreeIterator#getParentName
 * @brief Returns the identifier of the parent node in the annotationt ree of the current annotation.
 *
 * @signature TName getParentName(it);
 *
 * @param[in] it The AnnotationTreeIterator to get the name for.
 *
 * @return TName The name of the parent of the current annotation.  This is a reference to the corresponding value in
 *               @link FragmentStore::annotationNameStore @endlink.
 */

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationNameStore>::Type
getParentName(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TID					TID;

	TID parentID = it.store->annotationStore[it._id].parentID;
	if (parentID == TAnnotation::INVALID_ID) parentID = it._id;
	return getAnnoName(*it.store, parentID);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AnnotationTreeIterator#getType
 * @brief Returns the type name of the current annotation.
 *
 * @signature TSeq getType(it);
 *
 * @param[in] it The AnnotationTreeIterator to query for its type.
 *
 * @return The name of the current annotation, e.g. "exon" or "mRNA".  This is a reference to an entry in @link
 *         FragmentStore::annotationTypeStore @endlink.
 *
 * @see AnnotationTreeIterator#setType
 */

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationTypeStore>::Type
getType(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return getAnnoType(*it.store, getAnnotation(it).typeID);
}

/*!
 * @fn AnnotationTreeIterator#setType
 * @brief Sets the type name of the current annotation.
 *
 * @signature void setType(it, typeName);
 *
 * @param[in,out] it       The iterator to the annotation to set the type name for.
 * @param[in]     typeName The name of the type (e.g. "exon" or "mRNA").  Type: @link ContainerConcept @endlink.
 *
 * @see AnnotationTreeIterator#getType
 */

template <typename TFragmentStore, typename TSpec, typename TTypeName>
inline void
setType(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it, TTypeName & typeName)
{
	_storeAppendType(*it.store, getAnnotation(it).typeID, typeName);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AnnotationTreeIterator#getUniqueName
 * @brief Returns a unique name of the current annotation.
 *
 * @signature CharString getUniqueName(it);
 *
 * @param[in] it The AnnotationTreeIterator to query.
 *
 * @return CharString A unique name of the current annotation.
 *
 * Some annotation file formats do not require that every annotation has a non-empty name.  This function returns the
 * name if non-emtpy and otherwise generates one using the type an id.
 *
 * @see AnnotationTreeIterator#getName
 */

template <typename TFragmentStore, typename TSpec>
inline CharString
getUniqueName(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return getAnnoUniqueName(*it.store, it._id);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AnnotationTreeIterator#clearValues
 * @brief Clear all key-value pairs of a given annotation.
 *
 * @signature void clearValues(it);
 *
 * @param[in,out] it The AnnotationTreeIterator to clear all key-value pairs of the current annotation.
 *
 * @see AnnotationTreeIterator#assignValueByKey
 * @see AnnotationTreeIterator#getValueByKey
 */

template <typename TFragmentStore, typename TSpec>
inline void 
clearValues(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	clear(getAnnotation(it).values);
}

/*!
 * @fn AnnotationTreeIterator#assignValueByKey
 * @brief Add or update a key-value pair of the current annotation.
 *
 * @signature void assignValueByKey(it, key, value);
 *
 * @param[in,out] it    The AnnotationTreeIterator to
 * @param[in]     key   The key whose value should be changed.  Type: @link ContainerConcept @endlink.
 * @param[in]     value The value that should be assigned.  Type: @link ContainerConcept @endlink.
 *
 * @see AnnotationTreeIterator#clearValues
 * @see AnnotationTreeIterator#getValueByKey
 */

template <typename TFragmentStore, typename TSpec, typename TKey, typename TValue>
inline void 
assignValueByKey(
	Iter<TFragmentStore, AnnotationTree<TSpec> > & it,
	TKey const & key,
	TValue const & value)
{
	annotationAssignValueByKey(*it.store, getAnnotation(it), key, value);
}

/*!
 * @fn AnnotationTreeIterator#getValueByKey
 * @brief Retrieve a key's value in the current annotation.
 *
 * @signature bool getValueByKey(it, key);
 * @signature bool getValueByKey(value, it, key);
 *
 * @param[in]  it    The AnnotationTreeIterator for which to retrieve the key.
 * @param[in]  key   The key to get the value for.
 * @param[out] value The resulting value.
 *
 * @see AnnotationTreeIterator#clearValues
 * @see AnnotationTreeIterator#assignValueByKey
 */

template <typename TFragmentStore, typename TSpec, typename TKey, typename TValue>
inline bool 
getValueByKey(
    TValue & value,
    Iter<TFragmentStore, AnnotationTree<TSpec> > const & it,
	TKey const & key)
{
	return annotationGetValueByKey(*it.store, getAnnotation(it), key, value);
}

template <typename TFragmentStore, typename TSpec, typename TKey>
inline CharString
getValueByKey(
	Iter<TFragmentStore, AnnotationTree<TSpec> > const & it,
	TKey const & key)
{
	return annotationGetValueByKey(*it.store, getAnnotation(it), key);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline void
goBegin(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	it._id = 0;
}

template <typename TFragmentStore, typename TSpec>
inline void
goEnd(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

	it._id = TAnnotation::INVALID_ID;
}

/*!
 * @fn AnnotationTreeIterator#goTo
 * @brief Moves the iterator to an arbitrary node giving its annotationID.
 *
 * @signature void goTo(it, annoID);
 *
 * @param[in,out] it     The AnnotationTreeIterator to move.
 * @param[in]     annoID The id of the annotation to move to.
 */

template <typename TFragmentStore, typename TSpec, typename TID>
inline void
goTo(Iter<TFragmentStore, AnnotationTree<TSpec> > & it, TID _id)
{
	it._id = _id;
}

template <typename TFragmentStore, typename TSpec>
inline void
clear(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

	it._id = TAnnotation::INVALID_ID;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline bool
atBegin(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	return it._id == 0;
}

template <typename TFragmentStore, typename TSpec>
inline bool
atEnd(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

	return it._id == TAnnotation::INVALID_ID;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline void
goNext(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	// preorder dfs
	if (!goDown(it) && !goRight(it))
		while (goUp(it) && !goRight(it)) ;
	if (isRoot(it)) {
		clear(it);
		return;
	}
}

/*!
 * @fn AnnotationTreeIterator#goNextRight
 * @brief Go to the next node in preorder DFS skipping the current node's subtree.
 *
 * @signature void goNextUp(it);
 *
 * @param[in,out] it The AnnotationTreeIterator to
 */

template <typename TFragmentStore, typename TSpec>
inline void
goNextRight(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	// preorder dfs
	if (!goRight(it))
		while (goUp(it) && !goRight(it)) ;
	if (isRoot(it)) {
		clear(it);
		return;
	}
}

/*!
 * @fn AnnotationTreeIterator#goNextUp
 * @brief Go to the next node in preorder DFS skipping the subtrees of the current node and all of its siblings.
 *
 * @signature void goNextUp(it);
 *
 * @param[in,out] it The AnnotationTreeIterator to
 */

template <typename TFragmentStore, typename TSpec>
inline void
goNextUp(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	// preorder DFS
	while (goUp(it) && !goRight(it)) ;
	if (isRoot(it)) {
		clear(it);
		return;
	}
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AnnotationTreeIterator#goRoot
 * @brief Move the iterator down to the tree root.
 *
 * @signature void goRoot(it);
 *
 * @param[in,out] it The AnnotationTreeIterator to move.
 */

template <typename TFragmentStore, typename TSpec>
inline void
goRoot(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	it._id = 0;
}

/*!
 * @fn AnnotationTreeIterator#goUp
 * @brief Move the iterator down to the parent in the annotation tree.
 *
 * @signature bool goUp(it);
 *
 * @param[in,out] it The AnnotationTreeIterator to move.
 *
 * @return bool <tt>true</tt> if the iterator could be moved and <tt>false</tt> otherwise.
 *
 * @see AnnotationTreeIterator#goDown
 * @see AnnotationTreeIterator#goRight
 */

template <typename TFragmentStore, typename TSpec>
inline bool
goUp(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TID					TID;
	
	TID parentID = getAnnotation(it).parentID;
	if (parentID != TAnnotation::INVALID_ID)
	{
		it._id = parentID;
		return true;
	}
	return false;
}

/*!
 * @fn AnnotationTreeIterator#goDown
 * @brief Move the iterator down to the leftmost child in the annotation tree.
 *
 * @signature bool goDown(it);
 *
 * @param[in,out] it The AnnotationTreeIterator to move.
 *
 * @return bool <tt>true</tt> if the iterator could be moved and <tt>false</tt> otherwise.
 *
 * @see AnnotationTreeIterator#goUp
 * @see AnnotationTreeIterator#goRight
 */

template <typename TFragmentStore, typename TSpec>
inline bool
goDown(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TID					TID;
	
	TID lastChildID = getAnnotation(it).lastChildID;
	if (lastChildID != TAnnotation::INVALID_ID)
	{
		it._id = it.store->annotationStore[lastChildID].nextSiblingID;
		return true;
	}
	return false;
}

/*!
 * @fn AnnotationTreeIterator#goRight
 * @brief Move the iterator right to the next sibling in the annotation tree.
 *
 * @signature bool goRight(it);
 *
 * @param[in,out] it The AnnotationTreeIterator to move.
 *
 * @return bool <tt>true</tt> if the iterator could be moved and <tt>false</tt> otherwise.
 *
 * @see AnnotationTreeIterator#goUp
 * @see AnnotationTreeIterator#goDown
 */

template <typename TFragmentStore, typename TSpec>
inline bool
goRight(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TID					TID;
	
	TAnnotation const &anno = getAnnotation(it);
	TID nextSiblingID = anno.nextSiblingID;
	if (nextSiblingID != TAnnotation::INVALID_ID)
	{
		TID lastChildID = it.store->annotationStore[anno.parentID].lastChildID;
		if (it._id != lastChildID)
		{
			it._id = nextSiblingID;
			return true;
		}
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AnnotationTreeIterator#nodeUp
 * @brief Returns a new iterator to the parent node of the current annotation in the annotation tree.
 *
 * @signature TIter nodeUp(it);
 *
 * @param[in] it The AnnotationTreeIterator to query.
 *
 * @return TIter Iterator to the parent node.
 *
 * @see AnnotationTreeIterator#nodeDown
 * @see AnnotationTreeIterator#nodeRight
 */

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
nodeUp(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	Iter<TFragmentStore, AnnotationTree<TSpec> > tmp(it);
	goUp(tmp);
	return tmp;
}

/*!
 * @fn AnnotationTreeIterator#nodeDown
 * @brief Returns a new iterator to the first child node of the current annotation in the annotation tree.
 *
 * @signature TIter nodeDown(it);
 *
 * @param[in] it The AnnotationTreeIterator to query.
 *
 * @return TIter Iterator to the first child node of it.
 *
 * @see AnnotationTreeIterator#nodeUp
 * @see AnnotationTreeIterator#nodeRight
 */

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
nodeDown(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	Iter<TFragmentStore, AnnotationTree<TSpec> > tmp(it);
	goDown(tmp);
	return tmp;
}

/*!
 * @fn AnnotationTreeIterator#nodeRight
 * @brief Returns a new iterator to the right sibling of the current annotation in the annotation tree.
 *
 * @signature TIter nodeRight(it);
 *
 * @param[in] it The AnnotationTreeIterator to query.
 *
 * @return TIter Iterator to the right sibling of it.
 *
 * @see AnnotationTreeIterator#nodeUp
 * @see AnnotationTreeIterator#nodeDown
 */

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
nodeRight(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	Iter<TFragmentStore, AnnotationTree<TSpec> > tmp(it);
	goRight(tmp);
	return tmp;
}

//////////////////////////////////////////////////////////////////////////////

// insert a new _id into a cyclic list and returns new last child _id
template <typename TAnnotationStore, typename TID>
inline TID
_cyclicListFrontInsert(TAnnotationStore & annotationStore, TID newID, TID lastChildID)
{
	typedef typename Value<TAnnotationStore>::Type TAnnotation;

	TID nextID, newLastID;
	if (lastChildID != TAnnotation::INVALID_ID)
	{
		// get last node in the cycle
		TAnnotation &lastChild = annotationStore[lastChildID];
		// last child points to first child
		nextID = lastChild.nextSiblingID;
		// insert new node between last and first
		lastChild.nextSiblingID = newID;
		// last child remains the same
		newLastID = lastChildID;
	} else
		// cyclic list was empty
		newLastID = nextID = newID;
	
	// link new node to former first node
	annotationStore[newID].nextSiblingID = nextID;
	
	return newLastID;
}

// delete an _id from a cyclic list and returns new last child _id
template <typename TAnnotationStore, typename TID>
inline TID
_cyclicListSearchPrev(TAnnotationStore & annotationStore, TID _id, TID lastChildID)
{
	typedef typename Value<TAnnotationStore>::Type TAnnotation;

	if (lastChildID == TAnnotation::INVALID_ID)
		return TAnnotation::INVALID_ID;
	
	TID prevID, i = lastChildID;
	do {
		prevID = i;
		i = annotationStore[i].nextSiblingID;
		if (i == _id) break;
	} while (i != lastChildID);

	if (i == _id)
		return prevID;
	else
		return TAnnotation::INVALID_ID;
}

// delete an _id from a cyclic list and returns new last child _id
template <typename TAnnotationStore, typename TID>
inline TID
_cyclicListRemove(TAnnotationStore & annotationStore, TID _id, TID lastChildID)
{
	typedef typename Value<TAnnotationStore>::Type TAnnotation;

	TID prevID = _cyclicListSearchPrev(annotationStore, _id, lastChildID);
	
	if (prevID != TAnnotation::INVALID_ID)
	{
		annotationStore[prevID].nextSiblingID = annotationStore[_id].nextSiblingID;
		
		if (_id == lastChildID)
		{
			if (prevID != _id)
				return prevID;
			else
				return TAnnotation::INVALID_ID;
		} else
			return lastChildID;
	}
	return lastChildID;
}

/*!
 * @fn AnnotationTreeIterator#createLeftChild
 * @brief Create a new leftmost child of the current node and returns an iterator to it.
 *
 * @signature TIter createLeftChild(it);
 *
 * @param[in,out] it The AnnotationTreeIterator to create a left child for.
 *
 * @return TIter Iterator to the new child.
 *
 * @see AnnotationTreeIterator#createRightChild
 * @see AnnotationTreeIterator#createSibling
 */

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
createLeftChild(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TID					TID;
	
	appendValue(it.store->annotationStore, getAnnotation(it));
	TAnnotation &parentAnno = getAnnotation(it);

	TID childID = length(it.store->annotationStore) - 1;
	TAnnotation &childAnno = it.store->annotationStore[childID];
	
	parentAnno.lastChildID = _cyclicListFrontInsert(it.store->annotationStore, childID, parentAnno.lastChildID);
	childAnno.parentID = it._id;
	childAnno.lastChildID = TAnnotation::INVALID_ID;
	
	Iter<TFragmentStore, AnnotationTree<TSpec> > childIter(it);
	childIter._id = childID;
	return childIter;
}

/*!
 * @fn AnnotationTreeIterator#createRightChild
 * @brief Creates a new rightmost child of the current node and returns an iterator to it.
 *
 * @signature TIter createRightChild(it);
 *
 * @param[in,out] it The AnnotationTreeIterator that was just created.
 *
 * @return TIter Iterator to the new child.
 *
 * @see AnnotationTreeIterator#createLeftChild
 * @see AnnotationTreeIterator#createSibling
 */

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
createRightChild(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TID					TID;
	
	appendValue(it.store->annotationStore, getAnnotation(it));
	TAnnotation &parentAnno = getAnnotation(it);

	TID childID = length(it.store->annotationStore) - 1;
	TAnnotation &childAnno = it.store->annotationStore[childID];
	
	_cyclicListFrontInsert(it.store->annotationStore, childID, parentAnno.lastChildID);
	parentAnno.lastChildID = childID;
	childAnno.parentID = it._id;
	childAnno.lastChildID = TAnnotation::INVALID_ID;
	
	Iter<TFragmentStore, AnnotationTree<TSpec> > childIter(it);
	childIter._id = childID;
	return childIter;
}

/*!
 * @fn AnnotationTreeIterator#createSibling
 * @brief Creates a new right sibling of the current node and return an iterator to it.
 *
 * @signature TIter createSibling(it);
 *
 * @param[in] it The AnnotationTreeIterator to create the sibling for.
 *
 * @return TIter Iterator to the new sibling.
 *
 * @see AnnotationTreeIterator#createLeftChild
 * @see AnnotationTreeIterator#createRightChild
 */

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
createSibling(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TID					TID;

	appendValue(it.store->annotationStore, getAnnotation(it));
	TAnnotation &anno = getAnnotation(it);

	TID siblingID = length(it.store->annotationStore) - 1;

	TAnnotation &parentAnno = it.store->annotationStore[anno.parentID];
	if (parentAnno.lastChildID == it._id)
		parentAnno.lastChildID = siblingID;

	TAnnotation &siblingAnno = it.store->annotationStore[siblingID];
	siblingAnno.nextSiblingID = anno.nextSiblingID;
	siblingAnno.parentID = anno.parentID;
	siblingAnno.lastChildID = TAnnotation::INVALID_ID;
	anno.nextSiblingID = siblingID;
	
	Iter<TFragmentStore, AnnotationTree<TSpec> > siblingIter(it);
	siblingIter._id = siblingID;
	return siblingIter;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn AnnotationTreeIterator#isRoot
 * @brief Return a boolean indicating whether the annotation is the root.
 *
 * @signature bool isRoot(it);
 *
 * @param[in] it The AnnotationTreeIterator to query.
 */

template <typename TFragmentStore, typename TSpec>
inline bool
isRoot(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

//	if (it._id >= length(it.store->annotationStore)) return false;
	return it.store->annotationStore[it._id].parentID == TAnnotation::INVALID_ID;
}

/*!
 * @fn AnnotationTreeIterator#isLeaf
 * @brief Return a boolean indicating whether the annotation is a leaf.
 *
 * @signature bool isLeaf(it);
 *
 * @param[in] it The AnnotationTreeIterator to query.
 */

template <typename TFragmentStore, typename TSpec>
inline bool
isLeaf(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

//	if (it._id >= length(it.store->annotationStore)) return false;
	return it.store->annotationStore[it._id].lastChildID == TAnnotation::INVALID_ID;
}

/*!
 * @fn AnnotationTreeIterator#isLastChild
 * @brief Returns a boolean value that indicates whether the current node is the last child.
 *
 * @signature bool isLastChild(it);
 *
 * @param[in] it The AnnotationTreeIterator to query.
 */

template <typename TFragmentStore, typename TSpec>
inline bool
isLastChild(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TID					TID;

	TAnnotation const &anno = getAnnotation(it);
	TID nextSiblingID = anno.nextSiblingID;
	if (nextSiblingID != TAnnotation::INVALID_ID)
	{
		TID lastChildID = it.store->annotationStore[anno.parentID].lastChildID;
		return it._id == lastChildID;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAnnotationStore>
inline void
_storeClearAnnoBackLinks(TAnnotationStore & me)
{
	typedef typename Value<TAnnotationStore>::Type				TAnnotation;
	typedef typename Iterator<TAnnotationStore, Standard>::Type TAnnoIter;

	TAnnoIter it = begin(me, Standard());
	TAnnoIter itEnd = end(me, Standard());
	
	for (; it != itEnd; ++it)
	{
		(*it).lastChildID = TAnnotation::INVALID_ID;
		(*it).nextSiblingID = TAnnotation::INVALID_ID;
	}
}

template <typename TAnnotationStore>
inline void
_storeCreateAnnoBackLinks(TAnnotationStore & me)
{
	typedef typename Value<TAnnotationStore>::Type				TAnnotation;
	typedef typename TAnnotation::TID							TID;
	typedef typename Iterator<TAnnotationStore, Standard>::Type TAnnoIter;
	
	TAnnoIter itBegin = begin(me, Standard());
	TAnnoIter itEnd = end(me, Standard());
	TID _id = (itEnd - itBegin) - 1;
	TAnnoIter it = itBegin + _id;
	
	for (; itBegin <= it; --it, --_id)
	{
		if ((*it).parentID != TAnnotation::INVALID_ID)
		{
			TAnnoIter parent = itBegin + (*it).parentID;
			if ((*parent).lastChildID == TAnnotation::INVALID_ID)
			{
				(*parent).lastChildID = _id;
				(*it).nextSiblingID = _id;
			}

			if ((*it).nextSiblingID == TAnnotation::INVALID_ID)
			{
				TAnnoIter lastChild = itBegin + (*parent).lastChildID;
				(*it).nextSiblingID = (*lastChild).nextSiblingID;
				(*lastChild).nextSiblingID = _id;
			}
		}
		else 
			(*it).nextSiblingID = TAnnotation::INVALID_ID;
	}
}

template <typename TPos, typename TSpec>
inline std::ostream &
operator << (std::ostream & out, AnnotationStoreElement<TPos, TSpec> const & anno)
{
    out << "parentID:     \t" << anno.parentID << std::endl;
    out << "contigID:     \t" << anno.contigID << std::endl;
    out << "countID:      \t" << anno.countID << std::endl;
    out << "typeID:       \t" << anno.typeID << std::endl;
    out << "beginPos:     \t" << anno.beginPos << std::endl;
    out << "endPos:       \t" << anno.endPos << std::endl;
    out << "lastChildID:  \t" << anno.lastChildID << std::endl;
    out << "nextSiblingID:\t" << anno.nextSiblingID << std::endl;
    
    return out;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
