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

#ifndef SEQAN_HEADER_STORE_IO_GFF_H
#define SEQAN_HEADER_STORE_IO_GFF_H

namespace SEQAN_NAMESPACE_MAIN {

//////////////////////////////////////////////////////////////////////////////
// Read Gff
//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec = void>
struct IOContextGff_
{
    typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotation;
    typedef typename TAnnotation::TID                   TID;

    CharString contigName;
    CharString typeName;
    CharString annotationName;
    CharString parentKey;
    CharString parentName;

    CharString _key;
    CharString _value;
    StringSet<CharString> keys;
    StringSet<CharString> values;

    CharString gtfGeneID;
    CharString gtfGeneName;
    CharString gtfTranscriptName;       // transcipt_id is stored in parentName

    TID annotationID;
    TAnnotation annotation;
};

template <typename TFragmentStore, typename TSpec>
inline void clear(IOContextGff_<TFragmentStore, TSpec> & ctx)
{
    typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotation;

    clear(ctx.contigName);
    clear(ctx.typeName);
    clear(ctx.annotationName);
    clear(ctx.parentKey);
    clear(ctx.parentName);
    clear(ctx._key);
    clear(ctx._value);
    clear(ctx.gtfGeneID);
    clear(ctx.gtfGeneName);
    clear(ctx.gtfTranscriptName);
    clear(ctx.keys);
    clear(ctx.values);
    ctx.annotationID = TAnnotation::INVALID_ID;
    clear(ctx.annotation.values);
}

//////////////////////////////////////////////////////////////////////////////
// _readOneAnnotation
//
// reads in one annotation line from a Gff file

template <typename TFragmentStore, typename TSpec>
inline void
_readOneAnnotation(
    IOContextGff_<TFragmentStore, TSpec> & ctx,
    GffRecord const & record)
{
//IOREV _nodoc_ _hasCRef_
    typedef typename TFragmentStore::TContigPos         TContigPos;

    clear(ctx);

    // read column 1: contig name
    ctx.contigName = record.ref;

    // skip column 2
    // read column 3: type
    ctx.typeName = record.type;

    // read column 4 and 5: begin and endposition
    ctx.annotation.beginPos = record.beginPos;
    ctx.annotation.endPos = record.endPos;

    // skip column 6
    // read column 7: orientation
    if (record.strand == '-')
    {
        TContigPos tmp = ctx.annotation.beginPos;
        ctx.annotation.beginPos = ctx.annotation.endPos;
        ctx.annotation.endPos = tmp;
    }

    // skip column 8
    // read column 9: name
    for (unsigned i = 0; i < length(record.tagNames); ++i)
    {
        ctx._key = record.tagNames[i];
        ctx._value = record.tagValues[i];
        if (ctx._key == "ID")
        {
            ctx.annotationName = ctx._value;
        }
        else if (!empty(ctx._key) && !empty(ctx._value))
        {
            appendValue(ctx.keys, ctx._key);
            appendValue(ctx.values, ctx._value);
        }

        if (ctx._key == "Parent" || ctx._key == "ParentID" || ctx._key == "transcript_id")
        {
            ctx.parentKey = ctx._key;
            ctx.parentName = ctx._value;
        }
        else if (ctx._key == "transcript_name")
        {
            ctx.gtfTranscriptName = ctx._value;
        }
        else if (ctx._key == "gene_id")
        {
            ctx.gtfGeneID = ctx._value;
        }
        else if (ctx._key == "gene_name")
        {
            ctx.gtfGeneName = ctx._value;
        }

        clear(ctx._key);
        clear(ctx._value);
    }
}

template <typename TAnnotation>
inline void
_adjustParent(
    TAnnotation & parent,
    TAnnotation const & child)
{
    if (child.contigID == TAnnotation::INVALID_ID || child.beginPos == TAnnotation::INVALID_POS || child.endPos == TAnnotation::INVALID_POS)
        return;

    parent.contigID = child.contigID;

    // Has parent an invalid begin and end position?
    if ((parent.beginPos == TAnnotation::INVALID_POS) && (parent.endPos == TAnnotation::INVALID_POS))
    {
        parent.beginPos = child.beginPos;
        parent.endPos = child.endPos;
        return;
    }

    if ((parent.beginPos == TAnnotation::INVALID_POS) || (parent.endPos == TAnnotation::INVALID_POS))
        return;

    typename TAnnotation::TPos childBegin, childEnd;
    if (child.beginPos < child.endPos)
    {
        childBegin = child.beginPos;
        childEnd = child.endPos;
    }
    else
    {
        childBegin = child.endPos;
        childEnd = child.beginPos;
    }

    // Keep parent's orientation and maximize begin and end using child's boundaries.
    if (parent.beginPos < parent.endPos)
    {
        if (parent.beginPos == TAnnotation::INVALID_POS || parent.beginPos > childBegin)
            parent.beginPos = childBegin;
        if (parent.endPos == TAnnotation::INVALID_POS || parent.endPos < childEnd)
            parent.endPos = childEnd;
    }
    else
    {
        if (parent.endPos == TAnnotation::INVALID_POS || parent.endPos > childBegin)
            parent.endPos = childBegin;
        if (parent.beginPos == TAnnotation::INVALID_POS || parent.beginPos < childEnd)
            parent.beginPos = childEnd;
    }
}

template <typename TFragmentStore, typename TSpec>
inline void
_storeOneAnnotation(
    TFragmentStore & fragStore,
    IOContextGff_<TFragmentStore, TSpec> & ctx)
{
    typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotation;
    typedef typename TAnnotation::TID                   TID;

    SEQAN_ASSERT_EQ(length(fragStore.annotationStore), length(fragStore.annotationNameStore));

    // for lines in Gtf format get/add the parent gene first
    TID geneID = TAnnotation::INVALID_ID;
    if (!empty(ctx.gtfGeneID))
        _storeAppendAnnotationName(fragStore, geneID, ctx.gtfGeneID, (TID) TFragmentStore::ANNO_GENE);

    // if we have a parent transcript, get/add the parent transcript then
    if (!empty(ctx.parentName))
    {
// From now, we support gtf files with genes/transcripts having the same name.
//
//        // if gene and transcript names are equal (like in some strange gtf files)
//        // try to make the transcript name unique
//        if (ctx.gtfGeneID == ctx.parentName)
//            append(ctx.parentName, "_1");

        if (ctx.parentKey == "transcript_id")
            // type is implicitly given (mRNA)
            _storeAppendAnnotationName(fragStore, ctx.annotation.parentID, ctx.parentName, (TID) TFragmentStore::ANNO_MRNA);
        else
            // type is unknown
            _storeAppendAnnotationName(fragStore, ctx.annotation.parentID, ctx.parentName);
    }
    else
        ctx.annotation.parentID = 0;    // if we have no parent, we are a child of the root

    // add contig and type name
    _storeAppendContig(fragStore, ctx.annotation.contigID, ctx.contigName);
    _storeAppendType(fragStore, ctx.annotation.typeID, ctx.typeName);

    // add annotation name of the current line
    _storeAppendAnnotationName(fragStore, ctx.annotationID, ctx.annotationName, ctx.annotation.typeID);

    for (unsigned i = 0; i < length(ctx.keys); ++i)
    {
        // don't store gene_name as key/value pair unless it is a gene
        if (ctx.keys[i] == "gene_name" && ctx.annotation.typeID != TFragmentStore::ANNO_GENE)
            continue;

        // don't store transcript_name as key/value pair unless it is a transcript
        if (ctx.keys[i] == "transcript_name" && ctx.annotation.typeID != TFragmentStore::ANNO_MRNA)
            continue;

        // don't store Parent, transcript_id or gene_id as key/value pair (the are used to link annotations)
        if (ctx.keys[i] != ctx.parentKey && ctx.keys[i] != "gene_id")
            annotationAssignValueByKey(fragStore, ctx.annotation, ctx.keys[i], ctx.values[i]);
    }

    fragStore.annotationStore[ctx.annotationID] = ctx.annotation;

    TAnnotation & parent = fragStore.annotationStore[ctx.annotation.parentID];
    if (ctx.annotation.parentID != 0 && parent.parentID == TAnnotation::INVALID_ID)
        parent.parentID = 0;    // if our parent has no parent, it becomes a child of the root

    if (geneID != TAnnotation::INVALID_ID)
    {
        // link and adjust our gtf ancestors
        TAnnotation & gene = fragStore.annotationStore[geneID];
//        TAnnotation & transcript = fragStore.annotationStore[ctx.annotation.parentID];

        gene.parentID = 0;
        gene.typeID = TFragmentStore::ANNO_GENE;
        _adjustParent(gene, ctx.annotation);

        if (!empty(ctx.gtfGeneName))
            annotationAssignValueByKey(fragStore, gene, "gene_name", ctx.gtfGeneName);

        parent.parentID = geneID;
        parent.typeID = TFragmentStore::ANNO_MRNA;
        _adjustParent(parent, ctx.annotation);
        if (!empty(ctx.gtfTranscriptName))
            annotationAssignValueByKey(fragStore, parent, "transcript_name", ctx.gtfTranscriptName);
    }
}

template <typename TSpec, typename TConfig>
inline void
readRecords(FragmentStore<TSpec, TConfig> & fragStore,
            GffFileIn & gffFile)
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;

    if (atEnd(gffFile))
        return;

    refresh(fragStore.contigNameStoreCache);
    refresh(fragStore.annotationNameStoreCache);
    refresh(fragStore.annotationTypeStoreCache);

    GffRecord record;
    IOContextGff_<TFragmentStore> ctx;

    while (!atEnd(gffFile))
    {
        readRecord(record, gffFile);
        _readOneAnnotation(ctx, record);
        _storeOneAnnotation(fragStore, ctx);
    }
    _storeClearAnnoBackLinks(fragStore.annotationStore);
    _storeCreateAnnoBackLinks(fragStore.annotationStore);
    _storeRemoveTempAnnoNames(fragStore);
}

//////////////////////////////////////////////////////////////////////////////
// Write Gff
//////////////////////////////////////////////////////////////////////////////

// This function write the information that are equal for gff and gtf files.
template <typename TSpec, typename TConfig, typename TAnnotation, typename TID>
inline void
_writeCommonGffGtfInfo(
    GffRecord & record,
    FragmentStore<TSpec, TConfig> & store,
    TAnnotation & annotation,
    TID /*id*/)
{
    typedef FragmentStore<TSpec, TConfig>       TFragmentStore;
    typedef typename TFragmentStore::TContigPos TContigPos;

    clear(record);

    // write column 1: contig name
    if (annotation.contigID < length(store.contigNameStore))
        if (length(store.contigNameStore[annotation.contigID]) > 0u)
            record.ref = store.contigNameStore[annotation.contigID];

    // skip column 2: source
    record.source = ".";

    // write column 3: type
    if (annotation.typeID < length(store.annotationTypeStore))
        if (length(store.annotationTypeStore[annotation.typeID]) > 0u)
            record.type = store.annotationTypeStore[annotation.typeID];

    TContigPos beginPos = annotation.beginPos;
    TContigPos endPos = annotation.endPos;
    char orientation = '+';
    if (endPos < beginPos)
    {
        TContigPos tmp = beginPos;
        beginPos = endPos;
        endPos = tmp;
        orientation = '-';
    }

    // write column 4: begin position
    if (beginPos != TAnnotation::INVALID_POS)
        record.beginPos = beginPos;

    // write column 5: end position
    if (endPos != TAnnotation::INVALID_POS)
        record.endPos = endPos;

    // skip column 6: score

    // write column 7: orientation
    record.strand = orientation;
}

template <typename TRecord, typename TSpec, typename TConfig, typename TAnnotation, typename TID>
inline bool
_fillAnnotationRecord(
    TRecord & record,
    FragmentStore<TSpec, TConfig> & store,
    TAnnotation & annotation,
    TID id,
    Gff)
{
    if (id == 0)
        return false;

    _writeCommonGffGtfInfo(record, store, annotation, id);

    // write column 9: group
    // write column 9.1: annotation id
    if (id < length(store.annotationNameStore) && !empty(getAnnoName(store, id)))
    {
        appendValue(record.tagNames, "ID");
        appendValue(record.tagValues, getAnnoName(store, id));
    }
    else if (annotation.lastChildID != TAnnotation::INVALID_ID)
    {
        appendValue(record.tagNames, "ID");
        appendValue(record.tagValues, getAnnoUniqueName(store, id));
    }

    // write column 9.2: parent id
    if (store.annotationStore[annotation.parentID].typeID > 1)  // ignore root/deleted nodes
    {
        appendValue(record.tagNames, "Parent");
        appendValue(record.tagValues, getAnnoUniqueName(store, annotation.parentID));
    }

    // write column 9.3-...: key, value pairs
    for (unsigned keyID = 0; keyID < length(annotation.values); ++keyID)
        if (!empty(annotation.values[keyID]))
        {
            appendValue(record.tagNames, store.annotationKeyStore[keyID]);
            appendValue(record.tagValues, annotation.values[keyID]);
        }

    return true;
}

template <typename TRecord, typename TSpec, typename TConfig, typename TAnnotation, typename TID>
inline bool
_fillAnnotationRecord(
    TRecord & record,
    FragmentStore<TSpec, TConfig> & store,
    TAnnotation & annotation,
    TID id,
    Gtf)
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;

    if (annotation.typeID <= TFragmentStore::ANNO_MRNA)
        return false;

    _writeCommonGffGtfInfo(record, store, annotation, id);

    // write column 9: group

    // step up until we reach a transcript
    TID transcriptID = annotation.parentID;
    while (transcriptID < length(store.annotationStore) && store.annotationStore[transcriptID].typeID != TFragmentStore::ANNO_MRNA)
        transcriptID = store.annotationStore[transcriptID].parentID;

    // step up until we reach a gene
    TID geneID = transcriptID;
    while (geneID < length(store.annotationStore) && store.annotationStore[geneID].typeID != TFragmentStore::ANNO_GENE)
        geneID = store.annotationStore[geneID].parentID;

    typename ID<TAnnotation>::Type valueID;
    if (geneID < length(store.annotationStore) &&
        (valueID = annotationGetValueIDByKey(store, store.annotationStore[geneID], "gene_name")) != TAnnotation::INVALID_ID)
    {
        appendValue(record.tagNames, "gene_name");
        appendValue(record.tagValues, store.annotationStore[geneID].values[valueID]);
    }
    if (transcriptID < length(store.annotationStore) &&
        (valueID = annotationGetValueIDByKey(store, store.annotationStore[transcriptID], "transcript_name")) != TAnnotation::INVALID_ID)
    {
        appendValue(record.tagNames, "transcript_name");
        appendValue(record.tagValues, store.annotationStore[transcriptID].values[valueID]);
    }

    if (id < length(store.annotationNameStore) && !empty(getAnnoName(store, id)))
    {
        appendValue(record.tagNames, "ID");
        appendValue(record.tagValues, getAnnoName(store, id));
    }

    // write key, value pairs
    for (unsigned keyID = 0; keyID < length(annotation.values); ++keyID)
        if (!empty(annotation.values[keyID]))
        {
            appendValue(record.tagNames, store.annotationKeyStore[keyID]);
            appendValue(record.tagValues, annotation.values[keyID]);
        }

    // The GTF format version 2.2 requires the keys gene_id and transcript_id to be the last keys of line
    // read http://mblab.wustl.edu/GTF22.html and http://www.bioperl.org/wiki/GTF

    if (geneID < length(store.annotationStore))
    {
        appendValue(record.tagNames, "gene_id");
        appendValue(record.tagValues, getAnnoUniqueName(store, geneID));
    }

    if (transcriptID < length(store.annotationStore))
    {
        appendValue(record.tagNames, "transcript_id");
        appendValue(record.tagValues, getAnnoUniqueName(store, transcriptID));
    }
    return true;
}

// support for dynamically chosen file formats
template <typename TRecord, typename TSpec, typename TConfig, typename TAnnotation, typename TID>
inline bool
_fillAnnotationRecord(
    TRecord & /*record*/,
    FragmentStore<TSpec, TConfig> & /*store*/,
    TAnnotation & /*annotation*/,
    TID /*id*/,
    TagSelector<> const & /*format*/)
{
    SEQAN_FAIL("AnnotationStore: File format not specified.");
    return false;
}

template <typename TRecord, typename TSpec, typename TConfig, typename TAnnotation, typename TID, typename TTagList>
inline bool
_fillAnnotationRecord(
    TRecord & record,
    FragmentStore<TSpec, TConfig> & store,
    TAnnotation & annotation,
    TID id,
    TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        return _fillAnnotationRecord(record, store, annotation, id, TFormat());
    else
        return _fillAnnotationRecord(record, store, annotation, id, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}


template <typename TTargetStream, typename TSpec, typename TConfig, typename TFormat>
inline void
_writeGffGtf(
    TTargetStream & target,
    FragmentStore<TSpec, TConfig> & store,
    TFormat const &format)
{
    typedef FragmentStore<TSpec, TConfig>                           TFragmentStore;
    typedef typename TFragmentStore::TAnnotationStore               TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type                  TAnnotation;
    typedef typename Iterator<TAnnotationStore, Standard>::Type     TAnnoIter;
    typedef typename ID<TAnnotation>::Type                          TID;

    TAnnoIter it = begin(store.annotationStore, Standard());
    TAnnoIter itEnd = end(store.annotationStore, Standard());

    GffRecord record;
    typename DirectionIterator<TTargetStream, Output>::Type iter = directionIterator(target, Output());

    for (TID id = 0; it != itEnd; ++it, ++id)
    {
        if (_fillAnnotationRecord(record, store, *it, id, format))
            writeRecord(iter, record, format);
    }
}

template <typename TSpec, typename TFSSpec, typename TFSConfig>
inline void
writeRecords(FormattedFile<Gff, Output, TSpec> & gffFile,
             FragmentStore<TFSSpec, TFSConfig> & store)
{
    _writeGffGtf(gffFile, store, format(gffFile));
}

} // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
