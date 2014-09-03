// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Smart file for reading/writing files in GFF or GTF format.
// ==========================================================================

#ifndef SEQAN_GFF_IO_GFF_FILE_H_
#define SEQAN_GFF_IO_GFF_FILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Typedefs
// ============================================================================

typedef SmartFile<Gff, Input>   GffFileIn;
typedef SmartFile<Gff, Output>  GffFileOut;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction SmartFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct SmartFileContext<SmartFile<Gff, TDirection, TSpec>, TStorageSpec>
{
    typedef StringSet<CharString>                                   TNameStore;
    typedef NameStoreCache<TNameStore>                              TNameStoreCache;
    typedef GffIOContext<TNameStore, TNameStoreCache, TStorageSpec> Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormats<SmartFile<Gff, TDirection, TSpec> >
{
    typedef TagSelector<
                TagList<Gff,
                TagList<Gtf
                > >
            > Type;
};

// --------------------------------------------------------------------------
// _mapGffFormatToCompressionFormat
// --------------------------------------------------------------------------

inline Nothing
_mapFileFormatToCompressionFormat(Gff)
{
    return Nothing();
}

inline Nothing
_mapFileFormatToCompressionFormat(Gtf)
{
    return Nothing();
}

// ----------------------------------------------------------------------------
// Function read(); GffRecord
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec, typename TForwardIter>
inline void
readRecord(GffRecord & /* record */,
           GffIOContext<TNameStore, TNameStoreCache, TStorageSpec> & /* context */,
           TForwardIter & /* iter */,
           TagSelector<> const & /* format */)
{
    SEQAN_FAIL("GffFileIn: File format not specified.");
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec, typename TForwardIter, typename TTagList>
inline void
readRecord(GffRecord & record,
           GffIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        readRecord(record, context, iter, TFormat());
    else
        readRecord(record, context, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// convient GffFile variant
template <typename TSpec>
inline void
readRecord(GffRecord & record, SmartFile<Gff, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function write(); GffRecord
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
write(TTarget & /* target */,
      GffRecord & /* record */,
      GffIOContext<TNameStore, TNameStoreCache, TStorageSpec> & /* context */,
      TagSelector<> const & /* format */)
{
    SEQAN_FAIL("GffFileOut: File format not specified.");
}

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec, typename TTagList>
inline void
writeRecord(TTarget & target,
            GffRecord & record,
            GffIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
            TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        write(target, record, context, TFormat());
    else
        write(target, record, context, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

template <typename TSpec>
inline void
writeRecord(SmartFile<Gff, Output, TSpec> & file, GffRecord & record)
{
    write(file.iter, record, context(file), file.format);
}

}  // namespace seqan

#endif // SEQAN_GFF_IO_GFF_FILE_H_