// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Our own implementation of a streambuf_iterator. We could not use the STL's
// iterator as we need to have access to the underlying streambuf which is a
// private member of the STL iterator.
// ==========================================================================
// TODO(esiragusa): tests

#ifndef SEQAN_STREAM_ITER_H_
#define SEQAN_STREAM_ITER_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

template <typename TDirection>
struct StreamIterator {};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class StreamIterator
// ----------------------------------------------------------------------------

// Unfortunately some of the most useful members of basic_streambuf are
// protected, so we define a subclass to cast and access them
template <typename TValue, typename TTraits_ = std::char_traits<TValue> >
class StreamBuffer: public std::basic_streambuf<TValue, TTraits_>
{
public:
    typedef TTraits_ TTraits;
    typedef std::basic_streambuf<TValue, TTraits_> TBase;

    using typename TBase::pos_type;
    using typename TBase::off_type;

    using TBase::eback;
    using TBase::gptr;
    using TBase::egptr;

    using TBase::pbase;
    using TBase::pptr;
    using TBase::epptr;

    size_t chunkSize(Input)
    {
        return egptr() - gptr();
    }

    size_t chunkSize(Output)
    {
        return epptr() - pptr();
    }

    template <typename TOffset>
    void advanceChunk(TOffset ofs, Input)
    {
        this->gbump(ofs);
    }

    template <typename TOffset>
    void advanceChunk(TOffset ofs, Output)
    {
        this->pbump(ofs);
    }

    void reserveChunk(Input)
    {
        if (gptr() == egptr())
            this->underflow();
    }

    void reserveChunk(Output)
    {
        if (pptr() == epptr())
            this->overflow(EOF);
    }

    template <typename TOffset>
    std::streampos seekoff(TOffset ofs, std::ios_base::seekdir way, std::ios_base::openmode which)
    {
        return this->seekoff(ofs, way, which);
    }

    template <typename TOffset, typename TDirection>
    void goFurther(TOffset ofs, TDirection dir)
    {
        size_t left = chunkSize(dir);
        if (SEQAN_LIKELY((size_t)ofs <= left))
        {
            advanceChunk(ofs, dir);
            return;
        }

        while (true)
        {
            size_t adv = std::min((size_t)ofs, left);
            advanceChunk(adv, dir);
            ofs -= adv;
            if (ofs == 0)
                return;

            this->underflow();
            left = chunkSize(dir);

            if (SEQAN_UNLIKELY(left == 0))
            {
                // if chunking isn't available try to seek
                pos_type res = seekoff(ofs,
                                       std::ios_base::cur,
                                       (IsSameType<TDirection, Input>::VALUE)? std::ios_base::in: std::ios_base::out);

                // if seek doesn't work manually skip characters (when reading)
                if (IsSameType<TDirection, Input>::VALUE)
                    if (res == pos_type(off_type(-1)))
                        for (; ofs != 0; --ofs)
                            this->sbumpc();

                return;
            }
        }
    }

//    template <typename TOffset, typename TDirection>
//    void goFurther(TOffset ofs, TDirection dir)
//    {
//        typedef typename MakeUnsigned<TOffset>::Type TUOffset;
//        while (ofs != 0)
//        {
//            reserveChunk(dir);
//            size_t adv = chunkSize(dir);
//
//            if (SEQAN_UNLIKELY(adv == 0))
//            {
//                // if chunking isn't available try to seek
//                pos_type res = seekoff(ofs,
//                                       std::ios_base::cur,
//                                       (IsSameType<TDirection, Input>::VALUE)? std::ios_base::in: std::ios_base::out);
//
//                // if seek doesn't work manually skip characters (when reading)
//                if (IsSameType<TDirection, Input>::VALUE)
//                    if (res == pos_type(off_type(-1)))
//                        for (; ofs != 0; --ofs)
//                            this->sbumpc();
//
//                return;
//            }
//
//            if (adv > (TUOffset)ofs)
//                adv = ofs;
//            ofs -= adv;
//
//            advanceChunk(adv, dir);
//        }
//    }

};

template <typename TStream>
class Iter<TStream, StreamIterator<Input> >
{
public:
    typedef typename Value<TStream>::Type   TValue;
    typedef std::basic_istream<TValue>      TIStream;
    typedef std::basic_streambuf<TValue>    TBasicBuffer;
    typedef StreamBuffer<TValue>            TStreamBuffer;

    TStreamBuffer *streamBuf;

    Iter():
        streamBuf()
    {}

    Iter(TIStream& stream):
        streamBuf(static_cast<StreamBuffer<TValue> *>(stream.rdbuf()))
    {}

    Iter(TStreamBuffer *buf):
        streamBuf(static_cast<StreamBuffer<TValue> *>(buf))
    {}
};

template <typename TStream>
class Iter<TStream, StreamIterator<Output> >
{
public:
    typedef typename Value<TStream>::Type   TValue;
    typedef std::basic_ostream<TValue>      TOStream;
    typedef std::basic_streambuf<TValue>    TBasicBuffer;
    typedef StreamBuffer<TValue>            TStreamBuffer;

    TStreamBuffer *streamBuf;

    Iter():
        streamBuf()
    {}

    Iter(TOStream& stream):
        streamBuf(static_cast<StreamBuffer<TValue> *>(stream.rdbuf()))
    {}

    Iter(TBasicBuffer *buf):
        streamBuf(static_cast<StreamBuffer<TValue> *>(buf))
    {}

    template <typename TValue2>
    TValue2 & operator=(TValue2 &val)
    {
        setValue(*this, val);
        return val;
    }

    template <typename TValue2>
    TValue2 const & operator=(TValue2 const &val)
    {
        setValue(*this, val);
        return val;
    }
};


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TStream>
struct Reference<Iter<TStream, StreamIterator<Input> > >:
    Value<Iter<TStream, StreamIterator<Input> > > {};

template <typename TStream>
struct Reference<Iter<TStream, StreamIterator<Input> > const>:
    Value<Iter<TStream, StreamIterator<Input> > > {};

template <typename TStream>
struct Reference<Iter<TStream, StreamIterator<Output> > >
{
    typedef Iter<TStream, StreamIterator<Output> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename TStream>
struct GetValue<Iter<TStream, StreamIterator<Input> > >:
    Reference<Iter<TStream, StreamIterator<Input> > const> {};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TStream, typename TDirection>
struct Value<Iter<TStream, StreamIterator<TDirection> > > : Value<TStream> {};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TStream, typename TDirection>
struct Position<Iter<TStream, StreamIterator<TDirection> > > : Position<TStream> {};

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <typename TStream, typename TDirection>
struct Difference<Iter<TStream, StreamIterator<TDirection> > > : Difference<TStream> {};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TStream, typename TDirection>
struct Size<Iter<TStream, StreamIterator<TDirection> > > : Size<TStream> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TObject, typename TDirection>
struct DirectionIterator:
    If<Is<StreamConcept<TObject> >,
       Iter<TObject, StreamIterator<TDirection> >,
       typename Iterator<TObject, Rooted>::Type>
{};

// ----------------------------------------------------------------------------
// Function value() - Input
// ----------------------------------------------------------------------------

template <typename TStream>
inline typename Reference<Iter<TStream, StreamIterator<Input> > >::Type
value(Iter<TStream, StreamIterator<Input> > &iter)
{
    SEQAN_ASSERT(iter.streamBuf != NULL);
    return iter.streamBuf->sgetc();
}
template <typename TStream>
inline typename Reference<Iter<TStream, StreamIterator<Input> > const>::Type
value(Iter<TStream, StreamIterator<Input> > const &iter)
{
    SEQAN_ASSERT(iter.streamBuf != NULL);
    return iter.streamBuf->sgetc();
}

// ----------------------------------------------------------------------------
// Function value() - Ouput
// ----------------------------------------------------------------------------

template <typename TStream>
inline Iter<TStream, StreamIterator<Output> > &
value(Iter<TStream, StreamIterator<Output> > & iter)
{
    return iter;
}
template <typename TStream>
inline Iter<TStream, StreamIterator<Output> > const &
value(Iter<TStream, StreamIterator<Output> > const & iter)
{
    return iter;
}

template <typename TStream, typename TValue>
inline void
setValue(Iter<TStream, StreamIterator<Output> > & iter, TValue const &val)
{
    return setValue(const_cast<Iter<TStream, StreamIterator<Output> > const &>(iter), val);
}
template <typename TStream, typename TValue>
inline void
setValue(Iter<TStream, StreamIterator<Output> > const & iter, TValue const &val)
{
    SEQAN_ASSERT(iter.streamBuf != NULL);
    iter.streamBuf->sputc(val);
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

template <typename TStream>
inline void
goNext(Iter<TStream, StreamIterator<Input> > & iter)
{
    SEQAN_ASSERT(iter.streamBuf != NULL);
    iter.streamBuf->sbumpc();
}

template <typename TStream>
inline void
goNext(Iter<TStream, StreamIterator<Output> > &)
{
    // We do nothing here, as the stream is advanced by sputc whenever you assign
    // a value to the iterator with *iter= or setValue
}

// ----------------------------------------------------------------------------
// Function goFurther()
// ----------------------------------------------------------------------------

template <typename TStream, typename TOffset, typename TDirection>
inline void
goFurther(Iter<TStream, StreamIterator<TDirection> > &iter, TOffset ofs)
{
    SEQAN_ASSERT(iter.streamBuf != NULL);
    iter.streamBuf->goFurther(ofs, TDirection());
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

template <typename TStream, typename TDirection>
inline typename Position<Iter<TStream, StreamIterator<TDirection> > const>::Type
position(Iter<TStream, StreamIterator<TDirection> > const & iter)
{
    SEQAN_ASSERT(iter.streamBuf != NULL);
    iter.stream->seekpos(0, std::ios_base::cur, (IsSameType<TDirection, Input>::VALUE)? std::ios_base::in: std::ios_base::out);
}

// ----------------------------------------------------------------------------
// Function setPosition()
// ----------------------------------------------------------------------------

template <typename TStream, typename TDirection, typename TPosition>
inline void
setPosition(Iter<TStream, StreamIterator<TDirection> > const & iter, TPosition pos)
{
    SEQAN_ASSERT(iter.streamBuf != NULL);
    iter.stream->seekpos(pos, (IsSameType<TDirection, Input>::VALUE)? std::ios_base::in: std::ios_base::out);
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TStream>
inline bool
atEnd(Iter<TStream, StreamIterator<Input> > const & iter)
{
    typedef typename Value<Iter<TStream, StreamIterator<Input> > >::Type TValue;
    typedef StreamBuffer<TValue> TStreamBuffer;

    if (SEQAN_UNLIKELY(iter.streamBuf == NULL))
    {
        return true;
    }
    else
    {
        TStreamBuffer *buf = static_cast<TStreamBuffer*>(iter.streamBuf);
        if (SEQAN_LIKELY(buf->gptr() < buf->egptr()))
            return false;
        else
            return TStreamBuffer::TTraits::eq_int_type(buf->sgetc(), TStreamBuffer::TTraits::eof());
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_ITER_H_