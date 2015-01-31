#ifndef APPS_BS_TOOLS_BISAR_H_
#define APPS_BS_TOOLS_BISAR_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <fstream>
#include <seqan/file.h>
#include <seqan/store.h>
#include <seqan/bam_io.h>
#include <seqan/score.h>

using namespace seqan;

class Times
{
public:
    double time_all;
    double time_globalAlignment;
    double time_writeBsAlignment;

    static Times & instance()
    {
        static Times times;
        return times;
    }

private:
    Times() :
        time_all(0),
        time_globalAlignment(0),
        time_writeBsAlignment(0)
    {}
};


// Realign with 4-letter alphabet and set alignedReadStore entries
template<typename TReadGaps, typename TContigGaps, typename TFragmentStore, typename TID, typename TBsScoreCTLeft, typename TBsScoreCTRight, typename TBsScoreGALeft, typename TBsScoreGARight, typename TOptions>
inline int
reAlign4(TReadGaps &readGaps,
         TContigGaps &contigGaps,
         TFragmentStore &store,
         TID &id,
         TBsScoreCTLeft &scoringSchemeCTLeft,
         TBsScoreCTRight &scoringSchemeCTRight,
         TBsScoreGALeft &scoringSchemeGALeft,
         TBsScoreGARight &scoringSchemeGARight,
         TOptions &options)
{
    typedef typename Iterator<TReadGaps>::Type      TReadGapsIterator;
    typedef typename Iterator<TContigGaps>::Type    TContigGapsIterator;

    typedef double TValue;
    TValue scoreA;

#ifdef POST_PRO_PROFILE
    double timeStamp = sysTime();
#endif
    // Do alignment
    //std::cout << " Do alignment: " << std::endl;
    int band = 6;   // dep. on interOffset?
    if (getMateNo(store, store.alignedReadStore[id].readID) != -1 )     // If paired
    {
        if (getMateNo(store, store.alignedReadStore[id].readID) == 0 && store.alignedReadStore[id].beginPos < store.alignedReadStore[id].endPos)
        {
            scoreA = (double)globalAlignment(contigGaps, readGaps, scoringSchemeCTLeft, AlignConfig<true, false, false, true>(), -band, band)/10000.0;    // convert back from int to double
        }
        else if (getMateNo(store, store.alignedReadStore[id].readID) == 0 && store.alignedReadStore[id].beginPos > store.alignedReadStore[id].endPos)
        {
            scoreA = (double)globalAlignment(contigGaps, readGaps, scoringSchemeGALeft, AlignConfig<true, false, false, true>(), -band, band)/10000.0;
        }
        else if (getMateNo(store, store.alignedReadStore[id].readID) == 1 && store.alignedReadStore[id].beginPos > store.alignedReadStore[id].endPos)
        {
            scoreA = (double)globalAlignment(contigGaps, readGaps, scoringSchemeCTRight, AlignConfig<true, false, false, true>(), -band, band)/10000.0;
        }
        else //if (getMateNo(store, store.alignedReadStore[id].readID) == 1 && store.alignedReadStore[id].beginPos < store.alignedReadStore[id].endPos)
        {
            scoreA = (double)globalAlignment(contigGaps, readGaps, scoringSchemeGARight, AlignConfig<true, false, false, true>(), -band, band)/10000.0;
        }
    }
    else
    {
        if (store.alignedReadStore[id].beginPos < store.alignedReadStore[id].endPos)
        {
            scoreA = (double)globalAlignment(contigGaps, readGaps, scoringSchemeCTLeft, AlignConfig<true, false, false, true>(), -band, band)/10000.0;
        }
        else
        {
            scoreA = (double)globalAlignment(contigGaps, readGaps, scoringSchemeGALeft, AlignConfig<true, false, false, true>(), -band, band)/10000.0;
        }
    }
#ifdef POST_PRO_PROFILE
    Times::instance().time_globalAlignment += (sysTime() - timeStamp);
#endif

    // Get number of errors (mismatches, indels)
    // Remove gaps at beginning and end of read caused by inaccurate begin and end positions of contig infix
    TReadGapsIterator itR = begin(readGaps);
    TContigGapsIterator itC = begin(contigGaps);
    while (isGap(itR))
    {
        ++itR;
        ++itC;
    }
    setBeginPosition(readGaps, beginPosition(readGaps));

    if (itC != begin(contigGaps)) // If read has gaps at beginning -> itC != begin(contigGaps) -> set new clipped position; otherwise do nothing
        setClippedBeginPosition(contigGaps, position(itC));
    TReadGapsIterator itREnd = end(readGaps);
    unsigned countEndGaps = 0;
    while(isGap(--itREnd))    // start behind last position?
        ++countEndGaps;

    setEndPosition(readGaps, endPosition(readGaps)); // ? Set the end position of the clipped gapped sequence, given a source position
    if (countEndGaps > 0) setEndPosition(contigGaps, endPosition(contigGaps)-countEndGaps); // Only shorten if there are gaps at end of read, not if gaps at end of contig

    itREnd = end(readGaps);
    itR = begin(readGaps);
    itC = begin(contigGaps);
    unsigned gaps = 0;
    unsigned mismatches = 0;
    unsigned matches = 0;

    if (false)
    {
        std::cout << "align: (reAlign4 2) " << store.readNameStore[store.alignedReadStore[id].readID] << std::endl;
        std::cout << contigGaps << std::endl;
        std::cout << readGaps << std::endl;
    }

    // Set beginPos and endPos new
    if (store.alignedReadStore[id].beginPos < store.alignedReadStore[id].endPos)
    {
        if (store.alignedReadStore[id].beginPos < options.intervalOffset)
            store.alignedReadStore[id].beginPos = 0 + beginPosition(contigGaps);
        else
            store.alignedReadStore[id].beginPos = store.alignedReadStore[id].beginPos - options.intervalOffset + beginPosition(contigGaps);

        store.alignedReadStore[id].endPos = store.alignedReadStore[id].beginPos + endPosition(contigGaps) - beginPosition(contigGaps);
    }
    else
    {
        if (store.alignedReadStore[id].beginPos < options.intervalOffset)
            store.alignedReadStore[id].endPos = 0 + beginPosition(contigGaps);
        else
            store.alignedReadStore[id].endPos = store.alignedReadStore[id].endPos - options.intervalOffset + beginPosition(contigGaps);
        store.alignedReadStore[id].beginPos = store.alignedReadStore[id].endPos + endPosition(contigGaps) - beginPosition(contigGaps);
    }
    // Count errors
    if (getMateNo(store, store.alignedReadStore[id].readID) != -1 )     // If paired
    {
        if ((getMateNo(store, store.alignedReadStore[id].readID) == 0 && store.alignedReadStore[id].beginPos < store.alignedReadStore[id].endPos) ||    // Mapped against CT ref
            (getMateNo(store, store.alignedReadStore[id].readID) == 1 && store.alignedReadStore[id].beginPos > store.alignedReadStore[id].endPos) )
            for (; !atEnd(itR) && !atEnd(itC); ++itR, ++itC)
            {
                if (isGap(itR) || isGap(itC)) ++gaps;
                else if (*itR == *itC || (*itR == 'T' && *itC == 'C')) ++matches;
                else ++mismatches;
            }
        else    // Mapped against GA ref
            for (; !atEnd(itR) && !atEnd(itC); ++itR, ++itC)
            {
                if (isGap(itR) || isGap(itC)) ++gaps;
                else if (*itR == *itC || (*itR == 'A' && *itC == 'G') ) ++matches;
                else ++mismatches;
            }
    }
    else
    {
        if (store.alignedReadStore[id].beginPos < store.alignedReadStore[id].endPos)
            for (; !atEnd(itR) && !atEnd(itC); ++itR, ++itC)
            {
                if (isGap(itR) || isGap(itC)) ++gaps;
                else if (*itR == *itC || (*itR == 'T' && *itC == 'C') ) ++matches;
                else ++mismatches;
            }
        else
            for (; !atEnd(itR) && !atEnd(itC); ++itR, ++itC)
            {
                if (isGap(itR) || isGap(itC)) ++gaps;
                else if (*itR == *itC || (*itR == 'A' && *itC == 'G')) ++matches;
                else ++mismatches;
            }
    }
    // Rescale score and update store entries
    resize(store.alignQualityStore, length(store.alignedReadStore), Generous());
    store.alignQualityStore[id].errors = gaps + mismatches;
    //store.alignQualityStore[id].score = (scoreA - (double)(mismatches+matches+gaps)*1)*(-10); // scale back to original scaling and convert to phred scale
    store.alignQualityStore[id].score = scoreA*(-10);
    store.alignedReadStore[id].gaps = _dataAnchors(readGaps);

    return 0;
}


template <typename TFragmentStore, typename TContigGaps, typename TID, typename TScore, typename TOptions>
inline int
writeBsAlignment(seqan::BamFileOut & bamFileOut,
                 TFragmentStore &store,
                 TContigGaps &contigGaps,
                 TID &bestID,
                 TScore &mapq,
                 BamAlignmentRecord &record,
                 TOptions &/*options*/)
{
#ifdef POST_PRO_PROFILE
    double timeStamp = sysTime();
#endif
    typedef typename TFragmentStore::TAlignedReadStore                          TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type                             TAlignedRead;
    typedef typename TAlignedRead::TGapAnchors                                  TReadGapAnchors;
    typedef typename TFragmentStore::TReadSeq                                   TReadSeq;
    typedef Gaps<TReadSeq, AnchorGaps<TReadGapAnchors> >                        TReadGaps;

    // Create (mate independent) record entries
    record.qName = store.readNameStore[store.alignedReadStore[bestID].readID];
    if (store.alignedReadStore[bestID].beginPos < store.alignedReadStore[bestID].endPos)
        record.beginPos = store.alignedReadStore[bestID].beginPos;
    else
    {
        record.beginPos = store.alignedReadStore[bestID].endPos;
        record.flag |= 0x0010;
    }
	record.mapQ = (unsigned)mapq;

    CharString md;
    String<CigarElement<> > cigar;

    TReadSeq readSeq = store.readSeqStore[store.alignedReadStore[bestID].readID];

    if (store.alignedReadStore[bestID].beginPos < store.alignedReadStore[bestID].endPos)
    {
        record.seq = readSeq;
    }
    else
    {
        Dna5String tmpReadSeq = readSeq;
        record.seq = Dna5StringReverseComplement(tmpReadSeq);
    }
    TReadGaps readGaps(record.seq, store.alignedReadStore[bestID].gaps);

    getCigarString(cigar, contigGaps, readGaps);
    getMDString(md, contigGaps, readGaps);
    record.cigar = cigar;

    clear(record.qual);
    resize(record.qual, length(readSeq));
    unsigned avgQual = 0;
    if (store.alignedReadStore[bestID].beginPos < store.alignedReadStore[bestID].endPos)
    {
        unsigned i = 0;
        for (unsigned j = 0; j < length(readSeq); ++j, ++i)
        {
            record.qual[i] = (char)(getQualityValue(readSeq[j]) + 33);
            avgQual += getQualityValue(readSeq[j]);
        }
        avgQual = avgQual/(i+1);
    }
    else
    {
        unsigned i = 0;
        for (int j = length(readSeq) -1; j >= 0; --j, ++i)
        {
            record.qual[i] = (char)(getQualityValue(readSeq[j]) + 33);
            avgQual += getQualityValue(readSeq[j]);
        }
        avgQual = avgQual/(i+1);
    }

    if (store.alignedReadStore[bestID].contigID == TAlignedRead::INVALID_ID)
    {
        record.rID = seqan::BamAlignmentRecord::INVALID_REFID;
        record.beginPos = -1;
    }
    else
    {
        record.rID = store.alignedReadStore[bestID].contigID;
    }

    seqan::BamTagsDict tagsDict(record.tags);
    setTagValue(tagsDict, "NM", (int)store.alignQualityStore[bestID].errors);
    if (!empty(md))
        setTagValue(tagsDict, "MD", md);
    // for testing
    setTagValue(tagsDict, "DD", (int)avgQual);

#ifdef POST_PRO_PROFILE
    Times::instance().time_writeBsAlignment += (sysTime() - timeStamp);
#endif

    writeRecord(bamFileOut, record);

    return 0;
}

// TODO Does it make sense to count errors? (with score it would have to be exactly the same alignment then) threshold?

template<typename TFragmentStore, typename TID>
inline unsigned
countEqualHits(TFragmentStore &store, TID &id)
{
    // Count number of hits which have the same number of errors as second best hit
    unsigned count = 0;
    unsigned errors = store.alignQualityStore[id].errors;
    for (unsigned i = 0; i < length(store.alignQualityStore); ++i)
            if (store.alignQualityStore[i].errors == errors)
                ++count;
    return count;
}

// Compute pseud worst score for read
// Assuming errors are at best quality positions
// No gaps assumed
// Take max. allowed errors + 1
template<typename TScore, typename TReadSeq, typename TOptions>
inline void
computePseudoWorstScore(TScore &score, TReadSeq &readSeq, TOptions &options)
{
    // Get best qualitites

    String<int> bestQuals;
    String<int> lowQuals;
    unsigned max3Errors = floor(options.max3Error * length(readSeq)/ 100.0);
    resize(bestQuals, max3Errors + 1, 0);
    for (unsigned i = 0; i < length(readSeq); ++i)
    {
        // Quite inefficient, but in practice there will be only around 2 errors
        // So sorting will be fast
        if (getQualityValue(readSeq[i]) > back(bestQuals))
        {
            int j = max3Errors - 1;
            for(; j >= 0; --j)
            {
                if (getQualityValue(readSeq[i]) < bestQuals[j]) continue;
            }
            insert(bestQuals, j+1, getQualityValue(readSeq[i]));
            if (back(bestQuals) > 0 ) appendValue(lowQuals, back(bestQuals));             // move() ?
            eraseBack(bestQuals);
         }
        else
            appendValue(lowQuals, getQualityValue(readSeq[i]));
    }
    // Compute pseudo worst score
    score = 0;
    for (unsigned i = 0; i < length(bestQuals); ++i)
    {
        //long double e = pow(10, -(long double)bestQuals[i]/10.0);
        score += options.gapOpenScore; //std::log10(options.scoreMismatch*(1.0-e) + options.scoreMatch*e);  // Assume mismatch
        //score += std::log10((0.1/12.0)/(1.0/16.0)); //-1.0; //-4.0;  // Assume gap opening (worst score at the moment) TODO
        //std::cout << " mismatch: " << std::log10(10.0*(e/3.0)) << std::endl;
    }
    for (unsigned i = 0; i < length(lowQuals); ++i)
    {
        long double e = pow(10, -(long double)lowQuals[i]/10.0);
        if (e == 1) e = 0.9;    // e.g. in case of Ns: take as low probability into account
        score += std::log10((((options.seqIdentity-options.refNRate)/4.0)/(1.0/16.0))*options.pseudoMatchScale*(1.0 - e) );
    }
    SEQAN_ASSERT_EQ(length(readSeq), length(bestQuals)+length(lowQuals));
    score *= -10;
    //score += length(readSeq);
}
// Use avg. quality, faster to compute and shouldn't make such a big difference, if we just use -1 for mismatches and quality only for matches
template<typename TScore, typename TReadSeq, typename TOptions>
inline void
computePseudoWorstScore2(TScore &score, TReadSeq &readSeq, TOptions &options)
{
    double avgQual = 0;
    for (unsigned i = 0; i < length(readSeq); ++i)
    {
        avgQual += getQualityValue(readSeq[i]);
    }
    avgQual = avgQual/(double)length(readSeq);
    // Compute pseudo worst score
	unsigned max3Errors = (unsigned)(options.max3Error * length(readSeq) / 100.0);
    unsigned pseudoErrors = max3Errors + 1;
    score = 0;
    score += pseudoErrors * options.gapOpenScore; //std::log10((0.01/12.0)/(1.0/16.0)); //options.scoreMismatch; // ??? (-1.0);   // Assumed mismatches
    double e = pow(10, -(long double)avgQual/10.0);
    score += (length(readSeq) - pseudoErrors) * std::log10((((options.seqIdentity-options.refNRate)/4.0)/(1.0/16.0))*options.pseudoMatchScale*(1.0 - e) ); // std::log10(((0.99/4.0)/(1.0/16.0))*(1.0 - e) + (((0.01/12.0)/(1.0/16.0)))*e);   // Assumed matches
    //score -= length(readSeq);
    score *= -10;
}


struct VerifiedRead
{
    typedef FragmentStore<MyFragmentStoreConfig>::TMappingQuality          TScore;
    typedef FragmentStore<MyFragmentStoreConfig>::TAlignedReadStore        TAlignedReadStore;
    typedef Value<TAlignedReadStore>::Type                                 TAlignedReadStoreElement;
    typedef TAlignedReadStoreElement::TID                                  TID;

    TID alignedReadID;
    BamAlignmentRecord record;
    TScore mapq;

	VerifiedRead():
		alignedReadID(0),
		mapq(0) {}

	VerifiedRead(TID _alignedReadID, BamAlignmentRecord _record, TScore _mapq):
		alignedReadID(_alignedReadID),
		record(_record),
		mapq(_mapq) {}

};

template<typename TFragmentStore, typename TOptions>
inline VerifiedRead
verifyRead(TFragmentStore &store, TOptions &options)
{
    typedef typename TFragmentStore::TMappingQuality        TScore;

    // Find best and second best hit (match mate hit) -> lowest score
    TScore bestScore = store.alignQualityStore[0].score;
    unsigned bestID = 0;
    TScore secBestScore = 10000; //store.alignQualityStore[0].score;
    unsigned secBestID = 0;
    for (unsigned i = 1; i < length(store.alignedReadStore); ++i)
    {
        if (bestScore > store.alignQualityStore[i].score)
        {
            secBestScore = bestScore;
            secBestID = bestID;
            bestScore = store.alignQualityStore[i].score;
            bestID = i;
        }
        else if (secBestScore > store.alignQualityStore[i].score)
        {
            secBestScore = store.alignQualityStore[i].score;
            secBestID = i;
        }
    }

    unsigned countHits = length(store.alignedReadStore);
    TScore mapq;
    computeMapq(mapq, bestID, countHits, secBestID, store, options);

    BamAlignmentRecord record;
    record.flag = 0;
    record.rNextID = BamAlignmentRecord::INVALID_REFID;
    record.pNext = BamAlignmentRecord::INVALID_POS;
    record.tLen = 0;

    VerifiedRead verifiedRead(bestID, record, mapq);
    return verifiedRead;
}

template<typename TFragmentStore, typename TID>
inline unsigned
countEqualMateHits(TFragmentStore &store, TID &mateID)
{
    // Count number of hits of one mate which have the same number of errors
    unsigned count = 0;
    unsigned errors = store.alignQualityStore[mateID].errors;
    if (getMateNo(store, store.alignedReadStore[mateID].readID) == 0)
    {
        for (unsigned i = 0; i < length(store.alignQualityStore); ++i)
            if (getMateNo(store, store.alignedReadStore[i].readID) == 0 && store.alignQualityStore[i].errors == errors)
                ++count;
    }
    else if (getMateNo(store, store.alignedReadStore[mateID].readID) == 1)
    {
        for (unsigned i = 0; i < length(store.alignQualityStore); ++i)
            if (getMateNo(store, store.alignedReadStore[i].readID) == 1 && store.alignQualityStore[i].errors == errors)
                ++count;
    }
    return count;
}

template<typename TFragmentStore, typename TID>
inline unsigned
countEqualMateHits(TFragmentStore &store, TID &mateIDL, TID &mateIDR)
{
    typedef typename TFragmentStore::TAlignedReadStore      TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type         TAlignedReadStoreElement;
    // Count number of match mate hits which have the same number of errors
    unsigned count = 0;

    unsigned errors = store.alignQualityStore[mateIDL].errors + store.alignQualityStore[mateIDR].errors;

    for (unsigned i = 0; i < length(store.alignQualityStore); ++i)
        if (getMateNo(store, store.alignedReadStore[i].readID) == 0 && store.alignedReadStore[i].pairMatchID != TAlignedReadStoreElement::INVALID_ID) // Only if left mate and has match mate
            if (store.alignQualityStore[i].errors + store.alignQualityStore[store.alignedReadStore[i].pairMatchID].errors == errors)
                ++count;

    return count;
}


template<typename TScore, typename TID, typename TFragmentStore, typename TOptions>
inline void
computeMapq(TScore &mapq, TID &bestID, unsigned countHits, TID &secBestID, TFragmentStore &store, TOptions &options)
{
    typedef typename TFragmentStore::TReadSeq               TReadSeq;

    unsigned countBest;
    unsigned countSecBest;
    // If only one pair hit: use single hit counts
    // If multiple pair hits: use pair hit counts
    if (!empty(store.matePairStore))
    {
        countBest = countEqualMateHits(store, bestID);
        countSecBest = countEqualMateHits(store, secBestID);
    }
    else
    {
        countBest = countEqualHits(store, bestID);
        countSecBest = countEqualHits(store, secBestID);
    }
    if (countHits == 1)
    {
        TScore worstScore;
        TReadSeq readSeq = store.readSeqStore[store.alignedReadStore[bestID].readID];
        computePseudoWorstScore2(worstScore, readSeq, options);
        mapq = worstScore - store.alignQualityStore[bestID].score;
        if (false) // store.readNameStore[store.alignedReadStore[bestID].readID] == "simulated_1.42" || store.readNameStore[store.alignedReadStore[bestID].readID] == "simulated_1.88" || store.readNameStore[store.alignedReadStore[bestID].readID] == "simulated_1.1042" )
        {
            std::cout << " Only one hit" << std::endl;
            std::cout << "      readName" <<  store.readNameStore[store.alignedReadStore[bestID].readID] << std::endl;
            std::cout << "      ***: score :" << store.alignQualityStore[bestID].score << "  worstScore: " << worstScore << std::endl;
            std::cout << "      ***: mapq :" << mapq << std::endl;
            std::cout << "      ***: errors: " << static_cast<unsigned int>(store.alignQualityStore[bestID].errors) << std::endl;
        }
    }
    else if (countBest != 1) // Same edit distance
    {
        mapq = 0;
    }
    else
    {
        mapq = store.alignQualityStore[secBestID].score  - store.alignQualityStore[bestID].score - 10*std::log10((double)countSecBest);

        //std::cout << " Unique:  read:" << store.readNameStore[store.alignedReadStore[bestID].readID] << std::endl;
        //std::cout << "***: sec best score :" << store.alignQualityStore[secBestID].score  << " best score: " << store.alignQualityStore[bestID].score << "  count2: " << countSecBest << std::endl;
        //std::cout << "***: sec best id :" << secBestID  << " best id: " <<bestID  << std::endl;
        //std::cout << "*** mapq: " << mapq << std::endl;
        //  can be negative, since it's possible that there exist single hits, which are better than left or right mate match hit
    }
}

struct VerifiedMates
{
    typedef FragmentStore<MyFragmentStoreConfig>::TMappingQuality          TScore;
    typedef FragmentStore<MyFragmentStoreConfig>::TAlignedReadStore        TAlignedReadStore;
    typedef Value<TAlignedReadStore>::Type                                 TAlignedReadStoreElement;
    typedef TAlignedReadStoreElement::TID                                  TID;

    TID alignedReadIDL;
    TID alignedReadIDR;
    BamAlignmentRecord recordL;
    BamAlignmentRecord recordR;
    TScore mapqL;
    TScore mapqR;

	VerifiedMates():
		alignedReadIDL(0),
		alignedReadIDR(0),
		mapqL(0),
		mapqR(0) {}

	VerifiedMates(TID _alignedReadIDL,
                  TID _alignedReadIDR,
	              BamAlignmentRecord _recordL,
	              BamAlignmentRecord _recordR,
	              TScore _mapqL,
	              TScore _mapqR):
		alignedReadIDL(_alignedReadIDL),
		alignedReadIDR(_alignedReadIDR),
		recordL(_recordL),
		recordR(_recordR),
		mapqL(_mapqL),
		mapqR(_mapqR) {}
};


template<typename TFragmentStore, typename TOptions>
inline VerifiedMates
verifyMates(TFragmentStore &store, TOptions &options)
{
    typedef typename TFragmentStore::TMatePairStore         TMatePairStore;
    typedef typename Value<TMatePairStore>::Type            TMatePairStoreElement;
    typedef typename TFragmentStore::TAlignedReadStore      TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type         TAlignedReadStoreElement;
    typedef typename TAlignedReadStoreElement::TID          TID;
    typedef typename TFragmentStore::TContigPos             TContigPos;
    typedef typename TFragmentStore::TMappingQuality        TScore;

    // Find best and second best hit (match mate hit) -> lowest score
    TScore bestScoreL = store.alignQualityStore[0].score;
    TID bestIDL = 0;
    TScore secBestScoreL = store.alignQualityStore[0].score;
    TID secBestIDL = 0;
    TScore bestScoreR = store.alignQualityStore[0].score;
    TID bestIDR = 0;
    TScore secBestScoreR = store.alignQualityStore[0].score;
    TID secBestIDR = 0;

    unsigned countPairHits = 0;
    unsigned countHitsL = 0;
    unsigned countHitsR = 0;
    TScore bestMateScore = 1000000;   // Set to maximum
    TID bestMateID = 0;
    TScore secBestMateScore = 1000000;
    TID secBestMateID = 0;
    if (store.matePairStore[0].readID[0] != TMatePairStoreElement::INVALID_ID && store.matePairStore[0].readID[1] != TMatePairStoreElement::INVALID_ID) // If read is paired
    {
        for (TID i = 0; i < length(store.alignedReadStore); ++i)
        {
            // Get single left and right scores
            if (getMateNo(store, store.alignedReadStore[i].readID) == 0 )
            {
                ++countHitsL;
                if (bestScoreL >= store.alignQualityStore[i].score)
                {
                    secBestScoreL = bestScoreL;
                    secBestIDL = bestIDL;
                    bestScoreL = store.alignQualityStore[i].score;
                    bestIDL = i;
                }
                else if (secBestScoreL >= store.alignQualityStore[i].score)
                {
                    secBestScoreL = store.alignQualityStore[i].score;
                    secBestIDL = i;
                }
            }
            else if (getMateNo(store, store.alignedReadStore[i].readID) == 1 )
            {
                ++countHitsR;
                if (bestScoreR >= store.alignQualityStore[i].score)
                {
                    secBestScoreR = bestScoreR;
                    secBestIDR = bestIDR;
                    bestScoreR = store.alignQualityStore[i].score;
                    bestIDR = i;
                }
                else if (secBestScoreR >= store.alignQualityStore[i].score)
                {
                    secBestScoreR = store.alignQualityStore[i].score;
                    secBestIDR = i;
                }
            }
            // Get mate match scores
            if (getMateNo(store, store.alignedReadStore[i].readID) == 0 && store.alignedReadStore[i].pairMatchID != TAlignedReadStoreElement::INVALID_ID ) // Only if left mate and has match mate
            {
                ++countPairHits;

                double currMateScore = store.alignQualityStore[i].score + store.alignQualityStore[store.alignedReadStore[i].pairMatchID].score;
                if (currMateScore <= bestMateScore)
                {
                    secBestMateScore = bestMateScore;
                    secBestMateID = bestMateID;             // ID of left aligned read
                    bestMateScore = currMateScore;
                    bestMateID = i;
                }
                else if (currMateScore <= secBestMateScore)
                {
                    secBestMateScore = currMateScore;
                    secBestMateID = i;
                }
            }
        }
    }
    //std::cerr << "      countPairHits: " << countPairHits << "  countHitsL  " << countHitsL << "  countHitsR  " << countHitsR << "  length alignedReadStore: " << length(store.alignedReadStore) << std::endl;

    // Compute mapping qualities dependent on case
    unsigned count1 = 0;
    if (countPairHits >= 1)
        count1 = countEqualMateHits(store, bestMateID, store.alignedReadStore[bestMateID].pairMatchID);

    TScore mapqL;
    TScore mapqR;
    BamAlignmentRecord recordL;
    BamAlignmentRecord recordR;
    if (countPairHits == 1)  // Case1: Only one pair hit (single hits can be nonunique)
    {
        TID bestMateIDR = store.alignedReadStore[bestMateID].pairMatchID;
        computeMapq(mapqL, bestMateID, countHitsL, secBestIDL, store, options);
        computeMapq(mapqR, bestMateIDR, countHitsR, secBestIDR, store, options);

        mapqL = mapqL + mapqR;
        mapqR = mapqL;

        // Set record entries with matepair information
        recordL.flag = 0;
        recordR.flag = 0;
        recordL.flag |= 0x001;
        recordR.flag |= 0x001;
        recordL.flag |= 0x002;
        recordR.flag |= 0x002;
        recordL.flag |= 0x0040;  // Is left read
        recordR.flag |= 0x0080;  // Is right read
        recordL.rNextID = store.alignedReadStore[bestMateIDR].contigID;
        recordR.rNextID = store.alignedReadStore[bestMateID].contigID;

        if (store.alignedReadStore[bestMateID].beginPos < store.alignedReadStore[bestMateID].endPos)
        {
            recordR.pNext = store.alignedReadStore[bestMateID].beginPos;
        }
        else
        {
            recordR.flag |= 0x0020;
            recordR.pNext = store.alignedReadStore[bestMateID].endPos;
        }
        if (store.alignedReadStore[bestMateIDR].beginPos < store.alignedReadStore[bestMateIDR].endPos)
        {
            recordL.pNext = store.alignedReadStore[bestMateIDR].beginPos;
        }
        else
        {
            recordL.flag |= 0x0020; // Mate is reverse complemented
            recordL.pNext = store.alignedReadStore[bestMateIDR].endPos;
        }
        TContigPos beginPair = store.alignedReadStore[bestMateID].beginPos;
        TContigPos endPair = store.alignedReadStore[bestMateIDR].beginPos;
        if (beginPair < endPair)
        {
            recordL.tLen = endPair - beginPair;
            recordR.tLen = -(endPair - beginPair);
        }
        else
        {
            recordL.tLen = -(beginPair - endPair);
            recordR.tLen = beginPair - endPair;
        }

        //std::cout << "One hit:  mapqL: " << mapqL << "  errors: " << static_cast<unsigned>(store.alignQualityStore[bestMateID].errors) << "  score: " << store.alignQualityStore[bestMateID].score << std::endl;
        //std::cout << "  right:  mapqR: " << mapqR << "  errors: " << static_cast<unsigned>(store.alignQualityStore[bestMateIDR].errors) << "  score: " << store.alignQualityStore[bestMateIDR].score << std::endl;

        VerifiedMates verifiedMates(bestMateID, bestMateIDR, recordL, recordR, mapqL, mapqR);
        return verifiedMates;

    }
    else if (countPairHits > 1)   // Case 2: Multiple mate matches (no summing up of mapqs)
    {
            // Compare best mate match score against second best mate match (not single)
            // -> Ensure that mate match get higher mapping quality than it would get in not paired mode
            // -> Increases prob. that score of each single mate is higher that that of second best mate match
            // (is mate match is classified as best, doesn't mean there are no single reads which are better than single mates)
            TID bestMateIDR = store.alignedReadStore[bestMateID].pairMatchID;
            TID secBestMateIDR = store.alignedReadStore[secBestMateID].pairMatchID;

            computeMapq(mapqL, bestMateID, countHitsL, secBestMateID, store, options);
            computeMapq(mapqR, bestMateIDR, countHitsR, secBestMateIDR, store, options);

            // If pair hit is unique
            if (count1 == 1)
            {
                mapqL = mapqL + mapqR;
                mapqR = mapqL;
            }

            // Set record entries with matepair information
            recordL.flag = 0;
            recordR.flag = 0;
            recordL.flag |= 0x001;
            recordR.flag |= 0x001;
            recordL.flag |= 0x002;
            recordR.flag |= 0x002;
            recordL.flag |= 0x0040;  // Is left read
            recordR.flag |= 0x0080;  // Is right read
            recordL.rNextID = store.alignedReadStore[bestMateIDR].contigID;
            recordR.rNextID = store.alignedReadStore[bestMateID].contigID;

            if (store.alignedReadStore[bestMateID].beginPos < store.alignedReadStore[bestMateID].endPos)
            {
                recordR.pNext = store.alignedReadStore[bestMateID].beginPos;
            }
            else
            {
                recordR.flag |= 0x0020;
                recordR.pNext = store.alignedReadStore[bestMateID].endPos;
            }
            if (store.alignedReadStore[bestMateIDR].beginPos < store.alignedReadStore[bestMateIDR].endPos)
            {
                recordL.pNext = store.alignedReadStore[bestMateIDR].beginPos;
            }
            else
            {
                recordL.flag |= 0x0020; // Mate is reverse complemented
                recordL.pNext = store.alignedReadStore[bestMateIDR].endPos;
            }
            TContigPos beginPair = store.alignedReadStore[bestMateID].beginPos;
            TContigPos endPair = store.alignedReadStore[bestMateIDR].beginPos;
            if (beginPair < endPair)
            {
                recordL.tLen = endPair - beginPair;
                recordR.tLen = -(endPair - beginPair);
            }
            else
            {
                recordL.tLen = -(beginPair - endPair);
                recordR.tLen = beginPair - endPair;
            }

            VerifiedMates verifiedMates(bestMateID, bestMateIDR, recordL, recordR, mapqL, mapqR);
            return verifiedMates;
    }
    // else
    VerifiedMates verifiedMates;
    return verifiedMates;

    /*else if (options.outputSingleMates && countPairHits == 0)   // Case 3: Only single matches // input of single reads by mapper possible?
    {
        if (countHitsL > 0)
        {
            // Compute mapqs
            unsigned count1L = countEqualHits(store, bestIDL);
            unsigned count2L = countEqualHits(store, secBestIDL);
            if (countHitsL == 1)
                mapqL = 255 - bestScoreL; //-10*log10(1-pow(10, -bestScoreL/10.0));
            else if (count1L != 1)
                mapqL = 0;
            else
                mapqL = secBestScoreL - bestScoreL - 10*std::log10(count2L);

            // Set record entries with matepair information
            recordL.flag = 0;
            recordL.flag |= 0x008;   // Mate is unmapped (or not paired with this one)
            recordL.flag |= 0x0040;  // Is left read
            recordL.rNextID = '*';
            recordL.pNext = '*';
            recordL.tLen = 0;

        }
        if (countHitsR > 0)
        {
            // Compute mapqs
            unsigned count1R = countEqualHits(store, bestIDR);
            unsigned count2R = countEqualHits(store, secBestIDR);
            if (countHitsR == 1)
                mapqR = 255 - bestScoreR; //-10*log(1-pow(10, -bestScoreR/10.0));
            else if (count1R != 1)
                mapqR = 0;
            else
                mapqR = secBestScoreR - bestScoreR - 10*std::log10(count2R);

            // Set record entries with matepair information
            recordR.flag = 0;
            recordR.flag |= 0x008;   // Mate is unmapped (or not paired with this one)
            recordR.flag |= 0x0080;  // Is right read
            recordL.rNextID = '*';
            recordL.pNext = '*';
            recordR.tLen = 0;

        }
    }*/
}


template <typename TAlignedRead, typename TMInfo, typename TFragStore, typename TOptions>
inline int
_compareAlignedReadAndMateInfo2(TAlignedRead const &a, TMInfo const &b, TFragStore const &fragStore, TOptions &options)
{
    if ((__int32)a.contigID < b.contigID) return -1;
    if ((__int32)a.contigID > b.contigID) return 1;

    typename TFragStore::TContigPos posA = _min(a.beginPos, a.endPos);
    if (posA < b.beginPos - options.intervalOffset)
    {
        return -1;
    }
    if (posA > b.beginPos + options.intervalOffset) return 1;

    bool reversedA = (a.beginPos > a.endPos);
    if (!reversedA && b.reversed)
    {
        return -1;
    }
    if (reversedA && !b.reversed) return 1;

    typedef typename TFragStore::TMatePairStore     TMatePairStore;
    typedef typename Value<TMatePairStore>::Type    TMatePair;
    typename TMatePair::TID matePairIDA = TMatePair::INVALID_ID;

    if (a.readID < length(fragStore.readStore))
        matePairIDA = fragStore.readStore[a.readID].matePairID;

    if (matePairIDA < b.matePairID)
    {
        return -1;
    }
    if (matePairIDA > b.matePairID) return 1;
    return 0;
}

// Use pairMatchIDs to store alignID of other match mate, which can still be used to access this entry in alignedReadStore (if order not touched since building!)
template<typename TFragmentStore, typename TMatchMateInfos, typename TOptions>
inline void
_generatePseudoPairMatchIDs (
    TFragmentStore &store,
    TMatchMateInfos & matchMateInfos,
    TOptions &options)
{
    typedef typename TFragmentStore::TAlignedReadStore              TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type                 TAlignedRead;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type    TIter;
    typedef typename Iterator<TMatchMateInfos, Standard>::Type      TMIter;

    TIter it = begin(store.alignedReadStore, Standard());
    TIter itEnd = end(store.alignedReadStore, Standard());
    TMIter mit = begin(matchMateInfos, Standard());
    TMIter mitEnd = end(matchMateInfos, Standard());

    if (it == itEnd || mit == mitEnd) return;

    // sort the aligned read store by: begin position, contig name
    //std::sort(it,  itEnd,  AlignedMateLess_<TFragmentStore>(fragStore)); // Should be sorted by position for one read anyway
    std::sort(mit, mitEnd, MatchMateInfoLess_());

    while (true)
    {
        // skip already aligned reads
        while (it->pairMatchID != TAlignedRead::INVALID_ID)
            if (++it == itEnd) return;

        int cmp = _compareAlignedReadAndMateInfo2(*it, *mit, store, options);
        if (cmp == 0)   // both are equal -> link them
        {
            (*it).pairMatchID = (*mit).pairMatchID;
            store.alignedReadStore[(*it).pairMatchID].pairMatchID = (*it).id; // Only possible if order of alignedReadStore not touched, alignID == position
        }
        if (cmp >= 0)   // MateInfo is less or equal
        {
            if (++mit == mitEnd) return;
        }
        if (cmp <= 0)   // AlignedRead is less or equal
        {
            if (++it == itEnd) return;
        }
    }
}


template<typename TOptions, typename TModel>
inline bool
postProcessMain(TOptions &options, TModel const &)
{
    typedef FragmentStore<MyFragmentStoreConfig> TFragmentStore;

    typedef typename TFragmentStore::TReadNameStore             TReadNameStore;
    typedef NameStoreCache<TReadNameStore, CharString>          TReadNameStoreCache;

    typedef typename TFragmentStore::TContigStore                                   TContigStore;
    typedef typename Value<TContigStore>::Type                                      TContig;
    typedef typename TContig::TContigSeq                                            TContigSeq;
    typedef typename TFragmentStore::TReadSeq                                       TReadSeq;
    typedef typename TFragmentStore::TAlignedReadStore                              TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type                                 TAlignedRead;

    typedef typename TAlignedRead::TGapAnchors                                      TReadGapAnchors;
    typedef Gaps<TReadSeq, AnchorGaps<TReadGapAnchors> >                            TReadGaps;
    typedef typename TContig::TGapAnchors                                           TContigGapAnchors;
    typedef Gaps<TContigSeq, AnchorGaps<TContigGapAnchors> >                        TContigGaps;
    typedef typename TFragmentStore::TContigPos                                     TContigPos;

    typedef typename TAlignedRead::TID                                              TID;

    typedef StringSet<TContigGaps>                                                  TSetContigGaps;

    typedef typename Value<typename TFragmentStore::TMatePairStore>::Type           TMatePairStoreElement;
    typedef MatchMateInfo_<TID>                                                     TMatchMateInfo;
    typedef String<TMatchMateInfo>                                                  TMatchMateInfos;

    // Initialize aligment scores
    typedef double  TValue;
    typedef Score<int, BsTagList<BsCaseCT, TModel, Left> >           TBsScoreCTLeft;
    typedef Score<int, BsTagList<BsCaseCT, TModel, Right> >          TBsScoreCTRight;
    typedef Score<int, BsTagList<BsCaseGA, TModel, Left> >           TBsScoreGALeft;
    typedef Score<int, BsTagList<BsCaseGA, TModel, Right> >          TBsScoreGARight;

    BsSubstitutionMatrix<TValue, BsCaseCT, BsSimple> bsSubstitutionMatrixCT(options.globalMethRate, options.bsConversionRate, options.seqIdentity, options.refNRate);
    BsSubstitutionMatrix<TValue, BsCaseGA, BsSimple> bsSubstitutionMatrixGA(options.globalMethRate, options.bsConversionRate, options.seqIdentity, options.refNRate);
    TValue const * seqErrorFreqs;
    TValue const * insErrorFreqs;
    TValue const * delErrorFreqs;
    if (options.nonSimpleSubstErrors) seqErrorFreqs = SeqErrorFreqs<TValue, BsNonSimple>::getData();
    else seqErrorFreqs = SeqErrorFreqs<TValue, BsSimple>::getData();
    if (options.nonSimpleInsErrors) insErrorFreqs = InsErrorFreqs<TValue, BsNonSimple>::getData();
    else insErrorFreqs = InsErrorFreqs<TValue, BsSimple>::getData();
    if (options.nonSimpleDelErrors)
    {
        delErrorFreqs = DelErrorFreqs<TValue, BsNonSimple>::getData();
        options.scalingFactorDelErrors = 3.5;       // Compute automatically...
    }
    else delErrorFreqs = DelErrorFreqs<TValue, BsSimple>::getData();

    TBsScoreCTLeft  scoringSchemeCTLeft(options, bsSubstitutionMatrixCT, seqErrorFreqs, insErrorFreqs, delErrorFreqs);
    TBsScoreCTRight scoringSchemeCTRight(options, bsSubstitutionMatrixCT, seqErrorFreqs, insErrorFreqs, delErrorFreqs);
    TBsScoreGALeft  scoringSchemeGALeft(options, bsSubstitutionMatrixGA, seqErrorFreqs, insErrorFreqs, delErrorFreqs);
    TBsScoreGARight scoringSchemeGARight(options, bsSubstitutionMatrixGA, seqErrorFreqs, insErrorFreqs, delErrorFreqs);
    //std::cout << "Computation of alignment scores finished." << std::endl;

    // Load original reads into fragmentStore
    TFragmentStore store;
    if (options.readFileName2 == "")
        loadReadsCroppedID(store, options.readFileName);
    else
        loadReadsCroppedID(store, options.readFileName, options.readFileName2);

    TReadNameStoreCache readNameCache(store.readNameStore);
    refresh(store.readNameStoreCache);

    // Load reference into fragmentStore
    loadContigs(store, options.refFileName);
    refresh(store.contigNameStoreCache);

    // Parse SAM file ....
    seqan::BamFileOut bamFileOut;
    if (!open(bamFileOut, toCString(options.outputFileName)))
    {
        std::cerr << "ERROR: Could not open " << toCString(options.outputFileName) << " for writing.\n";
        return 1;
    }
    seqan::BamFileIn bamFileIn(bamFileOut);  // make sequence names etc. known to bamFileOut.
    if (!open(bamFileIn, toCString(options.samFileName)))
    {
        std::cerr << "ERROR: Could not open " << toCString(options.samFileName) << " for reading.\n";
        return 1;
    }

    BamHeader header;
    BamAlignmentRecord record;
    // Read header
    readHeader(header, bamFileIn);
    // Check that this fits to the contig names in the store.
    if (length(contigNames(context(bamFileIn))) != length(store.contigNameStore))
    {
        std::cerr << "ERROR: Different number of contigs in reference and SAM/BAM file!\n";
        return 1;
    }
    for (unsigned i = 0; i < length(contigNames(context(bamFileIn))); ++i)
        if (contigNames(context(bamFileIn))[i] != store.contigNameStore[i])
        {
            std::cerr << "ERROR: Different contig names/order in reference and SAM/BAM file!\n";
            return 1;
        }
    // Write out header again. Maybe add information about BS mapping ?
    writeHeader(bamFileOut, header);

    // Parse SAM file, verify and write to output
    // data structure to temporarily store information about match mates
    TMatchMateInfos matchMateInfos;

    CharString currReadName;
    TID readID;
    //TSetContigGapAnchors setContigGapAnchors;   // TODO only store anchorGaps and assign contig seq later for output
    TSetContigGaps setContigGaps;

    int helper = 0;
    while (!atEnd(bamFileIn))
    {
        ++helper;
        readRecord(record, bamFileIn);
        // TODO generalize
        if ((record.qName[length(record.qName)-2] == '/' && (!empty(currReadName) && prefix(record.qName, length(record.qName)-2) != prefix(currReadName, length(currReadName)-2)  )) ||
            (record.qName[length(record.qName)-2] != '/' &&  record.qName != currReadName) ||
            (empty(currReadName) && record.qName != currReadName) )   // First entry for read: Verify previous read and read in new read
        {
            if (!empty(store.alignedReadStore))   // After all entries for one read are loaded
            {
                if (!empty(store.matePairStore))
                {
                    // set the match mate IDs using the information stored in matchMateInfos
                    _generatePseudoPairMatchIDs(store, matchMateInfos, options);

                    VerifiedMates verifiedMates = verifyMates(store, options);  // For mates: Find best alignments and compute mapq, verify
                    //std::cout << "alignedReadIDL: " << verifiedMates.alignedReadIDL << "  alignedReadIDR: " << verifiedMates.alignedReadIDR << "  length(alignedReadStore) " <<  length(store.alignedReadStore) << std::endl;
                    //std::cout << "Contig gaps: "  << std::endl;
                    //std::cout << setContigGaps[verifiedMates.alignedReadIDL] << std::endl;
					unsigned max4Errors = (unsigned)(options.max4Error*length(store.readSeqStore[store.alignedReadStore[verifiedMates.alignedReadIDL].readID]) / 100.0);

                    if (verifiedMates.mapqL >= options.minMapq && store.alignQualityStore[verifiedMates.alignedReadIDL].errors <= max4Errors && store.alignQualityStore[verifiedMates.alignedReadIDL].score <= options.maxScore)
                        writeBsAlignment(bamFileOut, store, setContigGaps[verifiedMates.alignedReadIDL], verifiedMates.alignedReadIDL, verifiedMates.mapqL, verifiedMates.recordL, options);
					max4Errors = (unsigned)(options.max4Error*length(store.readSeqStore[store.alignedReadStore[verifiedMates.alignedReadIDR].readID]) / 100.0);
                    if (verifiedMates.mapqR >= options.minMapq && store.alignQualityStore[verifiedMates.alignedReadIDR].errors <= max4Errors && store.alignQualityStore[verifiedMates.alignedReadIDR].score <= options.maxScore)
                        writeBsAlignment(bamFileOut, store, setContigGaps[verifiedMates.alignedReadIDR], verifiedMates.alignedReadIDR, verifiedMates.mapqR, verifiedMates.recordR, options);
                }
                else
                {
                    VerifiedRead verifiedRead = verifyRead(store, options);   // Find best alignment and compute mapq, verify
					unsigned max4Errors = (unsigned)(options.max4Error*length(store.readSeqStore[store.alignedReadStore[verifiedRead.alignedReadID].readID]) / 100.0);
                    if (verifiedRead.mapq >= options.minMapq && static_cast<unsigned int>(store.alignQualityStore[verifiedRead.alignedReadID].errors) <= max4Errors && store.alignQualityStore[verifiedRead.alignedReadID].score <= options.maxScore)
                    {
                        writeBsAlignment(bamFileOut, store, setContigGaps[verifiedRead.alignedReadID], verifiedRead.alignedReadID, verifiedRead.mapq, verifiedRead.record, options);
                    }
               }
            }
            ////////////////////////////////////////////////////////////////////////////////////
            // Read entry for next read
            clear(store.alignedReadStore);
            clear(store.alignQualityStore);
            clear(matchMateInfos);
            clear(setContigGaps);

            currReadName = record.qName;
        }

        if (hasFlagUnmapped(record)) continue;     // Read is unmapped
        // Get readID (could be curr. read or mate) -> Only ID, seq. we will get the original from readSeqStore
        // If read name not found, skip entry
        if (!getIDByName(store.readNameStore, record.qName, readID, readNameCache)) continue;

        if (hasFlagMultiple(record))   //)    // If paired: Get readID for current mate
        {
            // if the read is in the store and paired
            // check the mate pair store if it is the same mate of the pair
            // assuming that only one flag 0x040 or 0x0080 is 1
            int inPair = 1 - ((record.flag & 0x40) >> 6);	// bit 7 is set => inPair = 0
                                                             // else inPair = 1 (even if bits 6 and 7 are not set)
            TID matePairID = store.readStore[readID].matePairID;
            if (matePairID != TMatePairStoreElement::INVALID_ID)
            {
                readID = store.matePairStore[matePairID].readID[inPair];
                if (readID == TMatePairStoreElement::INVALID_ID) continue;
            }
        }
        // Get positions
        unsigned len = 0;
        _getLengthInRef(len, record.cigar);

        TContigPos beginPos = record.beginPos;
        TContigPos endPos = beginPos + len;
        if (hasFlagRC(record)) // Reverse
        {
            TContigPos temp = beginPos;
            beginPos = endPos;
            endPos = temp;
        }
        // Create alignedReadStore entry
        TID id = length(store.alignedReadStore);
        TReadGapAnchors readGapAnchors;
        TAlignedRead alignedRead = TAlignedRead(id, readID, record.rID, beginPos, endPos, readGapAnchors);
        appendValue(store.alignedReadStore, alignedRead);

        TReadSeq readSeq =  store.readSeqStore[readID];
        TContigSeq contigInf;
        TReadGaps readGaps;
        // Set readGaps source and get contig infix
        TContigPos beginInf;
        TContigPos endInf;
        if (beginPos < endPos)
        {
            setSource(readGaps, readSeq);
            if (beginPos < options.intervalOffset) beginInf = 0;
            else beginInf = beginPos-options.intervalOffset;
            if (endPos > (long)length(store.contigStore[record.rID].seq) - options.intervalOffset) endInf = length(store.contigStore[record.rID].seq);
            else endInf = endPos+options.intervalOffset;
        }
        else
        {
            TReadSeq readSeq = store.readSeqStore[readID];
            reverseComplement(readSeq);
            for (int i = length(store.readSeqStore[readID]) - 1; i >= 0; --i)   // assign rev. compl qualities
            {
                assignQualityValue(readSeq[i], getQualityValue(store.readSeqStore[readID][length(store.readSeqStore[readID]) - 1 - i]));
            }
            assignSource(readGaps, readSeq);
            if (endPos < options.intervalOffset) beginInf = 0;
            else beginInf = endPos-options.intervalOffset;
            if (beginPos > (long)length(store.contigStore[record.rID].seq) - options.intervalOffset) endInf = length(store.contigStore[record.rID].seq);
            else endInf = beginPos+options.intervalOffset;
        }
        // Create contigGaps
        TContigGaps contigGaps;
        contigInf = infix(store.contigStore[record.rID].seq, beginInf, endInf);
        assignSource(contigGaps, contigInf);
        // Realign and set alignedReadStore entries
        reAlign4(readGaps, contigGaps, store, id, scoringSchemeCTLeft, scoringSchemeCTRight, scoringSchemeGALeft, scoringSchemeGARight, options);
        appendValue(setContigGaps, contigGaps, Generous());

        // Store information about mate
        if (getMateNo(store, readID) == 0 && !hasFlagNextUnmapped(record) )  // store info only if read is the first mate and if the second mate is mapped
        {
            typename TMatePairStoreElement::TID matePairID = store.readStore[readID].matePairID;

            TMatchMateInfo matchMateInfo;
            matchMateInfo.readID = readID;
            matchMateInfo.contigID = (TID)record.rID;
            matchMateInfo.pairMatchID = id;
            matchMateInfo.matePairID = matePairID;
            matchMateInfo.reversed = hasFlagNextRC(record);
            matchMateInfo.beginPos = record.pNext;
            appendValue(matchMateInfos, matchMateInfo);
            back(store.alignedReadStore).pairMatchID = id;  // pairMatchID == alignedRead id of first mate // abuse
        }
    }
    if(!empty(store.alignedReadStore))
    {
        // Deal with last read ...
        if (!empty(store.matePairStore))
        {
            // set the match mate IDs using the information stored in matchMateInfos
            _generatePseudoPairMatchIDs(store, matchMateInfos, options);
            VerifiedMates verifiedMates = verifyMates(store, options);  // For mates: Find best alignments and compute mapq, verify
            unsigned max4Errors = (unsigned)(options.max4Error*length(store.readSeqStore[store.alignedReadStore[verifiedMates.alignedReadIDL].readID])/100.0);
            if (verifiedMates.mapqL >= options.minMapq && store.alignQualityStore[verifiedMates.alignedReadIDL].errors <= max4Errors && store.alignQualityStore[verifiedMates.alignedReadIDL].score <= options.maxScore)
                writeBsAlignment(bamFileOut, store, setContigGaps[verifiedMates.alignedReadIDL], verifiedMates.alignedReadIDL, verifiedMates.mapqL, verifiedMates.recordL, options);
			max4Errors = (unsigned)(options.max4Error*length(store.readSeqStore[store.alignedReadStore[verifiedMates.alignedReadIDR].readID]) / 100.0);
            if (verifiedMates.mapqR >= options.minMapq && store.alignQualityStore[verifiedMates.alignedReadIDR].errors <= max4Errors && store.alignQualityStore[verifiedMates.alignedReadIDR].score <= options.maxScore)
                writeBsAlignment(bamFileOut, store, setContigGaps[verifiedMates.alignedReadIDR], verifiedMates.alignedReadIDR, verifiedMates.mapqR, verifiedMates.recordR, options);

        }
        else
        {
            VerifiedRead verifiedRead = verifyRead(store, options);   // Find best alignment and compute mapq, verify
			unsigned max4Errors = (unsigned)(options.max4Error*length(store.readSeqStore[store.alignedReadStore[verifiedRead.alignedReadID].readID]) / 100.0);
            if (verifiedRead.mapq >= options.minMapq && store.alignQualityStore[verifiedRead.alignedReadID].errors <= max4Errors && store.alignQualityStore[verifiedRead.alignedReadID].score <= options.maxScore)
                writeBsAlignment(bamFileOut, store, setContigGaps[verifiedRead.alignedReadID], verifiedRead.alignedReadID, verifiedRead.mapq, verifiedRead.record, options);
        }
    }
    return 0;
}

#endif

