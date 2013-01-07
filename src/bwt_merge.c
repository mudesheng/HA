#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include "bwt_gen64.h"
#include "utils.h"

/*int compareIndex(const void *p, const void *q)
{
	uint64_t *p1 = (uint64_t*)p, *q1 = (uint64_t*)q ;
	if(*p1 < *q1) return -1 ;
	else if(*p1 > *q1) return 1 ;
	else return 0 ;
} */
int64_t FindMostLargeOfLower(const uint64_t *array, const int64_t len, const uint64_t query)
{
	int64_t low = 0, high = len-1, mid ;
	int64_t base = 0;
	while(1)
	{
		if(low >= high)
		{
			if(query >= array[high]) { base = array[high]; break; } 
			else if(query >= array[high -1]) { base = array[high -1] ; break; }
			else { printf("[FindMostLargeOfLower] algorithm is error....exit\n");  exit(1); }
		}
		mid = (high + low) / 2 ;
		if(query >= array[mid]) low = mid + 1 ; 
		else high = mid - 1 ;
	}	

	return base ;
}

void DecodeSuperBlockCompressedBWT(const BWT64 *bwt64, const int64_t sBIndex, uint32_t *buf)
{
	// check arguments
	if(sizeof(buf) < SUPERBLOCK_SIZE/CHAR_PER_WORD)
	{
		printf("[DecodeCompressedBWT]sizeof(buf) < SUPERBLOCK_SIZE/CHAR_PER_WORD: %lu < %d\n", sizeof(buf), SUPERBLOCK_SIZE/CHAR_PER_WORD);
	}
	
	int index = 0;
	int64_t end = 0;	
	if(sBIndex >= bwt64->textLength / SUPERBLOCK_SIZE){ end = bwt64->compressed_size; }
	else end = bwt64->superBlockIndex[sBIndex + 1];
	for(int64_t i = bwt64->superBlockIndex[sBIndex] + (SUPERBLOCK_SIZE/BLOCK_SIZE) * 2; i < end; i++)
	{
		uint8_t count = bwt64->bwtCode[i] | 0x3f ;
		uint8_t c = bwt64->bwtCode[i] >> 6 | 0x3 ;
		for(int j = 0; j < count ; j++)
		{
			buf[index>>4] <<= 2 ;
			buf[index>>4] |= c  ;
			index++;
		}	
	}
	if(index % CHAR_PER_WORD != 0)
	{
		buf[index>>4] <<= ((CHAR_PER_WORD - index % CHAR_PER_WORD )<<1) ;
	}
}

BWT64 *BWT2BWT64(BWT *bwt, BWTInc64 *bwtInc64)
{
    BWT64 *bwt64 ;

    bwt64 = (BWT64*)calloc(1, sizeof(BWT64));
    
    bwt64->textLength = bwt->textLength ;
    bwt64->inverseSa0 = bwt->inverseSa0 ;
    bwt64->cumulativeFreq = (uint64_t*)calloc(ALPHABET_SIZE + 1, sizeof(uint64_t));
    for(int i = 0; i < ALPHABET_SIZE + 1; i++ )
    {    bwt64->cumulativeFreq[i] = bwt->cumulativeFreq[i]; }

    bwt64->decodeTable = (unsigned*)calloc(DNA_OCC_CNT_TABLE_SIZE_IN_WORD, sizeof(unsigned int));
    memcpy(bwt64->decodeTable, bwt->decodeTable, DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned int));
    
    bwt64->bwtSizeInWord = bwt->bwtSizeInWord ;
    bwt64->occMajorSizeInWord = bwt->occMajorSizeInWord ;
    bwt64->occSizeInWord = bwt->occSizeInWord ;
    bwtInc64->availableWord = (uint64_t)((bwt64->bwtSizeInWord + bwt64->occMajorSizeInWord * 2 + bwt64->occSizeInWord) * 1.2);
    bwtInc64->workingMemory = (unsigned*)calloc(bwtInc64->availableWord, BYTES_IN_WORD);
    bwt64->occValueMajor = (uint64_t*)bwtInc64->workingMemory ;
    bwt64->occValue = (uint32_t*)(bwtInc64->workingMemory + bwt64->occMajorSizeInWord * 2 + OCC_INTERVAL);
    bwt64->bwtCode = (uint32_t*)(bwtInc64->workingMemory + bwtInc64->availableWord - bwt64->bwtSizeInWord - BWT_ADDR_OFFSET );
    memcpy(bwt64->bwtCode, bwt->bwtCode, bwt64->bwtSizeInWord * sizeof(uint32_t));

    for(int64_t i = 0 ; i < bwt64->occMajorSizeInWord ; i++)
    { bwt64->occValueMajor[i] = bwt->occValueMajor[i]; }

    memcpy(bwt64->occValue, bwt->occValue, bwt64->occSizeInWord * sizeof(uint32_t));
    return bwt64 ;
}

BWTInc64 *BWTInc2BWTInc64(BWTInc *bwtInc)
{
    BWTInc64 *bwtInc64 ;
    uint32_t i ;

    bwtInc64 = (BWTInc64*)calloc(1, sizeof(BWTInc64));
    bwtInc64->numberOfIterationDone = 1 ;
    bwtInc64->firstCharInLastIteration = bwtInc->firstCharInLastIteration ;
    bwtInc64->bwt = BWT2BWT64(bwtInc->bwt, bwtInc64);
    /*// cut down firstCharInLastIteration in the cumulativeFreq
    bwtInc64->bwt->cumulativeFreq[1] -= (bwtInc64->firstCharInLastIteration <= 0);
    bwtInc64->bwt->cumulativeFreq[2] -= (bwtInc64->firstCharInLastIteration <= 1);
    bwtInc64->bwt->cumulativeFreq[3] -= (bwtInc64->firstCharInLastIteration <= 2);
    bwtInc64->bwt->cumulativeFreq[4] -= (bwtInc64->firstCharInLastIteration <= 3);*/

    bwtInc64->targetNBit = bwtInc->targetNBit ;
    bwtInc64->packedShift = (unsigned*)calloc(CHAR_PER_WORD, sizeof(unsigned int));
    for(i = 0; i < CHAR_PER_WORD ; i++)
    { bwtInc64->packedShift[i] = BITS_IN_WORD - (i+1)* BIT_PER_CHAR; }

    return bwtInc64 ;    
}

void BWTIncSetMergeSize64(BWTInc64 *bwtInc64, BWTInc * bwtInc)
{
    int64_t minWordSize ;
    minWordSize = bwtInc->bwt->textLength * (sizeof(uint64_t)/BYTES_IN_WORD) + 
        bwtInc64->bwt->bwtSizeInWord + bwtInc->bwt->bwtSizeInWord + (bwtInc64->bwt->occMajorSizeInWord + bwtInc->bwt->occMajorSizeInWord) * 2 + bwtInc64->bwt->occSizeInWord + bwtInc->bwt->occSizeInWord + 5 * OCC_INTERVAL_MAJOR;
    if(minWordSize > bwtInc64->availableWord)
    {
        bwtInc64->workingMemory = (unsigned*)realloc(bwtInc64->workingMemory, minWordSize * BYTES_IN_WORD);
        bwtInc64->bwt->occValueMajor = (uint64_t*)bwtInc64->workingMemory ;
        bwtInc64->bwt->occValue = (uint32_t*)(bwtInc64->workingMemory + bwtInc64->bwt->occMajorSizeInWord * 2 + OCC_INTERVAL);
        bwtInc64->bwt->bwtCode = (uint32_t*)(bwtInc64->workingMemory + bwtInc64->availableWord - bwtInc64->bwt->bwtSizeInWord - BWT_ADDR_OFFSET);
        bwtInc64->availableWord = minWordSize ;
    }
    bwtInc64->bwt->bwtCode =  (uint32_t*)memmove(bwtInc64->workingMemory + bwtInc64->availableWord - bwtInc64->bwt->bwtSizeInWord , bwtInc64->bwt->bwtCode, bwtInc64->bwt->bwtSizeInWord * sizeof(uint32_t));
    bwtInc64->bwt->occValue = (uint32_t*)memmove(bwtInc64->workingMemory + bwtInc->bwt->textLength * (sizeof(uint64_t) /BYTES_IN_WORD) + bwtInc64->bwt->occMajorSizeInWord * 2 + OCC_INTERVAL , bwtInc64->bwt->occValue, bwtInc64->bwt->occSizeInWord * BYTES_IN_WORD);
    bwtInc64->bwt->occValueMajor = (uint64_t*)memmove(bwtInc64->workingMemory + bwtInc->bwt->textLength * (sizeof(uint64_t)/BYTES_IN_WORD), bwtInc64->bwt->occValueMajor, bwtInc64->bwt->occMajorSizeInWord * sizeof(uint64_t));
}

uint64_t buildRelativeRank64(uint64_t *relativeRank, const BWT64 *bwt64, const BWT *bwt, const unsigned int *packedShift, const int64_t start, const int64_t end, const int64_t bwt_sa, const int64_t bwt64_sa)
{
    int64_t i , saIndex, rIndex ;
    int64_t last_rindex;
    const int divideCount = bwt64->divideNumber ;
    //uint32_t  s ;
    int64_t bp  = bwt_sa > bwt->inverseSa0 ? bwt_sa -1 : bwt_sa ;
    uint32_t c = (bwt->bwtCode[bp>>4] >> packedShift[bp & 0xf]) & 0x3 ; // the last char of partition in bwt code 
    saIndex = bwt64_sa; // the last char of divide partition order by divideCount 
    relativeRank[bp] = bwt64_sa ; // relative rank of last char in the oldBWT
    rIndex = bwt->cumulativeFreq[c] + BWTOccValue(bwt, bwt_sa, c) + 1; // the rIndex of insert BWT code 
    for(i = 0 ; i < (end - start) - 1 ; i++)
    {
        saIndex = bwt64->cumulativeFreq[c] + BWTOccValue64(bwt64, saIndex, c) + divideCount ; // +divideCount for special $ in the end string 
        if(rIndex > bwt->inverseSa0)
        {
            relativeRank[rIndex - 1] = saIndex ;
            c = (bwt->bwtCode[(rIndex - 1)>>4] >> packedShift[(rIndex - 1) & 0xf]) & 0x3 ;
        } else {
            relativeRank[rIndex] = saIndex ;
            c = (bwt->bwtCode[rIndex>>4] >> packedShift[rIndex & 0xf]) & 0x3 ;
        }
        last_rindex = rIndex;
        rIndex = bwt->cumulativeFreq[c] + BWTOccValue(bwt,rIndex, c) + 1 ;
        if(i > end - start - 100) fprintf(stderr, "%u", c);
    }
    fprintf(stderr, "\n");
#ifdef DEBUG
    fprintf(stderr, "[buildRelativeRank64] rIndex : %ld\n", rIndex);
    if(start == 0 && rIndex != bwt->inverseSa0) 
    {
        fprintf(stderr, "[buildRelativeRank64]rIndex != bwt->inverseSa0\n");
        assert(rIndex == bwt->inverseSa0);
    }
#endif
    saIndex = bwt64->cumulativeFreq[c] + BWTOccValue64(bwt64,saIndex,c) + divideCount ;
    fprintf(stderr, "[buildRelativeRank64] saIndex : %ld\n", saIndex);

    return saIndex + rIndex ;
}

int64_t buildRelativeRankFromOriginPac(const char *fileName, uint64_t *relativeRank, const BWT64 *bwt64, const BWT *bwt, const int64_t start, const int64_t end, const int64_t bwt64_sa)
{
    FILE *fp = xopen(fileName, "rb");
    int buf_size = 1000000;
    uint8_t *buf = (uint8_t*)xcalloc(buf_size, sizeof(uint8_t));
    int64_t cor_pos = (end + CHAR_PER_BYTE -1)/CHAR_PER_BYTE ;
    int last_bufSize = 0;
    if(cor_pos < buf_size) { fprintf(stderr, "[buildRelativeRankFromOriginPac] cor_pos < buf_size, exit..\n"); exit(1); }
    fseek(fp, cor_pos, SEEK_SET);
    bgint_t saIndex = bwt64_sa;
    uint32_t c ;
    for(int64_t i = end - 1; i >= start; i--)
    {
        relativeRank[i] = saIndex;
        if(i < cor_pos * CHAR_PER_BYTE)
        {
            cor_pos -= buf_size;
            if(cor_pos < 0) { cor_pos = 0; fseek(fp, 0, SEEK_SET);}
            else fseek(fp, -(buf_size + last_bufSize), SEEK_CUR);
            fread(buf, sizeof(uint8_t), buf_size, fp);
            last_bufSize = buf_size;
        }
        c = GET_CHAR_FROM_PACKED(buf, i - cor_pos * CHAR_PER_BYTE);
        saIndex = bwt64->cumulativeFreq[c] + BWTOccValue64(bwt64, saIndex, c) + bwt64->divideNumber;
    }

    // sort partition of relativeRank
    BWTQSort64(relativeRank + start, end - start);
    // free and clean work
    free(buf);
    fclose(fp);
    
    return saIndex + bwt->inverseSa0;
}

uint32_t *mergeBWTInc64_1209(const uint64_t* __restrict relativeRank, BWTInc64 *bwtInc64, const BWTInc *bwtInc)
{
    uint32_t *mergedBWT ;
    int64_t oIndex = 0, mIndex = 0 ;
    const int divideCount = bwtInc64->bwt->divideNumber ;
    int64_t oldSALen = bwtInc64->bwt->textLength  + divideCount, insertLen = bwtInc->bwt->textLength , oldSAIndex = 0 ;
    int sentinelCount = 0 ;
    uint32_t *packedShift = bwtInc64->packedShift ;
    uint32_t *buf = xcalloc(BWT_BUF_SIZE, sizeof(uint32_t));
    int bcount = 0 ;

    // set mergedBWT 
    bwtInc64->bwt->bwtSizeInWord = BWTResidentSizeInWord64(bwtInc64->bwt->textLength + bwtInc->bwt->textLength);
    //mergedBWT = bwtInc64->workingMemory + bwtInc64->availableWord - bwtInc64->bwt->bwtSizeInWord - BWT_ADDR_OFFSET ;
    mergedBWT = xcalloc(bwtInc64->bwt->bwtSizeInWord, sizeof(uint32_t));
    for(int64_t i = 0 ; i < insertLen ; i++)
    {
        while(oldSAIndex < oldSALen && oldSAIndex < relativeRank[i])
        {
            if(sentinelCount < divideCount && oldSAIndex == bwtInc64->bwt->sentinelPosition[sentinelCount])
            {
                if(mIndex + sentinelCount >= bwtInc64->bwt->inverseSa0)
                {
                    bwtInc64->bwt->sentinelPosition[sentinelCount] = mIndex + sentinelCount + 1 ;
                }else bwtInc64->bwt->sentinelPosition[sentinelCount] = mIndex + sentinelCount ;
                sentinelCount++ ;
            } else {
                if(bcount >= (BWT_BUF_SIZE<<4)) 
                {
                    memcpy(mergedBWT + ((mIndex-bcount)>>4), buf, BWT_BUF_SIZE * sizeof(uint32_t));
                    memset(buf , 0 , BWT_BUF_SIZE * sizeof(uint32_t));
                    bcount = 0 ;
                }
                buf[bcount>>4] |= (((bwtInc64->bwt->bwtCode[oIndex>>4] >> packedShift[oIndex & 0xf]) & 0x3) << packedShift[mIndex & 0xf]) ;
                mIndex++ ; oIndex++ ; bcount++;
            }
            oldSAIndex++;
        }
        // check
        if(i == bwtInc->bwt->inverseSa0)
        {
            if(mIndex + sentinelCount != bwtInc64->bwt->inverseSa0)
            {
                fprintf(stderr, "[mergeBWT64] i(%ld) + mBWTIndex(%ld) + sentinelCount(%d) != bwtInc64->bwt->inverseSa0(%ld), exit..\n", i, mIndex, sentinelCount, bwtInc64->bwt->inverseSa0);
                //exit(1);
            }
        }
        if(bcount >= (BWT_BUF_SIZE<<4)) 
        {
            memcpy(mergedBWT + ((mIndex-bcount)>>4), buf, BWT_BUF_SIZE * sizeof(uint32_t));
            memset(buf , 0 , BWT_BUF_SIZE * sizeof(uint32_t));
            bcount = 0 ;
        }
        buf[bcount>>4] |= (((bwtInc->bwt->bwtCode[i>>4] >> packedShift[i & 0xf]) & 0x3) << packedShift[mIndex & 0xf]);
        mIndex++ ; bcount++;
    }
    // copy unfinished the oldBWT 
    while(oldSAIndex < oldSALen)
    {
        if(sentinelCount < divideCount && oldSAIndex == bwtInc64->bwt->sentinelPosition[sentinelCount])
        {
            if(mIndex + sentinelCount >= bwtInc64->bwt->inverseSa0)
            {
                bwtInc64->bwt->sentinelPosition[sentinelCount] = mIndex + sentinelCount + 1 ;
            }else bwtInc64->bwt->sentinelPosition[sentinelCount] = mIndex + sentinelCount ;
            sentinelCount++ ;
        }else {
            if(bcount >= (BWT_BUF_SIZE<<4)) 
            {
                memcpy(mergedBWT + ((mIndex-bcount)>>4), buf, BWT_BUF_SIZE * sizeof(uint32_t));
                memset(buf , 0 , BWT_BUF_SIZE * sizeof(uint32_t));
                bcount = 0 ;
            }
            buf[bcount>>4] |= (((bwtInc64->bwt->bwtCode[oIndex>>4] >> packedShift[oIndex & 0xf]) & 0x3) << packedShift[mIndex & 0xf]);
            mIndex++ ; oIndex++ ; bcount++;
        }
        oldSAIndex++ ;
    }
    if(bcount > 0)
    {
        memcpy(mergedBWT + ((mIndex - bcount)>>4), buf, ((bcount + CHAR_PER_WORD - 1)>>4) * sizeof(uint32_t));
    }
    assert(mIndex == oldSALen + insertLen - divideCount );
    free(buf);
    
    return mergedBWT ;
}

uint32_t *mergeBWTInc64_single(const uint64_t* __restrict relativeRank, BWTInc64 *bwtInc64, const BWTInc *bwtInc)
{
    const uint32_t *bwtCode64 = bwtInc64->bwt->bwtCode , *bwtCode = bwtInc->bwt->bwtCode;
	uint32_t *mergedBWT = (uint32_t*)(bwtInc64->workingMemory + bwtInc64->availableWord - bwtInc64->bwt->bwtSizeInWord - BWT_ADDR_OFFSET) ;		
    int64_t iBWTIndex = 0 , oBWTIndex = 0, mBWTIndex = 0 ;
    uint32_t iWIndex , oWIndex , mWIndex  ;
    const int divideCount = bwtInc64->bwt->divideNumber ;
    uint32_t iWord, oWord , mWord ;
    int64_t oldSALen = bwtInc64->bwt->textLength + bwtInc64->bwt->divideNumber, insertBWTLen = bwtInc->bwt->textLength, oldSAIndex = 0 ;
    int sentinelCount = 0 ;
    uint32_t *packedShift = bwtInc64->packedShift ;
    // set mergedBWT 
    //uint32_t *mergedBWT = bwtInc64->workingMemory + bwtInc64->availableWord - bwtInc64->bwt->bwtSizeInWord - BWT_ADDR_OFFSET;
    //uint32_t *mergedBWT = (uint32_t*)xcalloc(bwtInc64->bwt->bwtSizeInWord, sizeof(uint32_t));

    iWord = bwtCode[iBWTIndex >> 4] ;  iWIndex = 0 ;
    oWord = bwtCode64[oBWTIndex >>4]; oWIndex = 0 ;
    mWord = 0 ;                     mWIndex = 0 ;
    for(int64_t i = 0; i < insertBWTLen; i++)
    {
        // copy from oldBWT
        int64_t beginIndex = oldSAIndex;
        while(oldSAIndex < relativeRank[i] && oldSAIndex < oldSALen)
        {
            if(oldSAIndex == bwtInc64->bwt->sentinelPosition[sentinelCount] && sentinelCount < divideCount)
            {
                if(mBWTIndex + (oldSAIndex - beginIndex) + sentinelCount >= bwtInc64->bwt->inverseSa0)
                {
                    bwtInc64->bwt->sentinelPosition[sentinelCount] = mBWTIndex + (oldSAIndex - beginIndex) + sentinelCount + 1;
                } else bwtInc64->bwt->sentinelPosition[sentinelCount] = mBWTIndex + (oldSAIndex - beginIndex) + sentinelCount;
                sentinelCount++;
                beginIndex++; // that SA is sentinel, not in BWT string, need skip this position
            } 
            /*// check
            if(bwtInc64->bwt->inverseSa0 - bwtInc->bwt->inverseSa0 == oldSAIndex)
            {
                if(mBWTIndex + sentinelCount + (oldSAIndex - beginIndex) != bwtInc64->bwt->inverseSa0 || i != bwtInc->bwt->inverseSa0)
                {
                    fprintf(stderr, "[mergeBWT64] mBWTIndex(%ld) + sentinelCount(%d) + (oldSAIndex(%ld) - beginIndex(%ld))!= bwtInc64->bwt->inverseSa0(%ld) or i(%d) != bwtInc->bwt->inverseSa0(%ld) exit..\n", mBWTIndex, sentinelCount, oldSAIndex, beginIndex, bwtInc64->bwt->inverseSa0, i, bwtInc->bwt->inverseSa0);
                    //exit(1);
                }
            } */
            oldSAIndex++;
        }
        int len = oldSAIndex - beginIndex;
        for(int j = 0; j < len;)
        {
            if(oWIndex >= CHAR_PER_WORD) { oWord = bwtCode64[oBWTIndex >>4] ; oWIndex = 0;  } 
            if(mWIndex >= CHAR_PER_WORD) { mergedBWT[(mBWTIndex - mWIndex)>>2] = mWord ; mWord = 0 ; mWIndex = 0; }
            int gap = (oWIndex <= mWIndex ? CHAR_PER_WORD - mWIndex : CHAR_PER_WORD - oWIndex);
            gap = (gap <= len - j ? gap : len -j);
            mWord <<= (gap<<1);
            uint32_t swap = (oWord >> ((CHAR_PER_WORD - (oWIndex + gap))<<1));
            mWord |= (swap & masklow32[gap]);
            mWIndex += gap; oWIndex += gap;
            mBWTIndex += gap; oBWTIndex += gap;
            j += gap;
        }

        // copy from insertBWT
        if(iWIndex >= CHAR_PER_WORD) { iWord = bwtCode[iBWTIndex >> 4] ;  iWIndex = 0;}
        if(mWIndex >= CHAR_PER_WORD) { mergedBWT[(mBWTIndex - mWIndex)>>4] = mWord ; mWord = 0 ; mWIndex = 0; }
        mWord <<= BIT_PER_CHAR ;
        mWord |= ((iWord >> packedShift[iWIndex & 0xf]) & 0x3) ;
        mWIndex++ ; iWIndex++;
        mBWTIndex++; iBWTIndex++;
    }
    // copy unfinished the oldBWT 
    while(oldSAIndex < oldSALen)
    {
        if(oldSAIndex == bwtInc64->bwt->sentinelPosition[sentinelCount] && sentinelCount < divideCount)
        {
            if(mBWTIndex + sentinelCount >= bwtInc64->bwt->inverseSa0)
            {
                bwtInc64->bwt->sentinelPosition[sentinelCount] = mBWTIndex + sentinelCount + 1;
            } else bwtInc64->bwt->sentinelPosition[sentinelCount] = mBWTIndex + sentinelCount;
            sentinelCount++;
        } else {
            if(oWIndex >= CHAR_PER_WORD) { oWord = bwtCode64[oBWTIndex >>4]; oWIndex = 0 ; }
            if(mWIndex >= CHAR_PER_WORD) { mergedBWT[(mBWTIndex - mWIndex)>>4] = mWord ; mWord = 0; mWIndex = 0; }
            mWord <<= BIT_PER_CHAR ;
            mWord |= ((oWord >> packedShift[oWIndex & 0xf]) & 0x3);
            mWIndex++ ; oWIndex++ ;
            mBWTIndex++ ; oBWTIndex++ ;
        }
        oldSAIndex++;
    }
    // clear mergedBWT buffer 
    if(mWIndex > 0)
    {
        mWord <<= ((CHAR_PER_WORD - mWIndex)<<1) ;
        mergedBWT[(mBWTIndex - mWIndex)>>4] = mWord ;
    }
    assert(mBWTIndex == oldSALen + insertBWTLen - divideCount);
    
    return mergedBWT ;
} 
                    
void mergeBWTInc64_part(const uint64_t *relativeRank, const int64_t start, const int64_t end, uint32_t *part_merged, const uint32_t *part_bwtCode64, const int64_t bwtCode64_start, const int64_t bwtCode64_end, BWTInc64 *bwtInc64, const BWTInc *bwtInc, const int64_t *sentinelPosition)
{
    const uint32_t *bwtCode = bwtInc->bwt->bwtCode; 
    int64_t iBWTIndex = start, oBWTIndex = bwtCode64_start, mBWTIndex = 0;
    uint32_t iWIndex, oWIndex, mWIndex;
    uint32_t iWord, oWord, mWord;
    uint32_t *packedShift = bwtInc64->packedShift;
    const int divideCount = bwtInc64->bwt->divideNumber;
    int sentinelCount;
    for(sentinelCount = 0; sentinelCount < bwtInc64->bwt->divideNumber; sentinelCount++)
    {
        if(sentinelPosition[sentinelCount] >= relativeRank[start]) break;
    }
    int64_t oldSAIndex = bwtCode64_start + sentinelCount ;
    int64_t bwt64SA_end;
    if(bwtCode64_end == bwtInc64->bwt->textLength)
    {
        bwt64SA_end = bwtInc64->bwt->textLength + bwtInc64->bwt->divideNumber;
    } else {
        bwt64SA_end = relativeRank[end];
    }
    int64_t bwt64Base = oBWTIndex / CHAR_PER_WORD * CHAR_PER_WORD;
    iWord = bwtCode[start>>4]; iWIndex = start % CHAR_PER_WORD;
    oWord = part_bwtCode64[(oBWTIndex - bwt64Base)>>4]; oWIndex = oBWTIndex % CHAR_PER_WORD;
    mWord = 0;  mWIndex = 0;
    for(int64_t i = start; i <  end; i++)
    {
        // copy from oldBWT
        int64_t beginIndex = oldSAIndex;
        while(oldSAIndex < relativeRank[i] && oldSAIndex < bwt64SA_end)
        {
            if(sentinelCount < divideCount && oldSAIndex == bwtInc64->bwt->sentinelPosition[sentinelCount])
            {
                if(oldSAIndex + i >= bwtInc64->bwt->inverseSa0)
                {
                    bwtInc64->bwt->sentinelPosition[sentinelCount] += (i + 1);
                } else bwtInc64->bwt->sentinelPosition[sentinelCount] += i;
                sentinelCount++;
                beginIndex++; // that SA is sentinel, not int BWT string, need skip this position 
            }
            oldSAIndex++;
        }
        int len = oldSAIndex - beginIndex;
        for(int j = 0; j < len;)
        {
            if(oWIndex >= CHAR_PER_WORD) {oWord = part_bwtCode64[(oBWTIndex - bwt64Base)>>4]; oWIndex = 0;}
            if(mWIndex >= CHAR_PER_WORD) { part_merged[(mBWTIndex - mWIndex)>>4] = mWord; mWord = 0; mWIndex = 0; }
            int gap = (oWIndex <= mWIndex ? CHAR_PER_WORD - mWIndex : CHAR_PER_WORD - oWIndex);
            gap = (gap <= len - j ? gap : len - j);
            mWord <<= (gap<<1);
            uint32_t swap = (oWord >> ((CHAR_PER_WORD - (oWIndex + gap))<<1));
            mWord |= (swap & masklow32[gap]);
            mWIndex += gap; oWIndex += gap;
            mBWTIndex += gap; oBWTIndex += gap;
            j += gap;
        }

        // copy from insertBWT
        if(iWIndex >= CHAR_PER_WORD) { iWord = bwtCode[iBWTIndex >> 4]; iWIndex = 0; }
        if(mWIndex >= CHAR_PER_WORD) { part_merged[(mBWTIndex - mWIndex)>>4] = mWord; mWord = 0; mWIndex = 0; }
        mWord <<= BIT_PER_CHAR;
        mWord |= ((iWord >> packedShift[iWIndex & 0xf]) & 0x3) ;
        mWIndex++ ; iWIndex++;
        mBWTIndex++; iBWTIndex++;
    }

    // copy unfinished the oldBWT
    while(oldSAIndex < bwt64SA_end)
    {
        if(oldSAIndex == bwtInc64->bwt->sentinelPosition[sentinelCount] && sentinelCount < divideCount)
        {
            if(oldSAIndex + end >= bwtInc64->bwt->inverseSa0)
            {
                bwtInc64->bwt->sentinelPosition[sentinelCount] += (end+ 1);
            } else bwtInc64->bwt->sentinelPosition[sentinelCount] += end;
            sentinelCount++;
        } else {
            if(oWIndex >= CHAR_PER_WORD) {oWord = part_bwtCode64[(oBWTIndex - bwt64Base)>>4]; oWIndex = 0;}
            if(mWIndex >= CHAR_PER_WORD) { part_merged[(mBWTIndex - mWIndex)>>4] = mWord; mWord = 0; mWIndex = 0; }
            mWord <<= BIT_PER_CHAR ;
            mWord |= ((oWord >> packedShift[oWIndex & 0xf]) & 0x3);
            mWIndex++ ; oWIndex++ ;
            mBWTIndex++ ; oBWTIndex++ ;
        }
        oldSAIndex++;
    }
    // clear mergedBWT buffer 
    if(mWIndex > 0)
    {
        mWord <<= ((CHAR_PER_WORD - mWIndex)<<1) ;
        part_merged[(mBWTIndex - mWIndex)>>4] = mWord ;
    }
}

uint32_t *mergeBWTInc64(const uint64_t* __restrict relativeRank, BWTInc64 *bwtInc64, const BWTInc *bwtInc, const int num_threads)
{
    uint32_t *mergedBWT;
    int error_flag = 0;
	bwtInc64->bwt->bwtSizeInWord = BWTResidentSizeInWord64(bwtInc64->bwt->textLength + bwtInc->bwt->textLength);
    if(num_threads <= 1)
    {
        mergedBWT = mergeBWTInc64_single(relativeRank, bwtInc64, bwtInc);
    } else {
        int64_t index_RR[num_threads+1], bwtCode64_index[num_threads+1];
        if((num_threads + 1) * BWT_SINGLE_OFFSET > BWT_ADDR_OFFSET)
        {
            error_flag = 1;
            fprintf(stderr, "[mergeBWTInc64] num_threads bigger than %d\n", BWT_ADDR_OFFSET/BWT_SINGLE_OFFSET - 1);
            //exit(1);
        } else {
            index_RR[0] = 0; bwtCode64_index[0] = 0;
            for(int i = 1; i < num_threads; i++)
            {
                index_RR[i] = (bwtInc->bwt->textLength / num_threads) * i;
                while(index_RR[i] < (bwtInc->bwt->textLength / num_threads) * (i + 0.1))
                {
                    while(relativeRank[index_RR[i]-1] == relativeRank[index_RR[i]])
                    {
                        index_RR[i]++;
                    }
                    int sc = 0;
                    for(int j = 0; j < bwtInc64->bwt->divideNumber; j++)
                    {
                        if(bwtInc64->bwt->sentinelPosition[j] < relativeRank[index_RR[i]]) sc++;
                        else break;
                    }
                    bwtCode64_index[i] = relativeRank[index_RR[i]] - sc;
                    if((index_RR[i]  + bwtCode64_index[i])%CHAR_PER_WORD == 0) break;
                    index_RR[i]++;
                }
                if(index_RR[i] >= (bwtInc->bwt->textLength / num_threads) * (i + 0.1))
                {
                    error_flag = 1;
                    break;
                    fprintf(stderr, "[mergeBWTInc64] divide index bigger than bwtInc->bwt->textLength\n");
                    //exit(1);
                }
            }
            index_RR[num_threads] = bwtInc->bwt->textLength ;
            bwtCode64_index[num_threads] = bwtInc64->bwt->textLength;
            if(error_flag == 1)
            {
                mergedBWT = mergeBWTInc64_single(relativeRank, bwtInc64, bwtInc);
            } else {
                uint32_t *part_bwtCode64[num_threads], *part_merged[num_threads];
                // set mergedBWT 
                mergedBWT = bwtInc64->workingMemory + bwtInc64->availableWord - bwtInc64->bwt->bwtSizeInWord - BWT_ADDR_OFFSET;
                int64_t cumulativeLen = 0;
                for(int i = 0; i < num_threads; i++)
                {
                    part_merged[i] = mergedBWT + cumulativeLen;
                    cumulativeLen = (index_RR[i+1] + bwtCode64_index[i+1] + CHAR_PER_WORD -1) / CHAR_PER_WORD + (i+1) * BWT_SINGLE_OFFSET;
                    part_bwtCode64[i] = mergedBWT + cumulativeLen - (bwtCode64_index[i+1] - bwtCode64_index[i] + CHAR_PER_WORD - 1)/CHAR_PER_WORD - 10;
                    memmove(part_bwtCode64[i], bwtInc64->bwt->bwtCode + bwtCode64_index[i]/CHAR_PER_WORD, ((bwtCode64_index[i+1] + CHAR_PER_WORD -1)/CHAR_PER_WORD - bwtCode64_index[i]/CHAR_PER_WORD)* BYTES_IN_WORD); 

                }
                int64_t sp[bwtInc64->bwt->divideNumber];
                memcpy(sp, bwtInc64->bwt->sentinelPosition, bwtInc64->bwt->divideNumber * sizeof(int64_t));
                #pragma omp parallel for schedule(dynamic)
                for(int i = 0; i < num_threads; i++)
                {
                    mergeBWTInc64_part(relativeRank, index_RR[i], index_RR[i+1], part_merged[i], part_bwtCode64[i], bwtCode64_index[i], bwtCode64_index[i+1], bwtInc64, bwtInc, sp);
                }

                // merging part BWT
                for(int i = 1; i < num_threads; i++)
                {
                    memmove(mergedBWT + (index_RR[i] + bwtCode64_index[i])/CHAR_PER_WORD, part_merged[i], (index_RR[i+1] + bwtCode64_index[i+1] - index_RR[i] - bwtCode64_index[i] + CHAR_PER_WORD -1)/CHAR_PER_WORD * BYTES_IN_WORD);
                }
            }
        }
    }

    return mergedBWT;
}


uint64_t smallOrCloseSearch(const uint64_t *relativeRank, const uint64_t bwtLocation, const uint64_t size)
{
    uint64_t low = 0 , high = size - 1 , mid ;
    while(low <= high)
    {
        mid = (low + high)/ 2 ;
        if(relativeRank[mid] == bwtLocation) return mid + 1;
        if(relativeRank[mid] > bwtLocation) high = mid - 1 ;
        else low = mid + 1 ;
    }
    if(relativeRank[high] < bwtLocation)
    {
        while(relativeRank[high] < bwtLocation) { mid = high ; high++ ;}
    } else {
        while(relativeRank[high] > bwtLocation) { mid = high ; high-- ;}
    }
    return mid + 1 ;
}
/*  transform the position of SA to the position in the BWT string ,
 *  if encounter a sentinelPosition, return -1, else return the position of BWT string
 */
static inline int64_t posSA2PosBWT64(BWT64 *bwt, const int64_t isa)
{
    int64_t ret_v = -2 ;
    int i ;
    for(i = 0; i < bwt->divideNumber; i++)
    {
        if(isa < bwt->sentinelPosition[i]) break ;
        else if(isa == bwt->sentinelPosition[i]) {
            fprintf(stderr, "[posSA2PosBWT64] isa == bwt->sentinelPosition[%d](%lu\n)", i, bwt->sentinelPosition[i]);
            ret_v = -1 ;
            break;
        }
    }
    if(ret_v == -2) ret_v = isa - i ; 
    return ret_v ;
}

void mergeSortRank64(const uint64_t *src1, const int64_t src1_len, const uint64_t *src2, const int64_t src2_len, uint64_t *dest)
{
    int64_t destIndex = 0, src2Index = 0;
    for(int64_t src1Index = 0; src1Index < src1_len; src1Index++)
    {
        while(src2Index < src2_len && src2[src2Index] <= src1[src1Index])
        {
            dest[destIndex++] = src2[src2Index++];
        }
        dest[destIndex++] = src1[src1Index];
    }
    while(src2Index < src2_len)
    {
        dest[destIndex++] = src2[src2Index++];
    }
}

BWTInc64 *mergeBWT64(BWTInc64 *bwtInc64, const BWTInc *bwtInc, const int divideCount, const char *fileName, const int num_threads)
{
    uint64_t *relativeRank ;
	int64_t cur_sa , bwt64_cur_sa ;
    int64_t last_pos = bwtInc->bwt->textLength ;
	FILE *fp = xopen(fileName, "rb");
    time_t timeval;
    clock_t t = clock();
    relativeRank = (uint64_t*)bwtInc64->workingMemory ;
    int64_t p[num_threads], end[num_threads], sarray[num_threads], bwt64_sarray[num_threads];
    sarray[0] = 0; bwt64_sarray[0] = divideCount;
    // the argument for parallel threads
    for(int i = 0; i < num_threads; i++)
    {
        int64_t pos = last_pos - (int64_t)bwtInc->bwt->textLength/num_threads - (int64_t)bwtInc->bwt->textLength % num_threads ;
        if(i < num_threads - 1)
        {
            int len = CHAR_PER_BYTE * 800;
            int flag = 0;
            int64_t cor_pos = pos - len;
            cor_pos = cor_pos/CHAR_PER_BYTE * CHAR_PER_BYTE ; 
            uint8_t str[len/CHAR_PER_BYTE + 3];
            int64_t low = -5 , high = -1, bwt64_low = -5, bwt64_high = -1 ;
            while(cor_pos > 0)
            {
                fseek(fp, (cor_pos/CHAR_PER_BYTE) * sizeof(uint8_t), SEEK_SET);
                int n = fread(str, sizeof(uint8_t), len/CHAR_PER_BYTE + 1 , fp);
                if(n <= len/CHAR_PER_BYTE) 
                {
                    fprintf(stderr, "[mergeBWT64] read sequence number less than expect, exit...\n");
                    exit(1);
                }
                for(int j = len - 1; j >= 0; j--)
                {
                    uint32_t c = GET_CHAR_FROM_PACKED(str, j);
                    if(low == -5 && high == -1)
                    {
                        low = bwtInc->bwt->cumulativeFreq[c] + 1;
                        high = bwtInc->bwt->cumulativeFreq[c+1] +1 - 1;
                        bwt64_low = bwtInc64->bwt->cumulativeFreq[c] + divideCount ;
                        bwt64_high = bwtInc64->bwt->cumulativeFreq[c+1] + divideCount - 1;
                    } else {
                        low = bwtInc->bwt->cumulativeFreq[c] + BWTOccValue(bwtInc->bwt, low, c) + 1 ; 
                        high = bwtInc->bwt->cumulativeFreq[c] + BWTOccValue(bwtInc->bwt, high, c) + 1 ;
                        bwt64_low = bwtInc64->bwt->cumulativeFreq[c] + BWTOccValue64(bwtInc64->bwt, bwt64_low, c) + divideCount ;
                        bwt64_high = bwtInc64->bwt->cumulativeFreq[c] + BWTOccValue64(bwtInc64->bwt, bwt64_high, c) + divideCount;
                    }
                    fprintf(stderr, "%u", c);
					/* debug 
					{
						int i ;
						for(i = 0 ; i < bwtInc64->bwt->divideNumber; i++)
						{
							if(bwt64_low <= bwtInc64->bwt->sentinelPosition[i]) break ;
						}
						int64_t isa = bwt64_low - i ;
						uint32_t tmp = (bwtInc64->bwt->bwtCode[isa>>5] >> (((~isa)& 0xf)<<1))&0x3 ;
						fprintf(stderr, "%u", tmp);
					} // end debug */

                    if(low + 1 == high && bwt64_low == bwt64_high)
                    {
                        pos = cor_pos + j ;
                        //if(low > bwtInc->bwt->inverseSa0) low-- ;
                        //bwt64_low = posSA2PosBWT64(bwtInc64->bwt, bwt64_low);
                        if(i > 0) 
                        {
                            bwt64_sarray[i] = bwt64_cur_sa;
                            sarray[i] = cur_sa;
                        }
                        bwt64_cur_sa = bwt64_low ;
                        cur_sa = low ;
                        flag = 1 ; break;
                    }
                }
                fprintf(stderr, "\n");
                if(flag == 1) break;
                else cor_pos -= len;
            }
            if(flag != 1) 
			{
				fprintf(stderr, "[mergeBWT64] flag != 1, error, program exit...\n"); exit(1); 
			}
        } else {
            pos = 0 ;
            if(i > 0)
            {
                sarray[i] = cur_sa ;
                bwt64_sarray[i] = bwt64_cur_sa ;
            }
            uint8_t str[30];
            fseek(fp, 0, SEEK_SET);
            fread(str, sizeof(uint8_t), 30, fp);
            for(int j = 100 - 1; j >= 0; j--)
            {
                int c = GET_CHAR_FROM_PACKED(str, j);
                fprintf(stderr, "%u", c);
            }
            fprintf(stderr, "\n");
        }
        end[i] = last_pos;
        p[i] = pos;
        last_pos = pos;
    }
    time(&timeval); fprintf(stderr, "begin build relative rank\ntime:\t%s\n", ctime(&timeval)); 
    //bwtInc64->bwt->divideBoundary[divideCount].bwtLocation = bwtInc64->bwt->inverseSa0 ;
    // buildRelativeRank64 function maybe need parallel ,if will be a bottleneck of time complexity
    omp_set_num_threads(num_threads);
    t = clock();
	//#pragma omp parallel for firstprivate(relativeRank, bwtInc64, bwtInc, p, end, sarray, bwt64_sarray) schedule(dynamic)
	#pragma omp parallel for schedule(dynamic) 
	for(int i = 0; i < num_threads; i++)
	{
        if(p[i] == 0 )
		{
            bwtInc64->bwt->inverseSa0 = buildRelativeRankFromOriginPac(fileName, relativeRank, bwtInc64->bwt, bwtInc->bwt, p[i], end[i], bwt64_sarray[i]);
        } else {
		    buildRelativeRankFromOriginPac(fileName, relativeRank, bwtInc64->bwt, bwtInc->bwt, p[i], end[i], bwt64_sarray[i]);
        }
    }
    //merge sort
    {
        int64_t processed_len = 0;
        for(int i = 1; i < num_threads; i++)
        {
            int64_t len = end[i] - p[i];
            uint64_t *swapbuf = (uint64_t*)xcalloc(len, sizeof(uint64_t));
            memcpy(swapbuf, relativeRank + p[i], len * sizeof(uint64_t));
            processed_len += (end[i-1] - p[i-1]);
            mergeSortRank64(swapbuf, len, relativeRank + p[i-1], processed_len, relativeRank + p[i]);

            // free and clean work
            free(swapbuf);
        }
    }
    
    fprintf(stderr, "[mergeBWT64]buildRelativeRankFromOriginPac used %.2f CPU secs\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    time(&timeval); fprintf(stderr, "finished building relative rank\ntime:\t%s\n", ctime(&timeval));
    t = clock();
    //uint32_t *mergeBWT1 = mergeBWTInc64_1209(relativeRank, bwtInc64, bwtInc);
    //bwtInc64->bwt->bwtCode = mergeBWTInc64_single(relativeRank, bwtInc64, bwtInc);
    bwtInc64->bwt->bwtCode = mergeBWTInc64(relativeRank, bwtInc64, bwtInc, num_threads);

    fprintf(stderr, "[mergeBWT64]mergeBWTInc64 used %.2f CPU secs\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    time(&timeval); fprintf(stderr, "finished merging BWTInc\ntime:\t%s\n", ctime(&timeval));
    // adjust SA of the sentinels 
    {
        int p, j ;
        for(j = 0 ; j < divideCount ; j++)
        {
            if(bwtInc64->bwt->inverseSa0 <= bwtInc64->bwt->sentinelPosition[j]) break ;
        }
        p = j ;
        for(j = divideCount -1 ; j >= p ; j--)
        {
            bwtInc64->bwt->sentinelPosition[j+1] = bwtInc64->bwt->sentinelPosition[j];
            bwtInc64->bwt->sentinelSA[j+1] = bwtInc64->bwt->sentinelSA[j];
        }
        bwtInc64->bwt->sentinelPosition[p] = bwtInc64->bwt->inverseSa0 ;
        bwtInc64->bwt->sentinelSA[p] = bwtInc->bwt->saIndexRange->startSaIndex ;
    }

    // build auxiliary structure and update info and pointers in BWT
    bwtInc64->bwt->textLength += bwtInc->bwt->textLength ;
    // bwtInc64->bwt->bwtSizeInWord has been set in mergeBWTInc64();
    bwtInc64->bwt->occSizeInWord = BWTOccValueMinorSizeInWord64(bwtInc64->bwt->textLength);
    bwtInc64->bwt->occMajorSizeInWord = BWTOccValueMajorSizeInWord64(bwtInc64->bwt->textLength);
    bwtInc64->bwt->occValueMajor = (uint64_t*)bwtInc64->workingMemory ;
    bwtInc64->bwt->occValue = (uint32_t*)(bwtInc64->workingMemory + bwtInc64->bwt->occMajorSizeInWord * 2 + OCC_INTERVAL);
    BWTGenerateOccValueFromBwt64(bwtInc64->bwt->bwtCode, bwtInc64->bwt->occValue, bwtInc64->bwt->occValueMajor , bwtInc64->bwt->textLength, bwtInc64->bwt->decodeTable);
    bwtInc64->firstCharInLastIteration = bwtInc->firstCharInLastIteration ;
    // update divideBoundary location in the bwt
    bwtInc64->bwt->divideBoundary[divideCount].textLocation = bwtInc->bwt->saIndexRange->endSaIndex  ;
    bwtInc64->bwt->divideBoundary[divideCount].bwtLocation = divideCount ;
    bwtInc64->numberOfIterationDone++ ;
    bwtInc64->bwt->cumulativeFreq[1] += bwtInc->bwt->cumulativeFreq[1] ;
    bwtInc64->bwt->cumulativeFreq[2] += bwtInc->bwt->cumulativeFreq[2] ;
    bwtInc64->bwt->cumulativeFreq[3] += bwtInc->bwt->cumulativeFreq[3] ;
    bwtInc64->bwt->cumulativeFreq[4] += bwtInc->bwt->cumulativeFreq[4] ;
	// clean and free work
    time(&timeval); fprintf(stderr, "finished merging proceed\ntime:\t%s\n", ctime(&timeval));
	fclose(fp);
    return bwtInc64 ;
}
