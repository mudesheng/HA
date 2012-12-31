/*

   BWTConstruct.h		BWT-Index Construction

   This module constructs BWT and auxiliary data structures.

   Copyright (C) 2004, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef BWT_GEN64_H
#define BWT_GEN64_H

#include"bwt_gen.h"
/*
#define ALPHABET_SIZE				4
#define BIT_PER_CHAR				2
#define CHAR_PER_WORD				16
#define CHAR_PER_BYTE				4

#define BITS_IN_WORD 32
#define BITS_IN_BYTE 8
#define BYTES_IN_WORD 4

#define DNA_OCC_CNT_TABLE_SIZE_IN_WORD	65536

#define BITS_PER_OCC_VALUE			16
#define OCC_VALUE_PER_WORD			2
#define OCC_INTERVAL				256
#define OCC_INTERVAL_MAJOR			65536

#define TRUE    1
#define FALSE   0

#define BWTINC_INSERT_SORT_NUM_ITEM 7

#define average(value1, value2)					( ((value1) & (value2)) + ((value1) ^ (value2)) / 2 )
#define min(value1, value2)						( ((value1) < (value2)) ? (value1) : (value2) )
#define max(value1, value2)						( ((value1) > (value2)) ? (value1) : (value2) )
#define med3(a, b, c)							( a<b ? (b<c ? b : a<c ? c : a) : (b>c ? b : a>c ? c : a))
#define swap(a, b, t);							t = a; a = b; b = t;
#define truncateLeft(value, offset)				( (value) << (offset) >> (offset) )
#define truncateRight(value, offset)			( (value) >> (offset) << (offset) )
#define DNA_OCC_SUM_EXCEPTION(sum)			((sum & 0xfefefeff) == 0)
*/

#define ALL_ONE_MASK64 0xFFFFFFFFFFFFFFFF
#define BWTINT uint64_t
#define OCC_INTERVAL_MAJOR64	    0X100000000 //4294967296
#define BWT_ADDR_OFFSET 80000
#define BWT_SINGLE_OFFSET 800
#define BWT_BUF_SIZE 0x1000
#define SUPERBLOCK_SIZE 65536
#define BLOCK_SIZE 4096
#define MAX_COMPRESS_BLOCK 63 // max compress block 2**6 -1, just 6 bit encode for counter  

#define GET_CHAR_FROM_PACKED(str, offset) ((str[(offset)/4] >> (((~(offset)) & 0X3)<<1)) & 0X3)

/*
typedef struct SaIndexRange64 {
	uint64_t startSaIndex;
	uint64_t endSaIndex;
} SaIndexRange64;
*/
/* every SAValue store two SA value, first value ((high << 32) | low)*/
typedef struct SAValue {
	uint32_t low[2];
	uint16_t high[2];
} SAValue;


typedef struct DivideBoundary {
    uint64_t textLocation ;
    uint64_t bwtLocation ;
} DivideBoundary ;

typedef struct BWT64 {
	uint64_t textLength;			// length of the text
	uint64_t saInterval;			// interval between two SA values stored explicitly
	uint64_t inverseSaInterval;		// interval between two inverse SA stored explicitly
	uint64_t inverseSa0;			// SA-1[0]
	uint64_t *cumulativeFreq;		// cumulative frequency
	uint64_t *superBlockIndex;		// super block index of bwtCode
	uint32_t *bwtCode;				// BWT code compression state 
	unsigned int *occValue;				// Occurrence values stored explicitly
	uint64_t *occValueMajor;		// Occurrence major values stored explicitly
	SAValue *saValue;				// SA values stored explicitly
	uint64_t *inverseSa;			// Inverse SA stored explicitly
	SaIndexRange *saIndexRange;			// SA index range
	int64_t saIndexRangeNumOfChar;			// Number of characters indexed in SA index range
	uint64_t *saValueOnBoundary;	// Pre-calculated frequently referred data
	unsigned int *decodeTable;			// For decoding BWT by table lookup
	unsigned int decodeTableGenerated;	// == TRUE if decode table is generated on load and will be freed
	uint64_t bwtSizeInWord;			// Temporary variable to hold the memory allocated
	uint64_t occSizeInWord;			// Temporary variable to hold the memory allocated
	uint64_t occMajorSizeInWord;	// Temporary variable to hold the memory allocated
	uint64_t saValueSize;			// Temporary variable to hold the memory allocated
	uint64_t inverseSaSize;			// Temporary variable to hold the memory allocated
	uint64_t saIndexRangeSize;		// Temporary variable to hold the memory allocated
    uint8_t isDivide ;              // Set if divided construction method
    DivideBoundary *divideBoundary ;// The boundary of divided part 
    uint32_t divideNumber ;         // The number of parts divided
    uint64_t *sentinelPosition ;    // The positions of sentinels '$' in the BWT code 
    uint64_t *sentinelSA ;          // The SA of the sentinels corresponding sentinelPosition
	uint64_t compressed_size ;		// total compressed size of the bwtCode 
} BWT64;

typedef struct BWTInc64 {
	BWT64 *bwt;
	unsigned int numberOfIterationDone;
	uint64_t *cumulativeCountInCurrentBuild;
	uint64_t availableWord;
	uint64_t targetTextLength;
	float targetNBit;
	unsigned int buildSize;
	unsigned int initialMaxBuildSize;
	unsigned int incMaxBuildSize;
	unsigned int firstCharInLastIteration;
	unsigned int *workingMemory;
	unsigned int *packedText;
	unsigned char *textBuffer;
	unsigned int *packedShift;
} BWTInc64;

uint64_t BWTOccValue64(const BWT64 *bwt, uint64_t index, const unsigned int character);
void BWTGenerateOccValueFromBwt64(const unsigned int*  bwt, unsigned int* __restrict occValue, uint64_t* __restrict occValueMajor, const uint64_t textLength, const unsigned int*  decodeTable) ;
uint64_t BWTOccValueMajorSizeInWord64(const uint64_t numChar);
uint64_t BWTOccValueMinorSizeInWord64(const uint64_t numChar);
uint64_t BWTResidentSizeInWord64(const uint64_t numChar) ;
BWTInc64 *BWTIncConstructFromPacked64(const char *inputFileName, const float targetNBit , const uint64_t initialMaxBuildSize, const uint64_t incMaxBuildSize) ;
BWTInc64 *BWTInc2BWTInc64(BWTInc *bwtInc) ;
void BWTIncSetMergeSize64(BWTInc64 *bwtInc64, BWTInc * bwtInc) ;
BWTInc64 *mergeBWT64(BWTInc64 *bwtInc64, const BWTInc *bwtInc, const int divideCount, const char *fileName, const int num_threads) ;
void BWTSaveBwtCodeAndOcc64(const BWT64 *bwt, const char *bwtFileName, const char *occValueFileName) ;
void BWTQSort64(bgint_t* __restrict key, const bgint_t numItem);
void BWTIncFree64(BWTInc64 *bwtInc) ;


#endif
