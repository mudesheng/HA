/*

   BWTConstruct.c		BWT-Index Construction

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "bwt_gen.h"
#include "bwt_gen64.h"
#include "QSufSort.h"

static uint64_t TextLengthFromBytePacked64(uint64_t bytePackedLength, unsigned int bitPerChar,
											 unsigned int lastByteLength)
{
	if (bytePackedLength > ALL_ONE_MASK64 / (BITS_IN_BYTE / bitPerChar)) {
		fprintf(stderr, "TextLengthFromBytePacked64(): text length > 2^64!\n");
		exit(1);
	}
	return (bytePackedLength - 1) * (BITS_IN_BYTE / bitPerChar) + lastByteLength;
}

static void initializeVAL64(uint64_t *startAddr, const uint64_t length, const int64_t initValue)
{
	uint64_t i;
	for (i=0; i<length; i++) startAddr[i] = initValue;
}

static void initializeVAL(unsigned int *startAddr, const uint64_t length, const unsigned int initValue)
{
	uint64_t i;
	for (i=0; i<length; i++) startAddr[i] = initValue;
}

static void GenerateDNAOccCountTable(unsigned int *dnaDecodeTable)
{
	unsigned int i, j, c, t;

	for (i=0; i<DNA_OCC_CNT_TABLE_SIZE_IN_WORD; i++) {
		dnaDecodeTable[i] = 0;
		c = i;
		for (j=0; j<8; j++) {
			t = c & 0x00000003;
			dnaDecodeTable[i] += 1 << (t * 8);
			c >>= 2;
		}
	}

}
// for BWTIncCreate()
uint64_t BWTOccValueMajorSizeInWord64(const uint64_t numChar)
{
	uint64_t numOfOccValue;
	uint64_t numOfOccIntervalPerMajor;
	numOfOccValue = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL + 1; // Value at both end for bi-directional encoding
	numOfOccIntervalPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
	return (numOfOccValue + numOfOccIntervalPerMajor - 1) / numOfOccIntervalPerMajor * ALPHABET_SIZE;
}
// for BWTIncCreate()
uint64_t BWTOccValueMinorSizeInWord64(const uint64_t numChar)
{
	uint64_t numOfOccValue;
	numOfOccValue = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;		// Value at both end for bi-directional encoding
	return (numOfOccValue + OCC_VALUE_PER_WORD - 1) / OCC_VALUE_PER_WORD * ALPHABET_SIZE;
}
// for BWTIncCreate()
uint64_t BWTResidentSizeInWord64(const uint64_t numChar) {

	uint64_t numCharRoundUpToOccInterval;

	// The $ in BWT at the position of inverseSa0 is not encoded
	numCharRoundUpToOccInterval = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL;

	return (numCharRoundUpToOccInterval + CHAR_PER_WORD - 1) / CHAR_PER_WORD;

}

static void BWTIncSetBuildSizeAndTextAddr64(BWTInc64 *bwtInc)
{
	uint64_t maxBuildSize;

	if (bwtInc->bwt->textLength == 0) {
		// initial build
		// Minus 2 because n+1 entries of seq and rank needed for n char
		maxBuildSize = (bwtInc->availableWord - 2 - OCC_INTERVAL / CHAR_PER_WORD)
							/ (2 * CHAR_PER_WORD + 1) * CHAR_PER_WORD;
		if (bwtInc->initialMaxBuildSize > 0) {
			bwtInc->buildSize = min(bwtInc->initialMaxBuildSize, maxBuildSize);
		} else {
			bwtInc->buildSize = maxBuildSize;
		}
	} else {
		// Minus 3 because n+1 entries of sorted rank, seq and rank needed for n char
		// Minus numberOfIterationDone because bwt slightly shift to left in each iteration
		maxBuildSize = (bwtInc->availableWord - bwtInc->bwt->bwtSizeInWord - bwtInc->bwt->occSizeInWord - 3
							 - bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR) 
							 / 4;
		if (maxBuildSize < CHAR_PER_WORD) {
			fprintf(stderr, "BWTIncSetBuildSizeAndTextAddr64(): Not enough space allocated to continue construction!\n");
			exit(1);
		}
		if (bwtInc->incMaxBuildSize > 0) {
            bwtInc->buildSize = min(bwtInc->incMaxBuildSize, maxBuildSize);
		} else {
			bwtInc->buildSize = maxBuildSize;
		}
		if (bwtInc->buildSize < CHAR_PER_WORD) {
			bwtInc->buildSize = CHAR_PER_WORD;
		}
	}

	if (bwtInc->buildSize < CHAR_PER_WORD) {
		fprintf(stderr, "BWTIncSetBuildSizeAndTextAddr64(): Not enough space allocated to continue construction!\n");
		exit(1);
	}

	bwtInc->buildSize = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;

	bwtInc->packedText = bwtInc->workingMemory + 3 * (bwtInc->buildSize + 1);
	bwtInc->textBuffer = (unsigned char*)(bwtInc->workingMemory + 2 * (bwtInc->buildSize + 1));

}

// for ceilLog2()
static unsigned int leadingZero(const unsigned int input)
{
	unsigned int l;
	const static unsigned int leadingZero8bit[256] = {8,7,6,6,5,5,5,5,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
											 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
											 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
											 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
											 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	if (input & 0xFFFF0000) {
		if (input & 0xFF000000) {
			l = leadingZero8bit[input >> 24];
		} else {
			l = 8 + leadingZero8bit[input >> 16];
		}
	} else {
		if (input & 0x0000FF00) {
			l = 16 + leadingZero8bit[input >> 8];
		} else {
			l = 24 + leadingZero8bit[input];
		}
	}
	return l;

}
// for BitPerBytePackedChar()
static unsigned int ceilLog2(const unsigned int input)
{
	if (input <= 1) return 0;
	return BITS_IN_WORD - leadingZero(input - 1);

}
// for ConvertBytePackedToWordPacked()
static unsigned int BitPerBytePackedChar(const unsigned int alphabetSize)
{
	unsigned int bitPerChar;
	bitPerChar = ceilLog2(alphabetSize);
	// Return the largest number of bit that does not affect packing efficiency
	if (BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar) > bitPerChar)
		bitPerChar = BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar);
	return bitPerChar;
}
// for ConvertBytePackedToWordPacked()
static unsigned int BitPerWordPackedChar(const unsigned int alphabetSize)
{
	return ceilLog2(alphabetSize);
}

static void ConvertBytePackedToWordPacked64(const unsigned char *input, unsigned int *output, const unsigned int alphabetSize,
										  const uint64_t textLength)
{
	uint64_t i, j, k;
	unsigned int c;
	unsigned int bitPerBytePackedChar;
	unsigned int bitPerWordPackedChar;
	unsigned int charPerWord;
	unsigned int charPerByte;
	unsigned int bytePerIteration;
	uint64_t byteProcessed = 0;
	uint64_t wordProcessed = 0;
	unsigned int mask, shift;
	
	unsigned int buffer[BITS_IN_WORD];

	bitPerBytePackedChar = BitPerBytePackedChar(alphabetSize);
	bitPerWordPackedChar = BitPerWordPackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerBytePackedChar;
	charPerWord = BITS_IN_WORD / bitPerWordPackedChar;

	bytePerIteration = charPerWord / charPerByte;
	mask = truncateRight(ALL_ONE_MASK, BITS_IN_WORD - bitPerWordPackedChar);
	shift = BITS_IN_WORD - BITS_IN_BYTE + bitPerBytePackedChar - bitPerWordPackedChar;

	while ((wordProcessed + 1) * charPerWord < textLength) {

		k = 0;
		for (i=0; i<bytePerIteration; i++) {
			c = (unsigned int)input[byteProcessed] << shift;
			for (j=0; j<charPerByte; j++) {
				buffer[k] = c & mask;
				c <<= bitPerBytePackedChar;
				k++;
			}
			byteProcessed++;
		}

		c = 0;
		for (i=0; i<charPerWord; i++) {
			c |= buffer[i] >> bitPerWordPackedChar * i;
		}
		output[wordProcessed] = c;
		wordProcessed++;

	}

	k = 0;
	for (i=0; i < (textLength - wordProcessed * charPerWord - 1) / charPerByte + 1; i++) {
		c = (unsigned int)input[byteProcessed] << shift;
		for (j=0; j<charPerByte; j++) {
			buffer[k] = c & mask;
			c <<= bitPerBytePackedChar;
			k++;
		}
		byteProcessed++;
	}

	c = 0;
	for (i=0; i<textLength - wordProcessed * charPerWord; i++) {
		c |= buffer[i] >> bitPerWordPackedChar * i;
	}
	output[wordProcessed] = c;
}
/*
static void dumpAndWriteStorage(const unsigned char *input,const int64_t start, const int64_t end, const char *fileName)
{
    int64_t textLength = end - start;
    uint8_t lastByteLen ;
    unsigned char *textBuffer;
    FILE *fp ;

    fp = xopen(fileName, "wb");
    textBuffer = (unsigned char*)calloc(textLength/ CHAR_PER_BYTE + 1, sizeof(unsigned char));
    for(uint64_t textLocation = 0 , i = start ; i < end ; i++, textLocation++)
    {
        textBuffer[textLocation>>2] |= ((input[i>>2] >> ((~i & 0x3)<<1)) & 0x3);
        (textLocation & 0x3 == 0x3 ) ? 0 : textBuffer[textLocation>>2] <<= BIT_PER_CHAR;
    }
    fwrite(textBuffer, sizeof(unsigned char), textLength / CHAR_PER_BYTE + 1 , fp);
    lastByteLen = textLength % CHAR_PER_BYTE ;
    fwrite(&lastByteLen, sizeof(uint8_t), 1, fp);
    fclose(fp);
    free(textBuffer);
} */

BWT64 *BWTCreate64(const uint64_t textLength, unsigned int *decodeTable)
{
	BWT64 *bwt;

	bwt = (BWT64*)calloc(1, sizeof(BWT64));

	bwt->textLength = 0;
	bwt->inverseSa = 0;

	bwt->cumulativeFreq = (uint64_t*)calloc((ALPHABET_SIZE + 1), sizeof(uint64_t));
	initializeVAL64(bwt->cumulativeFreq, ALPHABET_SIZE + 1, 0);

	bwt->bwtSizeInWord = 0;
	bwt->saValueOnBoundary = NULL;

	// Generate decode tables
	if (decodeTable == NULL) {
		bwt->decodeTable = (unsigned*)calloc(DNA_OCC_CNT_TABLE_SIZE_IN_WORD, sizeof(unsigned int));
		GenerateDNAOccCountTable(bwt->decodeTable);
	} else {
		bwt->decodeTable = decodeTable;
	}

	bwt->occMajorSizeInWord = BWTOccValueMajorSizeInWord64(textLength);
	bwt->occValueMajor = (uint64_t*)calloc(bwt->occMajorSizeInWord, sizeof(uint64_t));

	bwt->occSizeInWord = 0;
	bwt->occValue = NULL;

	bwt->saInterval = ALL_ONE_MASK;
	bwt->saValueSize = 0;
	bwt->saValue = NULL;

	bwt->inverseSaInterval = ALL_ONE_MASK;
	bwt->inverseSaSize = 0;
	bwt->inverseSa = NULL;

	return bwt;
}

BWTInc64 *BWTIncCreate64(const uint64_t textLength, const float targetNBit,
					 const uint64_t initialMaxBuildSize, const uint64_t incMaxBuildSize)
{
	BWTInc64 *bwtInc;
	unsigned int i;

	if (targetNBit == 0) {
		fprintf(stderr, "BWTIncCreate64() : targetNBit = 0!\n");
		exit(1);
	}
	
	bwtInc = (BWTInc64*)calloc(1, sizeof(BWTInc64));
	bwtInc->numberOfIterationDone = 0;
	bwtInc->bwt = BWTCreate64(textLength, NULL);
	bwtInc->initialMaxBuildSize = initialMaxBuildSize;
	bwtInc->incMaxBuildSize = incMaxBuildSize;
	bwtInc->targetNBit = targetNBit;
	bwtInc->cumulativeCountInCurrentBuild = (uint64_t*)calloc((ALPHABET_SIZE + 1), sizeof(uint64_t));  
	initializeVAL64(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

	// Build frequently accessed data
	bwtInc->packedShift = (unsigned*)calloc(CHAR_PER_WORD, sizeof(unsigned int));
	for (i=0; i<CHAR_PER_WORD; i++) {
		bwtInc->packedShift[i] = BITS_IN_WORD - (i+1) * BIT_PER_CHAR;
	}

	bwtInc->targetTextLength = textLength;
	bwtInc->availableWord = (uint64_t)((textLength + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL / BITS_IN_WORD * bwtInc->targetNBit);
	if (bwtInc->availableWord < BWTResidentSizeInWord64(textLength) + BWTOccValueMinorSizeInWord64(textLength)) {
		fprintf(stderr, "BWTIncCreate64() : targetNBit is too low!\n");
		exit(1);
	}
	bwtInc->workingMemory = (unsigned*)calloc(bwtInc->availableWord, BYTES_IN_WORD);

	return bwtInc;

}
// for BWTIncConstruct()
static void BWTIncPutPackedTextToRank(const unsigned int *packedText, unsigned int* __restrict rank,
									  uint64_t * __restrict cumulativeCount, const unsigned int numChar)
{
	unsigned int i, j;
	unsigned int c, t;
	unsigned int packedMask;
	unsigned int rankIndex;
	unsigned int lastWord, numCharInLastWord;

	lastWord = (numChar - 1) / CHAR_PER_WORD;
	numCharInLastWord = numChar - lastWord * CHAR_PER_WORD;

	packedMask = ALL_ONE_MASK >> (BITS_IN_WORD - BIT_PER_CHAR);
	rankIndex = numChar - 1;

	t = packedText[lastWord] >> (BITS_IN_WORD - numCharInLastWord * BIT_PER_CHAR);
	for (i=0; i<numCharInLastWord; i++) {
		c = t & packedMask;
		cumulativeCount[c+1]++;
		rank[rankIndex] = c;
		rankIndex--;
		t >>= BIT_PER_CHAR;
	}

	for (i=lastWord; i--;) {	// loop from lastWord - 1 to 0
		t = packedText[i];
		for (j=0; j<CHAR_PER_WORD; j++) {
			c = t & packedMask;
			cumulativeCount[c+1]++;
			rank[rankIndex] = c;
			rankIndex--;
			t >>= BIT_PER_CHAR;
		}
	}

	// Convert occurrence to cumulativeCount
	cumulativeCount[2] += cumulativeCount[1];
	cumulativeCount[3] += cumulativeCount[2];
	cumulativeCount[4] += cumulativeCount[3];
}


static void ForwardDNAAllOccCountNoLimit(const unsigned int*  dna, const unsigned int index,
										 uint64_t * __restrict occCount, const unsigned int*  dnaDecodeTable)
{
	static const unsigned int truncateRightMask[16] = { 0x00000000, 0xC0000000, 0xF0000000, 0xFC000000,
											   0xFF000000, 0xFFC00000, 0xFFF00000, 0xFFFC0000,
											   0xFFFF0000, 0xFFFFC000, 0xFFFFF000, 0xFFFFFC00,
											   0xFFFFFF00, 0xFFFFFFC0, 0xFFFFFFF0, 0xFFFFFFFC };

	unsigned int iteration, wordToCount, charToCount;
	unsigned int i, j, c;
	unsigned int sum;

	occCount[0] = 0;
	occCount[1] = 0;
	occCount[2] = 0;
	occCount[3] = 0;

	iteration = index / 256;
	wordToCount = (index - iteration * 256) / 16;
	charToCount = index - iteration * 256 - wordToCount * 16;

	for (i=0; i<iteration; i++) {

		sum = 0;
		for (j=0; j<16; j++) {
			sum += dnaDecodeTable[*dna >> 16];
			sum += dnaDecodeTable[*dna & 0x0000FFFF];
			dna++;
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
			occCount[0] += sum & 0x000000FF;	sum >>= 8;
			occCount[1] += sum & 0x000000FF;	sum >>= 8;
			occCount[2] += sum & 0x000000FF;	sum >>= 8;
			occCount[3] += sum;
		} else {
			// only some or all of the 3 bits are on
			// in reality, only one of the four cases are possible
			if (sum == 0x00000100) {
				occCount[0] += 256;
			} else if (sum == 0x00010000) {
				occCount[1] += 256;
			} else if (sum == 0x01000000) {
				occCount[2] += 256;
			} else if (sum == 0x00000000) {
				occCount[3] += 256;
			} else {
				fprintf(stderr, "ForwardDNAAllOccCountNoLimit(): DNA occ sum exception!\n");
				exit(1);
			}
		}

	}

	sum = 0;
	for (j=0; j<wordToCount; j++) {
		sum += dnaDecodeTable[*dna >> 16];
		sum += dnaDecodeTable[*dna & 0x0000FFFF];
		dna++;
	}

	if (charToCount > 0) {
		c = *dna & truncateRightMask[charToCount];	// increase count of 'a' by 16 - c;
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0xFFFF];
		sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
	}

	occCount[0] += sum & 0x000000FF;	sum >>= 8;
	occCount[1] += sum & 0x000000FF;	sum >>= 8;
	occCount[2] += sum & 0x000000FF;	sum >>= 8;
	occCount[3] += sum;
}

static void BWTIncBuildPackedBwt(const unsigned int *relativeRank, unsigned int* __restrict bwt, const unsigned int numChar,
								 const uint64_t *cumulativeCount, const unsigned int *packedShift) {

	unsigned int i, c, r;
	unsigned int previousRank, currentRank;
	unsigned int wordIndex, charIndex;
	unsigned int inverseSa0;

	inverseSa0 = previousRank = relativeRank[0];

	for (i=1; i<=numChar; i++) {
		currentRank = relativeRank[i];
		// previousRank > cumulativeCount[c] because $ is one of the char
		c = (previousRank > cumulativeCount[1]) + (previousRank > cumulativeCount[2]) 
											    + (previousRank > cumulativeCount[3]);
		// set bwt for currentRank
		if (c > 0) {
			// c <> 'a'
			r = currentRank;
			if (r > inverseSa0) {
				// - 1 because $ at inverseSa0 is not encoded			
				r--;
			}
			wordIndex = r / CHAR_PER_WORD;
			charIndex = r - wordIndex * CHAR_PER_WORD;
			bwt[wordIndex] |= c << packedShift[charIndex];
		}
		previousRank = currentRank;
	}
}

static inline uint64_t BWTOccValueExplicit(const BWT64 *bwt, const uint64_t occIndexExplicit,
											   const unsigned int character)
{
	uint64_t occIndexMajor;

	occIndexMajor = occIndexExplicit * OCC_INTERVAL / OCC_INTERVAL_MAJOR;

	if (occIndexExplicit % OCC_VALUE_PER_WORD == 0) {
		return bwt->occValueMajor[occIndexMajor * ALPHABET_SIZE + character] +
			   (bwt->occValue[occIndexExplicit / OCC_VALUE_PER_WORD * ALPHABET_SIZE + character] >> 16);

	} else {
		return bwt->occValueMajor[occIndexMajor * ALPHABET_SIZE + character] +
			   (bwt->occValue[occIndexExplicit / OCC_VALUE_PER_WORD * ALPHABET_SIZE + character] & 0x0000FFFF);
	}
}


static unsigned int ForwardDNAOccCount(const unsigned int*  dna, const unsigned int index, const unsigned int character,
									   const unsigned int*  dnaDecodeTable)
{
	static const unsigned int truncateRightMask[16] = { 0x00000000, 0xC0000000, 0xF0000000, 0xFC000000,
											   0xFF000000, 0xFFC00000, 0xFFF00000, 0xFFFC0000,
											   0xFFFF0000, 0xFFFFC000, 0xFFFFF000, 0xFFFFFC00,
											   0xFFFFFF00, 0xFFFFFFC0, 0xFFFFFFF0, 0xFFFFFFFC };

	unsigned int wordToCount, charToCount;
	unsigned int i, c;
	unsigned int sum = 0;

	wordToCount = index / 16;
	charToCount = index - wordToCount * 16;

	for (i=0; i<wordToCount; i++) {
		sum += dnaDecodeTable[dna[i] >> 16];
		sum += dnaDecodeTable[dna[i] & 0x0000FFFF];
	}

	if (charToCount > 0) {
		c = dna[i] & truncateRightMask[charToCount];	// increase count of 'a' by 16 - c;
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0xFFFF];
		sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
	}

	return (sum >> (character * 8)) & 0x000000FF;

}

static unsigned int BackwardDNAOccCount(const unsigned int*  dna, const unsigned int index, const unsigned int character,
										const unsigned int*  dnaDecodeTable)
{
	static const unsigned int truncateLeftMask[16] =  { 0x00000000, 0x00000003, 0x0000000F, 0x0000003F,
											   0x000000FF, 0x000003FF, 0x00000FFF, 0x00003FFF,
											   0x0000FFFF, 0x0003FFFF, 0x000FFFFF, 0x003FFFFF,
											   0x00FFFFFF, 0x03FFFFFF, 0x0FFFFFFF, 0x3FFFFFFF };

	unsigned int wordToCount, charToCount;
	unsigned int i, c;
	unsigned int sum = 0;

	wordToCount = index / 16;
	charToCount = index - wordToCount * 16;

	dna -= wordToCount + 1;

	if (charToCount > 0) {
		c = *dna & truncateLeftMask[charToCount];	// increase count of 'a' by 16 - c;
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0xFFFF];
		sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
	}
	
	for (i=0; i<wordToCount; i++) {
		dna++;
		sum += dnaDecodeTable[*dna >> 16];
		sum += dnaDecodeTable[*dna & 0x0000FFFF];
	}

	return (sum >> (character * 8)) & 0x000000FF;

}

uint64_t BWTOccValue64(const BWT64 *bwt, uint64_t index, const unsigned int character) 
{
	uint64_t occValue;
	uint64_t occExplicitIndex, occIndex;
    int i = 0 ;
	// $ is supposed to be positioned at inverseSa0 but it is not encoded
	// therefore index is subtracted by 0~divideCount for adjustment
    for(i = 0 ; i < bwt->divideNumber; i++)
    {
        if(index <= bwt->sentinelPosition[i]) break ;
    }
    index -= i ;

	occExplicitIndex = (index + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;	// Bidirectional encoding
	occIndex = occExplicitIndex * OCC_INTERVAL;
	occValue = BWTOccValueExplicit(bwt, occExplicitIndex, character);

	if (occIndex == index) {
		return occValue;
	}

	if (occIndex < index) {
		return occValue + ForwardDNAOccCount(bwt->bwtCode + occIndex / CHAR_PER_WORD, index - occIndex, character, bwt->decodeTable);
	} else {
		return occValue - BackwardDNAOccCount(bwt->bwtCode + occIndex / CHAR_PER_WORD, occIndex - index, character, bwt->decodeTable);
	}

}

static unsigned int BWTIncGetAbsoluteRank64(BWT64 *bwt, uint64_t* __restrict absoluteRank, unsigned int* __restrict seq,
										  const unsigned int *packedText, const unsigned int numChar,
										  const uint64_t *cumulativeCount, const unsigned int firstCharInLastIteration)
{
	uint64_t saIndex;
	unsigned int lastWord;
	unsigned int packedMask;
	unsigned int i, j;
	unsigned int c, t;
	unsigned int rankIndex;
	unsigned int shift;
	unsigned int seqIndexFromStart[ALPHABET_SIZE];
	unsigned int seqIndexFromEnd[ALPHABET_SIZE];

	for (i=0; i<ALPHABET_SIZE; i++) {
		seqIndexFromStart[i] = cumulativeCount[i];
		seqIndexFromEnd[i] = cumulativeCount[i+1] - 1;
	}

	shift = BITS_IN_WORD - BIT_PER_CHAR;
	packedMask = ALL_ONE_MASK >> shift;
	saIndex = bwt->inverseSa0;
	rankIndex = numChar - 1;

	lastWord = numChar / CHAR_PER_WORD;
	for (i=lastWord; i--;) {	// loop from lastWord - 1 to 0
		t = packedText[i];
		for (j=0; j<CHAR_PER_WORD; j++) {
			c = t & packedMask;
			saIndex = bwt->cumulativeFreq[c] + BWTOccValue64(bwt, saIndex, c) + 1;
			// A counting sort using the first character of suffix is done here
			// If rank > inverseSa0 -> fill seq from end, otherwise fill seq from start -> to leave the right entry for inverseSa0
			if (saIndex > bwt->inverseSa0) {
				seq[seqIndexFromEnd[c]] = rankIndex;
				absoluteRank[seqIndexFromEnd[c]] = saIndex;
				seqIndexFromEnd[c]--;
			} else {
				seq[seqIndexFromStart[c]] = rankIndex;
				absoluteRank[seqIndexFromStart[c]] = saIndex;
				seqIndexFromStart[c]++;
			}
			rankIndex--;
			t >>= BIT_PER_CHAR;
		}
	}

	absoluteRank[seqIndexFromStart[firstCharInLastIteration]] = bwt->inverseSa0;	// representing the substring of all preceding characters
	seq[seqIndexFromStart[firstCharInLastIteration]] = numChar;

	return seqIndexFromStart[firstCharInLastIteration];
}

//void BWTQSort64(bgint_t* __restrict key, const bgint_t numItem)
void BWTQSort64(bgint_t* key, const bgint_t numItem)
{
	#define EQUAL_KEY_THRESHOLD	4	// Partition for equal key if data array size / the number of data with equal value with pivot < EQUAL_KEY_THRESHOLD

	int64_t lowIndex, highIndex, midIndex;
	int64_t lowPartitionIndex, highPartitionIndex;
    int max_iteration = 32;
	int64_t lowStack[max_iteration], highStack[max_iteration];
	int stackDepth;
	int64_t i, j;
	//bgint_t tempSeq, tempKey;
	bgint_t tempKey;
	int64_t numberOfEqualKey;

	if (numItem < 2) return;

	stackDepth = 0;

    lowIndex = 0;
    highIndex = numItem - 1;

	for (;;) {

		for (;;) {

			// Sort small array of data
			if (highIndex - lowIndex < BWTINC_INSERT_SORT_NUM_ITEM) {	 // Insertion sort on smallest arrays
				for (i=lowIndex+1; i<=highIndex; i++) {
					tempKey = key[i];
					for (j = i; j > lowIndex && key[j-1] > tempKey; j--) {
						key[j] = key[j-1];
					}
					if (j != i) {
						key[j] = tempKey;
					}
				}
				break;
			}

			// Choose pivot as median of the lowest, middle, and highest data; sort the three data

			midIndex = average(lowIndex, highIndex);
			if (key[lowIndex] > key[midIndex]) {
				tempKey = key[lowIndex];
				key[lowIndex] = key[midIndex];
				key[midIndex] = tempKey;
			}
			if (key[lowIndex] > key[highIndex]) {
				tempKey = key[lowIndex];
				key[lowIndex] = key[highIndex];
				key[highIndex] = tempKey;
			}
			if (key[midIndex] > key[highIndex]) {
				tempKey = key[midIndex];
				key[midIndex] = key[highIndex];
				key[highIndex] = tempKey;
			}

			// Partition data

			numberOfEqualKey = 0;

			lowPartitionIndex = lowIndex + 1;
			highPartitionIndex = highIndex - 1;

			for (;;) {
				while (lowPartitionIndex <= highPartitionIndex && key[lowPartitionIndex] <= key[midIndex]) {
					numberOfEqualKey += (key[lowPartitionIndex] == key[midIndex]);
					lowPartitionIndex++;
				}
				while (lowPartitionIndex < highPartitionIndex) {
					if (key[midIndex] >= key[highPartitionIndex]) {
						numberOfEqualKey += (key[midIndex] == key[highPartitionIndex]);
						break;
					}
					highPartitionIndex--;
				}
				if (lowPartitionIndex >= highPartitionIndex) {
					break;
				}
				tempKey = key[lowPartitionIndex];
				key[lowPartitionIndex] = key[highPartitionIndex];
				key[highPartitionIndex] = tempKey;
				if (highPartitionIndex == midIndex) {
					// partition key has been moved
					midIndex = lowPartitionIndex;
				}
				lowPartitionIndex++;
				highPartitionIndex--;
			}

			// Adjust the partition index
			highPartitionIndex = lowPartitionIndex;
			lowPartitionIndex--;

			// move the partition key to end of low partition
			tempKey = key[midIndex];
			key[midIndex] = key[lowPartitionIndex];
			key[lowPartitionIndex] = tempKey;

			if (highIndex - lowIndex + BWTINC_INSERT_SORT_NUM_ITEM <= EQUAL_KEY_THRESHOLD * numberOfEqualKey) {

				// Many keys = partition key; separate the equal key data from the lower partition
		
				midIndex = lowIndex;

				for (;;) {
					while (midIndex < lowPartitionIndex && key[midIndex] < key[lowPartitionIndex]) {
						midIndex++;
					}
					while (midIndex < lowPartitionIndex && key[lowPartitionIndex] == key[lowPartitionIndex - 1]) {
						lowPartitionIndex--;
					}
					if (midIndex >= lowPartitionIndex) {
						break;
					}
					tempKey = key[midIndex];
					key[midIndex] = key[lowPartitionIndex - 1];
					key[lowPartitionIndex - 1] = tempKey;
					midIndex++;
					lowPartitionIndex--;
				}

			}

			if (lowPartitionIndex - lowIndex > highIndex - highPartitionIndex) {
				// put the larger partition to stack
				lowStack[stackDepth] = lowIndex;
				highStack[stackDepth] = lowPartitionIndex - 1;
				stackDepth++;
				// sort the smaller partition first
				lowIndex = highPartitionIndex;
			} else {
				// put the larger partition to stack
				lowStack[stackDepth] = highPartitionIndex;
				highStack[stackDepth] = highIndex;
				stackDepth++;
				// sort the smaller partition first
				if (lowPartitionIndex > lowIndex) {
					highIndex = lowPartitionIndex - 1;
				} else {
					// all keys in the partition equals to the partition key
					break;
				}
			}
            if(stackDepth >= max_iteration)
            {
                fprintf(stderr, "[BWTQSort] stackDepth > max_iteration, exit....\n");
                exit(1);
            }
			continue;
		}

		// Pop a range from stack
		if (stackDepth > 0) {
			stackDepth--;
			lowIndex = lowStack[stackDepth];
			highIndex = highStack[stackDepth];
			continue;
		} else return;
	}
}

static void BWTIncSortKey64(uint64_t* __restrict key, unsigned int* __restrict seq, const unsigned int numItem)
{
	#define EQUAL_KEY_THRESHOLD	4	// Partition for equal key if data array size / the number of data with equal value with pivot < EQUAL_KEY_THRESHOLD

	int lowIndex, highIndex, midIndex;
	int lowPartitionIndex, highPartitionIndex;
	int lowStack[32], highStack[32];
	int stackDepth;
	int i, j;
	unsigned int tempSeq;
    uint64_t tempKey;
	int numberOfEqualKey;

	if (numItem < 2) return;

	stackDepth = 0;

    lowIndex = 0;
    highIndex = numItem - 1;

	for (;;) {

		for (;;) {

			// Sort small array of data
			if (highIndex - lowIndex < BWTINC_INSERT_SORT_NUM_ITEM) {	 // Insertion sort on smallest arrays
				for (i=lowIndex+1; i<=highIndex; i++) {
					tempSeq = seq[i];
					tempKey = key[i];
					for (j = i; j > lowIndex && key[j-1] > tempKey; j--) {
						seq[j] = seq[j-1];
						key[j] = key[j-1];
					}
					if (j != i) {
						seq[j] = tempSeq;
						key[j] = tempKey;
					}
				}
				break;
			}

			// Choose pivot as median of the lowest, middle, and highest data; sort the three data

			midIndex = average(lowIndex, highIndex);
			if (key[lowIndex] > key[midIndex]) {
				tempSeq = seq[lowIndex];
				tempKey = key[lowIndex];
				seq[lowIndex] = seq[midIndex];
				key[lowIndex] = key[midIndex];
				seq[midIndex] = tempSeq;
				key[midIndex] = tempKey;
			}
			if (key[lowIndex] > key[highIndex]) {
				tempSeq = seq[lowIndex];
				tempKey = key[lowIndex];
				seq[lowIndex] = seq[highIndex];
				key[lowIndex] = key[highIndex];
				seq[highIndex] = tempSeq;
				key[highIndex] = tempKey;
			}
			if (key[midIndex] > key[highIndex]) {
				tempSeq = seq[midIndex];
				tempKey = key[midIndex];
				seq[midIndex] = seq[highIndex];
				key[midIndex] = key[highIndex];
				seq[highIndex] = tempSeq;
				key[highIndex] = tempKey;
			}

			// Partition data

			numberOfEqualKey = 0;

			lowPartitionIndex = lowIndex + 1;
			highPartitionIndex = highIndex - 1;

			for (;;) {
				while (lowPartitionIndex <= highPartitionIndex && key[lowPartitionIndex] <= key[midIndex]) {
					numberOfEqualKey += (key[lowPartitionIndex] == key[midIndex]);
					lowPartitionIndex++;
				}
				while (lowPartitionIndex < highPartitionIndex) {
					if (key[midIndex] >= key[highPartitionIndex]) {
						numberOfEqualKey += (key[midIndex] == key[highPartitionIndex]);
						break;
					}
					highPartitionIndex--;
				}
				if (lowPartitionIndex >= highPartitionIndex) {
					break;
				}
				tempSeq = seq[lowPartitionIndex];
				tempKey = key[lowPartitionIndex];
				seq[lowPartitionIndex] = seq[highPartitionIndex];
				key[lowPartitionIndex] = key[highPartitionIndex];
				seq[highPartitionIndex] = tempSeq;
				key[highPartitionIndex] = tempKey;
				if (highPartitionIndex == midIndex) {
					// partition key has been moved
					midIndex = lowPartitionIndex;
				}
				lowPartitionIndex++;
				highPartitionIndex--;
			}

			// Adjust the partition index
			highPartitionIndex = lowPartitionIndex;
			lowPartitionIndex--;

			// move the partition key to end of low partition
			tempSeq = seq[midIndex];
			tempKey = key[midIndex];
			seq[midIndex] = seq[lowPartitionIndex];
			key[midIndex] = key[lowPartitionIndex];
			seq[lowPartitionIndex] = tempSeq;
			key[lowPartitionIndex] = tempKey;

			if (highIndex - lowIndex + BWTINC_INSERT_SORT_NUM_ITEM <= EQUAL_KEY_THRESHOLD * numberOfEqualKey) {

				// Many keys = partition key; separate the equal key data from the lower partition
		
				midIndex = lowIndex;

				for (;;) {
					while (midIndex < lowPartitionIndex && key[midIndex] < key[lowPartitionIndex]) {
						midIndex++;
					}
					while (midIndex < lowPartitionIndex && key[lowPartitionIndex] == key[lowPartitionIndex - 1]) {
						lowPartitionIndex--;
					}
					if (midIndex >= lowPartitionIndex) {
						break;
					}
					tempSeq = seq[midIndex];
					tempKey = key[midIndex];
					seq[midIndex] = seq[lowPartitionIndex - 1];
					key[midIndex] = key[lowPartitionIndex - 1];
					seq[lowPartitionIndex - 1] = tempSeq;
					key[lowPartitionIndex - 1] = tempKey;
					midIndex++;
					lowPartitionIndex--;
				}

			}

			if (lowPartitionIndex - lowIndex > highIndex - highPartitionIndex) {
				// put the larger partition to stack
				lowStack[stackDepth] = lowIndex;
				highStack[stackDepth] = lowPartitionIndex - 1;
				stackDepth++;
				// sort the smaller partition first
				lowIndex = highPartitionIndex;
			} else {
				// put the larger partition to stack
				lowStack[stackDepth] = highPartitionIndex;
				highStack[stackDepth] = highIndex;
				stackDepth++;
				// sort the smaller partition first
				if (lowPartitionIndex > lowIndex) {
					highIndex = lowPartitionIndex - 1;
				} else {
					// all keys in the partition equals to the partition key
					break;
				}
			}
			continue;
		}

		// Pop a range from stack
		if (stackDepth > 0) {
			stackDepth--;
			lowIndex = lowStack[stackDepth];
			highIndex = highStack[stackDepth];
			continue;
		} else return;
	}
}


static void BWTIncBuildRelativeRank64(uint64_t* __restrict sortedRank, unsigned int* __restrict seq,
									unsigned int* __restrict relativeRank, const unsigned int numItem,
									uint64_t oldInverseSa0, const uint64_t *cumulativeCount)
{
	unsigned int i, c;
	unsigned int s;
    uint64_t r;
    uint64_t lastRank ;
	unsigned int  lastIndex;
	unsigned int oldInverseSa0RelativeRank = 0;
	unsigned int freq;

	lastIndex = numItem;
	lastRank = sortedRank[numItem];
	if (lastRank > oldInverseSa0) {
		sortedRank[numItem]--;	// to prepare for merging; $ is not encoded in bwt
	}
	s = seq[numItem];
	relativeRank[s] = numItem;
	if (lastRank == oldInverseSa0) {
		oldInverseSa0RelativeRank = numItem;
		oldInverseSa0++;	// so that this segment of code is not run again
		lastRank++;			// so that oldInverseSa0 become a sorted group with 1 item
	}

	c = ALPHABET_SIZE - 1;
	freq = cumulativeCount[c];

	for (i=numItem; i--;) {	// from numItem - 1 to 0
		r = sortedRank[i];
		if (r > oldInverseSa0) {
			sortedRank[i]--;	// to prepare for merging; $ is not encoded in bwt
		}
		s = seq[i];
		if (i < freq) {
			if (lastIndex >= freq) {
				lastRank++;	// to trigger the group across alphabet boundary to be split
			}
			c--;
			freq = cumulativeCount[c];
		}
		if (r == lastRank) {
			relativeRank[s] = lastIndex;
		} else {
			if (i == lastIndex - 1) {
				if (lastIndex < numItem && (int)seq[lastIndex + 1] < 0) {
					seq[lastIndex] = seq[lastIndex + 1] - 1;
				} else {
					seq[lastIndex] = (unsigned int)-1;
				}
			}
			lastIndex = i;
			lastRank = r;
			relativeRank[s] = i;
			if (r == oldInverseSa0) {
				oldInverseSa0RelativeRank = i;
				oldInverseSa0++;	// so that this segment of code is not run again
				lastRank++;			// so that oldInverseSa0 become a sorted group with 1 item
			}
		}
	}

}

static void BWTIncBuildBwt(unsigned int*  seq, const unsigned int *relativeRank, const unsigned int numChar,
						   const uint64_t *cumulativeCount)
{
	unsigned int i, c;
	unsigned int previousRank, currentRank;

	previousRank = relativeRank[0];

	for (i=1; i<=numChar; i++) {
		currentRank = relativeRank[i];
		c = (previousRank >= cumulativeCount[1]) + (previousRank >= cumulativeCount[2])
											  	 + (previousRank >= cumulativeCount[3]);
		seq[currentRank] = c;
		previousRank = currentRank;
	}
}

static void BWTIncMergeBwt64(const uint64_t *sortedRank, const unsigned int* oldBwt, const unsigned int *insertBwt,
						   unsigned int* __restrict mergedBwt, const uint64_t numOldBwt, const unsigned int numInsertBwt)
{
	unsigned int bitsInWordMinusBitPerChar;
	uint64_t leftShift, rightShift;
	uint64_t o;
	uint64_t oIndex, iIndex, mIndex;
	uint64_t mWord, mChar, oWord, oChar;
	uint64_t numInsert;

	bitsInWordMinusBitPerChar = BITS_IN_WORD - BIT_PER_CHAR;

	oIndex = 0;
	iIndex = 0;
	mIndex = 0;

	mWord = 0;
	mChar = 0;

	mergedBwt[0] = 0;	// this can be cleared as merged Bwt slightly shift to the left in each iteration

	while (oIndex < numOldBwt) {

		// copy from insertBwt
		while (iIndex <= numInsertBwt && sortedRank[iIndex] <= oIndex) {
			if (sortedRank[iIndex] != 0) {	// special value to indicate that this is for new inverseSa0
				mergedBwt[mWord] |= insertBwt[iIndex] << (BITS_IN_WORD - (mChar + 1) * BIT_PER_CHAR);
				mIndex++;
				mChar++;
				if (mChar == CHAR_PER_WORD) {
					mChar = 0;
					mWord++;
					mergedBwt[mWord] = 0;	// no need to worry about crossing mergedBwt boundary
				}
			}
			iIndex++;
		}

		// Copy from oldBwt to mergedBwt
		if (iIndex <= numInsertBwt) {
			o = sortedRank[iIndex];
		} else {
			o = numOldBwt;
		}
		numInsert = o - oIndex;

		oWord = oIndex / CHAR_PER_WORD;
		oChar = oIndex - oWord * CHAR_PER_WORD;
		if (oChar > mChar) {
			leftShift = (oChar - mChar) * BIT_PER_CHAR;
			rightShift = (CHAR_PER_WORD + mChar - oChar) * BIT_PER_CHAR;
			mergedBwt[mWord] = mergedBwt[mWord]
								| (oldBwt[oWord] << (oChar * BIT_PER_CHAR) >> (mChar * BIT_PER_CHAR))
								| (oldBwt[oWord+1] >> rightShift);
			oIndex += min(numInsert, CHAR_PER_WORD - mChar);
			while (o > oIndex) {
				oWord++;
				mWord++;
				mergedBwt[mWord] = (oldBwt[oWord] << leftShift) | (oldBwt[oWord+1] >> rightShift);
				oIndex += CHAR_PER_WORD;
			}
		} else if (oChar < mChar) {
			rightShift = (mChar - oChar) * BIT_PER_CHAR;
			leftShift = (CHAR_PER_WORD + oChar - mChar) * BIT_PER_CHAR;
			mergedBwt[mWord] = mergedBwt[mWord] 
								| (oldBwt[oWord] << (oChar * BIT_PER_CHAR) >> (mChar * BIT_PER_CHAR));
			oIndex += min(numInsert, CHAR_PER_WORD - mChar);
			while (o > oIndex) {
				oWord++;
				mWord++;
				mergedBwt[mWord] = (oldBwt[oWord-1] << leftShift) | (oldBwt[oWord] >> rightShift);
				oIndex += CHAR_PER_WORD;
			}
		} else { // oChar == mChar
			mergedBwt[mWord] = mergedBwt[mWord] | truncateLeft(oldBwt[oWord], mChar * BIT_PER_CHAR);
			oIndex += min(numInsert, CHAR_PER_WORD - mChar);
			while (o > oIndex) {
				oWord++;
				mWord++;
				mergedBwt[mWord] = oldBwt[oWord];
				oIndex += CHAR_PER_WORD;
			}
		}
		oIndex = o;
		mIndex += numInsert;

		// Clear the trailing garbage in mergedBwt
		mWord = mIndex / CHAR_PER_WORD;
		mChar = mIndex - mWord * CHAR_PER_WORD;
		if (mChar == 0) {
			mergedBwt[mWord] = 0;
		} else {
			mergedBwt[mWord] = truncateRight(mergedBwt[mWord], (BITS_IN_WORD - mChar * BIT_PER_CHAR));
		}

	}

	// copy from insertBwt
	while (iIndex <= numInsertBwt) {
		if (sortedRank[iIndex] != 0) {
			mergedBwt[mWord] |= insertBwt[iIndex] << (BITS_IN_WORD - (mChar + 1) * BIT_PER_CHAR);
			mIndex++;
			mChar++;
			if (mChar == CHAR_PER_WORD) {
				mChar = 0;
				mWord++;
				mergedBwt[mWord] = 0;	// no need to worry about crossing mergedBwt boundary
			}
		}
		iIndex++;
	}
}

void BWTClearTrailingBwtCode64(BWT64 *bwt)
{
	uint64_t bwtResidentSizeInWord;
	uint64_t wordIndex, offset;
	uint64_t i;

	bwtResidentSizeInWord = BWTResidentSizeInWord64(bwt->textLength);

	wordIndex = bwt->textLength / CHAR_PER_WORD;
	offset = (bwt->textLength - wordIndex * CHAR_PER_WORD) * BIT_PER_CHAR;
	if (offset > 0) {
		bwt->bwtCode[wordIndex] = truncateRight(bwt->bwtCode[wordIndex], BITS_IN_WORD - offset);
	} else {
		if (wordIndex < bwtResidentSizeInWord) {
			bwt->bwtCode[wordIndex] = 0;
		}
	}

	for (i=wordIndex+1; i<bwtResidentSizeInWord; i++) {
		bwt->bwtCode[i] = 0;
	}
}


void BWTGenerateOccValueFromBwt64(const unsigned int*  bwt, unsigned int* __restrict occValue,
								uint64_t* __restrict occValueMajor,
								const uint64_t textLength, const unsigned int*  decodeTable)
{
	uint64_t numberOfOccValueMajor, numberOfOccValue;
	unsigned int wordBetweenOccValue;
	unsigned int numberOfOccIntervalPerMajor;
	unsigned int c;
	uint64_t i, j;
	uint64_t occMajorIndex;
	uint64_t occIndex, bwtIndex;
	uint64_t sum;
	uint64_t tempOccValue0[ALPHABET_SIZE], tempOccValue1[ALPHABET_SIZE];

	wordBetweenOccValue = OCC_INTERVAL / CHAR_PER_WORD;

	// Calculate occValue
	// [lh3] by default: OCC_INTERVAL_MAJOR=65536, OCC_INTERVAL=256
	numberOfOccValue = (textLength + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;				// Value at both end for bi-directional encoding
	numberOfOccIntervalPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
	numberOfOccValueMajor = (numberOfOccValue + numberOfOccIntervalPerMajor - 1) / numberOfOccIntervalPerMajor;

	tempOccValue0[0] = 0;
	tempOccValue0[1] = 0;
	tempOccValue0[2] = 0;
	tempOccValue0[3] = 0;
	occValueMajor[0] = 0;
	occValueMajor[1] = 0;
	occValueMajor[2] = 0;
	occValueMajor[3] = 0;

	occIndex = 0;
	bwtIndex = 0;
	for (occMajorIndex=1; occMajorIndex<numberOfOccValueMajor; occMajorIndex++) {

		for (i=0; i<numberOfOccIntervalPerMajor/2; i++) {

			sum = 0;
			tempOccValue1[0] = tempOccValue0[0];
			tempOccValue1[1] = tempOccValue0[1];
			tempOccValue1[2] = tempOccValue0[2];
			tempOccValue1[3] = tempOccValue0[3];

			for (j=0; j<wordBetweenOccValue; j++) {
				c = bwt[bwtIndex];
				sum += decodeTable[c >> 16];
				sum += decodeTable[c & 0x0000FFFF];
				bwtIndex++;
			}
			if (!DNA_OCC_SUM_EXCEPTION(sum)) {
				tempOccValue1[0] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue1[1] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue1[2] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue1[3] += sum;
			} else {
				if (sum == 0x00000100) {
					tempOccValue1[0] += 256;
				} else if (sum == 0x00010000) {
					tempOccValue1[1] += 256;
				} else if (sum == 0x01000000) {
					tempOccValue1[2] += 256;
				} else {
					tempOccValue1[3] += 256;
				}
			}
			occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
			occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
			occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
			occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];
			tempOccValue0[0] = tempOccValue1[0];
			tempOccValue0[1] = tempOccValue1[1];
			tempOccValue0[2] = tempOccValue1[2];
			tempOccValue0[3] = tempOccValue1[3];
			sum = 0;

			occIndex++;

			for (j=0; j<wordBetweenOccValue; j++) {
				c = bwt[bwtIndex];
				sum += decodeTable[c >> 16];
				sum += decodeTable[c & 0x0000FFFF];
				bwtIndex++;
			}
			if (!DNA_OCC_SUM_EXCEPTION(sum)) {
				tempOccValue0[0] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue0[1] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue0[2] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue0[3] += sum;
			} else {
				if (sum == 0x00000100) {
					tempOccValue0[0] += 256;
				} else if (sum == 0x00010000) {
					tempOccValue0[1] += 256;
				} else if (sum == 0x01000000) {
					tempOccValue0[2] += 256;
				} else {
					tempOccValue0[3] += 256;
				}
			}
		}

		occValueMajor[occMajorIndex * 4 + 0] = occValueMajor[(occMajorIndex - 1) * 4 + 0] + tempOccValue0[0];
		occValueMajor[occMajorIndex * 4 + 1] = occValueMajor[(occMajorIndex - 1) * 4 + 1] + tempOccValue0[1];
		occValueMajor[occMajorIndex * 4 + 2] = occValueMajor[(occMajorIndex - 1) * 4 + 2] + tempOccValue0[2];
		occValueMajor[occMajorIndex * 4 + 3] = occValueMajor[(occMajorIndex - 1) * 4 + 3] + tempOccValue0[3];
		tempOccValue0[0] = 0;
		tempOccValue0[1] = 0;
		tempOccValue0[2] = 0;
		tempOccValue0[3] = 0;

	}

	while (occIndex < (numberOfOccValue-1)/2) {
		sum = 0;
		tempOccValue1[0] = tempOccValue0[0];
		tempOccValue1[1] = tempOccValue0[1];
		tempOccValue1[2] = tempOccValue0[2];
		tempOccValue1[3] = tempOccValue0[3];
		for (j=0; j<wordBetweenOccValue; j++) {
			c = bwt[bwtIndex];
			sum += decodeTable[c >> 16];
			sum += decodeTable[c & 0x0000FFFF];
			bwtIndex++;
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
			tempOccValue1[0] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[1] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[2] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[3] += sum;
		} else {
			if (sum == 0x00000100) {
				tempOccValue1[0] += 256;
			} else if (sum == 0x00010000) {
				tempOccValue1[1] += 256;
			} else if (sum == 0x01000000) {
				tempOccValue1[2] += 256;
			} else {
				tempOccValue1[3] += 256;
			}
		}
		occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
		occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
		occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
		occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];
		tempOccValue0[0] = tempOccValue1[0];
		tempOccValue0[1] = tempOccValue1[1];
		tempOccValue0[2] = tempOccValue1[2];
		tempOccValue0[3] = tempOccValue1[3];
		sum = 0;
		occIndex++;

		for (j=0; j<wordBetweenOccValue; j++) {
			c = bwt[bwtIndex];
			sum += decodeTable[c >> 16];
			sum += decodeTable[c & 0x0000FFFF];
			bwtIndex++;
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
			tempOccValue0[0] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue0[1] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue0[2] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue0[3] += sum;
		} else {
			if (sum == 0x00000100) {
				tempOccValue0[0] += 256;
			} else if (sum == 0x00010000) {
				tempOccValue0[1] += 256;
			} else if (sum == 0x01000000) {
				tempOccValue0[2] += 256;
			} else {
				tempOccValue0[3] += 256;
			}
		}
	}

	sum = 0;
	tempOccValue1[0] = tempOccValue0[0];
	tempOccValue1[1] = tempOccValue0[1];
	tempOccValue1[2] = tempOccValue0[2];
	tempOccValue1[3] = tempOccValue0[3];

	if (occIndex * 2 < numberOfOccValue - 1) {
		for (j=0; j<wordBetweenOccValue; j++) {
			c = bwt[bwtIndex];
			sum += decodeTable[c >> 16];
			sum += decodeTable[c & 0x0000FFFF];
			bwtIndex++;
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
			tempOccValue1[0] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[1] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[2] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[3] += sum;
		} else {
			if (sum == 0x00000100) {
				tempOccValue1[0] += 256;
			} else if (sum == 0x00010000) {
				tempOccValue1[1] += 256;
			} else if (sum == 0x01000000) {
				tempOccValue1[2] += 256;
			} else {
				tempOccValue1[3] += 256;
			}
		}
	}

	occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
	occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
	occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
	occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];

}

static void BWTIncConstruct64(BWTInc64 *bwtInc, const uint64_t numChar)
{
	uint64_t i;
	uint64_t mergedBwtSizeInWord, mergedOccSizeInWord;
	unsigned int firstCharInThisIteration;

	uint64_t *sortedRank;
    unsigned int *seq, *relativeRank, *insertBwt, *mergedBwt;
	uint64_t newInverseSa0RelativeRank, oldInverseSa0RelativeRank, newInverseSa0;

	#ifdef DEBUG
	if (numChar > bwtInc->buildSize) {
		fprintf(stderr, "BWTIncConstruct64(): numChar > buildSize!\n");
		exit(1);
	}
	#endif

	mergedBwtSizeInWord = BWTResidentSizeInWord64(bwtInc->bwt->textLength + numChar);
	mergedOccSizeInWord = BWTOccValueMinorSizeInWord64(bwtInc->bwt->textLength + numChar);

	initializeVAL64(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

	if (bwtInc->bwt->textLength == 0) {		// Initial build

		// Set address
		seq = bwtInc->workingMemory;
		relativeRank = bwtInc->workingMemory + bwtInc->buildSize + 1;
		mergedBwt = insertBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord;	// build in place

		BWTIncPutPackedTextToRank(bwtInc->packedText, relativeRank, bwtInc->cumulativeCountInCurrentBuild, numChar);

		firstCharInThisIteration = relativeRank[0];
		relativeRank[numChar] = 0;

		// Sort suffix
		QSufSortSuffixSort((int*)relativeRank, (int*)seq, (int)numChar, (int)ALPHABET_SIZE - 1, 0, FALSE);
		newInverseSa0 = relativeRank[0];

		// Clear BWT area
		initializeVAL(insertBwt, mergedBwtSizeInWord, 0);

		// Build BWT
		BWTIncBuildPackedBwt(relativeRank, insertBwt, numChar, bwtInc->cumulativeCountInCurrentBuild, bwtInc->packedShift);

		// so that the cumulativeCount is not deducted
		bwtInc->firstCharInLastIteration = ALPHABET_SIZE;

	} else {		// Incremental build
		// Set address
		sortedRank = (uint64_t*)bwtInc->workingMemory;
		seq = bwtInc->workingMemory + 2 * (bwtInc->buildSize + 1);
		insertBwt = seq;
		relativeRank = seq + bwtInc->buildSize + 1;

		// Store the first character of this iteration
		firstCharInThisIteration = bwtInc->packedText[0] >> (BITS_IN_WORD - BIT_PER_CHAR);

		// Count occurrence of input text
		ForwardDNAAllOccCountNoLimit(bwtInc->packedText, numChar, bwtInc->cumulativeCountInCurrentBuild + 1, bwtInc->bwt->decodeTable);
		// Add the first character of the previous iteration to represent the inverseSa0 of the previous iteration
		bwtInc->cumulativeCountInCurrentBuild[bwtInc->firstCharInLastIteration + 1]++;
		bwtInc->cumulativeCountInCurrentBuild[2] += bwtInc->cumulativeCountInCurrentBuild[1];
		bwtInc->cumulativeCountInCurrentBuild[3] += bwtInc->cumulativeCountInCurrentBuild[2];
		bwtInc->cumulativeCountInCurrentBuild[4] += bwtInc->cumulativeCountInCurrentBuild[3];

		// Get rank of new suffix among processed suffix
		// The seq array is built into ALPHABET_SIZE + 2 groups; ALPHABET_SIZE groups + 1 group divided into 2 by inverseSa0 + inverseSa0 as 1 group
		oldInverseSa0RelativeRank = BWTIncGetAbsoluteRank64(bwtInc->bwt, sortedRank, seq, bwtInc->packedText, 
														  numChar, bwtInc->cumulativeCountInCurrentBuild, bwtInc->firstCharInLastIteration);

		// Sort rank by ALPHABET_SIZE + 2 groups (or ALPHABET_SIZE + 1 groups when inverseSa0 sit on the border of a group)
		for (i=0; i<ALPHABET_SIZE; i++) {
			if (bwtInc->cumulativeCountInCurrentBuild[i] > oldInverseSa0RelativeRank ||
				bwtInc->cumulativeCountInCurrentBuild[i+1] <= oldInverseSa0RelativeRank) {
				BWTIncSortKey64(sortedRank + bwtInc->cumulativeCountInCurrentBuild[i], seq + bwtInc->cumulativeCountInCurrentBuild[i], bwtInc->cumulativeCountInCurrentBuild[i+1] - bwtInc->cumulativeCountInCurrentBuild[i]);
			} else {
				if (bwtInc->cumulativeCountInCurrentBuild[i] < oldInverseSa0RelativeRank) {
					BWTIncSortKey64(sortedRank + bwtInc->cumulativeCountInCurrentBuild[i], seq + bwtInc->cumulativeCountInCurrentBuild[i], oldInverseSa0RelativeRank - bwtInc->cumulativeCountInCurrentBuild[i]);
				}
				if (bwtInc->cumulativeCountInCurrentBuild[i+1] > oldInverseSa0RelativeRank + 1) {
					BWTIncSortKey64(sortedRank + oldInverseSa0RelativeRank + 1, seq + oldInverseSa0RelativeRank + 1, bwtInc->cumulativeCountInCurrentBuild[i+1] - oldInverseSa0RelativeRank - 1);
				}
			}
		}

		// build relative rank; sortedRank is updated for merging to cater for the fact that $ is not encoded in bwt
		// the cumulative freq information is used to make sure that inverseSa0 and suffix beginning with different characters are kept in different unsorted groups)
		BWTIncBuildRelativeRank64(sortedRank, seq, relativeRank, numChar, bwtInc->bwt->inverseSa0, bwtInc->cumulativeCountInCurrentBuild);
#ifdef DEBUG
		if (relativeRank[numChar] != oldInverseSa0RelativeRank) {
			fprintf(stderr, "BWTIncConstruct(): relativeRank[numChar] != oldInverseSa0RelativeRank!\n");
			exit(1);
		}
#endif

		// Sort suffix
		QSufSortSuffixSort((int*)relativeRank, (int*)seq, (int)numChar, (int)numChar, 1, TRUE);

		newInverseSa0RelativeRank = relativeRank[0];
		newInverseSa0 = sortedRank[newInverseSa0RelativeRank] + newInverseSa0RelativeRank;

		sortedRank[newInverseSa0RelativeRank] = 0;	// a special value so that this is skipped in the merged bwt

		// Build BWT
		BWTIncBuildBwt(seq, relativeRank, numChar, bwtInc->cumulativeCountInCurrentBuild);

		// Merge BWT
		mergedBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord 
				    - bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR;
					// minus numberOfIteration * occInterval to create a buffer for merging
		BWTIncMergeBwt64(sortedRank, bwtInc->bwt->bwtCode, insertBwt, mergedBwt, bwtInc->bwt->textLength, numChar);

	}

	// Build auxiliary structure and update info and pointers in BWT
	bwtInc->bwt->textLength += numChar;
	bwtInc->bwt->bwtCode = mergedBwt;
	bwtInc->bwt->bwtSizeInWord = mergedBwtSizeInWord;
	bwtInc->bwt->occSizeInWord = mergedOccSizeInWord;
	if (mergedBwt < bwtInc->workingMemory + mergedOccSizeInWord) {
		fprintf(stderr, "BWTIncConstruct() : Not enough memory allocated!\n");
		exit(1);
	}

	bwtInc->bwt->occValue = mergedBwt - mergedOccSizeInWord;

	BWTClearTrailingBwtCode64(bwtInc->bwt);
	BWTGenerateOccValueFromBwt64(bwtInc->bwt->bwtCode, bwtInc->bwt->occValue, bwtInc->bwt->occValueMajor,
							   bwtInc->bwt->textLength, bwtInc->bwt->decodeTable);

	bwtInc->bwt->inverseSa0 = newInverseSa0;
	
	bwtInc->bwt->cumulativeFreq[1] += bwtInc->cumulativeCountInCurrentBuild[1] - (bwtInc->firstCharInLastIteration <= 0);
	bwtInc->bwt->cumulativeFreq[2] += bwtInc->cumulativeCountInCurrentBuild[2] - (bwtInc->firstCharInLastIteration <= 1);
	bwtInc->bwt->cumulativeFreq[3] += bwtInc->cumulativeCountInCurrentBuild[3] - (bwtInc->firstCharInLastIteration <= 2);
	bwtInc->bwt->cumulativeFreq[4] += bwtInc->cumulativeCountInCurrentBuild[4] - (bwtInc->firstCharInLastIteration <= 3);

	bwtInc->firstCharInLastIteration = firstCharInThisIteration;

	// Set build size and text address for the next build
	BWTIncSetBuildSizeAndTextAddr64(bwtInc);
	bwtInc->numberOfIterationDone++;

}

BWTInc64 *BWTIncConstructFromPacked64(const char *inputFileName, const float targetNBit,
								  const uint64_t initialMaxBuildSize, const uint64_t incMaxBuildSize)
{

	FILE *packedFile;
	int64_t packedFileLen;
	uint64_t totalTextLength;
	uint64_t textToLoad, textSizeInByte;
	uint64_t processedTextLength;
	unsigned char lastByteLength;

	BWTInc64 *bwtInc;

	packedFile = (FILE*)fopen(inputFileName, "rb");

	if (packedFile == NULL) {
		fprintf(stderr, "BWTIncConstructFromPacked64() : Cannot open inputFileName!\n");
		exit(1);
	}

	fseek(packedFile, -1, SEEK_END);
	packedFileLen = ftell(packedFile);
	if ((int64_t)packedFileLen < 0) {
		fprintf(stderr, "BWTIncConstructFromPacked64: Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, packedFile);
	totalTextLength = TextLengthFromBytePacked64(packedFileLen, BIT_PER_CHAR, lastByteLength);

	bwtInc = BWTIncCreate64(totalTextLength, targetNBit, initialMaxBuildSize, incMaxBuildSize);

	BWTIncSetBuildSizeAndTextAddr64(bwtInc);

	if (bwtInc->buildSize > totalTextLength) {
		textToLoad = totalTextLength;
	} else {
		textToLoad = totalTextLength - ((totalTextLength - bwtInc->buildSize + CHAR_PER_WORD - 1) / CHAR_PER_WORD * CHAR_PER_WORD);
	}
	textSizeInByte = textToLoad / CHAR_PER_BYTE;	// excluded the odd byte

	fseek(packedFile, -2, SEEK_CUR);
	fseek(packedFile, -((int64_t)textSizeInByte), SEEK_CUR);
	fread(bwtInc->textBuffer, sizeof(unsigned char), textSizeInByte + 1, packedFile);
	fseek(packedFile, -((int64_t)textSizeInByte + 1), SEEK_CUR);

	ConvertBytePackedToWordPacked64(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
	BWTIncConstruct64(bwtInc, textToLoad);

	processedTextLength = textToLoad;

	while (processedTextLength < totalTextLength) {
		textToLoad = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;
		if (textToLoad > totalTextLength - processedTextLength) {
			textToLoad = totalTextLength - processedTextLength;
		}
		textSizeInByte = textToLoad / CHAR_PER_BYTE;
		fseek(packedFile, -((int64_t)textSizeInByte), SEEK_CUR);
		fread(bwtInc->textBuffer, sizeof(unsigned char), textSizeInByte, packedFile);
		fseek(packedFile, -((int64_t)textSizeInByte), SEEK_CUR);
		ConvertBytePackedToWordPacked64(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
		BWTIncConstruct64(bwtInc, textToLoad);
		processedTextLength += textToLoad;
		if (bwtInc->numberOfIterationDone % 10 == 0) {
			printf("[BWTIncConstructFromPacked] %u iterations done. %lu characters processed.\n",
				   bwtInc->numberOfIterationDone, processedTextLength);
		}
	}
	return bwtInc;
}

void BWTFree64(BWT64 *bwt)
{
	if (bwt == 0) return;
	free(bwt->cumulativeFreq);
	free(bwt->bwtCode);
	free(bwt->occValue);
	free(bwt->occValueMajor);
	free(bwt->saValue);
	free(bwt->inverseSa);
	free(bwt->decodeTable);
	free(bwt->saIndexRange);
	free(bwt->saValueOnBoundary);
    free(bwt->divideBoundary);
    free(bwt->sentinelPosition);
	free(bwt);
}

void BWTIncFree64(BWTInc64 *bwtInc)
{
	if (bwtInc == 0) return;
    free(bwtInc->bwt->divideBoundary);
    free(bwtInc->bwt->sentinelPosition);
	free(bwtInc->bwt->superBlockIndex) ;
    free(bwtInc->bwt);
    //free(bwtInc->cumulativeCountInCurrentBuild);
    //free(bwtInc->packedText);
    //free(bwtInc->textBuffer);
    //free(bwtInc->packedShift);
	free(bwtInc->workingMemory);
	free(bwtInc);
}

static uint64_t BWTFileSizeInWord64(const uint64_t numChar)
{
	// The $ in BWT at the position of inverseSa0 is not encoded
	return (numChar + CHAR_PER_WORD - 1) / CHAR_PER_WORD;
}

void BWTSaveBwtCodeAndOcc64(const BWT64 *bwt, const char *bwtFileName, const char *occValueFileName)
{
	gzFile bwtFile;
/*	FILE *occValueFile; */
	uint64_t bwtLength;

	bwtFile = xzopen(bwtFileName, "w");
	if (bwtFile == NULL) {
		fprintf(stderr, "BWTSaveBwtCodeAndOcc64(): Cannot open BWT code file!\n");
		exit(1);
	}

	gzwrite(bwtFile, &bwt->inverseSa0, sizeof(uint64_t)* 1);
	gzwrite(bwtFile, bwt->cumulativeFreq + 1, sizeof(uint64_t) * ALPHABET_SIZE);
    gzwrite(bwtFile, &bwt->isDivide, sizeof(uint8_t) * 1);
    gzwrite(bwtFile, &bwt->divideNumber, sizeof(uint32_t) * 1);
    if(bwt->isDivide == 1)
    {
        gzwrite(bwtFile, bwt->divideBoundary, sizeof(DivideBoundary) * bwt->divideNumber);
    }
    gzwrite(bwtFile, bwt->sentinelPosition, sizeof(uint64_t) * bwt->divideNumber);
    gzwrite(bwtFile, bwt->sentinelSA, sizeof(uint64_t) * bwt->divideNumber);
	bwtLength = BWTFileSizeInWord64(bwt->textLength);
	gzwrite(bwtFile, bwt->bwtCode , sizeof(unsigned int) * bwtLength);
	gzclose(bwtFile);
/*
	occValueFile = (FILE*)fopen(occValueFileName, "wb");
	if (occValueFile == NULL) {
		fprintf(stderr, "BWTSaveBwtCodeAndOcc(): Cannot open occ value file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, occValueFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, occValueFile);
	fwrite(bwt->occValue, sizeof(unsigned int), bwt->occSizeInWord, occValueFile);
	fwrite(bwt->occValueMajor, sizeof(unsigned int), bwt->occMajorSizeInWord, occValueFile);
	fclose(occValueFile);
*/
}

/*
int bwt_bwtgen_main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: bwtgen <in.pac> <out.bwt>\n");
		return 1;
	}
	bwt_bwtgen64(argv[1], argv[2]);
	return 0;
}

#ifdef MAIN_BWT_GEN

int main(int argc, char *argv[])
{
	return bwt_bwtgen_main(argc, argv);
}

#endif
*/
