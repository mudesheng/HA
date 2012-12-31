#ifndef BWTINDEX_H
#define BWTINDEX_H
   
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <zlib.h>
#include <omp.h>
#include <assert.h>
#include <stdint.h>
#include "seqIO.h"
#include "hwgsa_para.h"
#include "bwt.h"
#include "utils.h"
#include "bwt_gen64.h"

#define NUM_LIB 4 // initial number of library
#define INI_SIZE 0x1000 // initial array size 4K
//#define DEBUG 1 // debug flag
#define BWTINT uint64_t 
#define MIN_KMERFREQ 4
#define MIN_ZERO_GAP_NUM 6 // must bigger than 4 
#define NUM_FILE 4 // the number of files in the one library 
//#define QUAL_BASE 0x7F // the transform baseline of base quality
#define DEFAULT_DEPTH 40  // default sequence coverage or depth
#define AVERAGE_CAPACITY 10000 // the average number kmer of every block  
#define DNA_SYMBOLS_NUMBER 4
#define DEBUG 1

#define SET_LCP_FLAG(lcp, k) ((lcp)[(k)>>3] |= SET_BIT[(k) & 0X7])
#define GET_LCP_FLAG(lcp, k) ((lcp)[(k)>>3] & SET_BIT[(k) & 0X7])
#define GET_CHAR_FROM_UINT64(pac, k) (((pac) >> (((~(k)) & 0X1F)<<1)) & 0X3)
#define GET_CHAR_FROM_PAC(pac, k) (((pac)[(k)>>2] >> (((~(k)) & 0X3)<<1)) & 0X3)

extern const uint8_t SET_BIT[];
extern const uint8_t RESET_BIT[];
// global variable
int32_t KMER ;
int32_t NCPU ;


/*
typedef struct SaIndexRange {
    BWTINT startSaIndex ;
    BWTINT endSaIndex ;
} SaIndexRange ;
*/

/*
typedef struct BWT {
    BWTINT textLength ;             // length of the text 
    BWTINT saInterval ;             // interval between two SA values stored explicitly
    BWTINT inverseSaInterval ;      // interval between two inverse SA stored explicitly
    BWTINT inverseSa0 ;             // SA-1[0]
    BWTINT *cumulativeFreq ;        // cumulative frequency
    BWTINT *bwtCode ;               // BWT code
    BWTINT *occValue ;              // Occurrence values stored explicitly
    BWTINT *occValueMajor ;         // Occurrence values stored explicitly
    BWTINT *saValue ;               // SA values stored explicitly
    BWTINT *inverseSa ;             // Inverse SA stored explicitly
    SaIndexRange *saIndexRange ;    // SA index range
    int64_t saIndexRangeNumOfChar ; // Number of characters indexed in SA index range
    BWTINT *saValueOnBoundary ;     // Pre-calculated frequently referred data
    BWTINT *decodeTable ;           // For decoding BWT by table lookup
    BWTINT decodeTableGenerated ;   // == TRUE if decode table is generated on load and will be freed
    BWTINT bwtSizeInWord ;          // Temporary variable to hold the memory allocated
    BWTINT occSizeInWord ;          // Temporary variable to hold the memory allocated
    BWTINT occMajorSizeInWord ;     // Temporary variable to hold the memory allocated
    BWTINT saValueSize ;            // Temporary variable to hold the memory allocated
    BWTINT inverseSaSize ;          // Temporary variable to hold the memory allocated
    BWTINT saIndexRangeSize ;       // Temporary variable to hold the memory allocated
} BWT ;


typedef struct BWTInc {
    BWT *bwt ;
    unsigned int numberOfIterationDone ;
    BWTINT *cumulativeCountIncurrentBuild ;
    unsigned int availableWord ;
    BWTINT targetTextLength ;
    float targetNBit ;
    BWTINT buildSize ;
    BWTINT initialMaxBuildSize ;
    BWTINT incMaxBuildSize ;
    BWTINT firstCharInLastIteration ;
    BWTINT *workingMemory ;
    BWTINT *packedText ;
    BWTINT *textBuffer ;
    BWTINT *packedShift ;
} BWTInc ;
*/

typedef struct LCPBound {
    int64_t low ; // low bound 
    int64_t high ; // high bound
    int16_t depth; // the len of backward cycle
    //uint8_t base:2 ; // the base char
    //uint8_t flag:1 ; // the flag of base char
    //int offset ; // offset of low
    //int offset_h ; // offset of high
} LCPBound ;

typedef struct LongZeroSA {
    uint64_t ISA[2] ;
    uint64_t SA[2] ;
} LongZeroSA ;

typedef struct KmerFreq {
    uint64_t sa:56 ;
    uint64_t freq:8 ;
} KmerFreq ;

typedef struct KmerInfo {
    uint64_t kmer[(MAX_KMER + 1 + 32 -1) / 32]; // len = KMER + 1 
    //uint64_t rkmer[(MAX_KMER_LEN + 1 + 32 -1) / 32]; // reverse complementary kmer
    uint8_t freq ;
} KmerInfo ;

typedef struct KmerFreqCurve {
    uint8_t peak ;
    uint8_t valley ;
} KmerFreqCurve;

#ifdef __cplusplus
extern "C" {
#endif

int bwtIndex(lib_info *libIndex_info, BntWriteArgs *bntWriteArgs,Arguments *arguments) ;
int hwgsa_index(int argc , char *argv[]);
KmerFreqCurve parseKmerFreqFile(const char *fn, Arguments *arguments);
extern int64_t posSA2PosBWT(const bwt_t *bwt, const int64_t isa);

#ifdef __cplusplus
}
#endif

#endif



