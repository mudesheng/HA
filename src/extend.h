#ifndef EXTEND_UNITIG_H
#define ExTEND_UNITIG_H

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include <omp.h>


#include "bwtindex.h"

#define STEP_LEN 5  // the step of kmer bin skipped
#define MAX_READ_LEN 500 
#define SA_SET_ACCESS_FLAG  0X0100000000000000 
#define SA_ACCESS_MASK_FLAG 0X00FFFFFFFFFFFFFF
#define MAX_UNMAP_REGION_ALLOW 20
#define MAX_INDEL_NUM 1
#define MAX_INDEL_LEN 2
#define MAX_MISMACTH 2 
#define MIN_FLANK 10 
#define MAX_MAPFILE_LINE_LENGTH 2028

typedef struct {
    int64_t low;
    int64_t high;
} Region ;

typedef struct {
    int chance_num ;
    int len ;
    uint8_t *s[3];
} Seq_region;

typedef struct {
    int64_t isa ; // the read base inverse SA
    int distance ; // the distance to the processing position
} BaseISA ;

typedef struct {
    BaseISA *baseIndex;
    int count, max ;
} ISABase ;

typedef struct {
    int64_t pair_head_sa;
    int64_t read_head_sa ;
    int64_t isa ;
    int read_len ;
	int ref_position;
    uint32_t indel_num:4 ;
    uint32_t mismatch_num:4 ; 
    uint32_t read_type:4 ; // 1 note Illumina/solexa, 2 note 454 
    uint32_t score:4 ; // the macth score, no mismatch is note 0, mismatch mapped note 1~ 14, 0XF note unmapped
    uint32_t deleted:4 ;
    uint32_t order:1 ; // 0 note the first read of pair end reads, 1 note the second read of pair end reads
    uint32_t strand:1 ; // 0 denote '+', 1 denote '-'
    uint32_t cigar_count:4 ;
    char cigar[MAX_UNMAP_REGION_ALLOW];
    char *ref_name;
} read_info ;

typedef struct {
    read_info *ri ;
    int count, max ;
} MapBase ;

typedef struct {
    int64_t sa ; // the SA of calculateSA
    int64_t isa ; // the last isa
    int len, max; // the base length
    uint8_t *s ; // the base sequence
    int mbIndex; // the index of MapBase
} SA_REV ;

typedef struct {
    int8_t max_indel_num; // max proportion of gap number allow, e.g define 2 , and read length 100, max allow two gaps
    int8_t max_indel_len; // max gap extention length allow
    int8_t max_mismatch; // max proportion of mismatch number allow
    int8_t min_flank; // min flank length of read two end not allow indels
    int8_t min_step; // the minumun distance between indels
} AlignLimit ;

typedef struct {
	char mapLine[MAX_MAPFILE_LINE_LENGTH];
} MapLineInfo;

typedef struct {
	MapLineInfo *mlinfo;
	int count, max;
} MapInfo;

typedef struct {
    uint64_t pos:47 ;
    uint64_t flag:1 ;
    uint64_t read_len:16 ;
    int unitigID;
    int ref_pos;
} ReadPos;

typedef struct {
    int unitigID[2];
    int pair_num;
} Unitig_pair;

typedef struct {
    Unitig_pair *unitig_pairs;
    int64_t count, max ;
} UnitigPInfo;

typedef struct {
    UnitigPInfo *unitigPInfo;
    int count, max;
} UnitigGroup;

typedef struct {
    int *insert_size;
    int count, max ;
} Insert_info;

typedef struct {
    int unitigID ;
    int64_t offset ;
} UnitigLoc;

typedef struct {
    UnitigLoc *ul;
    int max, count;
} UnitigLocTable;

typedef struct {
    int ugID; // unitigGroup ID
    UnitigPInfo upi; 
    faSeq_t fs;
    MapBase mb; 
} UGInfo; // unitigs group information

typedef struct {
    int unitigID[2];
	int len[2]; // the length of unitig pair
    int gap_len;
    int SD;
    double score;
} UnitigFlatten;

typedef struct {
    UnitigFlatten *uf;
    int count, max;
} UnitigScaff;

typedef struct {
    int gap_len;
    int SD;
} EstimateGap;

typedef struct{
    EstimateGap *eg;
    int count, max;
    int unitigID[2];
} UnitigPairGapInfo;

typedef struct {
    int unitigID ;
	int len; // unitig length
    int order;
    int gap_len;
    int SD ;
} OrderElem;

typedef struct{
    OrderElem *oe;
    int count, max;
} OrderUnitigInfo;

typedef struct {
    uint64_t pair_head_sa:56 ;
    uint64_t order:1 ;
    uint64_t flag:2 ; //  0 is initial state, 1 denote pair end read, 2 denote this is need map read , 3 denote this is deleted read 
    kstring_t seq, qual;
} ReadSeq;

typedef struct {
    ReadSeq *rs;
    int count, max;
} ReadPool;

int hwgsa_extend(int argc, char *argv[]);

#endif
