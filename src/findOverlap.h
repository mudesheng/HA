#include "bwtindex.h"
#include "hwgsa_para.h"

#define AVERAGE_READ_LEN 150


typedef struct {
    int64_t *db;
    int count, max;
} ISABase;

typedef struct {
    uint8_t *bnt;
    int count, max;
    int64_t isa, sa;
	ReadLocation rl;
} ReadInfo;

typedef struct {
    ReadInfo *ri;
    int count, max;
} ReadBase; 

typedef struct {
    int consis_len;
    int max_diff;
    uint8_t flag:1; // denote has been  mapped
} OverlapInfo;

int hwgsa_findOverlap(int argc, char *argv[]); 

