#include <stdio.h>

#define LEN_WIDTH 20 

typedef struct SaPos {
	uint32_t vID ; // the ID of vertex
	int loc ; // loc of SA
	int64_t sa ;  // note the begin of read in the .pac
	int retainLen ; // the retain length of upmapped
	ReadLocation rl ;
	char *cigar ; // map info string
	int size ; // the size of cigar
	char strand // strand of read map
} SaPos ;

typedef struct SPArr {
	int size, count ;
	SaPos *sp ;
} SPArr ;


typedef struct ISACache {
	int64_t isa;
	int loc ; // loc of kmer coordinate
} ISACache ;

typedef struct ReadMapInfo {
	char *cigar ; // store the mapping information , format : "readID\tstrand_orientation(+/-)\tReadLength\tMapPath(StartPosition:[vertexID:]+endPosition)\n"
	int size ;
} ReadMapInfo ;

typedef struct ISARegion {
	int64_t low ;
	int64_t high ;
} ISARegion ;

typedef struct ReadMapArg {
	uint32_t vID ;
	int startLoc ;
	IsaRegion isaRegion ;
	SPArr *spArr;
	IsaCache *isaCache ;
} ReadMapArg ;

typedef struct ReadMapStack {
	ReadMapArg *rma ;
	int size, count;
} ReadMapStack ;

typedef struct MapInfo {
	uint32_t strand:1 ; // 0 note + strand, 1 note - strand
	uint32_t len:LEN_WIDTH ; // the length of read mapped
  	char *cigar ;	
} MapInfo ;

typedef struct MapIndex {
	MapInfo *arr ;
	int64_t size, count;
} MapIndex ;

typedef struct ReadLoc {
	uint64_t ID:34 ; // the readID
	uint64_t loc:30 ; // the loc position of edges
} ReadLoc ;

typedef struct MapEdgeInfo {	
	ReadLoc *readLocArr ;
	int size, count ;	
} MapEdgeInfo ;

