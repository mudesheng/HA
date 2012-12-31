#ifndef HWGSA_PARA_H
#define HWGSA_PARA_H

#include <stdint.h>

#define MIN_KMER 33 
#define MAX_KMER 255
#define DEF_KMER 63
#define MIN_KMER_EXTENSION_LENGTH 10
#define LOW_QUAL_CUTOFF 10 // low quality cutoff
#define READ_SIZE uint16_t  // define the variable length of reads
#define FORWARD 3
#define BACKWARD 1
#define BIDIRECTION 2 
#define LEN_INTERVAL_MAJOR 512
#define PROGRAM_NAME "HWGSA"
#define GET_CHAR_FROM_PAC(pac, pos) (((pac)[(pos)/CHAR_PER_BYTE] >> (((~(pos)) & 0x3)<<1)) & 0x3)
#ifndef PATH_LEN 
#define PATH_LEN 1024 
#endif

typedef struct Arguments {
    char conffile[PATH_LEN] ; // configure file 
    char prefix[PATH_LEN] ; // the prefix of the output file (contain outpufDir)
    int NCPU ;
    int maxMem ; // suffix is the Giga(G)
    uint32_t K ;  // kmer len
    char outputDir[PATH_LEN] ;
    FILE *logfp ; // the file pointer of log
    int qual ; //minimum  quality requirement
    int G ;  // genome size ,suffix is the Mega(M)
    uint8_t H ; // heterozygosity of diploid genome 
    uint64_t seq_len ;
    int min_kmerfreq ; //minimum kmer frequency for construct De bruijn graph, must >=2 
    int depth ; // the length of kmer Graph stretch
    int min_unitig_len ; // minimum length of unitig for output to *.unitig file
	int maxReadLen ; 
	int minReadLen ;
} Arguments ;

typedef struct {
    uint32_t insert_size ;
    uint32_t insert_SD ;
    //uint8_t reverse: 1 ;
} PE ; // Pair End library  

typedef struct {
    char *name ;
    uint64_t offset ;
    uint32_t number_rd ;
    READ_SIZE read_len ;
    READ_SIZE *length ;
    uint64_t *lenMajor; // the index length of Major between LEN_INTERVAL_MAJOR
    uint8_t diverse:1 ;
    uint8_t paired:1 ;
    uint8_t asm_flag:2; // denote which assembly phase used, note 1 used for all step of assembly pipeline, note 2 used for scaffold phase, 3 used for filling gap 
    uint8_t seq_profile:3; // denote the data origin 
    uint8_t qual_benchmark; // the benchmark of quality score, strength suggest use phred+33
    PE  pe ;
} LIB ; // library 
typedef struct {
    int64_t len_pac ; // length of package 
    int64_t num_seqs ; // number of reads 
    int32_t num_lib; // number of library by libPE, libSE
    int max_rd_len; // maximum read length
    LIB *lib ;
} lib_info ; // library information  of package

typedef struct {
    uint8_t paired ;
    char f1name[PATH_LEN] ;
    char f2name[PATH_LEN] ;
} ReadName ;

typedef struct {
    ReadName *readName ;
    int max ;
    int count ;
    uint8_t *bntBuf ;
    uint8_t *qualBuf ;
    uint64_t seq_len ;
    //uint64_t unprocessed_rds ;
} BntWriteArgs ;

// parseConffile return struct 
typedef struct {
    lib_info *libIndex_info ;
    BntWriteArgs *bntWriteArgs ;
} PConfReturn ;

typedef struct {
    int64_t bound;
    int length[2];
	int paired:1; // if reads in pair end library, set paired = 1 else paired = 0
	int serial:2; // serial note the number of readID
	int size: 24; // size note the length of certain read length
	int64_t readID ;
	PE pe;
} ReadLocation;

char *KStrip(char *s) ;
void preProcessArgs(Arguments *arguments) ;
PConfReturn *parseConffile(char *conffile) ;
void free_Arguments(Arguments * arguments) ;
ReadLocation readBoundaryLookup(const lib_info *libIndex_info,const int64_t position) ;
int checkSA(const ReadLocation rl, const int64_t a, const int fixed_len) ;
void writeArgsToFile(const Arguments *arguments, FILE *fp);
int checkGenGraphArgs(Arguments *arguments, FILE *fp);
void writeInfo2Ann(const char *fn, lib_info *libIndex_info);
lib_info *lib_info_restore(const char *filename);
void lib_infoFree(lib_info *libIndex_info);


#endif

