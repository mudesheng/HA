#ifndef GENERATE_GRAPH_H
#define GENERATE_GRAPH_H

#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include "bwtindex.h"
#include "utils.h"
#include "seqIO.h"
#include "bwt_gen.h"

// global variable
extern int32_t KMER ;
int32_t DEPTH ;
        
#define PEAK_SD 0.3
#define MIN_TIMES 5
#define LINK_BIGGER_THAN_ONE(link) (link != 0 && link != 1 && link != 2 && link != 4 && link != 8)
#define LINK_EQUL_ONE(link) (link == 1 || link == 2 || link == 4 || link == 8)
#define BRANCH_TWO_SET_FLAG  0X80000000
#define VERTEX_ID_GET_FLAG 0X7FFFFFFF
#define IN 1 // the degree  flag
#define OUT 2 

//#define FORWARD 0 
//#define BACKWARD 1 


typedef struct kmer_info {
    uint32_t freq:8 ;
    uint32_t lLink:4 ;
    uint32_t rLink:4 ;
    uint32_t base:2 ;
    uint32_t SNP:2 ; // SNP 
    uint32_t variant:1 ; // flag whether has a SNP 
    uint32_t direction:1 ; // the strand of double strand
    uint32_t deleted:1 ; // tag whether deleteed
    uint32_t curprocess:1 ; // in the processing state 
    uint32_t processed:1 ; // denote if has been processed 
	uint32_t state:2 ; // note the state of kmer_t
} kmer_info ;

typedef union Ku{
    kmer_info ki;
    uint32_t hi; // lookup kmer_info struct as a single whole
} Ku;

typedef struct kmer_t {
    uint64_t kmer[(MAX_KMER + 32 - 1) / 32];
    Ku ku;
} kmer_t ;

typedef struct KmerToVer {
	kmer_t k ;
	int ID ; // the vertex ID 
} KmerToVer ;

typedef struct VerticesIndex {
	KmerToVer *skArr ; // start kmer&ID array
	KmerToVer *ekArr ; // end kmer&ID array
	long size ;
	long count ;
} VerticesIndex ;

typedef struct SNPinfo {
    int position ;
    uint8_t SNP:2 ; 
} SNPinfo ;

typedef struct SNPset {
    SNPinfo *info ;
    int size ;
    int count ;
} SNPset ;

typedef struct HashTable {
    kmer_t *table;
    uint64_t count ;
    uint64_t size ;
    double loadFactor ;
} HashTable ;

typedef struct KmerAddr{
    int8_t hID ;  // hash ID, -1 denote have crash before, not need rewrite again, -2 denote rewrite as virgin 
    int64_t addr ;
    kmer_t k ;
} KmerAddr ;

typedef struct KmerStack {
    KmerAddr *kmerAddr ;
    int size ;
    int count ;
} KmerStack ;

typedef struct KmerArr {
	kmer_t *arr ;
	int size, count;
} KmerArr ;

typedef struct kmer_graph {
    //uint64_t  *edge_contig ; // have been change to pointer 
	char *seq ; // the sequence string
    uint8_t direction ;
    int  stackOffset ; // the offset of kmerStack 
    int count ;
    int size ; 
    int depth ; 
    int origin_graph_index ;
    uint8_t freq ;
    uint8_t processed:1 ; // marking if have been processed
    int state ; // the state of terminal position, same as retV of ThreadArgs
	int edgeID ;
	int vertexID ; // if kmer has two side branchs,  -1 note left branch, 1 note rigth branch
    int highRepeatKmerNum ; // the kmer_t freq >= 255 will be count 
    int rGraph[4]; // the start position of edge_contig right link
    int lGraph[4]; // the start position of edge_contig left link
    //int erGraph_link[4]; // the end position of edge_contig right link
    //int elGraph_link[4]; // the end position of edge_contig left link
} kmer_graph ;

typedef struct ThreadArgs {
    KmerStack kmerStack ;
    kmer_t initialKmer ; // initial kmer to process 
    uint8_t flag ; // thread operation flag (0 is wait initial, 1 is wait for processing, 2 is on the processing state , 3 is finished ,4 is put kmer to hash table,  5 is exit from thread, )
    int kmer_graph_index ; 
    int origin_graph_index;
    uint8_t direction ;
    int depth ;
    uint8_t retV ; // 0 note exit normally,1 note encounter a node that have processed just now in this kmer graph, 2 note run up to a have been processed node, 3 note encounter a branch node, 4 note encounter a tip terminal
} ThreadArgs ;

typedef struct BranchKmer {
    kmer_t initialKmer ;
    uint8_t direction ; // FORWARD and BACKWARD 
    int depth ;
    int kmer_graph_index ;
    int origin_graph_index;
} BranchKmer ;

typedef struct BranchQueue {
    BranchKmer *branchKmers ;
    int size ;
    int start ;
    int end ;
} BranchQueue ;

typedef struct BranchStack {
    int *indexStack ;
    int size;
    int cur_pos ;
} BranchStack;


typedef struct Edge {
	int ID ; // same  as vertex ID
	int rvertex[4] ;
	int lvertex[4] ;
} Edge ;

typedef struct EdgeArr {
	long size, count ;
	Edge *arr ;
} EdgeArr ;  

typedef struct Vertex {
	int ID ;
	char *seq ;
	uint8_t *freq ; // the kmer freq array, freq array length equal to s_size - KMER + 1 
	int s_size ; // the seq length 
	uint8_t indegree:2 ; // the in-degree
	uint8_t outdegree:2 ; // the out-degree
	int rvertex[4];
	int lvertex[4];
} Vertex ;

typedef struct VertexArr {
	long size ;
	int count ;
	Vertex *arr ;
} VertexArr ;



/*
 * hashlittle2: return 2 32-bit hash values
 *
 * This is identical to hashlittle(), except it returns two 32-bit hash
 * values instead of just one.  This is good enough for hash table
 * lookup with 2^^64 buckets, or if you want a second hash if you're not
 * happy with the first, or if you want a probably-unique 64-bit ID for
 * the key.  *pc is better mixed than *pb, so use *pc first.  If you want
 * a 64-bit value do something like "*pc + (((uint64_t)*pb)<<32)".
 */
void hashlittle2( 
  const void *key,       /* the key to hash */
  size_t      length,    /* length of the key */
  uint32_t   *pc,        /* IN: primary initval, OUT: primary hash */
  uint32_t   *pb)   ;     /* IN: secondary initval, OUT: secondary hash */


/*
#define addToQueue(graphIndex, ibase, idirection, iDepth, ikmer) {                                                      \
    if((idirection) == FORWARD)                                                                                         \
    {                                                                                                                   \
        initialGraph[graphIndex].rGraph[(ibase)] = gCount;                                                              \
        getNextKmer(branchQueue.branchKmers[branchQueue.end].initialKmer.kmer, ikmer, ibase, KMER);                    \
    } else {                                                                                                            \
        initialGraph[graphIndex].lGraph[(ibase)] = gCount;                                                              \
        getPreviousKmer(branchQueue.branchKmers[branchQueue.end].initialKmer.kmer, ikmer, ibase, KMER);                 \
    }                                                                                                                   \
    branchQueue.branchKmers[branchQueue.end].direction = idirection ;                                                   \
    branchQueue.branchKmers[branchQueue.end].kmer_graph_index = gCount ;                                                \
    branchQueue.branchKmers[branchQueue.end].origin_graph_index = graphIndex;                                           \
    gCount++ ;                                                                                                          \
    branchQueue.branchKmers[branchQueue.end].depth = iDepth;                                                            \
    branchQueue.end++ ;                                                                                                 \
    // realloc memory if overflow                                                                                     \
    if(gCount >= graphNum)                                                                                              \
    {                                                                                                                   \
        graphNum <<= 1 ;                                                                                                \
        initialGraph = (kmer_graph*)xrecalloc(initialGraph, graphNum * sizeof(kmer_graph));                             \
    }                                                                                                                   \
    if(branchQueue.end >= branchQueue.size)                                                                             \
    {                                                                                                                   \
        branchQueue.size <<= 1 ;                                                                                        \
        branchQueue.branchKmers = (BranchKmer*)xrecalloc(branchQueue.branchKmers, branchQueue.size * sizeof(BranchKmer)); \
    }                                                                                                                   \
}                                                                                         

#define addBranch(tIndex, curKmer, preKmer, idirection)                                                                                     \
{                                                                                                                               \
    uint8_t base ;                                                                                                              \
    kmer_t tk ;                                                                                                                 \
    int len = (KMER + 32 -1)/32;                                                                                                \
    kmer_info ki = kmerStack.kmerAddr[initialGraph[tIndex].stackOffset + initialGraph[tIndex].count -1].k.ku.ki;               \
    if(ki.direction == 1)                                                                                                       \
    {                                                                                                                           \
        for(uint8_t z = 0; z < 4; z++)                                                                                          \
        {                                                                                                                       \
            if((ki.rLink & (1<<z)) > 0)                                                                                         \
            {                                                                                                                   \
                base = (~z) & 0x3;                                                                                              \
                getPreviousKmer(tk.kmer, curKmer.kmer, base, KMER);                                                             \
                if(memcmp(tk.kmer, preKmer.kmer, sizeof(uint64_t) * len) != 0)                                                  \
                {                                                                                                               \
                    addToQueue(tIndex, base, BACKWARD, initialGraph[tIndex].depth, curKmer.kmer);                               \
                }                                                                                                               \
            }                                                                                                                   \
            if((ki.lLink & (1<<z)) > 0)                                                                                         \
            {                                                                                                                   \
                base = (~z) & 0x3;                                                                                              \
                getNextKmer(tk.kmer, curKmer.kmer, base, KMER);                                                                 \
                if(memcmp(tk.kmer, preKmer.kmer, sizeof(uint64_t) * len) != 0)                                                  \
                {                                                                                                               \
                    addToQueue(tIndex, base, FORWARD, initialGraph[tIndex].depth, curKmer.kmer);                                \
                }                                                                                                               \
            }                                                                                                                   \
        }                                                                                                                       \
    } else {                                                                                                                    \
        for(uint8_t z = 0; z < 4; z++)                                                                                          \
        {                                                                                                                       \
            if((ki.rLink & (1<<z)) > 0)                                                                                         \
            {                                                                                                                   \
                base = z ;                                                                                                      \
                getNextKmer(tk.kmer, curKmer.kmer, base, KMER);                                                                 \
                if(memcmp(tk.kmer, preKmer.kmer, sizeof(uint64_t) * len) != 0)                                                  \
                {                                                                                                               \
                    addToQueue(tIndex, base, FORWARD, initialGraph[tIndex].depth, curKmer.kmer);                                \
                }                                                                                                               \
            }                                                                                                                   \
            if((ki.lLink & (1<<z)) >0)                                                                                          \
            {                                                                                                                   \
                base = z;                                                                                                       \
                getPreviousKmer(tk.kmer, curKmer.kmer, base, KMER);                                                             \
                if(memcmp(tk.kmer, preKmer.kmer, sizeof(uint64_t) * len) != 0)                                                  \
                {                                                                                                               \
                    addToQueue(tIndex, base, BACKWARD, initialGraph[tIndex].depth, curKmer.kmer);                               \
                }                                                                                                               \
            }                                                                                                                   \
        }                                                                                                                       \
    }                                                                                                                           \
}                                                                                                   
*/


void generateKmerGraph(Arguments *arguments) ;
void getPreviousKmer(uint64_t *k, uint64_t *input, uint8_t base , uint32_t len);
void getNextKmer(uint64_t *k, uint64_t *input, uint8_t base , uint32_t len) ;
VertexArr *DBGverticesRestore(char *name) ;
EdgeArr *DBGEdgeRestore(char *name);
void EdgeArrFree(EdgeArr *edgeArr);


int hwgsa_genGraph(int argc, char *argv[]);


#endif




