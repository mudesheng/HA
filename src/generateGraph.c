#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "bwt.h"
#include "generateGraph.h"

// global variable
extern int32_t KMER ;
extern int32_t NCPU ;
extern int32_t DEPTH ;

uint8_t MaskLink[] = {0xFE, 0xFD, 0xFB, 0xF7};

static inline uint64_t get_time_hash(const uint64_t *key, const int kmerLenByWord)
{
    uint64_t hash = 0 ;
    for(int i = 0; i < kmerLenByWord; i++)
    {
        hash += key[i];
        hash += (hash<<10);
        hash ^= (hash>>6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash ;
}


inline void getNextKmer(uint64_t *k ,uint64_t *input, uint8_t base , uint32_t len)
{
    int lw = (len + 32 - 1) / 32 ; 
    int r = len % 32 ;
    int i ;
    for( i = 0 ; i < lw - 2 ; i++)
    {
        k[i] = (input[i] << BIT_PER_CHAR);
        k[i] |= (input[i+1] >> (31<<1));
    }
    k[i] = (input[i] << BIT_PER_CHAR);
    k[i] |= ((input[i+1] >> ((r - 1)<<1)) & 0x3);
    k[i+1] = ((input[i+1] << ((32 - r + 1)<<1)) >> ((32 - r )<<1));
    k[i+1] |= base ;
}

inline void getPreviousKmer(uint64_t *k, uint64_t *input, uint8_t base , uint32_t len)
{
    int lw = (len + 32 - 1)/ 32 ;
    int r = len % 32, i ;
    uint64_t b ;
    k[lw - 1] = (input[lw-1] >> BIT_PER_CHAR);
    b = input[lw -2] & 0x3 ;
    k[lw - 1] |= (b << ((r -1)<<1));
    for(i = lw - 2 ; i >0 ; i--)
    {
        k[i] = (input[i]>> BIT_PER_CHAR);
        b = input[i-1] & 0x3 ;
        k[i] |= (b << (31<<1));
    }
    k[i] = (input[i] >> BIT_PER_CHAR);
    b = base ;
    k[i] |= (b << (31<<1));
}


static inline int is_prime_kh(uint64_t num) 
{
    uint64_t i, max ;
    if(num < 4) return 1 ;
    if(num % 2 == 0) return 0 ;
    max = (uint64_t) sqrt((float)num);
    for(i = 3; i < max ; i+=2) { if(num % i == 0) return 0 ; }
    return 1 ;
}

static inline uint64_t find_next_prime_kh(uint64_t num) 
{
    if(num % 2 == 0) num++;
    while(1) { if(is_prime_kh(num)) return num ; num += 2 ; }
}

// flag "1" denote  get first kmer, "2" denote get second kmer
static inline void processToKmer(KmerInfo *kmerBuf, kmer_t *kmer, const int  len, int flag)
{
    int kmerLenByWord = (len + 32 -1)/32;
    if(flag == 1)
    {
        for(int j = 0 ; j < kmerLenByWord ; j++) { kmer->kmer[j] = kmerBuf->kmer[j]; }
        kmer->ku.ki.rLink = (1<< (kmerBuf->kmer[kmerLenByWord - 1] & 0x3)) ;
        kmer->kmer[kmerLenByWord - 1] >>= BIT_PER_CHAR ;
    } else if(flag == 2) {
        kmer->ku.ki.lLink = (1<< ((kmerBuf->kmer[0]>>(31<<1)) & 0x3)) ;
        for(int j = 0 ; j < kmerLenByWord - 1 ; j++) 
        {
            kmer->kmer[j] = kmerBuf->kmer[j] << BIT_PER_CHAR ;
            if(j + 1 < kmerLenByWord - 1) kmer->kmer[j] |= ((kmerBuf->kmer[j+1] >> (31<<1)) & 0x3);
            else kmer->kmer[j] |= ((kmerBuf->kmer[j+1] >> ((len%32)<<1)) & 0x3);
            //kmer->kmer[j] |= ((j < kmerLenByWord - 2 ? (kmerBuf->kmer[j+1]>>(31<<1)) : (kmerBuf->kmer[j+1]>> ((KMER%32 - 1)<<1))) & 0x3) ;
        }
        kmer->kmer[kmerLenByWord -1] = (kmerBuf->kmer[kmerLenByWord - 1] << ((32 - len%32)<<1)) >> ((32 - len%32)<<1); 
    } else {
        fprintf(stderr, "[processToKmer] flag set error....exit\n");
        exit(1);
    }
}

static inline void pushToBuf(HashTable *buf, kmer_t *kmer, kmer_t *rkmer, const int kmerLenByWord, const uint32_t freq , const int hashNumber, int64_t *delay )
{
    int flag  = 0, n ;
    for(int j = 0 ; j < kmerLenByWord ; j++)
    {
        if(kmer->kmer[j] > rkmer->kmer[j]) { flag = 1; break ; }
        else if(kmer->kmer[j] < rkmer->kmer[j] ) break;
    }
    //n = (kmer->kmer[0] ^ rkmer->kmer[0]) % hashNumber ;
    if(flag == 1)
    {
        uint8_t link_rev[] = { 0, 8, 4, 0, 2, 0, 0, 0, 1};
        if(kmer->ku.ki.rLink > 0) rkmer->ku.ki.lLink = link_rev[kmer->ku.ki.rLink];
        if(kmer->ku.ki.lLink > 0) rkmer->ku.ki.rLink = link_rev[kmer->ku.ki.lLink];
        n = get_time_hash(rkmer->kmer, kmerLenByWord) % hashNumber;
        // if buf is full , wait for process
        while(buf[n].loadFactor != 0)  { usleep(10); *delay += 10; }
        memcpy(&buf[n].table[buf[n].count], rkmer, sizeof(kmer_t));
    } else {
        n = get_time_hash(kmer->kmer, kmerLenByWord) % hashNumber;
        // if buf is full , wait for process
        while(buf[n].loadFactor != 0)  { usleep(10);  *delay += 10; }
        memcpy(&buf[n].table[buf[n].count], kmer, sizeof(kmer_t));
    }
    buf[n].table[buf[n].count].ku.ki.freq = freq ;
    buf[n].count++ ;
    // check if buf overflow 
    if(buf[n].count >= buf[n].size)
    {
        // set overflow flag 
        ////#pragma omp atomic
        buf[n].loadFactor = 1 ;
    }
}

void assignKmerToHash(HashTable *hashTable, HashTable *buffer, const int ID, int64_t *delay)
{
    HashTable hashBuf ;
    int kmerLenByWord = (KMER + 32 - 1) / 32 ;
    // initial hashBuf 
    hashBuf.size = buffer->size ;
    hashBuf.count = 0 ;
    hashBuf.table = (kmer_t*)xcalloc(hashBuf.size, sizeof(kmer_t));
    int rep = 0 ;
    while(1)
    {
        if(buffer->loadFactor == 1)
        {
            uint64_t hashAdd ;
            uint32_t hashF, hashS ; // the fisrt part and second part of hash address
            memcpy(hashBuf.table, buffer->table, buffer->count * sizeof(kmer_t));
            hashBuf.count = buffer->count ;
            buffer->count = 0 ;
            ////#pragma omp atomic
            buffer->loadFactor = 0 ;

            for(int i = 0 ; i < hashBuf.count ; i++)
            {
                hashF = 0, hashS = 0 ;
                hashlittle2((void*)hashBuf.table[i].kmer, kmerLenByWord * sizeof(uint64_t), &hashF, &hashS ) ;
                hashAdd =  hashF + (((uint64_t)hashS)<<32) ;
                hashAdd %= hashTable->size ;
                while(1)
                {
                    if(hashTable->table[hashAdd].ku.ki.freq == 0 )
                    {
                        memcpy(&hashTable->table[hashAdd], &hashBuf.table[i], sizeof(kmer_t)); 	
                        hashTable->count++; break;
                    } else if(memcmp(hashTable->table[hashAdd].kmer, hashBuf.table[i].kmer, kmerLenByWord * sizeof(uint64_t)) == 0)
                    {
                        rep++;
                        hashTable->table[hashAdd].ku.ki.lLink |= hashBuf.table[i].ku.ki.lLink ;
                        hashTable->table[hashAdd].ku.ki.rLink |= hashBuf.table[i].ku.ki.rLink ;
                        if(hashTable->table[hashAdd].ku.ki.freq + hashBuf.table[i].ku.ki.freq > 255)
                            hashTable->table[hashAdd].ku.ki.freq = 255 ;
                        else hashTable->table[hashAdd].ku.ki.freq += hashBuf.table[i].ku.ki.freq ;
                        break;
                    } else {
                        hashAdd++;
                        if(hashAdd >= hashTable->size) hashAdd = 0;
                    }
                }
            }
            if((double)hashTable->count / hashTable->size > 0.75)
            {
                fprintf(stderr, "repeat count: %d\n", rep);
                fprintf(stderr, "the thread : %d hash load factor bigger than 0.75, program exit...\n", ID);
                exit(1);
            } else if((double)hashTable->count / hashTable->size > hashTable->loadFactor) {
                fprintf(stderr, "repeat count: %d\n", rep);
                fprintf(stderr, "the thread : %d hash load factor bigger than %f\n", ID, hashTable->loadFactor);
            }
        } else if(buffer->loadFactor == 2) {
            break ; // exit to cycle and finish  function 
        } else { // buffer->loadFactor == 0 
            usleep(10);
            __sync_add_and_fetch(delay,10 );
        }
    }
    free(hashBuf.table);
}
// if kmer has been processed , return -1 . else return hash address 
uint64_t getHashAddr(const HashTable *t , const uint64_t *kmer, const uint32_t len  )
{
    uint64_t hashAdd ;
    uint32_t hashF, hashS ;

    hashF = 0, hashS = 0;
    hashlittle2((void*)kmer, len * sizeof(uint64_t), &hashF, &hashS);
    hashAdd = hashF + (((uint64_t)hashS)<<32);
    hashAdd %= t->size ;
    while(1)
    {
        if(memcmp(t->table[hashAdd].kmer, kmer, len * sizeof(uint64_t)) == 0) break ;
        else if(t->table[hashAdd].ku.ki.freq == 0 && DEBUG == 1) {
            fprintf(stderr, "getHashAddr() not found kmer in the hashTable , Please check\n");
            exit(1); 
        } else hashAdd = ((hashAdd >= (t->size -1)) ? 0 : (hashAdd+1)) ;
    }
    return hashAdd ;
}
/*// if is a branch or tip return 1 , else if a transform backward direction kmer return -1, else return 0 
inline int findBranchOrTip(int direction , kmer_t *k )
{
    int lb, rb ;
    for(int i = 0 , rb = 0 ; i < 4 ; i++)
    {
        if((k->rLink >> i) & 0x1 == 1) rb++ ;
    }
    for(int i = 0 , lb = 0 ; i < 4 ; i++)
    {
        if((k->lLink >> i) & 0x1 == 1) lb++ ;
    }
    // find different direction is also branch 
    if(direction == FORWARD) 
    {
        if(lb > 1) return -1 ;
        else if(rb != 1) return 1 ; 
    } else if(direction == BACKWARD) {
        if(rb > 1) return -1 ;
        else if(lb != 1) return 1 ;
    }      
    return 0 ;
} */

void basememcat(uint8_t *dest,const size_t dest_size,const uint8_t *src,const size_t src_size)
{
    int r = dest_size % CHAR_PER_BYTE ;
    if(src_size <= (CHAR_PER_BYTE - r)) 
    {
        dest[dest_size/CHAR_PER_BYTE] <<= (src_size<<1);
        dest[dest_size/CHAR_PER_BYTE] |= src[0];
    } else {
        int i, j;
        for(i = 0, j = dest_size/CHAR_PER_BYTE ; i < src_size / CHAR_PER_BYTE ; i++)
        {
            dest[j+ i] <<= ((CHAR_PER_BYTE - r)<<1);
            dest[j + i] |= (src[i] >> (r<<1));
        }
        if(src_size % CHAR_PER_BYTE > (CHAR_PER_BYTE - r) )
        {
            int sr = src_size % CHAR_PER_BYTE ;
            dest[j + i] <<= ((CHAR_PER_BYTE - r)<<1);
            dest[j + i] |= (src[i] >> ((sr - (CHAR_PER_BYTE - r))<<1));
            dest[j + i + 1] = ((src[i] << ((CHAR_PER_BYTE - (sr - (CHAR_PER_BYTE -r)))<<1)) >> ((CHAR_PER_BYTE - (sr - (CHAR_PER_BYTE - r)))<<1));
        } else {
            int sr = src_size % CHAR_PER_BYTE ;
            dest[j + i] <<= (sr<<1);
            dest[j + i] |= src[i] ;
        }
    }
}

kmer_graph *initialGraphStruct(const int initial_size, KmerStack *kmerStack, BranchQueue *branchQueue)
{
	kmer_graph *initialGraph = (kmer_graph*)xcalloc(initial_size, sizeof(kmer_graph));
	kmerStack->size = INITIAL_STACK_SIZE ;
	kmerStack->kmerAddr = (KmerAddr*)xcalloc(kmerStack->size, sizeof(KmerAddr));
	kmerStack->count = 0;

	memset(branchQueue, 0, sizeof(BranchQueue));
	branchQueue->size = initial_size ;
	branchQueue->branchKmers = (BranchKmer*)xcalloc(branchQueue->size, sizeof(BranchKmer));
	

	return initialGraph ;
}

/* return the number of tips removed */
int clipTipsInGraph(kmer_graph *initialGraph, KmerStack *kmerStack, const KmerFreqCurve curve, const int graphCount)
{
    //int size = (DEPTH + CHAR_PER_BYTE - 1)/CHAR_PER_BYTE ;
    //uint8_t *edge_contig0 = (uint8_t*)xcalloc(size , sizeof(uint8_t));
    //int count = 0 ;
    //int lStackIndex = -1, rStackIndex = -1 ;
    int tip_rm = 0 ;

    // remove tips and remove error bubbles
    {
        for(int graphIndex = 0; graphIndex < graphCount ; graphIndex++)
        {
            if(initialGraph[graphIndex].processed == 0)
            {
                int tip_n = 0 ;
                // test right graph and traverse right link and remove tips
                for(int i = 0; i < 4; i++)
                {
					int gI = initialGraph[graphIndex].rGraph[i];
                    if(gI > 0 && initialGraph[gI].processed == 0 ) 
                    {
                        if((initialGraph[gI].freq <= 2 && initialGraph[gI].count < 2 * KMER) || 
                            (initialGraph[gI].count <= KMER && initialGraph[gI].state == 4))
						{
							// set kmer deleted flag 
							for(int j = 0; j < initialGraph[gI].count ; j++)
							{
								if((j == initialGraph[gI].count - 1) && (initialGraph[gI].state != 4)) continue ;
								kmerStack->kmerAddr[initialGraph[gI].stackOffset + j].k.ku.ki.deleted = 1 ;
							}
							if(kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.direction ==1)
							{
								uint8_t base = (~i) & 0x3;
								kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.lLink &= MaskLink[base];
							} else {
								kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.rLink &= MaskLink[i];
							}
							initialGraph[graphIndex].rGraph[i] = -1 ;
							initialGraph[gI].processed = 1;
							if(initialGraph[gI].state != 4) kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].hID = -1;
							tip_rm++ ;
						}
                    }
                }
                
				// test left graph and traverse left link
                tip_n = 0;  
                for(int i = 0; i < 4; i++)
                {
					int gI = initialGraph[graphIndex].lGraph[i];
                    if(gI > 0 && initialGraph[gI].processed == 0)
                    {
                        if((initialGraph[gI].freq <= 2 && initialGraph[gI].count < 2 * KMER) ||
                            (initialGraph[gI].count <= KMER  && initialGraph[gI].state == 4))
						{
							// set kmer deleted flag
							for(int j = 0; j < initialGraph[gI].count ; j++)
							{
								if((j == initialGraph[gI].count - 1) && (initialGraph[gI].state != 4)) continue ;
								kmerStack->kmerAddr[initialGraph[gI].stackOffset + j].k.ku.ki.deleted = 1 ;
							}
							if(kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.direction ==1)
							{
								uint8_t base = (~i) & 0x3;
								kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.rLink &= MaskLink[base];
							} else {
								kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.lLink &= MaskLink[i];
							}
							initialGraph[graphIndex].lGraph[i] = -1 ;
							initialGraph[gI].processed = 1;
							if(initialGraph[gI].state != 4) kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].hID = -1;
							tip_rm++ ;
						}
                    }
                }
				// set local kmer graph deleted flag
				if((initialGraph[graphIndex].freq <= 2 && initialGraph[graphIndex].count < 2 * KMER) ||
						(initialGraph[graphIndex].count <= KMER && initialGraph[graphIndex].state == 4))
				{
					// set kmer deleted flag
					for(int j = 0; j < initialGraph[graphIndex].count ; j++)
					{
						if((j == initialGraph[graphIndex].count - 1) && (initialGraph[graphIndex].state != 4)) continue ;
						kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + j].k.ku.ki.deleted = 1 ;
					}
					initialGraph[graphIndex].processed = 1;
					if(initialGraph[graphIndex].state == 2 || initialGraph[graphIndex].state == 1 ) kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].hID = -1;
					tip_rm++ ;
				}
                if(initialGraph[graphIndex].state == 2 || initialGraph[graphIndex].state == 1)
                {
                    kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].hID = -1 ;
                }
                initialGraph[graphIndex].processed = 1;
            }
        }
    }
	return tip_rm ;
}

int ReduceGraph(kmer_graph *initialGraph, KmerStack *kmerStack, const KmerFreqCurve curve)
{
	int tip_rm = 0 ;

	BranchStack branchStack ;
	branchStack.size = 10; branchStack.cur_pos = 0;
	branchStack.indexStack = (int*)xcalloc(branchStack.size, sizeof(int)) ;
	branchStack.indexStack[branchStack.cur_pos++] = 0 ;
	// remove tips and remove error bubbles
	{
		while(branchStack.cur_pos > 0)
		{
			int graphIndex = branchStack.indexStack[--branchStack.cur_pos];
			if(initialGraph[graphIndex].processed == 0)
			{
				int rbranch_n = 0 , lbranch_n = 0 ;
				// test rigth graph and traverse rigth link
				for(int i = 0; i < 4; i++)
				{
					int gI = initialGraph[graphIndex].rGraph[i];
                    if(gI > 0 && initialGraph[gI].processed == 0 ) 
                    {
                        if((initialGraph[gI].freq <= 2 && initialGraph[gI].count < 2 * KMER) || 
                            (initialGraph[gI].count <= KMER && initialGraph[gI].state == 4))
						{
							// set kmer deleted flag 
							for(int j = 0; j < initialGraph[gI].count ; j++)
							{
								if((j == initialGraph[gI].count - 1) && (initialGraph[gI].state != 4)) continue ;
								kmerStack->kmerAddr[initialGraph[gI].stackOffset + j].k.ku.ki.deleted = 1 ;
							}
							if(kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.direction ==1)
							{
								uint8_t base = (~i) & 0x3;
								kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.lLink &= MaskLink[base];
							} else {
								kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.rLink &= MaskLink[i];
							}
							initialGraph[graphIndex].rGraph[i] = -1 ;
							initialGraph[gI].processed = 1;
							if(initialGraph[gI].state != 4) kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].hID = -1;
							tip_rm++ ;
						}
                    } else { rbranch_n++ ; }
				}

				// test left graph and traverse left link
				for(int i = 0; i < 4; i++)
				{
					int gI = initialGraph[graphIndex].lGraph[i] ;
					if(gI > 0 && initialGraph[gI].processed == 0)
					{
						if((initialGraph[gI].freq <= 2 && initialGraph[gI].count < 2 * KMER) ||
							(initialGraph[gI].count <= KMER && initialGraph[gI].state == 4))
						{
							// set kmer deleted flag
							for(int j = 0; j < initialGraph[gI].count; j++)
							{
								if((j == initialGraph[gI].count - 1) && (initialGraph[gI].state != 4)) continue ;
								kmerStack->kmerAddr[initialGraph[gI].stackOffset + j].k.ku.ki.deleted = 1 ;
							}
							if(kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count - 1].k.ku.ki.direction == 1)
							{
								uint8_t base = (~i) & 0x3 ;
								kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.rLink &= MaskLink[base] ;
							} else {
								kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.lLink &= MaskLink[i] ;
							}
							initialGraph[graphIndex].lGraph[i] = -1 ;
							initialGraph[gI].processed = 1;
							if(initialGraph[gI].state != 4) kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].hID = -1 ;
							tip_rm++;
						}
					} else { lbranch_n++ ;}
				}
				
				// remove right error bubble
				if(rbranch_n > 1)
				{
					int error_branch = -1 ;
					int correct_branch = -1 ;
					for(int i = 0; i < 4; i++)
					{
						int gI = initialGraph[graphIndex].rGraph[i] ;
						if(gI > 0)
						{
							if(initialGraph[gI].freq <= 2 && (initialGraph[gI].state == 1 || initialGraph[gI].state == 3)) { error_branch = gI ; }
						   	else if(initialGraph[gI].freq > curve.valley && (initialGraph[gI].state == 1 || initialGraph[gI].state ==3)) {correct_branch = gI; }
						}	
					}
					
					if(error_branch > 0 && correct_branch > 0 &&
					  memcmp(&kmerStack->kmerAddr[initialGraph[error_branch].stackOffset + initialGraph[error_branch].count -1].k,
					  &kmerStack->kmerAddr[initialGraph[correct_branch].stackOffset + initialGraph[correct_branch].count -1].k, 
					  sizeof(kmer_t)) == 0 && initialGraph[error_branch].freq <= initialGraph[correct_branch].freq/5)
					{
						// set kmer deleted flag
						for(int j = 0; j < initialGraph[error_branch].count ; j++)
						{
							if((j == initialGraph[error_branch].count -1) && (initialGraph[error_branch].state != 4)) continue ;
							kmerStack->kmerAddr[initialGraph[error_branch].stackOffset + j].k.ku.ki.deleted = 1 ;
						}
						initialGraph[error_branch].processed = 1;
						kmerStack->kmerAddr[initialGraph[error_branch].stackOffset + initialGraph[error_branch].count - 1].hID = -1;
						tip_rm++;
					}		
				}

				// remove left error bubble
				if(lbranch_n > 1)
				{
					int error_branch = -1 ;
					int correct_branch = -1 ;
					for(int i = 0; i < 4; i++)
					{
						int gI = initialGraph[graphIndex].lGraph[i] ;
						if(gI > 0)
						{
							if(initialGraph[gI].freq <= 2 && (initialGraph[gI].state == 1 || initialGraph[gI].state == 3)) { error_branch = gI ; }
						   	else if(initialGraph[gI].freq > curve.valley && (initialGraph[gI].state == 1 || initialGraph[gI].state ==3)) {correct_branch = gI; }
						}	
					}
					
					if(error_branch > 0 && correct_branch > 0 &&
					  memcmp(&kmerStack->kmerAddr[initialGraph[error_branch].stackOffset + initialGraph[error_branch].count -1].k,
					  &kmerStack->kmerAddr[initialGraph[correct_branch].stackOffset + initialGraph[correct_branch].count -1].k,
					  sizeof(kmer_t)) == 0 && initialGraph[error_branch].freq <= initialGraph[correct_branch].freq/5)
					{
						// set kmer deleted flag
						for(int j = 0; j < initialGraph[error_branch].count ; j++)
						{
							if((j == initialGraph[error_branch].count -1) && (initialGraph[error_branch].state != 4)) continue ;
							kmerStack->kmerAddr[initialGraph[error_branch].stackOffset + j].k.ku.ki.deleted = 1 ;
						}
						initialGraph[error_branch].processed = 1;
						kmerStack->kmerAddr[initialGraph[error_branch].stackOffset + initialGraph[error_branch].count - 1].hID = -1;
						tip_rm++;
					}		
				}
			}
		}
	}

	// clean work
	free(branchStack.indexStack);
	return tip_rm  ;
}

/*
void clipTipsInGraph(kmer_graph *initialGraph, KmerStack *kmerStack, const KmerFreqCurve curve)
{
    //int size = (DEPTH + CHAR_PER_BYTE - 1)/CHAR_PER_BYTE ;
    //uint8_t *edge_contig0 = (uint8_t*)xcalloc(size , sizeof(uint8_t));
    //int count = 0 ;
    //int lStackIndex = -1, rStackIndex = -1 ;
    int graphIndex ;

    BranchStack branchStack;
    branchStack.size = 10 ; branchStack.cur_pos = 0 ;
    branchStack.indexStack = (int*)xcalloc(branchStack.size, sizeof(int)); 
    branchStack.indexStack[branchStack.cur_pos++] = 0 ;
    // remove tips and remove error bubbles
    {
        while(branchStack.cur_pos > 0)
        {
            graphIndex = branchStack.indexStack[--branchStack.cur_pos];
            if(initialGraph[graphIndex].processed == 0)
            {
                int tip_n = 0 , rbranch_n = 0, lbranch_n = 0, bubble_n = 0, normalB = 0, normalIndex ;
                // test right graph and traverse right link
                for(int i = 0; i < 4; i++)
                {
                    if(initialGraph[graphIndex].rGraph[i] > 0 ) 
                    {
                        int gI = initialGraph[graphIndex].rGraph[i];
                        if(initialGraph[gI].freq > curve.valley && (initialGraph[gI].state != 4 || 
                        initialGraph[gI].count > 2 * KMER)) 
                        {
                            normalB++;
                            normalIndex = gI;
                        } else if(((initialGraph[gI].freq <= curve.valley && initialGraph[gI].count < 2 * KMER) || 
                            initialGraph[gI].count <= KMER) && initialGraph[gI].state == 4)  tip_n++;
                        else if(graphIndex > 0 && initialGraph[graphIndex].direction == BACKWARD && 
                            initialGraph[gI].state == 2 && initialGraph[gI].freq <= curve.valley &&
                            MIN_TIMES * initialGraph[gI].freq < initialGraph[graphIndex].freq &&
                            initialGraph[gI].count < 2 * KMER) bubble_n++;
                        rbranch_n++;
                    }
                }
                // remove tips
                if((tip_n >= 1 ) || (graphIndex > 0 && initialGraph[graphIndex].direction == BACKWARD && bubble_n == 1))
                {
                    for(int i = 0; i < 4; i++)
                    {
                        if(initialGraph[graphIndex].rGraph[i] > 0)
                        {
                            int gI = initialGraph[graphIndex].rGraph[i];
                            if( (initialGraph[gI].freq <= curve.valley && initialGraph[gI].count < 2 * KMER && 
                                initialGraph[gI].state == 4) || (initialGraph[gI].count <= KMER  && 
                                initialGraph[gI].state == 4) || (graphIndex > 0 && 
                                initialGraph[graphIndex].direction == BACKWARD && initialGraph[gI].state == 2 && 
                                initialGraph[gI].freq <= curve.valley && initialGraph[gI].count < 2 * KMER &&
                                MIN_TIMES * initialGraph[gI].freq < initialGraph[graphIndex].freq) )
                            {
                                for(int j = 0; j < initialGraph[gI].count ; j++)
                                {
                                    kmerStack->kmerAddr[initialGraph[gI].stackOffset + j].k.ku.ki.deleted = 1 ;  
                                }
                                if(kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.direction ==1)
                                {
                                    uint8_t base = (~i) & 0x3;
                                    kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.lLink &= MaskLink[base];
                                } else {
                                    kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.rLink &= MaskLink[i];
                                }
                                initialGraph[graphIndex].rGraph[i] = 0 ;
                                initialGraph[gI].processed = 1;
                                if(initialGraph[gI].state == 2) kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].hID = -1;
                            }
                        }
                    }
                    if(normalB >= 1)
                    {
                        if(branchStack.cur_pos + normalB >= branchStack.size)
                        {
                            branchStack.size <<= 1;
                            branchStack.indexStack = (int*)realloc(branchStack.indexStack, branchStack.size * sizeof(int));
                        }
                        for(int j = 0; j < 4; j++)
                        {
                            int tI = initialGraph[graphIndex].rGraph[j];
                            if(tI > 0 && initialGraph[tI].freq > curve.valley && initialGraph[tI].state == 3 )
                            {
                                branchStack.indexStack[branchStack.cur_pos++] = tI;
                            } else if(initialGraph[tI].freq > curve.valley && initialGraph[tI].state == 2) {
                                kmerStack->kmerAddr[initialGraph[tI].stackOffset + initialGraph[tI].count -1].hID = -1;
                                initialGraph[tI].processed = 1;
                            }
                        }
                    }
                } else if(rbranch_n == 2 && normalB == 1) { // remove bubble and chimeric connection
                    int remove_flag = 0, bI;
                    for(int i = 0 ; i < 4; i++)
                    {
                        if(initialGraph[graphIndex].rGraph[i] > 0)
                        {
                            int gI = initialGraph[graphIndex].rGraph[i] ;
                            if(initialGraph[gI].freq <= curve.valley && 
                                MIN_TIMES * initialGraph[gI].freq < initialGraph[normalIndex].freq  && 
                                MIN_TIMES * initialGraph[gI].freq < initialGraph[graphIndex].freq &&
                                initialGraph[gI].count < 2 * KMER) 
                            {
                                remove_flag = 1 ; 
                                bI = gI ;
                            }
                        }
                    }
                    if(branchStack.cur_pos >= branchStack.size)
                    {
                        branchStack.size <<= 1;
                        branchStack.indexStack = (int*)realloc(branchStack.indexStack, branchStack.size * sizeof(int));
                    }
                    if(initialGraph[normalIndex].state == 3)  branchStack.indexStack[branchStack.cur_pos++] = normalIndex;
                    else if(initialGraph[normalIndex].freq > curve.valley && initialGraph[normalIndex].state == 2) {
                        kmerStack->kmerAddr[initialGraph[normalIndex].stackOffset + initialGraph[normalIndex].count -1].hID = -1;
                        initialGraph[normalIndex].processed = 1;
                    }
                    if(remove_flag == 1 && initialGraph[bI].state == 3)// remove chimeric connection or bubble
                    {
                        int rI, lI, rn = 0, ln = 0;
                        for(int i = 0; i < 4; i++)
                        {
                            if(initialGraph[bI].rGraph[i] > 0) { rI = initialGraph[bI].rGraph[i]; rn++; }
                            if(initialGraph[bI].lGraph[i] > 0) { lI = initialGraph[bI].lGraph[i]; ln++; }
                        }
                        if( rn == 1 && ln == 1)
                        {
                            // remove edge
                            if(MIN_TIMES * initialGraph[bI].freq < initialGraph[lI].freq && 
                                MIN_TIMES * initialGraph[bI].freq < initialGraph[rI].freq)
                            {
                                int j ;
                                for(j = 0; j < initialGraph[bI].count - 1; j++)
                                {
                                    kmerStack->kmerAddr[initialGraph[bI].stackOffset + j].k.ku.ki.deleted = 1;
                                }
                                // remove link from beginning of branch edge 
                                for(j = 0; j < 4; j++)
                                {
                                    if(initialGraph[graphIndex].rGraph[j] == bI) break; 
                                }
                                if(kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.direction ==1)
                                {
                                    uint8_t base = (~j) & 0x3;
                                    kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.lLink &= MaskLink[base];
                                } else {
                                    kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.rLink &= MaskLink[j];
                                }
                                initialGraph[graphIndex].rGraph[j] = 0;
                                initialGraph[bI].processed = 1;
                                // remove link from end of branch edge
                                for(j = 0; j < 4; j++)
                                {
                                    if(initialGraph[bI].lGraph[j] == lI ) break;
                                }
                                if(kmerStack->kmerAddr[initialGraph[bI].stackOffset + initialGraph[bI].count -1].k.ku.ki.direction == 1)
                                {
                                    uint8_t base = (~j) & 0x3;
                                    kmerStack->kmerAddr[initialGraph[bI].stackOffset + initialGraph[bI].count -1].k.ku.ki.rLink = (1 << base);
                                } else {
                                    kmerStack->kmerAddr[initialGraph[bI].stackOffset + initialGraph[bI].count -1].k.ku.ki.lLink = (1 << j);
                                }
                                if(branchStack.cur_pos + 2 >= branchStack.size)
                                {
                                    branchStack.size <<= 1;
                                    branchStack.indexStack = (int*)realloc(branchStack.indexStack, branchStack.size * sizeof(int));
                                }
                                if(initialGraph[lI].state == 3)  branchStack.indexStack[branchStack.cur_pos++] = lI;
                                if(initialGraph[rI].state == 3)  branchStack.indexStack[branchStack.cur_pos++] = rI;
                            }
                        } else {
                            if(branchStack.cur_pos + 8 >= branchStack.size)
                            {
                                branchStack.size <<= 1;
                                branchStack.indexStack = (int*)realloc(branchStack.indexStack, branchStack.size * sizeof(int));
                            }
                            for(int i = 0; i < 4; i++)
                            {
                                if(initialGraph[bI].rGraph[i] > 0 && initialGraph[initialGraph[bI].rGraph[i]].state == 3) 
                                {  branchStack.indexStack[branchStack.cur_pos++] = initialGraph[bI].rGraph[i]; }
                                if(initialGraph[bI].lGraph[i] > 0 && initialGraph[initialGraph[bI].lGraph[i]].state == 3) 
                                {  branchStack.indexStack[branchStack.cur_pos++] = initialGraph[bI].lGraph[i]; }
                            }
                        }
                    } else if(remove_flag == 1 && initialGraph[bI].state == 2 ) {// remove bubble that backward trace
                        int j;
                        for(j = 0; j < initialGraph[bI].count -1; j++)
                        {
                            kmerStack->kmerAddr[initialGraph[bI].stackOffset + j].k.ku.ki.deleted = 1;
                        }
                        // delete last kmer
                        kmerStack->kmerAddr[initialGraph[bI].stackOffset + initialGraph[bI].count -1].hID = -1;
                        // removing graphIndex right link to bI
                        for(j = 0; j < 4; j++)
                        {
                            if(initialGraph[graphIndex].rGraph[j] == bI) break;
                        }
                        if(kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.direction == 1)
                        {
                            uint8_t base = (~j) & 0x3;
                            kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.lLink &= MaskLink[base];
                        } else {
                            kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.rLink &= MaskLink[j];
                        }
                        initialGraph[graphIndex].rGraph[j] = 0;
                        initialGraph[bI].processed = 1;
                    }
                } else {
                    if(branchStack.cur_pos + 4 >= branchStack.size)
                    {
                        branchStack.size <<= 1;
                        branchStack.indexStack = (int*)realloc(branchStack.indexStack, branchStack.size * sizeof(int));
                    }
                    for(int i = 0; i < 4; i++)
                    {
                        int gI = initialGraph[graphIndex].rGraph[i];
                        if(gI > 0 && initialGraph[gI].freq > curve.valley && initialGraph[gI].state == 3 )
                        {
                            branchStack.indexStack[branchStack.cur_pos++] = gI;
                        } else if(initialGraph[gI].freq > curve.valley && initialGraph[gI].state == 2) {
                            kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].hID = -1;
                            initialGraph[gI].processed = 1;
                        }
                    }
                }
                // test left graph and traverse left link
                tip_n = 0; bubble_n = 0; normalB = 0; 
                for(int i = 0; i < 4; i++)
                {
                    if(initialGraph[graphIndex].lGraph[i] > 0)
                    {
                        int gI = initialGraph[graphIndex].lGraph[i];
                        if(initialGraph[gI].freq > curve.valley && (initialGraph[gI].state != 4 || 
                            initialGraph[gI].count > 2 * KMER))
                        {
                            normalB++;
                            normalIndex = gI;
                        } else if(((initialGraph[gI].freq <= curve.valley && initialGraph[gI].count < 2 * KMER) ||
                            (initialGraph[gI].count <= KMER )) && initialGraph[gI].state == 4)  tip_n++;
                        else if(graphIndex > 0 && initialGraph[graphIndex].direction == FORWARD &&
                            initialGraph[gI].state == 2 && initialGraph[gI].freq <= curve.valley &&
                            MIN_TIMES * initialGraph[gI].freq < initialGraph[graphIndex].freq &&
                            initialGraph[gI].count < 2 * KMER) bubble_n++;
                        lbranch_n++;
                    }
                }
                // remove tips
                if((tip_n >= 1) || ( graphIndex > 0 && (initialGraph[graphIndex].direction == FORWARD && bubble_n == 1)))
                {
                    for(int i = 0; i < 4; i++)
                    {
                        if(initialGraph[graphIndex].lGraph[i] > 0)
                        {
                            int gI = initialGraph[graphIndex].lGraph[i];
                            if( (initialGraph[gI].freq <= curve.valley && initialGraph[gI].count < 2 * KMER &&
                                initialGraph[gI].state == 4) || (initialGraph[gI].count <= KMER && 
                                initialGraph[gI].state == 4) || (graphIndex > 0 && initialGraph[gI].state == 2 && 
                                initialGraph[graphIndex].direction == FORWARD && initialGraph[gI].freq <= curve.valley &&
                                initialGraph[gI].count < 2 * KMER && 
                                MIN_TIMES * initialGraph[gI].freq < initialGraph[graphIndex].freq) ) 
                            {
                                for(int j = 0; j < initialGraph[gI].count; j++)
                                {
                                    kmerStack->kmerAddr[initialGraph[gI].stackOffset + j].k.ku.ki.deleted = 1;
                                }
                                if(kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.direction == 1)
                                {
                                    uint8_t base = (~i) & 0x3;
                                    kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.rLink &= MaskLink[base];
                                } else {
                                    kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.lLink &= MaskLink[i];
                                }
                                initialGraph[graphIndex].lGraph[i] = 0;
                                initialGraph[gI].processed = 1;
                                if(bubble_n == 1) kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].hID = -1;
                            }
                        }
                    }
                    if(normalB >= 1)
                    {
                        if(branchStack.cur_pos + normalB >= branchStack.size)
                        {
                            branchStack.size <<= 1;
                            branchStack.indexStack = (int*)realloc(branchStack.indexStack, branchStack.size * sizeof(int));
                        }
                        for(int j = 0; j < 4; j++)
                        {
                            int tI = initialGraph[graphIndex].lGraph[j];
                            if(tI > 0 && initialGraph[tI].freq > curve.valley && initialGraph[tI].state == 3)
                            {
                                branchStack.indexStack[branchStack.cur_pos++] = tI;
                            } else if(initialGraph[tI].freq > curve.valley && initialGraph[tI].state == 2) {
                                kmerStack->kmerAddr[initialGraph[tI].stackOffset + initialGraph[tI].count -1].hID = -1;
                                initialGraph[tI].processed = 1;
                            }
                        }
                    }
                } else if(lbranch_n == 2 && normalB == 1) { // remove bubble and chimeric connection
                    int remove_flag = 0, bI; 
                    for(int i = 0; i < 4; i++)
                    {
                        if(initialGraph[graphIndex].lGraph[i] > 0)
                        {
                            int gI = initialGraph[graphIndex].lGraph[i];
                            if(initialGraph[gI].freq < curve.valley && 
                                MIN_TIMES * initialGraph[gI].freq < initialGraph[normalIndex].freq &&
                                MIN_TIMES * initialGraph[gI].freq < initialGraph[graphIndex].freq &&
                                initialGraph[gI].count < 2 * KMER)
                            {
                                remove_flag = 1;
                                bI = gI;
                            }
                        }
                    }
                    if(branchStack.cur_pos >= branchStack.size)
                    {
                        branchStack.size <<= 1;
                        branchStack.indexStack = (int*)realloc(branchStack.indexStack, branchStack.size * sizeof(int));
                    }
                    if(initialGraph[normalIndex].state == 3) branchStack.indexStack[branchStack.cur_pos++] = normalIndex;
                    else if(initialGraph[normalIndex].freq > curve.valley && initialGraph[normalIndex].state == 2) {
                        kmerStack->kmerAddr[initialGraph[normalIndex].stackOffset + initialGraph[normalIndex].count -1].hID = -1;
                        initialGraph[normalIndex].processed = 1;
                    }
                    if(remove_flag == 1 && initialGraph[bI].state ==3) // remove chimeric connection or bubble
                    {
                        int rI, lI, rn = 0, ln = 0;
                        for(int i = 0; i < 4; i++)
                        {
                            if(initialGraph[bI].lGraph[i] > 0) { lI = initialGraph[bI].lGraph[i]; ln++; }
                            if(initialGraph[bI].rGraph[i] > 0) { rI = initialGraph[bI].rGraph[i]; rn++; }
                        }
                        if(rn == 1 && ln == 1)
                        {
                            // remove edge 
                            if(MIN_TIMES * initialGraph[bI].freq < initialGraph[rI].freq &&
                                MIN_TIMES * initialGraph[bI].freq < initialGraph[lI].freq)
                            {
                                int j;
                                for(j = 0; j < initialGraph[bI].count - 1; j++)
                                {
                                    kmerStack->kmerAddr[initialGraph[bI].stackOffset + j].k.ku.ki.deleted = 1;
                                }
                                // remove link from beginning of branch edge 
                                for(j = 0; j < 4; j++)
                                {
                                    if(initialGraph[graphIndex].lGraph[j] == bI) break;
                                }
                                if(kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.direction == 1)
                                {
                                    uint8_t base = (~j) & 0x3;
                                    kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.rLink &= MaskLink[base];
                                } else {
                                    kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.lLink &= MaskLink[j];
                                }
                                initialGraph[graphIndex].lGraph[j] = 0;
                                initialGraph[bI].processed = 1;
                                // remove link from end of branch edge 
                                for(j = 0; j < 4; j++)
                                {
                                    if(initialGraph[bI].rGraph[j] == rI) break;
                                }
                                if(kmerStack->kmerAddr[initialGraph[bI].stackOffset + initialGraph[bI].count -1].k.ku.ki.direction == 1)
                                {
                                    uint8_t base = (~j) & 0x3;
                                    kmerStack->kmerAddr[initialGraph[bI].stackOffset + initialGraph[bI].count - 1].k.ku.ki.lLink = (1 << base);
                                } else {
                                    kmerStack->kmerAddr[initialGraph[bI].stackOffset + initialGraph[bI].count - 1].k.ku.ki.rLink = (1 << j);
                                }
                                if(branchStack.cur_pos + 2 >= branchStack.size)
                                {
                                    branchStack.size <<= 1;
                                    branchStack.indexStack = (int*)realloc(branchStack.indexStack, branchStack.size * sizeof(int));
                                }
                                if(initialGraph[rI].state == 3) branchStack.indexStack[branchStack.cur_pos++] = rI;
                                if(initialGraph[lI].state == 3) branchStack.indexStack[branchStack.cur_pos++] = lI;
                            }
                        } else {
                            if(branchStack.cur_pos + 8 >= branchStack.size)
                            {
                                branchStack.size <<= 1;
                                branchStack.indexStack = (int*)realloc(branchStack.indexStack, branchStack.size * sizeof(int));
                            }
                            for(int i = 0; i < 4; i++)
                            {
                                if(initialGraph[bI].lGraph[i] > 0 && initialGraph[initialGraph[bI].lGraph[i]].state == 3)
                                {   branchStack.indexStack[branchStack.cur_pos++] = initialGraph[bI].lGraph[i];   }
                                if(initialGraph[bI].rGraph[i] > 0 && initialGraph[initialGraph[bI].rGraph[i]].state == 3)
                                {   branchStack.indexStack[branchStack.cur_pos++] = initialGraph[bI].rGraph[i];   }
                            }
                        }
                    } else if(remove_flag == 1 && initialGraph[bI].state == 2) {// remove bubble that backward trace
                        int j;
                        for(j = 0; j < initialGraph[bI].count -1; j++)
                        {
                            kmerStack->kmerAddr[initialGraph[bI].stackOffset + j].k.ku.ki.deleted = 1;
                        }
                        // delete last kmer
                        kmerStack->kmerAddr[initialGraph[bI].stackOffset + initialGraph[bI].count -1].hID = -1;
                        // removing graphIndex left link to bI
                        for(j = 0; j < 4; j++)
                        {
                            if(initialGraph[graphIndex].lGraph[j] == bI) break;
                        }
                        if(kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.direction == 1)
                        {
                            uint8_t base = (~j) & 0x3;
                            kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.rLink &= MaskLink[base];
                        } else {
                            kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.lLink &= MaskLink[j];
                        }
                        initialGraph[graphIndex].lGraph[j] = 0;
                        initialGraph[bI].processed = 1;
                    }
                } else {
                    if(branchStack.cur_pos + 4 >= branchStack.size)
                    {
                        branchStack.size <<= 1;
                        branchStack.indexStack = (int*)realloc(branchStack.indexStack, branchStack.size * sizeof(int));
                    }
                    for(int i = 0; i < 4; i++)
                    {
                        int gI = initialGraph[graphIndex].lGraph[i];
                        if(gI > 0 && initialGraph[gI].freq > curve.valley && initialGraph[gI].state == 3)
                        {
                            branchStack.indexStack[branchStack.cur_pos++] = gI;
                        } else if(initialGraph[gI].freq > curve.valley && initialGraph[gI].state == 2) {
                            kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].hID = -1;
                            initialGraph[gI].processed = 1;
                        }
                    }
                }
                // check whether a secondary tip
                {
                    int branch_n = 0 ;
                    for(int i = 0 ; i < 4; i++)
                    {
                        if(initialGraph[graphIndex].rGraph[i] > 0) branch_n++;
                        if(initialGraph[graphIndex].lGraph[i] > 0) branch_n++;
                    }
                    if(graphIndex > 0 && branch_n == 0) // remove edge 
                    {
                        int gI = initialGraph[graphIndex].origin_graph_index ;
                        int oI = 0 , total_b = 0;
                        for(int i = 0; i < 4; i++)
                        {
                            if(initialGraph[gI].rGraph[i] > 0)
                            {
                                if(initialGraph[gI].rGraph[i] != graphIndex) oI = initialGraph[gI].rGraph[i];
                                total_b++;
                            }
                            if(initialGraph[gI].lGraph[i] > 0)
                            {
                                if(initialGraph[gI].lGraph[i] != graphIndex) oI = initialGraph[gI].lGraph[i];
                                total_b++;
                            }
                        }
                        if(gI > 0 && total_b == 2 && initialGraph[gI].direction == initialGraph[oI].direction &&
                            2 * initialGraph[graphIndex].freq <= initialGraph[gI].freq &&
                            2 * initialGraph[graphIndex].freq <= initialGraph[oI].freq && 
                            initialGraph[graphIndex].count < KMER)
                        {
                            for(int i = 0; i < initialGraph[graphIndex].count; i++)
                            {
                                kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + i].k.ku.ki.deleted = 1;
                            }
                            // remove link from beginning of branch edge 
                            for(int i = 0; i < 4; i++)
                            {
                                if(initialGraph[gI].rGraph[i] == graphIndex)
                                {
                                    if(kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].k.ku.ki.direction == 1)
                                    {
                                        uint8_t base = (~i) & 0x3;
                                        kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].k.ku.ki.lLink &= MaskLink[base];
                                    } else {
                                        kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count - 1].k.ku.ki.rLink &= MaskLink[i];
                                    }
                                    initialGraph[gI].rGraph[i] = 0;
                                    break;
                                }
                            }
                            for(int i = 0; i < 4; i++)
                            {
                                if(initialGraph[gI].lGraph[i] == graphIndex)
                                {
                                    if(kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].k.ku.ki.direction == 1)
                                    {
                                        uint8_t base = (~i) & 0x3;
                                        kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].k.ku.ki.rLink &= MaskLink[base];
                                    } else {
                                        kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count - 1].k.ku.ki.lLink &= MaskLink[i];
                                    }
                                    initialGraph[gI].lGraph[i] = 0 ;
                                    break;
                                }
                            }
                        }
                    }
                }
                if(initialGraph[graphIndex].state == 2)
                {
                    kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].hID = -1 ;
                }
                if(initialGraph[graphIndex].state == 3 && (rbranch_n + lbranch_n == 0))
                { // denote as virgin
                    kmerStack->kmerAddr[initialGraph[graphIndex].stackOffset + initialGraph[graphIndex].count -1].k.ku.ki.processed = 0 ;
                }
                initialGraph[graphIndex].processed = 1;
            }
        }
   }
} */ 


inline void writeUnitig2File(const uint8_t *contig, const SNPset *set, const int len , FILE *fp, const int freq, const int contigID )
{
    char base[] = {'A','C','G','T'};
    char *s = (char*)xcalloc(len + set->count * 6 + 100 , sizeof(char));
    int s_count = 0 ;
    for(int i = 0 , snp_count = 0 ; i < len ; i++)
    {
        if(set->info[snp_count].position == i)
        {
            s[s_count++] = '{' ; s[s_count++] = base[(contig[i/CHAR_PER_BYTE] >> ((CHAR_PER_BYTE - i%CHAR_PER_BYTE - 1)<<1)) & 0x3] ; s[s_count++] = ','; s[s_count++] = base[set->info[snp_count].SNP] ; s[s_count++] = '}' ;
            if(snp_count < set->count - 1) snp_count++ ;
        } else {
            s[s_count++] = base[(contig[i/CHAR_PER_BYTE] >> ((CHAR_PER_BYTE - i%CHAR_PER_BYTE - 1)<<1) & 0x3)];
        }
    }
    // write to file 
    fprintf(fp, ">%d\t%d\n", contigID, freq);
    fprintf(fp, "%s\n", s);

    // clean and free work 
    free(s);
}

int comparePos(const void *p, const void *q)
{
    SNPinfo *p1 = (SNPinfo*)p, *q1 = (SNPinfo*)q ;
    if(p1->position < q1->position) return -1 ;
    else if(p1->position > q1->position) return 1 ;
    else return 0 ;
}

int writeKmerToHashTable(KmerStack *kmerStack, HashTable *hashTable)
{
    // set kmer process flag
    for(int i = 0; i < kmerStack->count; i++)
    {
        if(kmerStack->kmerAddr[i].hID != -1)
        {
            KmerAddr ka = kmerStack->kmerAddr[i];
            ka.k.ku.ki.direction = 0 ;
            if(ka.k.ku.ki.curprocess == 1 || ka.k.ku.ki.processed == 1) ka.k.ku.ki.processed = 1 ;
            ka.k.ku.ki.curprocess = 0; 
            hashTable[ka.hID].table[ka.addr] = ka.k;
        }
    }
    return 0 ;
}

void resetHashTable(HashTable *hashTable, const int hashNumber, const int index)
{
	HashTable ht = hashTable[index];
    for(int64_t i = 0; i < ht.size; i++)
    {
		kmer_t nk = ht.table[i];
        if(nk.ku.ki.freq > 0 && nk.ku.ki.deleted == 0)
        {
            nk.ku.ki.processed = 0;
			if(LINK_BIGGER_THAN_ONE(nk.ku.ki.rLink) || LINK_BIGGER_THAN_ONE(nk.ku.ki.lLink))
			{
				kmer_t tk = ResetLinkInfo(nk, hashTable, hashNumber);
				if(memcmp(&tk, &nk, sizeof(kmer_t)) != 0) { nk = tk; }
			}
			hashTable[index].table[i] = nk ;
        }
    }
}

void addToQueue(kmer_graph *initialGraph, const int graphIndex,const int graphCount, BranchQueue *branchQueue, const uint8_t base, int direction, uint64_t *kmerSeq)
{
	if(branchQueue->end >= branchQueue->size)
	{
		branchQueue->size <<= 1 ;
		branchQueue->branchKmers = (BranchKmer*)xrecalloc(branchQueue->branchKmers, branchQueue->size * sizeof(BranchKmer));
	}
	if(direction == FORWARD)
	{
		initialGraph[graphIndex].rGraph[base] = graphCount;
		getNextKmer(branchQueue->branchKmers[branchQueue->end].initialKmer.kmer, kmerSeq, base, KMER);	
	} else { // BACKWARD
		initialGraph[graphIndex].lGraph[base] = graphCount;
		getPreviousKmer(branchQueue->branchKmers[branchQueue->end].initialKmer.kmer, kmerSeq, base, KMER);
	}	
	branchQueue->branchKmers[branchQueue->end].direction = direction ;
	branchQueue->branchKmers[branchQueue->end].kmer_graph_index = graphCount ;
	branchQueue->branchKmers[branchQueue->end].origin_graph_index = graphIndex;
	branchQueue->branchKmers[branchQueue->end].depth = initialGraph[graphIndex].depth ;
	branchQueue->end++;
}


// return graphCount ;
int addBranch(kmer_graph *initialGraph, const int graphIndex, int graphCount, BranchQueue *branchQueue, const int direction, kmer_t curkmer, kmer_t prekmer)
{
	kmer_t tk;
	if(prekmer.ku.ki.direction == 1)
	{
		getRevKmer(prekmer.kmer, tk.kmer, KMER);
		prekmer = tk;
	}

	int len = (KMER + 32 -1)/32 ;
	kmer_info ki = curkmer.ku.ki;
	
	if(ki.direction == 1)
	{
		getRevKmer(curkmer.kmer, tk.kmer, KMER);
		curkmer = tk ;
		for(uint8_t z = 0; z < 4; z++)
		{
			if(((ki.rLink & (1<<z)) > 0) && ((direction == BACKWARD) || (direction == BIDIRECTION)))
			{
				uint8_t base = (~z) & 0x3 ;
				getPreviousKmer(tk.kmer, curkmer.kmer, base, KMER);
				if(memcmp(tk.kmer, prekmer.kmer, sizeof(uint64_t) * len) != 0)
				{
					addToQueue(initialGraph,graphIndex, graphCount, branchQueue, base, BACKWARD, curkmer.kmer); 
					graphCount++ ;
				}
			}

			if(((ki.lLink & (1<<z)) > 0) && ((direction == FORWARD || direction == BIDIRECTION)))
			{
				uint8_t base = (~z) & 0x3 ;
				getNextKmer(tk.kmer, curkmer.kmer, base, KMER);
				if(memcmp(tk.kmer, prekmer.kmer, sizeof(uint64_t) * len) != 0)
				{
					addToQueue(initialGraph, graphIndex, graphCount, branchQueue, base, FORWARD, curkmer.kmer);
					graphCount++ ;
				}
			}
		}
	} else {
		for(uint8_t z = 0; z < 4; z++)
		{
			if(((ki.rLink & (1<<z)) > 0) && ((direction == FORWARD || direction == BIDIRECTION)))
			{
				uint8_t base = z ;
				getNextKmer(tk.kmer, curkmer.kmer, base, KMER);
				if(memcmp(tk.kmer, prekmer.kmer, sizeof(uint64_t) * len) != 0)
				{
					addToQueue(initialGraph, graphIndex, graphCount, branchQueue, base, FORWARD, curkmer.kmer);
					graphCount++;
				}
			}

			if(((ki.lLink & (1<<z)) > 0) && ((direction == BACKWARD) || (direction == BIDIRECTION)))
			{
				uint8_t base = z ;
				getPreviousKmer(tk.kmer, curkmer.kmer, base, KMER);
				if(memcmp(tk.kmer, prekmer.kmer, sizeof(uint64_t) * len) != 0)
				{
					addToQueue(initialGraph, graphIndex, graphCount, branchQueue, base, BACKWARD, curkmer.kmer);
					graphCount++;
				}
			}
		}	
	}

	return graphCount ;
}

/*
int extractUnitig(kmer_graph *initialGraph, KmerStack *kmerStack, const KmerFreqCurve *curve, const int unitigID, gzFile fz, const int min_unitig_len)
{
    int rlen = 5000, llen = 5000, rcount = 0, lcount = 0 ;
    char *runitig  = (char*)xcalloc(rlen, sizeof(char)), *lunitig = (char*)xcalloc(llen, sizeof(char));
    int gI = 0, totalFreq = 0;
    // extract right link unitig
    while(1)
    {
        int rn = 0, ln = 0;
		int lIndex, rIndex  ;
        for(int i = 0 ; i < 4; i++)
        {
            if(initialGraph[gI].rGraph[i] > 0) { rn++; rIndex = initialGraph[gI].rGraph[i]; }
			if(initialGraph[gI].lGraph[i] > 0) { ln++; lIndex = initialGraph[gI].lGraph[i]; }
        }
		// process reverse direction branch 
		if(gI > 0 && ln >= 1) 
		{
			//initialGraph[gI].processed = 1 ;
			if(ln == 1 && rn == 1 &&  initialGraph[lIndex].count < 2 * KMER && MIN_TIMES * initialGraph[lIndex].freq < initialGraph[gI].freq &&
				abs(initialGraph[rIndex].freq - curve->peak) < curve->peak * PEAK_SD)
			{ 
				// add gI edge_contig to runitig 
				if(rcount + initialGraph[rIndex].count + 10 > rlen)
				{
					int old = rlen ;
					rlen <<= 1; rlen += initialGraph[rIndex].count ;
					runitig = (char*)realloc(runitig, rlen * sizeof(char));
					memset(runitig + old, 0 , (rlen - old) * sizeof(char));
				}
				bit64ToChar(runitig + rcount, initialGraph[rIndex].edge_contig, initialGraph[rIndex].count);
				rcount += initialGraph[rIndex].count;
				totalFreq += initialGraph[rIndex].count * initialGraph[rIndex].freq ;
				initialGraph[lIndex].processed = 1;	initialGraph[rIndex].processed = 1; 
				if(initialGraph[rIndex].state == 3) { gI = rIndex; continue; }
				else break;
			} else { break; }
		}
        if(rn == 2 ) 
        {  
            int g1 = 0, g2 = 0 ;
            for(int i = 0; i < 4; i++)
            {
                if(initialGraph[gI].rGraph[i] > 0) 
                { 
                    if(g1 == 0) g1 = initialGraph[gI].rGraph[i]; 
                    else g2 = initialGraph[gI].rGraph[i];
                }
            }
            if(rcount + initialGraph[g1].count + initialGraph[g2].count + 10 > rlen)
            {
				int old = rlen ;
                rlen <<= 1; rlen += initialGraph[g1].count + initialGraph[g2].count;
                runitig = (char*)realloc(runitig, rlen * sizeof(char));
                memset(runitig + old, 0, (rlen - old) * sizeof(char));
            }
            // make sure freq of g1 is the smaller than freq of g2
            if(initialGraph[g1].freq > initialGraph[g2].freq) { int tmp  = g1 ;  g1 = g2 ; g2 = tmp ; }
            if(initialGraph[g1].freq <= curve->valley && initialGraph[g1].count < 2 * KMER &&
                abs(initialGraph[g2].freq - curve->peak) < curve->peak * PEAK_SD &&
                MIN_TIMES * initialGraph[g1].freq < initialGraph[g2].freq)
            {   // process tips and chimeric link or bubble that cause by middle sequence error of the reads
                bit64ToChar(runitig + rcount, initialGraph[g2].edge_contig,initialGraph[g2].count);
                rcount += initialGraph[g2].count;
                totalFreq += initialGraph[g2].count * initialGraph[g2].freq ;
				initialGraph[g1].processed = 1; initialGraph[g2].processed = 1; 
                // update new graph index(gI)
                {
                    int n = 0 ;
					if((initialGraph[g1].state == 3 && initialGraph[g2].state == 3) || 
						(initialGraph[g1].state == 4 && initialGraph[g2].state ==3))
					{
						gI = g2 ;
					} else if((initialGraph[g1].state == 1 && initialGraph[g2].state == 3) || 
							(initialGraph[g2].state == 1 && initialGraph[g1].state == 3)){
						for(int j = 0; j < 4; j++)
						{
							if(initialGraph[g1].rGraph[j] > 0) { gI = initialGraph[g1].rGraph[j]; n++; }
							if(initialGraph[g2].rGraph[j] > 0) { gI = initialGraph[g2].rGraph[j]; n++; }
						}
						if(n == 1 && abs(initialGraph[gI].freq - curve->peak) < curve->peak * PEAK_SD)
						{
							bit64ToChar(runitig + rcount, initialGraph[gI].edge_contig, initialGraph[gI].count);
							rcount += initialGraph[gI].count;
							totalFreq += initialGraph[gI].count * initialGraph[gI].freq;
						} else break;						
					} else break;
                }
            } else if(initialGraph[g1].freq > curve->peak * (1 - PEAK_SD) && 
                        initialGraph[g2].freq > curve->peak * (1 - PEAK_SD)) {
                // process self-cycle branch
                // make sure g1 is the self-cycle branch
                if(initialGraph[g2].state == 1) {int tmp = g1; g1 = g2; g2 = tmp; }
                if(initialGraph[g1].state == 1 && memcmp(&kmerStack->kmerAddr[initialGraph[g1].stackOffset + initialGraph[g1].count -1], &kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1], sizeof(KmerAddr)) == 0 && initialGraph[g1].count < KMER && initialGraph[g2].freq < curve->peak * (1 + PEAK_SD))
                {
                    float cov = (float)initialGraph[g1].freq / initialGraph[g2].freq  ;
                    if(rcount + initialGraph[g1].count * ((int)cov + 1) + 10 + initialGraph[g2].count > rlen)
                    {
						int old = rlen;
                        rlen <<= 1 ; rlen += initialGraph[g2].count;
                        runitig = (char*)realloc(runitig, rlen * sizeof(char));
                        memset(runitig + old, 0, (rlen - old) * sizeof(char));
                    }
                    if((int)(cov + PEAK_SD) == (int)(cov + 2 * PEAK_SD) && initialGraph[g1].freq < 255)
                    {
                        int t = (int)(cov + PEAK_SD);
                        for(int j = 0; j < t; j++)
                        {
                            bit64ToChar(runitig + rcount, initialGraph[g1].edge_contig, initialGraph[g1].count);
                            rcount += initialGraph[g1].count;
                        }
                    }else {
						if(cov < 1) cov = 1;
                        for(int j = 0; j < (int)cov; j++)
                        {
                            bit64ToChar(runitig + rcount, initialGraph[g1].edge_contig, initialGraph[g1].count);
                            rcount += initialGraph[g1].count;
                        }
                        // add uncertain position
                        runitig[rcount++] = '{' ; runitig[rcount++] = ',';
                        bit64ToChar(runitig + rcount, initialGraph[g1].edge_contig, initialGraph[g1].count);
                        rcount += initialGraph[g1].count ;
						if(initialGraph[g1].freq == 255)
						{
							runitig[rcount++] = ',';
							bit64ToChar(runitig + rcount, initialGraph[g1].edge_contig, initialGraph[g1].count);
							rcount += initialGraph[g1].count ;
						}
                        runitig[rcount++] = '}' ;
                    }
					totalFreq += initialGraph[g1].count * initialGraph[g1].freq ;
                    bit64ToChar(runitig + rcount, initialGraph[g2].edge_contig, initialGraph[g2].count);
                    rcount += initialGraph[g2].count ;
                    totalFreq += initialGraph[g2].count * initialGraph[g2].freq;
					initialGraph[g1].processed = 1;	initialGraph[g2].processed = 1; 
                    // update new graphIndex
					gI = g2 ;
                } else break; // can't process
            } else if(abs(initialGraph[g1].freq - curve->peak/2) < (curve->peak/2) * PEAK_SD &&
                        abs(initialGraph[g2].freq - curve->peak/2) < (curve->peak/2) * PEAK_SD &&
                        initialGraph[g1].count < 2 * KMER && initialGraph[g2].count < 2 * KMER && 
                        memcmp(&kmerStack->kmerAddr[initialGraph[g1].stackOffset + initialGraph[g1].count -1], &kmerStack->kmerAddr[initialGraph[g2].stackOffset + initialGraph[g2].count -1], sizeof(KmerAddr)) == 0) {
                // process the SNP or polymorph by the sample
                int max_diff = 3 ,diff = 0, i , j, equal_flag, c = 0, merge_flag = 0 ;
                char s[256], unitig_g1[256], unitig_g2[256] ; 
                memset(s, 0 , 256); memset(unitig_g1, 0 , 256); memset(unitig_g2, 0, 256);
                bit64ToChar(unitig_g1, initialGraph[g1].edge_contig, initialGraph[g1].count);
                bit64ToChar(unitig_g2, initialGraph[g2].edge_contig, initialGraph[g2].count);
                if(initialGraph[g1].count == initialGraph[g2].count) equal_flag = 1 ;
                else equal_flag = 0 ;
                for(i = 0 , j = 0; i < initialGraph[g1].count && j < initialGraph[g2].count; i++, j++)
                {
                    if(unitig_g1[i] == unitig_g2[j])
                    {
                        s[c++] = unitig_g1[i];
                    } else {
                        if( i + 4 < initialGraph[g1].count && j + 4 < initialGraph[g2].count && 
                            unitig_g1[i+1] == unitig_g2[j+1] &&  unitig_g1[i+2] == unitig_g2[j+2] && unitig_g1[i+3] == unitig_g2[j+3] &&
							unitig_g1[i+4] == unitig_g2[j+4] && diff < max_diff)
                        {
                            s[c++] = '{' ; s[c++] = unitig_g1[i]; s[c++] = ','; 
                            s[c++] = unitig_g2[j]; s[c++] = '}';
                            diff++;
                        } else {
                            if(equal_flag == 0)
                            {
                                char last_base ;
                                if(c == 0) last_base = runitig[rcount-1];
                                else last_base = s[c-1];
                                if(initialGraph[g1].count > initialGraph[g2].count && unitig_g1[i] == last_base)
                                {
                                    int z = i + 1 ;
                                    while(z < initialGraph[g1].count)
                                    {
                                        if(unitig_g1[z] == last_base && diff < max_diff) { z++; diff++; }
                                        else if(unitig_g1[z] == unitig_g2[j]) break;
                                        else { merge_flag = 1; break; }
                                    }
                                    if(merge_flag == 0)
                                    {
                                        s[c++] = '{';  s[c++] = ',';
                                        for(int y = i ; y < z; y++)
                                        {
                                            s[c++] = unitig_g1[y];
                                        }
                                        s[c++] = '}';
                                    } else { merge_flag = 1; break; }
                                }else if(initialGraph[g1].count < initialGraph[g2].count && unitig_g2[j] == last_base){
                                    int z = j + 1;
                                    while(z < initialGraph[g2].count)
                                    {
                                        if(unitig_g2[z] == last_base && diff < max_diff) { z++; diff++; }
                                        else if(unitig_g2[z] == unitig_g1[i]) break;
                                        else { merge_flag = 1; break; }
                                    }
                                    if(merge_flag == 0)
                                    {
                                        s[c++] = '{'; s[c++] = ',';
                                        for(int y = j; y < z; y++)
                                        {
                                            s[c++] = unitig_g2[y];
                                        }
                                        s[c++] = '}';
                                    } else { merge_flag = 1;  break; }
                                } else { merge_flag = 1; break; }
                            } else { merge_flag = 1; break; }
                        }
                    }
                }
                // output sequence
                if(merge_flag == 0)
                {
                    strcat(runitig, s);
                    rcount += c ;
                    totalFreq += c * (initialGraph[g1].freq + initialGraph[g2].freq);
                }else {
                    runitig[rcount++] = '{' ;
                    strcat(runitig, unitig_g1); rcount += initialGraph[g1].count ;
                    runitig[rcount++] = ',';
                    strcat(runitig, unitig_g2); rcount += initialGraph[g2].count ;
                    runitig[rcount++] = '}';
                    totalFreq += initialGraph[g1].count * initialGraph[g1].freq + initialGraph[g2].count * initialGraph[g2].freq;
                }
				initialGraph[g1].processed = 1; initialGraph[g2].processed = 1; 
                // update new graphIndex 
                {
                    int n = 0 , newg ;
                    for(int j = 0; j < 4; j++)
                    {
                        if(initialGraph[g1].rGraph[j] > 0) { n++; newg = initialGraph[g1].rGraph[j]; }
                        if(initialGraph[g2].rGraph[j] > 0) { n++; newg = initialGraph[g2].rGraph[j]; }
                    }
                    if(n == 1 && initialGraph[newg].freq > curve->valley && 
                        initialGraph[newg].freq  < curve->peak * (1 + PEAK_SD))
                    {
						bit64ToChar(runitig + rcount, initialGraph[newg].edge_contig, initialGraph[newg].count);
						rcount += initialGraph[newg].count;
						totalFreq += initialGraph[newg].count * initialGraph[newg].freq;
						gI = newg ;
						initialGraph[newg].processed = 1; 
                    }else {  break; }              
                }
            } else break; // can't process
        } else if(rn == 1) {
            if(rcount + initialGraph[rIndex].count + 10 > rlen)
            {
				int old = rlen ;
                rlen <<= 1; rlen += initialGraph[rIndex].count ;
                runitig = (char*)realloc(runitig, rlen * sizeof(char));
                memset(runitig + old, 0 , (rlen - old) * sizeof(char));
            }
            if(abs(initialGraph[rIndex].freq - curve->peak) < curve->peak * PEAK_SD)
            {
                bit64ToChar(runitig + rcount, initialGraph[rIndex].edge_contig, initialGraph[rIndex].count);
                rcount += initialGraph[rIndex].count;
                totalFreq += initialGraph[rIndex].count * initialGraph[rIndex].freq ;
				initialGraph[rIndex].processed = 1; 
                if(initialGraph[rIndex].state == 3) gI = rIndex;
                else break;
            } else { break; }
        } else {  break; } // n == 0 or n > 2 
    }
    // extract left link unitig
    gI = 0 ;
    while(1)
    {
        int rn = 0, ln = 0 ;
		int lIndex, rIndex;
        for(int i = 0; i < 4; i++)
        {
            if(initialGraph[gI].lGraph[i] > 0) { ln++; lIndex = initialGraph[gI].lGraph[i]; }
			if(initialGraph[gI].rGraph[i] > 0) { rn++; rIndex = initialGraph[gI].rGraph[i]; }
        }
		// process reverse direction branch
		if(gI > 0 && rn >= 1)
		{
			if(ln == 1 && rn == 1 && initialGraph[rIndex].count < 2 * KMER && MIN_TIMES * initialGraph[rIndex].freq < initialGraph[gI].freq &&
				abs(initialGraph[lIndex].freq - curve->peak) < curve->peak * PEAK_SD)
			{
				//add gI edge_contig to lunitig
				if(lcount + initialGraph[lIndex].count + 10 > llen)
				{
					int old = llen;
					llen <<= 1; llen += initialGraph[lIndex].count;
					lunitig = (char*)realloc(lunitig, llen * sizeof(char));
					memmove(lunitig + llen - lcount, lunitig + old - lcount, lcount * sizeof(char));
				}
				bit64ToChar(lunitig + llen - lcount - initialGraph[lIndex].count , initialGraph[lIndex].edge_contig, initialGraph[lIndex].count);
				lcount += initialGraph[lIndex].count;
				totalFreq += initialGraph[lIndex].count * initialGraph[lIndex].freq ;
				initialGraph[rIndex].processed = 1; initialGraph[lIndex].processed = 1;
				if(initialGraph[lIndex].state == 3) { gI = lIndex; continue; }
				else break;
			} else break;
		}
        if(ln == 2)
        {
            int g1 = 0 , g2 = 0;
            for(int i = 0; i < 4; i++)
            {
                if(initialGraph[gI].lGraph[i] > 0)
                {
                    if(g1 == 0) g1 = initialGraph[gI].lGraph[i];
                    else g2 = initialGraph[gI].lGraph[i];
                }
            }
            if(lcount + initialGraph[g1].count + initialGraph[g2].count + 10 > llen)
            {
				int old = llen;
                llen <<= 1; llen += initialGraph[g1].count + initialGraph[g2].count;
                lunitig = (char*)realloc(lunitig, llen * sizeof(char));
                memmove(lunitig + llen - lcount, lunitig + old - lcount, lcount );
            }
            // make sure g1 is the smaller than g2 
            if(initialGraph[g1].freq > initialGraph[g2].freq) { int tmp = g1; g1 = g2; g2 = tmp; }
            if(initialGraph[g1].freq <= curve->valley && initialGraph[g1].count < 2 * KMER &&
                abs(initialGraph[g2].freq - curve->peak) < curve->peak * PEAK_SD &&
                MIN_TIMES * initialGraph[g1].freq < initialGraph[g2].freq)
            {   // process tips and chimeric link or bubble that cause by middle sequence error of the reads
                
                bit64ToChar(lunitig + llen - lcount - initialGraph[g2].count, initialGraph[g2].edge_contig, initialGraph[g2].count); 
                lcount += initialGraph[g2].count;
                totalFreq += initialGraph[g2].count * initialGraph[g2].freq ;
				initialGraph[g1].processed = 1;	initialGraph[g2].processed = 1;
				// update new graph index (gI)
                {
                    int n = 0 ;
					if((initialGraph[g1].state == 3 && initialGraph[g2].state == 3) || 
						(initialGraph[g1].state == 4 && initialGraph[g2].state ==3))
					{
						gI = g2 ;
					} else if((initialGraph[g1].state == 1 && initialGraph[g2].state == 3) || 
							(initialGraph[g2].state == 1 && initialGraph[g1].state == 3)){
						for(int j = 0; j < 4; j++)
						{
							if(initialGraph[g1].lGraph[j] > 0) { gI = initialGraph[g1].lGraph[j]; n++; }
							if(initialGraph[g2].lGraph[j] > 0) { gI = initialGraph[g2].lGraph[j]; n++; }
						}
						if(n == 1 && abs(initialGraph[gI].freq - curve->peak) < curve->peak * PEAK_SD)
						{
							bit64ToChar(lunitig + llen - lcount - initialGraph[gI].count, initialGraph[gI].edge_contig, initialGraph[gI].count);
							lcount += initialGraph[gI].count;
							totalFreq += initialGraph[gI].count * initialGraph[gI].freq ;
						} else break;
					} else break;
                }
            } else if(initialGraph[g1].freq > curve->peak * (1 - PEAK_SD) &&
                        initialGraph[g2].freq > curve->peak * (1 - PEAK_SD)) {
                // process seft-cycle branch
				// make sure g1 is the self-cycle branch
				if(initialGraph[g2].state == 1) { int tmp = g1; g1 = g2; g2 = tmp; }
				if(initialGraph[g1].state == 1 && memcmp(&kmerStack->kmerAddr[initialGraph[g1].stackOffset + initialGraph[g1].count -1], &kmerStack->kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1], sizeof(KmerAddr)) == 0 && initialGraph[g1].count < KMER && initialGraph[g2].freq < curve->peak * (1 + PEAK_SD))
				{
					float cov = (float)initialGraph[g1].freq / initialGraph[g2].freq ;
					if(lcount + initialGraph[g1].count * ((int)cov + 1) + 10 + initialGraph[g2].count > llen)
					{
						int old = llen;
						llen <<= 1; llen += initialGraph[g2].count;
						lunitig = (char*)realloc(lunitig, llen * sizeof(char));
						memmove(lunitig + llen - lcount , lunitig + old - lcount, lcount * sizeof(char));
					}
					if((int)(cov + PEAK_SD) == (int)(cov + 2 * PEAK_SD) && initialGraph[g1].freq < 255)
					{
						int t = (int)(cov + PEAK_SD);
						for(int j =0; j < t; j++)
						{
							bit64ToChar(lunitig + llen - lcount - initialGraph[g1].count, initialGraph[g1].edge_contig, initialGraph[g1].count);
							lcount += initialGraph[g1].count;
						}
				 	} else {
						// add uncertain position
						lunitig[llen - lcount - initialGraph[g1].count - 3] = '{'; lunitig[llen - lcount - initialGraph[g1].count -2] = ',';
						bit64ToChar(lunitig + llen - lcount - initialGraph[g1].count - 1, initialGraph[g1].edge_contig, initialGraph[g1].count);
						lunitig[llen - lcount - 1] = '}';
						lcount += initialGraph[g1].count + 3 ;
						if(cov < 1) cov = 1;
						for(int j = 0; j < (int)cov; j++)
						{
							bit64ToChar(lunitig + llen - lcount - initialGraph[g1].count, initialGraph[g1].edge_contig, initialGraph[g1].count);
							lcount += initialGraph[g1].count;
						}
					}
					bit64ToChar(lunitig + llen - lcount - initialGraph[g2].count, initialGraph[g2].edge_contig, initialGraph[g2].count);
					lcount += initialGraph[g2].count ;
                    totalFreq += initialGraph[g1].count * initialGraph[g1].freq ;
                    totalFreq += initialGraph[g2].count * initialGraph[g2].freq ;
					initialGraph[g1].processed = 1;	initialGraph[g2].processed = 1;
					// update new graphIndex[gI]
					gI = g2;
				} else break; // can't process
            } else if(abs(initialGraph[g1].freq - curve->peak/2) < (curve->peak/2) * PEAK_SD &&
                        abs(initialGraph[g2].freq - curve->peak/2) < (curve->peak/2) * PEAK_SD &&
                        initialGraph[g1].count < 2 * KMER && initialGraph[g2].count < 2 * KMER &&
                        memcmp(&kmerStack->kmerAddr[initialGraph[g1].stackOffset + initialGraph[g1].count -1], &kmerStack->kmerAddr[initialGraph[g2].stackOffset + initialGraph[g2].count -1], sizeof(KmerAddr)) == 0) {
                // process the SNP or polymorph by the sample
				int max_diff = 3, diff = 0, i , j, equal_flag, c = 0, merge_flag = 0;
				int len = 256;
				char s[len],unitig_g1[len], unitig_g2[len];
				memset(s, 0, len); memset(unitig_g1, 0, len); memset(unitig_g2, 0, len);
				bit64ToChar(unitig_g1, initialGraph[g1].edge_contig, initialGraph[g1].count);
				bit64ToChar(unitig_g2, initialGraph[g2].edge_contig, initialGraph[g2].count);
				if(initialGraph[g1].count == initialGraph[g2].count) equal_flag = 1;
				else equal_flag = 0;
				for(i = initialGraph[g1].count -1, j = initialGraph[g2].count -1; i >= 0 && j >= 0; i--, j--)
				{
					if(unitig_g1[i] == unitig_g2[j])
					{
						s[len - (++c)] = unitig_g1[i];
					} else {
						if(i - 4 >= 0 && i - 4 >= 0 && unitig_g1[i-1] == unitig_g2[j-1] && unitig_g1[i-2] == unitig_g2[j -2] && 
							unitig_g1[i-3] == unitig_g2[j-3] && unitig_g1[i-4] == unitig_g2[j-4] && diff < max_diff)
						{
							s[len - c - 5] = '{'; s[len - c - 4] = unitig_g1[i]; s[len - c - 3] = ',' ;
							s[len - c - 2] = unitig_g2[j]; s[len - c -1] = '}';
							c += 5  ;
							diff++;
						} else {
							if(equal_flag == 0)
							{
								char last_base;
								if(c == 0) last_base = lunitig[llen - lcount];
								else last_base = s[len - c];
								if(initialGraph[g1].count > initialGraph[g2].count && unitig_g1[i] == last_base)
								{
									int z = i - 1;
									while(z > 0)
									{
										if(unitig_g1[z] == last_base && diff < max_diff) { z--; diff++; }
                                        else if(unitig_g1[z] == unitig_g2[j]) break;
                                        else {merge_flag = 1; break; }
									}
                                    if(merge_flag == 0)
                                    {
                                        s[len - (++c)] = '}' ;
                                        for(int y = i ; y > z; y--)
                                        {
                                            s[len - (++c)] = unitig_g1[y];
                                        }
                                        s[len - (++c)] = ','; s[len - (++c)] = '{';
                                    } else { merge_flag = 1; break;  }
								} else if(initialGraph[g1].count < initialGraph[g2].count && unitig_g2[j] == last_base) {
                                    int z = j - 1 ;
                                    while(z  > 0)
                                    {
                                        if(unitig_g2[z] == last_base && diff < max_diff) { z--; diff++; }
                                        else if(unitig_g2[z] == unitig_g1[i]) break;
                                        else {merge_flag = 1; break;  }
                                    }
                                    if(merge_flag == 0)
                                    {
                                        s[len - (++c)] = '}';
                                        for(int y = j; y > z; y--)
                                        {
                                            s[len - (++c)] = unitig_g2[y];
                                        }
                                        s[len - (++c)] = ','; s[len - (++c)] = '{';
                                    } else { merge_flag = 1; break; }
                                } else { merge_flag = 1; break; }
							} else {merge_flag = 1; break; }
						}
					}
				}
                // output sequence
                if(merge_flag == 0)
                {
                    memcpy(lunitig + llen - lcount - c, s + len - c, c * sizeof(char));
                    lcount += c ;
                    totalFreq +=  c * (initialGraph[g1].freq + initialGraph[g2].freq);
                } else {
                    lunitig[llen - lcount++] = '}' ;
                    memcpy(lunitig + llen - lcount - initialGraph[g1].count, unitig_g1, initialGraph[g1].count);
                    lcount += initialGraph[g1].count; lunitig[llen - lcount++] = ',';
                    memcpy(lunitig + llen - lcount - initialGraph[g2].count, unitig_g2, initialGraph[g2].count);
                    lcount += initialGraph[g2].count; lunitig[llen - lcount++] = '{';
                    totalFreq += initialGraph[g1].count * initialGraph[g1].freq + initialGraph[g2].count * initialGraph[g2].freq;
                }
				initialGraph[g1].processed = 1;	initialGraph[g2].processed = 1;
				// update new graphIndex
                {
                    int n = 0, newg;
                    for(int j = 0; j < 4; j++)
                    {
                        if(initialGraph[g1].lGraph[j] > 0) { n++; newg = initialGraph[g1].lGraph[j]; }
                        if(initialGraph[g2].lGraph[j] > 0) { n++; newg = initialGraph[g2].lGraph[j]; }
                    }
                    if(n == 1 && initialGraph[newg].freq > curve->valley && 
                        initialGraph[newg].freq < curve->peak * (1 + PEAK_SD))
                    {
						bit64ToChar(lunitig + llen - lcount - initialGraph[newg].count, initialGraph[newg].edge_contig, initialGraph[newg].count);
						lcount += initialGraph[newg].count;
						totalFreq += initialGraph[newg].count * initialGraph[newg].freq;
						gI = newg ;
						initialGraph[newg].processed = 1;
                    } else {break; }
                } 
            } else break; // can't process
        } else if(ln ==1) {
            if(lcount + initialGraph[lIndex].count + 10 > llen)
            {
				int old = llen;
                llen <<= 1; llen += initialGraph[lIndex].count;
                lunitig = (char*)realloc(lunitig, llen * sizeof(char));
                memmove(lunitig + llen - lcount, lunitig + old - lcount, lcount);
            }
            if(abs(initialGraph[lIndex].freq - curve->peak) < curve->peak * PEAK_SD)
            {
                bit64ToChar(lunitig + llen - lcount - initialGraph[lIndex].count, initialGraph[lIndex].edge_contig, initialGraph[lIndex].count);
                lcount += initialGraph[lIndex].count;
                totalFreq += initialGraph[lIndex].count * initialGraph[lIndex].freq;
				initialGraph[lIndex].processed = 1;
                if(initialGraph[lIndex].state == 3) gI = lIndex;
                else break;
            } else { break; }
        } else { break; } // n == 0 or n > 2
    }
	initialGraph[0].processed = 1;
    // write to the gz file
    int write_flag = 0 ;
    {
        if(lcount + KMER + rcount  > min_unitig_len)
        {
            char s_head[256];
            char *s = (char*)xcalloc(lcount  + KMER + rcount + 20, sizeof(char));
            int c = 0 ;
            memset(s_head, 0 , 256);
            if(lcount > 1) { memcpy(s + c, lunitig + llen - lcount + 1, lcount -1); c += lcount - 1 ; }
            bit64ToChar(s + c, initialGraph[0].edge_contig, KMER);
            c += KMER;
            if(rcount > 1) { memcpy(s + c, runitig, rcount -1 ); c += rcount - 1; }
            s[c++] = '\n';
            sprintf(s_head, ">unitig%d\t%d\t%d\n", unitigID, c, totalFreq/c);
            gzwrite(fz, s_head, strlen(s_head));
            gzwrite(fz, s, strlen(s));
			free(s);
        } else { write_flag = 1; }
    }
	// free and clean work
	free(runitig); free(lunitig);
    return write_flag;
} */

int IsDeleted(const uint64_t *kmer, const HashTable *hashTable, const int hashNumber)
{
	int len = (KMER + 32 -1)/32 ;
	uint64_t rkmer[len];
	int flag = 0, n ;
	uint64_t hashAdd ;
	
	getRevKmer(kmer, rkmer, KMER) ;
	for(int j = 0; j < len; j++)
	{
		if(kmer[j] > rkmer[j]) { flag = 1 ; break; }
		else if(kmer[j] < rkmer[j]) break ;
	}

	if(flag == 1)
	{
		n = get_time_hash(rkmer, len) % hashNumber ;
		hashAdd = getHashAddr(&hashTable[n], rkmer, len) ;
	} else {
		n = get_time_hash(kmer, len) % hashNumber ;
		hashAdd = getHashAddr(&hashTable[n], kmer, len) ;
	}

	if(hashTable[n].table[hashAdd].ku.ki.deleted == 1) return 1 ;
	else return 0 ;
}

kmer_t  ResetLinkInfo(kmer_t k, HashTable *hashTable, const int hashNumber)
{
	kmer_info ki = k.ku.ki ;
	
	for(uint8_t z = 0; z < 4; z++)
	{
		if((ki.rLink & (1<<z)) > 0) 
		{
			kmer_t tk ;
			getNextKmer(tk.kmer, k.kmer, z, KMER);
		    if(IsDeleted(tk.kmer, hashTable, hashNumber))
			{
				k.ku.ki.rLink &= MaskLink[z] ;	
			} 	
		}
	}

	for(uint8_t z = 0; z < 4; z++)
	{
		if((ki.lLink & (1<<z)) > 0)
		{
			kmer_t tk ;
			getPreviousKmer(tk.kmer, k.kmer, z, KMER) ;
			if(IsDeleted(tk.kmer, hashTable, hashNumber))
			{
				k.ku.ki.lLink &= MaskLink[z];
			} 
		}
	}
	
	return k ;
}

// s_flag note step of removetips or simplify step 
static inline int getGraphEdge(BranchQueue *branchQueue, kmer_graph *initialGraph, KmerStack *kmerStack, HashTable *hashTable, const int hashNumber, const int s_flag)
{
    BranchKmer bk = branchQueue->branchKmers[branchQueue->start];
	branchQueue->start++;
    kmer_t k = bk.initialKmer;
    int len = (KMER + 32 -1)/32, i;
    // set depth
    int d  ;
	if(bk.depth < KMER) d = KMER ;
	else d = bk.depth ;
    //if(bk.depth > DEPTH) d = bk.depth;
    //else d = DEPTH;
    // set initialGraph
    initialGraph[bk.kmer_graph_index].direction = bk.direction ;
    initialGraph[bk.kmer_graph_index].stackOffset = kmerStack->count;
    initialGraph[bk.kmer_graph_index].origin_graph_index = bk.origin_graph_index;
    initialGraph[bk.kmer_graph_index].depth = bk.depth;
    
	for(i = 0; i < d; i++)
    {
        int flag = 0 , n ;
        uint64_t hashAdd ;
        uint8_t base ;
        uint64_t rkmer[len];
        int vFlag = 0 ; // record the count of link number 
        getRevKmer(k.kmer, rkmer, KMER);
        //getRevKmer(rkmer, tmp, KMER);
        for(int j = 0 ; j < len ; j++ )
        {
            if(k.kmer[j] > rkmer[j]) { flag = 1 ; break ; }
            else if(k.kmer[j] < rkmer[j]) break;
        }
        //n = (k.kmer[0] ^ rkmer[0]) % hashNumber ;
        if(flag == 1)
        {
            n = get_time_hash(rkmer, len ) % hashNumber;
            hashAdd = getHashAddr(&hashTable[n], rkmer, len);    
        }else {
            n = get_time_hash(k.kmer, len) % hashNumber;
            hashAdd = getHashAddr(&hashTable[n], k.kmer, len);
        }
        // set processed flag
        {
            Ku ku;
            ku.ki = hashTable[n].table[hashAdd].ku.ki;
            if(ku.ki.deleted == 1) {initialGraph[bk.kmer_graph_index].state = 4; break; }
            if(s_flag == 0)
            {
				uint32_t oldv, newv, curv;
                ku.ki.processed = 0; oldv = ku.hi;
                ku.ki.processed = 1; newv = ku.hi;
                if((curv = __sync_val_compare_and_swap(&(hashTable[n].table[hashAdd].ku.hi), oldv, newv)) != oldv)
                {
                    initialGraph[bk.kmer_graph_index].state = 2 ;   
                }
            } else if(s_flag == 1) {
                if(ku.ki.processed == 1) { initialGraph[bk.kmer_graph_index].state = 2; break; }
                else if(ku.ki.curprocess == 1) {initialGraph[bk.kmer_graph_index].state = 1; }
				else { hashTable[n].table[hashAdd].ku.ki.curprocess = 1 ; }
            } else { fprintf(stderr, "[getGraphEdge] s_flag set error, exit...\n"); exit(1); }
        }
        //memcpy(&(kmerStack->kmerAddr[kmerStack->count].k), &(hashTable[n].table[hashAdd]), sizeof(kmer_t));
        // realloc kmerAddr
        if(kmerStack->count >= kmerStack->size)
        {
            kmerStack->size <<= 1 ;
            kmerStack->kmerAddr = (KmerAddr*)xrecalloc(kmerStack->kmerAddr, kmerStack->size * sizeof(KmerAddr));
        }
        kmerStack->kmerAddr[kmerStack->count].k = hashTable[n].table[hashAdd];
        kmerStack->kmerAddr[kmerStack->count].hID = n ;
        kmerStack->kmerAddr[kmerStack->count].addr = hashAdd ;
        if(bk.direction == FORWARD)
        {
            if(flag == 1)
            {
                for(uint8_t j = 0; j <4; j++)
                {
                    if(((kmerStack->kmerAddr[kmerStack->count].k.ku.ki.lLink>>j) & 0x1) == 1)
                    {   base = (~j) & 0x3 ; vFlag++; }
                }
            } else {
                for(uint8_t j = 0 ; j < 4; j++)
                {
                    if(((kmerStack->kmerAddr[kmerStack->count].k.ku.ki.rLink>>j) & 0x1) == 1)
                    {   base = j ; vFlag++; }
                }
            }
            memcpy(rkmer, k.kmer, sizeof(uint64_t) * len);
            getNextKmer(k.kmer ,rkmer, base , KMER);
        } else {
            if(flag == 1)
            {
                for(uint8_t j = 0; j < 4; j++)
                {
                    if(((kmerStack->kmerAddr[kmerStack->count].k.ku.ki.rLink>>j) & 0x1) == 1)
                    {   base = (~j) & 0x3; vFlag++;   }
                }
            } else {
                for(uint8_t j = 0; j < 4; j++)
                {
                    if(((kmerStack->kmerAddr[kmerStack->count].k.ku.ki.lLink>>j) & 0x1) == 1)
                    {   base = j ; vFlag++;   }
                }
            }
            memcpy(rkmer, k.kmer, sizeof(uint64_t) * len);
            getPreviousKmer(k.kmer, rkmer, base, KMER);
        }
        kmerStack->kmerAddr[kmerStack->count].k.ku.ki.direction = flag ;
        kmerStack->count++ ;
        //  check return value
        if( (s_flag == 0 && initialGraph[bk.kmer_graph_index].state == 2) || 
            (s_flag == 1 && initialGraph[bk.kmer_graph_index].state == 1)) break;
        else if(vFlag == 0) { initialGraph[bk.kmer_graph_index].state = 4 ; break; }
        else if( vFlag > 1) { initialGraph[bk.kmer_graph_index].state = 3 ; break; } 
        else { // vFlag == 1
            int c = 0;
            for(uint8_t j = 0; j < 4; j++)
            {
                if(((kmerStack->kmerAddr[kmerStack->count -1].k.ku.ki.rLink>>j) & 0x1) == 1) c++;
                if(((kmerStack->kmerAddr[kmerStack->count -1].k.ku.ki.lLink>>j) & 0x1) == 1) c++;
            }
            if(c > 2) {  initialGraph[bk.kmer_graph_index].state = 3 ; break; }
            else if(c == 1 ) { initialGraph[bk.kmer_graph_index].state = 4; break; }
        }
    }
    initialGraph[bk.kmer_graph_index].depth -= i ;
    initialGraph[bk.kmer_graph_index].count = kmerStack->count - initialGraph[bk.kmer_graph_index].stackOffset ;
	// calculate kmer freq
	{
		int total_freq = 0 ;
		for(int i = initialGraph[bk.kmer_graph_index].stackOffset; i < initialGraph[bk.kmer_graph_index].stackOffset + initialGraph[bk.kmer_graph_index].count ; i++)
		{
			total_freq += kmerStack->kmerAddr[i].k.ku.ki.freq ;
		}
		if(initialGraph[bk.kmer_graph_index].count > 0)
		{
			initialGraph[bk.kmer_graph_index].freq = total_freq / initialGraph[bk.kmer_graph_index].count ;
		}
	}
    return 0  ;
}


void simplifyGraph(HashTable *hashTable, const int hashNumber, const int id, const KmerFreqCurve curve)
{
	int tip_rm = 0 ;
    for(int64_t i = hashTable[id].size -1 ; i >= 0; i--)
    {
        kmer_t k = hashTable[id].table[i];
        kmer_info ki = k.ku.ki;
        if(ki.deleted == 0 && ki.processed == 0 && ki.freq > curve.valley && ki.freq < 2 * curve.peak && 
            (LINK_BIGGER_THAN_ONE(ki.rLink) || LINK_BIGGER_THAN_ONE(ki.lLink)))
        {
			hashTable[id].table[i].ku.ki.curprocess = 1 ;
			int graphNum = 10 , gCount = 0 ;
			BranchQueue branchQueue ;
			KmerStack kmerStack ;
			kmer_graph *initialGraph = initialGraphStruct(graphNum, &kmerStack, &branchQueue) ;
			memcpy(&kmerStack.kmerAddr[kmerStack.count].k, &k, sizeof(kmer_t));
			kmerStack.kmerAddr[kmerStack.count].hID = id; kmerStack.kmerAddr[kmerStack.count].addr = i ;
			kmerStack.kmerAddr[kmerStack.count].k.ku.ki.direction = 0; kmerStack.kmerAddr[kmerStack.count].k.ku.ki.processed = 1 ;
			initialGraph[gCount].size = KMER ;
			initialGraph[gCount].stackOffset = kmerStack.count;
			initialGraph[gCount].origin_graph_index = -1 ;
			initialGraph[gCount].count = 1 ;
			initialGraph[gCount].freq = kmerStack.kmerAddr[0].k.ku.ki.freq ;
			initialGraph[gCount].depth = KMER * 4 ; 
			kmerStack.count++;
			gCount++;
			/* realloc memory if overflow */
			if(gCount + 8 >= graphNum)
			{
				graphNum <<= 1 ;
				initialGraph = (kmer_graph*)xrecalloc(initialGraph, graphNum * sizeof(kmer_graph));
			}
			
			//	int depth = DEPTH * 4 ;
			for(uint8_t z = 0; z < 4; z++)
			{
				if((ki.rLink & (1<<z)) > 0) { addToQueue(initialGraph, 0, gCount, &branchQueue, z, FORWARD, k.kmer); gCount++; }
				if((ki.lLink & (1<<z)) > 0) { addToQueue(initialGraph, 0, gCount, &branchQueue, z, BACKWARD, k.kmer); gCount++; }
			}
			while(branchQueue.start < branchQueue.end) // branchQueue is not NULL
			{
				int index = branchQueue.branchKmers[branchQueue.start].kmer_graph_index;
				getGraphEdge(&branchQueue, initialGraph, &kmerStack, hashTable, hashNumber, 1);
				// add new branch edges to the queue
				{
					if(initialGraph[index].depth > 0 && initialGraph[index].state == 3)
					{
						kmer_t ck, pk;
						if(initialGraph[index].count == 1)
						{
							int pi = initialGraph[index].origin_graph_index;
							pk = kmerStack.kmerAddr[initialGraph[pi].stackOffset + initialGraph[pi].count -1].k;
						} else if(initialGraph[index].count > 1) {
							pk = kmerStack.kmerAddr[initialGraph[index].stackOffset + initialGraph[index].count -2].k;
						} else {fprintf(stderr, "[getGraphEdge] initialGraph[index].count error, exit..\n"); exit(1);}
						ck = kmerStack.kmerAddr[initialGraph[index].stackOffset + initialGraph[index].count -1].k ;
						/* realloc memory if overflow */
						if(gCount + 8 >= graphNum)
						{
							graphNum <<= 1 ;
							initialGraph = (kmer_graph*)xrecalloc(initialGraph, graphNum * sizeof(kmer_graph));
						}
						gCount = addBranch(initialGraph, index, gCount, &branchQueue, BIDIRECTION, ck, pk);
					}
				}
			}
			// simplify kmer graph 
			tip_rm += ReduceGraph(initialGraph, &kmerStack, curve);

			// write kmer  info to the Hash Table
			writeKmerToHashTable(&kmerStack, hashTable);

			// free and clean work
			free(initialGraph); free(kmerStack.kmerAddr); free(branchQueue.branchKmers);
		}
	}
	fprintf(stderr, "removing tips: %d\n", tip_rm);
}

int CheckBranch(const kmer_t k, const int direction)
{
	kmer_info ki = k.ku.ki ;
	int b = 0;
	// other kmer branch must equal to 1
	if((ki.direction == 1 && direction == FORWARD) ||
	 (ki.direction == 0 && direction == BACKWARD))
	{
		for(uint8_t z = 0; z < 4; z++)
		{
			if((ki.lLink & (1<<z)) > 0) b++;
		}
	} else {
		for(uint8_t z = 0; z < 4; z++)
		{
			if((ki.rLink & (1<<z)) > 0) b++;
		}
	}
	return b ;
}

int CheckAndMergeSeq(const char *s1, const char *s2, char *s)
{
	int p1 = strlen(s1), p2 = strlen(s2);
	if(p1 >= KMER && p2 >= KMER) return 0 ;
	for(--p1, --p2; p1 >= 0 && p2 >= 0; p1--, p2--)
	{
		if(s1[p1] != s2[p2]) break;
	}
	
	if(p1 <= 6 && p2 <= 6) 
	{
		strcat(s, "{"); strncat(s, s1, p1+1); strcat(s, ","); strncat(s, s2, p2+1); strcat(s, "}");
		strcat(s, s1+p1+1);
		return 1 ;
	} else { return 0; }

}

char *getEdgeSeq(const kmer_graph *kg, const KmerStack *ks)
{
	char *s = xcalloc(kg->count + 2, sizeof(char));
	int c = 0;
	
	for(int i = kg->stackOffset ; i < kg->stackOffset + kg->count - 1; i++)
	{
		kmer_info ki = ks->kmerAddr[i].k.ku.ki ;
		int n = 0;
		if(kg->direction == FORWARD)
		{
			for(uint8_t z = 0; z < 4; z++)
			{
				if(ki.direction == 1)
				{
					if((ki.lLink & (1<<z)) > 0) 
					{
						uint8_t base = (~z) & 0x3 ;
						s[c++] = BIT_NT_CHAR[base] ;
						n++ ;
					}
				} else { // ki.direction == 0
					if((ki.rLink & (1<<z)) > 0)
					{
						s[c++] = BIT_NT_CHAR[z];
						n++;
					}
				}
			}
		} else { // BACKWARD
			for(uint8_t z = 0; z < 4; z++)
			{
				if(ki.direction == 1)
				{
					if((ki.rLink & (1<<z)) > 0)
					{
						uint8_t base = (~z) & 0x3 ;
						s[c++] = BIT_NT_CHAR[base];
						n++;
					} 
				} else { // ki.direction == 0
					if((ki.lLink & (1<<z)) > 0)
					{
						s[c++] = BIT_NT_CHAR[z];
						n++ ;
					}
				}
			}	
		}
	}

	return s ;
}

void extractUnitig(const Arguments *arguments, HashTable *hashTable, const int hashNumber, const KmerFreqCurve curve)
{
	char unitig_name[PATH_LEN];
	gzFile unitigF;
	int unitigID = 1 ;
	strcpy(unitig_name, arguments->prefix); strcat(unitig_name, ".unitig.gz");
	unitigF = xzopen(unitig_name, "w");
	for(int i = 0; i < hashNumber; i++)
	{
		for(int64_t j = 0; j < hashTable[i].size; j++)
		{
			kmer_t k = hashTable[i].table[j];
			kmer_info ki = k.ku.ki;
			if(ki.processed == 0 && ki.deleted == 0 && abs(ki.freq - curve.peak) < curve.peak * PEAK_SD && 
				LINK_EQUL_ONE(ki.rLink) && LINK_EQUL_ONE(ki.lLink))
			{
				int graphNum = 10, gCount = 0;
				BranchQueue branchQueue;
				KmerStack kmerStack;
				kmer_graph *initialGraph = initialGraphStruct(graphNum, &kmerStack, &branchQueue);
				hashTable[i].table[j].ku.ki.curprocess = 1; k.ku.ki.curprocess = 1 ; 
				memcpy(&kmerStack.kmerAddr[kmerStack.count].k, &k, sizeof(kmer_t));
				kmerStack.kmerAddr[kmerStack.count].k.ku.ki.direction = 0;
				kmerStack.kmerAddr[kmerStack.count].hID = i; kmerStack.kmerAddr[kmerStack.count].addr = j;
				initialGraph[gCount].size = KMER ;
				//initialGraph[0].edge_contig = (uint64_t*)xcalloc((initialGraph[0].size + 32 -1)/32, sizeof(uint64_t));
				int r_size = INITIAL_STACK_SIZE, l_size = INITIAL_STACK_SIZE, r_len = 0, l_len = 0;
				char *redge_contig = (char*)xcalloc(r_size, sizeof(char));
				char *ledge_contig = (char*)xcalloc(l_size, sizeof(char));
				initialGraph[gCount].stackOffset = kmerStack.count;
				//memcpy(initialGraph[0].edge_contig, kmerStack.kmerAddr[0].k.kmer, (KMER + 32 -1)/32);
				initialGraph[gCount].count = 1 ;
				initialGraph[gCount].origin_graph_index = -1 ;
				initialGraph[gCount].depth = INT_LEAST32_MAX - 1 ;
				kmerStack.count++;
				gCount++ ;

				//	extend right unitig
				{
					for(uint8_t z = 0; z < 4; z++)
					{
						if((ki.rLink & (1<<z)) > 0) { addToQueue(initialGraph, 0, gCount, &branchQueue, z, FORWARD, k.kmer); gCount++; }
						//if((ki.lLink & (1<<z)) > 0) { addToQueue(0, z, BACKWARD, INT_LEAST32_MAX, k.kmer); }
					}
					while(branchQueue.start < branchQueue.end) // branchQueue is not NULL
					{
						if(branchQueue.end - branchQueue.start == 1)
						{
							int gI = branchQueue.branchKmers[branchQueue.start].kmer_graph_index ;
							getGraphEdge(&branchQueue, initialGraph, &kmerStack, hashTable, hashNumber, 1);
							char *s = getEdgeSeq(&initialGraph[gI],&kmerStack);
							if(r_size <= r_len + strlen(s))
							{
								r_size = (r_len + strlen(s)) * 2 ;
								redge_contig = realloc(redge_contig, r_size);
								memset(redge_contig + r_len, 0, sizeof(char) * (r_size - r_len));	
							}
							strcat(redge_contig, s); r_len += strlen(s);
							free(s);
							// check right branchs and add new branch edges to the queue
							kmer_t ck = kmerStack.kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].k ;
							if(initialGraph[gI].state == 3 && CheckBranch(ck, BACKWARD) == 1 && CheckBranch(ck, FORWARD) == 2)
							{
								kmer_t pk;
								if(initialGraph[gI].count == 1)
								{
									int pi = initialGraph[gI].origin_graph_index;
									pk = kmerStack.kmerAddr[initialGraph[pi].stackOffset + initialGraph[pi].count -1].k;
								} else if(initialGraph[gI].count > 1) {
									pk = kmerStack.kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -2].k;
								} else {fprintf(stderr, "[getGraphEdge] initialGraph[index].count error, exit..\n"); exit(1);}
								
								/* realloc memory if overflow */
								if(gCount + 2 >= graphNum)
								{
									graphNum <<= 1 ;
									initialGraph = (kmer_graph*)xrecalloc(initialGraph, graphNum * sizeof(kmer_graph));
								}
								gCount = addBranch(initialGraph, gI, gCount, &branchQueue, FORWARD, ck, pk); 
							}
						} else if(branchQueue.end - branchQueue.start  == 2) {
							int gI1 = branchQueue.branchKmers[branchQueue.start].kmer_graph_index, gI2, tmp = -1 ;	
							getGraphEdge(&branchQueue, initialGraph, &kmerStack, hashTable, hashNumber, 1);
							gI2 = branchQueue.branchKmers[branchQueue.start].kmer_graph_index ;
							getGraphEdge(&branchQueue, initialGraph, &kmerStack, hashTable, hashNumber, 1);
							
							if(initialGraph[gI2].state == 3) { tmp = gI1; gI1 = gI2; gI2 = tmp; }
							kmer_t tk = kmerStack.kmerAddr[initialGraph[gI1].stackOffset + initialGraph[gI1].count -1].k ;
							if(initialGraph[gI1].state == 3 && initialGraph[gI2].state == 1 &&
								memcmp(&tk, &kmerStack.kmerAddr[initialGraph[gI2].stackOffset + initialGraph[gI2].count -1].k,
								sizeof(kmer_t)) == 0 && CheckBranch(tk, BACKWARD) == 2 && CheckBranch(tk, FORWARD == 1))
							{
								char *s1 = getEdgeSeq(&initialGraph[gI1], &kmerStack);
								char *s2 = getEdgeSeq(&initialGraph[gI2], &kmerStack);
								char *s = xcalloc(strlen(s1) + strlen(s2) + 5, sizeof(char));
								if(CheckAndMergeSeq(s1, s2, s))
								{
									if(r_size <= r_len + strlen(s))
									{
										r_size = (r_len + strlen(s)) * 2 ;
										redge_contig = realloc(redge_contig, r_size);
										memset(redge_contig + r_len, 0, sizeof(char) * (r_size - r_len));	
									}
									strcat(redge_contig, s); r_len += strlen(s);
									{
										kmer_t pk;
										if(initialGraph[gI1].count == 1)
										{
											int pi = initialGraph[gI1].origin_graph_index ;
											pk = kmerStack.kmerAddr[initialGraph[pi].stackOffset + initialGraph[pi].count -1].k ;
										} else if(initialGraph[gI1].count > 1) {
											pk = kmerStack.kmerAddr[initialGraph[gI1].stackOffset + initialGraph[gI1].count -2].k ;
										} else { fprintf(stderr, "[extractUnitig] initialGraph[gI1].count error, exit...\n"); exit(1); }
										/* realloc memory if overflow */
										if(gCount + 8 >= graphNum)
										{
											graphNum <<= 1 ;
											initialGraph = (kmer_graph*)xrecalloc(initialGraph, graphNum * sizeof(kmer_graph));
										}
										gCount = addBranch(initialGraph, gI1, gCount, &branchQueue, FORWARD, tk, pk);
									}
								}

								free(s1); free(s2); free(s);
							}

						} else if(branchQueue.end - branchQueue.start > 2) {
							// not defined 
							break;
						} else {
							fprintf(stderr, "[extractUnitig] branch number is error, ....exit!\n");
							exit(1);
						}
					}
				}


				//	extend left unitig
				{
					for(uint8_t z = 0; z < 4; z++)
					{
						if((ki.lLink & (1<<z)) > 0) { addToQueue(initialGraph, 0, gCount, &branchQueue, z, BACKWARD, k.kmer); gCount++; }
					}
					while(branchQueue.start < branchQueue.end) // branchQueue is not NULL
					{
						if(branchQueue.end - branchQueue.start == 1)
						{
							int gI = branchQueue.branchKmers[branchQueue.start].kmer_graph_index ;
							getGraphEdge(&branchQueue, initialGraph, &kmerStack, hashTable, hashNumber, 1);
							// branchQueue.start++;
							char *s = getEdgeSeq(&initialGraph[gI],&kmerStack);
							if(l_size <= l_len + strlen(s))
							{
								l_size = (l_len + strlen(s)) * 2 ;
								ledge_contig = realloc(ledge_contig, l_size);
								memset(ledge_contig + l_len, 0, sizeof(char) * (l_size - l_len));	
							}
							strcat(ledge_contig, s); l_len += strlen(s);
							free(s);
							// check right branchs and add new branch edges to the queue
							kmer_t ck = kmerStack.kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -1].k ;
							if(initialGraph[gI].state == 3 && CheckBranch(ck, BACKWARD) == 2 && CheckBranch(ck, FORWARD) == 1)
							{
								kmer_t pk;
								if(initialGraph[gI].count == 1)
								{
									int pi = initialGraph[gI].origin_graph_index;
									pk = kmerStack.kmerAddr[initialGraph[pi].stackOffset + initialGraph[pi].count -1].k;
								} else if(initialGraph[gI].count > 1) {
									pk = kmerStack.kmerAddr[initialGraph[gI].stackOffset + initialGraph[gI].count -2].k;
								} else {fprintf(stderr, "[getGraphEdge] initialGraph[index].count error, exit..\n"); exit(1);}
								
								/* realloc memory if overflow */
								if(gCount + 2 >= graphNum)
								{
									graphNum <<= 1 ;
									initialGraph = (kmer_graph*)xrecalloc(initialGraph, graphNum * sizeof(kmer_graph));
								}
								gCount = addBranch(initialGraph, gI, gCount, &branchQueue, BACKWARD, ck, pk); 
							}
						} else if(branchQueue.end - branchQueue.start  == 2) {
							int gI1 = branchQueue.branchKmers[branchQueue.start].kmer_graph_index, gI2, tmp = -1 ;	
							getGraphEdge(&branchQueue, initialGraph, &kmerStack, hashTable, hashNumber, 1);
							branchQueue.start++;
							gI2 = branchQueue.branchKmers[branchQueue.start].kmer_graph_index ;
							getGraphEdge(&branchQueue, initialGraph, &kmerStack, hashTable, hashNumber, 1);
							branchQueue.start++;
							if(initialGraph[gI2].state == 3) { tmp = gI1; gI1 = gI2; gI2 = tmp; }
							kmer_t tk = kmerStack.kmerAddr[initialGraph[gI1].stackOffset + initialGraph[gI1].count -1].k ;
							if(initialGraph[gI1].state == 3 && initialGraph[gI2].state == 1 &&
								memcmp(&tk, &kmerStack.kmerAddr[initialGraph[gI2].stackOffset + initialGraph[gI2].count -1],
								sizeof(kmer_t)) == 0 && CheckBranch(tk, BACKWARD) == 2 && CheckBranch(tk, FORWARD == 1))
							{
								char *s1 = getEdgeSeq(&initialGraph[gI1], &kmerStack);
								char *s2 = getEdgeSeq(&initialGraph[gI2], &kmerStack);
								char *s = xcalloc(strlen(s1) + strlen(s2) + 5, sizeof(char));
								if(CheckAndMergeSeq(s1, s2, s))
								{
									if(l_size <= l_len + strlen(s))
									{
										l_size = (l_len + strlen(s)) * 2 ;
										ledge_contig = realloc(ledge_contig, l_size);
										memset(ledge_contig + l_len, 0, sizeof(char) * (l_size - l_len));	
									}
									strcat(redge_contig, s); l_len += strlen(s);
									{
										kmer_t pk;
										if(initialGraph[gI1].count == 1)
										{
											int pi = initialGraph[gI1].origin_graph_index ;
											pk = kmerStack.kmerAddr[initialGraph[pi].stackOffset + initialGraph[pi].count -1].k ;
										} else if(initialGraph[gI1].count > 1) {
											pk = kmerStack.kmerAddr[initialGraph[gI1].stackOffset + initialGraph[gI1].count -2].k ;
										} else { fprintf(stderr, "[extractUnitig] initialGraph[gI1].count error, exit...\n"); exit(1); }
										/* realloc memory if overflow */
										if(gCount + 8 >= graphNum)
										{
											graphNum <<= 1 ;
											initialGraph = (kmer_graph*)xrecalloc(initialGraph, graphNum * sizeof(kmer_graph));
										}
										gCount = addBranch(initialGraph, gI1, gCount, &branchQueue, BACKWARD, tk, pk);
									}
								}

								free(s1); free(s2); free(s);
							}

						} else if(branchQueue.end - branchQueue.start > 2) {
							// not defined 
							break;
						} else {
							fprintf(stderr, "[extractUnitig] branch number is error, ....exit!\n");
							exit(1);
						}
					}
				}

				// write to unitig file
				{
					char *unitig = xcalloc(strlen(redge_contig) + strlen(ledge_contig) + KMER + 5, sizeof(char)) ;
					char tmp[PATH_LEN];
					memset(tmp, 0, PATH_LEN * sizeof(char));
					for(int i = strlen(ledge_contig)-1, j = 0; i >= 0; i--, j++)
					{
						if(ledge_contig[i] == '{') unitig[j] = '}';
						else if(ledge_contig[i] == '}') unitig[j] = '{' ;
						else unitig[j] = ledge_contig[i];
					}

					char *s = xcalloc(KMER + 2, sizeof(char));
					bit64ToChar(s, kmerStack.kmerAddr[0].k.kmer, KMER);
					strcat(unitig, s); strcat(unitig, redge_contig);
				
					sprintf(tmp, ">%d\t%lu\n", unitigID, strlen(unitig));	
					gzwrite(unitigF, tmp, strlen(tmp));
					strcat(unitig, "\n");
					gzwrite(unitigF, unitig, strlen(unitig));

					unitigID++;

					free(unitig);
				   	free(s);
				}

				// write kmer info to the Hash Table
				writeKmerToHashTable(&kmerStack, hashTable);

				// free and clean work
				free(redge_contig); free(ledge_contig);
				free(branchQueue.branchKmers); 
				free(kmerStack.kmerAddr); 
				free(initialGraph);
			}
		}
	}
	// free and clean work
	gzclose(unitigF);
}

// return the number of kmer_t added 
int AddToKmerStack(KmerStack *kmerStack, const kmer_t k)
{
	kmer_info ki = k.ku.ki ;

	int count = kmerStack->count ;

	// realloc kmerStack
	if(kmerStack->count + 8 > kmerStack->size)
	{
		kmerStack->size <<= 1 ;
		kmerStack->kmerAddr = (KmerAddr*)xrecalloc(kmerStack->kmerAddr, kmerStack->size * sizeof(KmerAddr));
	}
	if(ki.direction == 0)
	{
		for(uint8_t z = 0; z < 4 ; z++)
		{
			if((ki.rLink & (1<<z)) > 0) 
			{ 
				getNextKmer(kmerStack->kmerAddr[kmerStack->count].k.kmer, k.kmer, z, KMER) ;
				kmerStack->kmerAddr[kmerStack->count].addr = FORWARD ;
				kmerStack->count++;	
			} 
			if((ki.lLink & (1<<z)) > 0)
			{
				getPreviousKmer(kmerStack->kmerAddr[kmerStack->count].k.kmer, k.kmer, z, KMER);
				kmerStack->kmerAddr[kmerStack->count].addr = BACKWARD ;
				kmerStack->count++ ;
			}
		}
	} else { // ki.direction == 1
		kmer_t ck ;
		getRevKmer(k.kmer, ck.kmer, KMER);
		for(uint8_t z = 0; z < 4; z++)
		{
			if((ki.rLink & (1<<z)) > 0)
			{
				uint8_t base = (~z) & 0x3 ;
				getPreviousKmer(kmerStack->kmerAdrr[kmerStack->count].k.kmer, ck.kmer, base, KMER);
				kmerStack->kmerAddr[kmerStack->count].addr = BACKWARD ;
				kmerStack->count++;
			} 
			if((ki.lLink & (1<<z)) > 0)
			{
				uint8_t base = (~z) & 0x3 ;
				getNextKmer(kmerStack->kmerAddr[kmerStack->count].k.kmer, ck.kmer, base, KMER);
				kmerStack->kmerAddr[kmerStack->count].addr = FORWARD ;
				kmerStack->count++ ;
			}
		}
	}

	return kmerStack->count - count ;
	
}

kmer_t setEdgeDirectionFlag(HashTable *hashTable, const int hashNumber, const KmerAddr ka)
{
	kmer_t k = ka.k, ret_k ;
	memset(&ret_k, 0, sizeof(kmer_t));
	int len = (KMER + 32 -1)/32 ;
	while(1)
	{
		uint64_t hashAdd ;
		uint64_t rkmer[len] ;
		int flag = 0, n = -1 ;
		getRevKmer(k.kmer, rkmer, KMER);
		for(int i = 0; i < len; i++)
		{
			if(k.kmer[i] > rkmer[i]) { flag = 1 ; break; }
			else if(k.kmer[i] < rkmer[i]) break;
		}

		if(flag == 1)
		{
			n = get_time_hash(rkmer, len) % hashNumber ;
			hashAdd = getHashAddr(&hashTable[n], rkmer, len);
		} else {
			n = get_time_hash(k.kmer, len) % hashNumber ;
			hashAdd = getHashAddr(&hashTable[n], k.kmer, len);
		}
		kmer_info ki = hashTable[n].table[hashAdd].ku.ki ;
		if(ki.deleted == 1 || ki.processed == 1) { break; }
		else {
			// set kmer_t direction flag
			hashTable[n].table[hashAdd].ku.ki.direction = flag ;
			hashTable[n].table[hashAdd].ku.ki.processed = 1 ;
			if(LINK_BIGGER_THAN_ONE(ki.rLink) || LINK_BIGGER_THAN_THAN_ONE(ki.lLink)) 
			{
				ret_k = hashTable[n].table[hashAdd] ;
				break;
			} else { 
				if(ka.addr == FORWARD)
				{
					if(flag == 1)
					{
						if(ki.lLink == 0) { break; }
						else {
							for(uint8_t z = 0; z < 4; z++)
							{
								if((ki.lLink & (1<<z)) > 0)
								{
									uint8_t base = (~z) & 0x3 ;
									getNextKmer(k.kmer, rkmer, base, KMER);
								}
							}
						}
					} else {
						if(ki.rLink == 0) { break; }
						else {
							for(uint8_t z = 0; z < 4; z++)
							{
								if((ki.rLink) & (1<<z) > 0)
								{
									kmer_t tk = k ;
									getNextKmer(k.kmer, tk.kmer, z, KMER);
								}
							}
						}
					}
				} else { // ka.addr == BACKWARD
					if(flag == 1)
					{
						if(ki.rLink == 0) { break; }
						else {
							for(uint8_t z = 0; z < 4; z++)
							{
								if((ki.rLink & (1<<z)) > 0)
								{
									uint8_t base = (~z) & 0x3 ;
									getPreviousKmer(k.kmer, rkmer, base, KMER);
								}
							}
						}
					} else {
						if(ki.lLink == 0) { break; }
						else {
							for(uint8_t z = 0; z < 4; z++)
							{
								if((ki.lLink) & (1<<z) > 0)
								{
									kmer_t tk = k ;
									getPreviousKmer(k.kmer, tk.kmer, z, KMER);
								}
							}
						}
					}
				}
			}	
		}
	}
	
	return ret_k ;	

}


	// set kmer_info direction flag		
void SetDirectionFlag(HashTable hashTable, const int hashNumber)
{
	for(int i = 0; i < hashNumber; i++)
	{
		for(int64_t j = 0; j < hashTablep[i].size; j++)
		{
			kmer_t k = hashTable[i].table[j] ;
			kmer_info ki = k.ku.ki ;
			if( ki.processed == 0 && ki.deleted == 0 && 
				(LINK_BIGGER_THAN_ONE(ki.rLink) || LINK_BIGGER_THAN_ONE(ki.lLink) || ki.rLink == 0 || ki.lLink == 0))
			{
				hashTable[i].table[j].ku.ki.process = 1 ;
				KmerStack kmerStack ; // kmerStack.kmerAddr[0].addr as direction of kmer graph edge 
				kmerStack.size = 10 ; kmerStack.count = 0 ;
				kmerStack.kmerAddr = (KmerAddr*)xcalloc(kmerStack.size, sizeof(KmerAddr));
				AddToKmerStack(&kmerStack, k) ;

				while(kmerStack.count > 0)
				{
					kmerStack.count-- ;
					KmerAddr ka = kmerStack.kmerAddr[kmerStack.count] ;
					kmer_t k = setEdgeDirectionFlag(hashTable, hashNumber, ka);
					if(k.ku.ki.rLink != 0 || k.ku.ki.Llink != 0)
					{
						AddToKmerStack(&kmerStack, k);
					}
				}
				
				// free work 
				free(kmerStack.kmerAddr);	
			}
		}
	}
}

int compare_kmer(const void *p, const void *q)
{
	kmer_t *p1 = (kmer_t*)p, *q1 = (kmer_t*)q ;
	int len = (KMER + 32 -1) / 32 ;
	int retv = 0;
	for(int i = 0; i < len; i++)
	{
		if(p1->kmer[i] < q1->kmer[i]) { retv = -1; break;}
		else if(p1->kmer[i] > q1->kmer[i]) { retv = 1 ; break; }
	}
	
	return retv ;
}

int compare_KmerToVer(const void *p, const void *q)
{
	KmerToVer *p1 = (KmerToVer*)p, *q1 = (KmerToVer*)q ;
	int len = (KMER + 32 -1)/32 ;
	int retv = 0;
	for(int i = 0; i < len; i++)
	{
		if(p1->k.kmer[i] < q1->k.kmer[i]) { retv = -1; break; }
		else if(p1->k.kmer[i] > q1->k.kmer[i]) { retv = 1; break; }
	}

	return retv ;
}

inline kmer_t GetKmer(const kmer_t tk, const HashTable *hashTable, const int hashNumber)
{
	uint64_t hashAdd ;
	int len = (KMER + 32 -1)/32 ;
	uint64_t rkmer[len];
	int flag = 0, n = -1 ;
	getRevKmer(tk.kmer, rkmer, KMER);
	for(int i = 0; i < len; i++)
	{
		if(tk.kmer[i] > rkmer[i]) { flag = 1; break; }
		else if(tk.kmer[i] < rkmer[i]) { break; }
	}

	if(flag == 1)
	{
		n = get_time_hash(rkmer, len) % hashNumber ;
		hashAdd = getHashAddr(&hashTable[n], rkmer, len);
	} else {
		n = get_time_hash(tk.kmer, len) % hashNumber ;
		hashAdd = getHashAddr(&hashTable[n], tk.kmer, len);
	}
	
	return hashTable[n].table[hashAdd];	
}

int IsProcessed(const kmer_t tk, const  HashTable *hashTable, const int hashNumber)
{
	kmer_t ck = GetKmer(tk, hashTable, hashNumber);
	if(ck.processed == 1) return 1 ;
	else return 0 ;
}
inline char *GetKmerArrSeq(const KmerArr *ka)
{
	char *s = (char*)xcalloc(ka->count + KMER , sizeof(char));
	kmer_t tk ;
	if(ka->arr[0].ku.ki.direction == 0) { tk = ka->arr[0]; }
	else { getRevKmer(ka->arr[0].kmer, tk.kmer, KMER); }
	char *ks = transformKmerToString(tk.kmer, 0, KMER);
	strcpy(s, ks);
	for(int i = 0; i < ka->count-1; i++)
	{
		kmer_t ck = ka->arr[i];
		int degree = 0;
		if(ck.ku.ki.direction == 0)
		{
			for(uint8_t z = 0; z < 4; z++)
			{
				if((ck.ku.ki.rLink & (1<<z)) > 0)
				{
					degree++;
					char c = BIT_NT_CHAR[z];
					strncat(s, &c, 1);
				}
			}	
		} else {
			for(uint8_t z = 0; z < 4; z++)
			{
				if((ck.ku.ki.lLink & (1<<z)) > 0)
				{
					degree++ ;
					uint32_t base = (~z) & 0x3;
					char c = BIT_NT_CHAR[base];
					strncat(s, &c, 1);
				}
			}
		}
		if(degree != 1) { fprintf(stderr, "[GetKmerArrSeq] kmer degree set error, exit.....\n"); exit(1); }
		
	}

	return s ;	
}

inline uint8_t *GetKmerArrFreq(const KmerArr *ka)
{
	uint8_t *f = (uint8_t*)xcalloc(ka->count, sizeof(uint8_t));
	
	for(int i = 0; i < ka->count; i++)
	{
		f[i] = ka->arr[i].ku.ki.freq ;
	}
		
	return f ;
}

inline int GetVertexDegree(const kmer_t k, const int flag)
{
	uint32_t link ;
	int degree = 0;
	kmer_info ki = k.ku.ki;
	if(ki.direction == 0)
	{
		if(flag == IN) { link = ki.lLink; }
		else { link = ki.rLink; }
	} else {
		if(flag == IN) { link = ki.rLink; }
		else { link = ki.lLink ; }
	}
	for(int i = 0; i < 4 ; i++)
	{
		if((link & (1<<i)) > 0) { degree++; }
	}

	return degree ;
}

inline char *GetVertexString(const Vertex v)
{
	char *s = (char*)xcalloc(v.s_size * 7 + 20, sizeof(char));
	sprintf(s, ">%d\t%d\t%u\t%u\n%s\n+\n", v.ID, v.s_size, v.indegree, v.outdegree, v.seq);
	for(int i = 0; i < v.s_size - KMER + 1, i++)
	{
		sprintf(s + strlen(s), "%d ", v.freq[i]);	
	}
	sprintf(s + strlen(s), "\n");
	return s ;
}

void WriteVertexToFile(const uint32_t vID, const KmerArr *ke, gzFile vertexsgz)
{
	Vertex v ;
	if(vID > INT_LEAST32_MAX) 
	{	
		fprintf(stderr, "[saveDebruijinGraph] the number of vertices bigger than INT_LEAST32_MAX, exit.....\n");
		exit(1);
	}
	v.ID = vID ;
	v.seq = GetKmerArrSeq(ke);
	v.s_size = strlen(v.seq) ;
	v.freq = GetKmerArrFreq(ke);
	// set vertex indegree and outdegree 
	v.indegree = GetVertexDegree(ke->arr[0], IN);
	v.outdegree = GetVertexDegree(ke->arr[ke->count-1], OUT);

	char *s = GetVertexString(v);
	gzputs(vertexsgz, s);
	// final and free work
	free(s);
	free(v.freq);
	free(v.seq);
}


VerticesIndex *BulidVerticesIndexAndWriteFile(const HashTable *hashTable, const int hashNumber, gzFile vertexsgz)
{
	// construct vertices index of graph 
	uint32_t vID = 1 ;		// the minimum vertex ID
	VerticesIndex *vIdx = xcalloc(1, sizeof(VerticesIndex));
	vIdx->size = INITIAL_STACK_SIZE, vIdx->count = 0;
	vIdx->skArr = (KmerToVer*)xcalloc(vIdx->size, sizeof(KmerToVer));
	vIdx->ekArr = (KmerToVer*)xcalloc(vIdx->size, sizeof(KmerToVer));
	// get vertices kmer_t array
	for(int i = 0; i < hashNumber; i++)
	{
		for(int64_t j = 0; j < hashTable[i].size; j++)
		{
			kmer_t k = hashTable[i].table[j];
			kmer_info ki = k.ku.ki ;
			if(ki.deleted == 0 && ki.processed == 0)
			{
				int flag = 0 ;
				if(ki.direction == 0 && (ki.lLink == 0 || LINK_BIGGER_THAN_ONE(ki.lLink))) { flag = 1 ; }
				else if(ki.direction == 1 && (ki.rLink == 0 || LINK_BIGGER_THAN_ONE(ki.rLink))) { flag = 1; } 
				if(flag == 1)
				{
					KmerArr ka;
					ka.size = 10; ka.count = 0;
					ka.arr = (kmer_t*)xcalloc(ka.size, sizeof(kmer_t));
					memcpy(&ka.arr[ka.count].k, &k, sizeof(kmer_t));
					ka.count++ ;
					while(ka.count > 0)
					{
						ka.count--;
						kmer_t tk = ka.arr[ka.count];
						if(IsProcessed(tk, hashTable, hashNumber)) continue;
						KmerArr *ke = getForwardKmer(tk, hashTable, hashNumber);
						tk = ke->arr[ke->count-1] ;
						kmer_info tki = tk.ku.ki ;
						int rawCount = ke->count ;
						// realloc ka
						if(ka->count + 4 >= ka->size) 
						{
							ka->size <<= 1;
							ka->arr = (kmer_t*)xrecalloc(ka->arr, ka->size * sizeof(kmer_t));
						}
						if(vIdx->count >= vIdx->size)
						{
							vIdx->size <<= 1 ;
							vIdx->skArr = (KmerToVer*)xrecalloc(vIdx->skArr, vIdx->size * sizeof(KmerToVer));
							vIdx->ekArr = (KmerToVer*)xrecalloc(vIdx->ekArr, vIdx->size * sizeof(KmerToVer));
						}
						if(tki.direction == 1)
						{
							if(ke->count > 1 && LINK_BIGGER_THAN_ONE(tki.rLink)) 
							{
								ke->count-- ;
								tk = ke->arr[ke->count-1];
								tki = tk.ku.ki ;
							}
							// add to ka stack
							if(tki.processed == 0 && tki.lLink > 0)
							{
								kmer_t rk ;
								getRevKmer(tk.kmer, rk.kmer, KMER);
								for(uint8_t z = 0; z < 4; z++)
								{
									if((tki.lLink & (1<<z)) > 0)
									{
										uint8_t base = (~z) & 0x3 ;
										getNextKmer(tk.kmer, rk.kmer, base , KMER);
										ka->arr[ka->count] = tk ;
										ka->count++ ;
									} 
								}
							}
						}else { // tki.direction == 0
							if(ke->count > 1 && LINK_BIGGER_THAN_ONE(tki.lLink)) 
							{
								ke->count-- ;
								tk = ke->arr[ke->count-1];
								tki = tk.ku.ki ;
							}
							// add to ka stack
							if(tki.processed == 0 && tki.rLink > 0)
							{
								for(uint8_t z = 0; z < 4; z++)
								{
									if((tki.rLink & (1<<z)) > 0)
									{
										kmer_t ck ;
										getNextKmer(ck.kmer, tk.kmer, z, KMER);
										ka->arr[ka->count] = ck ;
										ka->count++ ;
									}
								}	
							}	
						}
						
						memcpy(&vIdx->skArr[vIdx->count].k, &ke->arr[0], sizeof(kmer_t));
						memcpy(&vIdx->ekArr[vIdx->count].k, &ke->arr[ke->count-1], sizeof(kmer_t));
						vIdx->skArr[vIdx->count].ID = vID ;
						vIdx->ekArr[vIdx->count].ID = vID ;
						WriteVertexToFile(vID, ke, vertexsgz);
						vID++ ;
						vIdx->count++;
						
						// final and free work
						free(ke->arr); free(ke);
					}	

					// final and free work
					free(ka.kmerAddr);
				}
			}
		}
	}

	// realloc the vIdx 
	vIdx->skArr = realloc(vIdx->skArr, vIdx->count * sizeof(KmerToVer));
	vIdx->ekArr = realloc(vIdx->ekArr, vIdx->count * sizeof(KmerToVer));
	// sort vertices array  by kmer sequence
	qsort(vIdx->skArr, vIdx->count, sizeof(KmerToVer), compare_KmerToVer);
	qsort(vIdx->ekArr, vIdx->count, sizeof(KmerToVer), compare_KmerToVer);

	return vIdx ;
}

KmerArr *getForwardKmer(const kmer_t tk, HashTable *hashTable, const int hashNumber)
{
	kmer_t k = tk ;
	KmerArr *ka = (KmerArr*)xcalloc(1, sizeof(KmerArr));
	ka->size = 10;
	ka->arr = (kmer_t*)xcalloc(ka->size, sizeof(kmer_t));
	int len = (KMER + 32 -1)/32 ;
	while(1)
	{
		uint64_t hashAdd ;
		uint64_t rkmer[len] ;
		int flag = 0, n = -1 ;
		getRevKmer(k.kmer. rkmer, KMER);
		for(int i = 0; i < len; i++)
		{
			if(k.kmer[i] > rkmer[i]) { flag = 1 ; break; }
			else if(k.kmer[i] < rkmer[i]) { break; }
		}
		
		if(flag == 1)
		{
			n = get_time_hash(rkmer, len) % hashNumber ;
			hashAdd = getHashAddr(&hashTable[n], rkmer, len);
		} else {
			n = get_time_hash(k.kmer, len) % hashNumber ;
			hashAdd = getHashAddr(&hashTable[n], k.kmer, len);
		}
		kmer_t ck = hashTable[n].table[hashAdd];
		kmer_info cki = ck.ku.ki ;
		if(ka->count >= ka->size)
		{
			ka->size <<= 1 ;
			ka->arr = (kmer_t*)xrecalloc(ka->arr, ka->size * sizeof(kmer_t));
		}
		memcpy(ka->arr[ka->count], &ck, sizeof(kmer_t));
		ka->count++ ;
		if(cki.deleted == 1 || cki.processed == 1) { break; }
		else if(LINK_BIGGER_THAN_ONE(cki.rLink) || LINK_BIGGER_THAN_ONE(cki.lLink)) {break; }
		else if((cki.direction == 0 && cki.rLink == 0) || (cki.direction == 1 && cki.lLink == 0)) { 
			hashTable[n].table[hashAdd].ku.ki.processed = 1 ;
		   	break; 
		} else {
			if(cki.direction == 0)
			{
				for(uint8_t z = 0; z < 4; z++)
				{
					if((cki.rLink & (1<<z)) > 0)
					{
						getNextKmer(k.kmer, ck.kmer, z, KMER);
					}	
				}
			} else { // cki.direction == 1
				for(uint8_t z = 0; z < 4; z++)
				{
					if((cki.lLink & (1<<z)) > 0)
					{
						uint8_t base = (~z) & 0x3 ;
						getNextKmer(k.kmer, rkmer, base, KMER);
					}
				}
			}
			hashTable[n].table[hashAdd].ku.ki.processed = 1 ;
		}
	}

	return ka ;
}


					
inline kmer_t GetMinKmer(const kmer_t k)
{
	kmer_t rk ;
	getRevKmer(k.kmer, rk.kmer, KMER);
	int flag = 0, len = (KMER + 32 -1)/32 ;
	for(int i = 0; i < len; i++)
	{
		if(k.kmer[i] > rk.kmer[i]) { flag = 1; break; }
		else if(k.kmer[i] < rk.kmer[i]) { break; }
	}
	
	if(flag == 1) return rk ;
	else return k;
}

void WriteEdgesTogz(const Edge *edge, const int size, gzFile edgesgz)
{
	for(int i = 1; i < size; i++)
	{
		char *s = (char*)xcalloc(PATH_LEN, sizeof(char));
		sprintf(s, "%d\t", i);
		for(int j = 0; j < 4; j++) 
		{
			if(edge[i].lvertex[j] > 0) { sprintf(s + strlen(s), "%d:", edge[i].lvertex[j]); }
		}
		sprintf(s + strlen(s), "\t");

		for(int j = 0; j < 4; j++)
		{
			if(edge[i].rvertex[j] > 0) { sprintf(s + strlen(s), "%d:", edge[i].rvertex[j]); }
		}
		sprintf(s + strlen(s), "\n");

		gzputs(edgesgz, s);

		// final and free work
		free(s);
	}
}
	
void WriteDBGedgesgz(const VerticesIndex *vIdx, gzFile edgesgz)
{
	int size = vIdx->count + 1 ;
	Edge *edge = (Edge*)xcalloc(size, sizeof(Edge));
	int len = (KMER + 32 -1)/32 ;

	for(int i = 0; i < vIdx->count; i++)
	{
		kmer_t k = vIdx->ekArr[i].k ;
		kmer_info ki = k.ku.ki ;
		if(ki.direction == 1)
		{
			kmer_t rk;
			getRevKmer(k.kmer, rk.kmer, KMER);
			for(uint8_t z = 0; z < 4; z++)
			{
				if((ki.lLink (1<<z)) > 0)
				{
					uint8_t base = (~z) & 0x3 ;
					kmer_t ck ;
					getNextKmer(ck.kmer, rk.kmer, base, KMER);
					KmerToVer ktv;
					ktv.k = GetMinKmer(ck);
					KmerToVer *kv = bsearch(ktv, vIdx->skArr, vIdx->count, sizeof(KmerToVer), compare_KmerToVer);
					if(kv == NULL) {fprintf(stderr, "[WriteDBGedgesgz] the start kmer of vertex set error, exit..\n"); exit(1); }
					edge[vIdx->ekArr[i].ID].rvertex[base] = kv->ID ;
				}
			}
		} else {
			for(uint8_t z = 0; z < 4; z++)
			{
				if((ki.rLink (1<<z)) > 0)
				{
					kmer_t ck ;
					getNextKmer(ck.kmer, k.kmer, z, KMER);
					KmerToVer ktv ;
					ktv.k = GetMinKmer(ck);
					KmerToVer *kv = bsearch(ktv, vIdx->skArr, vIdx->count, sizeof(KmerToVer), compare_KmerToVer);
					if(kv == NULL) { fprintf(stderr, "[WriteDBGedgesgz] the start kmer of vertex set error, exit...\n"); exit(1); }
					edge[vIdx->ekArr[i].ID].rvertex[z] = kv->ID ;
				}
			}
		}
	}

	// fill lvertex part
	for(int i = 1; i < size; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			int rID = edge[i].rvertex[j]; 
			if(rID > 0)
			{
				for(int z = 0; z < 4; z++)
				{
					if(edge[rID].lvertex[z] == 0) { edge[rID].lvertex[z] = i; }
				}
			}
		}
	}

	WriteEdgesTogz(edge, size, edgesgz);
	// final and free work
	free(edge);
}

inline Vertex SetVertexHead(char *s)
{
	Vertex v ; memset(&v, 0, sizeof(Vertex));
	char *t = strtok(s, " \t\n");	v.ID = atoi(t);
	t = strtok(NULL, " \t\n");		v.s_size = atoi(t):
	t = strtok(NULL, " \t\n");		v.indegree = atoi(t):
	t = strtok(NULL, " \t\n");		v.outdegree = atoi(t);

	return v ;
}

void SetVertexFreq(Vertex *v, const char *s)
{
	int size = v->s_size - KMER + 1;
	v->freq = (uint8_t*)xcalloc(size, sizeof(uint8_t));
	int i = 0;
	char *t = strtok(s, " \t\n");
	while(t != NULL)
	{
		if(i >= size) { fprintf(stderr, "[SetVertexFreq] size set error, exit...\n");  exit(1); }
		v->freq[i] = atoi(t);
		i++;
		t = strtok(s, " \t\n");
	}
}

VertexArr *DBGVerticesRestore(const char *gzfn)
{
	gzFile gzfp = xzopen(name, "r");
	VertexArr *vArr = (VertexArr*)xcalloc(1, sizeof(VertexArr));
	vArr->size = 100000 ;
	vArr->arr = (Vertex*)xcalloc(vArr->size, sizeof(Vertex));
	char buf[PATH_LEN];
	while(gzgets(gzfp, buf, PATH_LEN) != EOF)
	{
		Vertex v = SetVertexHead(buf);
		if(vArr->count <= v.ID) { vArr->count = v.ID + 1; }
		if(v.ID >= vArr->size)
		{
			vArr->size <<= 1;
			vArr->arr = (Vertex*)xrecalloc(vArr->arr, vArr->size * sizeof(Vertex));
		}
		vArr->arr[v.ID] = v ;
		int len = v.s_size * 5 + 10;
		char *s = (char*)xcalloc(len, sizeof(char));
		gzgets(gzfp, s, len);
		vArr->arr[v.ID].seq = strdup(strtok(s, " \t\n"));
		gzgets(gzfp, s, len); // ignore third identify ID line
		gzgets(gzfp, s, len);
		SetVertexFreq(&vArr->arr[v.ID], s);
		
		// final and free work
		free(s);
	}

	// realloc vArr->arr
	vArr->arr = (Vertex*)realloc(vArr->arr, vArr->count * sizeof(Vertex));
	// final and free work
	gzclose(gzfp);

	return vArr;
}

void DBGEdgesRestore(VertexArr *vArr, const char *gzfn)
{
	gzFile gzfp = xzopen(gzfn, "r");
	char *s = (char*)xcalloc(PATH_LEN, sizeof(char));

	while(gzgets(gzfp, s, PATH_LEN) != EOF)
	{
		char *t = strtok(s, " \t\n");
		int ID = atoi(t);
		t = strtok(NULL, " \t\n"); 
		char *t1 = strtok(NULL, " \t\n");
		int i = 0;
		char *e = strtok(t, " \t\n:");
		while(e != NULL)
		{
			vArr->arr[ID].lvertex[i] = atoi(e);
			i++;
			e = strtok(NULL, " \t\n:");
		}
		// set rigth vertex ID
		i = 0;
		e = strtok(t1, " \t\n:");
		while(e != NULL)
		{
			vArr->arr[ID].rvertex[i] = atoi(e);
			i++;
			e = strtok(NULL, " \t\n:");
		}
	}
	// final and free work
	free(s);
	gzclose(gzfp);
}

void VertexArrFree(const VertexArr *verArr)
{
	for(int i = 0; i < verArr->count; i++)
	{
		if(verArr->arr[i].seq != NULL) free(verArr->arr[i].seq);
		if(verArr->arr[i].freq != NULL) free(verArr->arr[i].freq);
	}

	free(verArr->arr);
	free(verArr);
}


void saveDebruijnGraph(const Arguments *arguments, HashTable *hashTable, const int hashNumber, const KmerFreqCurve curve)
{
	char fileName[PATH_LEN] ;
	gzFile vertexsgz, edgesgz;
	strcpy(fileName, arguments->prefix); strcat(fileName, ".DBGvertices.gz");
	vertexsgz = xzopen(fileName, "w");
	strcpy(flieName, arguments->prefix); strcat(flieName, ".DBGedges.gz");
	edgesgz = xzopen(fileName, "w");

	VerticesIndex *vIdx = BulidVerticesIndexAndWriteFile(hashTable,hashNumber, vertexsgz) ;

	// output logfile the summary 
	fprintf(arguments->logfp, "[saveDebruijinGraph] the kmer graph total have %d vertices\n", vIdx->count);

	WriteDBGedgesgz(vIdx, edgesgz);

	// final and free work
	free(vIdx->skArr); free(vIdx->ekArr); free(vIdx);
	gzclose(edgesgz);
	gzclose(vertexsgz);
}

void removingGraphTips(HashTable *hashTable, const int hashNumber, const int id, const KmerFreqCurve curve)
{
	int tip_rm = 0, depth = 0 ;
	if(DEPTH == 0) depth = 4 * KMER ;
	else depth = DEPTH ;
    for(int64_t i = 0; i < hashTable[id].size; i++)
    {
        kmer_t k = hashTable[id].table[i];
        kmer_info ki = k.ku.ki;
        if(ki.processed == 0 && (LINK_BIGGER_THAN_ONE(ki.rLink) || LINK_BIGGER_THAN_ONE(ki.lLink)))
        {
            uint32_t oldv, newv, curv;
            Ku ku;
            ku.ki = ki;
            ku.ki.processed = 0 ; oldv = ku.hi;
            ku.ki.processed = 1 ; newv = ku.hi;
            if((curv = __sync_val_compare_and_swap(&(hashTable[id].table[i].ku.hi), oldv, newv)) == oldv)
            {
                int graphNum = 10, gCount = 0 ;
                BranchQueue branchQueue;
                KmerStack kmerStack;
                kmer_graph *initialGraph = initialGraphStruct(initial_size, &kmerStack, &branchQueue);
                memcpy(&kmerStack.kmerAddr[kmerStack.count].k, &k, sizeof(kmer_t));
                kmerStack.kmerAddr[kmerStack.count].hID = id; kmerStack.kmerAddr[kmerStack.count].addr = i ;
                kmerStack.kmerAddr[kmerStack.count].k.ku.ki.direction = 0; kmerStack.kmerAddr[kmerStack.count].k.ku.ki.processed = 1 ;
                initialGraph[gCount].size = KMER ;
				initialGraph[gCount].depth = depth ;
                //initialGraph[0].edge_contig = (uint64_t*)xcalloc((initialGraph[0].size + 32 -1)/32, sizeof(uint64_t));
                initialGraph[gCount].stackOffset = kmerStack.count;
                initialGraph[gCount].origin_graph_index = -1 ;
                //memcpy(initialGraph[0].edge_contig, kmerStack.kmerAddr[0].k.kmer, (KMER + 32 -1)/32 * sizeof(uint64_t));
                initialGraph[gCount].count = 1 ;
				initialGraph[gCount].freq = kmerStack.kmerAddr[kmerStack.count].k.ku.ki.freq ;
                kmerStack.count++;
				gCount++;
				/* realloc memory if overflow */
				if(gCount + 8 >= graphNum)
				{
					graphNum <<= 1 ;
					initialGraph = (kmer_graph*)xrecalloc(initialGraph, graphNum * sizeof(kmer_graph));
				}

                for(uint8_t z = 0; z < 4; z++)
                {
                   if((ki.rLink & (1<<z)) > 0) { addToQueue(initialGraph, 0, gCount, &branchQueue, z, FORWARD, k.kmer); gCount++; }
                   if((ki.lLink & (1<<z)) > 0) { addToQueue(initialGraph, 0, gCount, &branchQueue, z, BACKWARD, k.kmer); gCount++; }
                }
                while(branchQueue.start < branchQueue.end) // branchQueue is not NULL
                {
					int index = branchQueue.branchKmers[branchQueue.start].kmer_graph_index;
                    getGraphEdge(&branchQueue, initialGraph, &kmerStack, hashTable, hashNumber, 0);
                    // add new branch edges to the queue
                    {
                        //int totalFreq = 0;
                        /*for(int j = 0; j < initialGraph[index].count; j++)
                        {
                            totalFreq += kmerStack.kmerAddr[initialGraph[index].stackOffset + j].k.ku.ki.freq ;
                        }
                        if(initialGraph[index].count > 0) initialGraph[index].freq = totalFreq / initialGraph[index].count;
                        else initialGraph[index].freq = 0; */
                        if(initialGraph[index].depth > 0 && initialGraph[index].state == 3)
                        {
                            kmer_t ck, pk;
                            if(initialGraph[index].count == 1)
                            {
                                int pi = initialGraph[index].origin_graph_index;
                                pk = kmerStack.kmerAddr[initialGraph[pi].stackOffset + initialGraph[pi].count -1].k;
                            } else if(initialGraph[index].count > 1) {
                                pk = kmerStack.kmerAddr[initialGraph[index].stackOffset + initialGraph[index].count -2].k;
                            } else {fprintf(stderr, "[getGraphEdge] initialGraph[index].count error, exit..\n"); exit(1);}
                            ck = kmerStack.kmerAddr[initialGraph[index].stackOffset + initialGraph[index].count -1].k ;
							/* realloc memory if overflow */
							if(gCount + 8 >= graphNum)
							{
								graphNum <<= 1 ;
								initialGraph = (kmer_graph*)xrecalloc(initialGraph, graphNum * sizeof(kmer_graph));
							}
                            gCount =  addBranch(initialGraph, index, gCount, &branchQueue, BIDIRECTION, ck, pk);
                        }
                    }
                }
                // clip tips 
                tip_rm += clipTipsInGraph(initialGraph, &kmerStack, curve, gCount);

                // write kmer info to the Hash Table
                writeKmerToHashTable(&kmerStack, hashTable);

                // free and clean work
                free(initialGraph); free(kmerStack.kmerAddr); free(branchQueue.branchKmers);
            }
        }
    }
	fprintf(stderr, "removing tips: %d\n", tip_rm) ;
}

void setHashInfo(HashTable *hashTable)
{
    for(int64_t i = 0; i < hashTable->size; i++)
    {
        if(hashTable->table[i].ku.ki.freq > 4)
        {
            // the index step has been index two times
            hashTable->table[i].ku.ki.freq = (hashTable->table[i].ku.ki.freq + 1)/4 ;
        } else if(hashTable->table[i].ku.ki.freq > 0) {
            hashTable->table[i].ku.ki.freq = 1 ;
        }
    }
}

// construct de bruijn graph by kmer file
void generateKmerGraph(Arguments *arguments)
{
    char kmerFile_name[PATH_LEN];
    FILE *kmerfp ;
    time_t timeval;
    char *prefix = arguments->prefix ;
    HashTable *hashTable , *buf;
    int64_t hashSize , kmerNum, bufSize = BUF_SIZE;
    int hashNumber ;
    KmerFreqCurve curve;
    char kmerFreq_name[PATH_LEN];
    strcpy(kmerFreq_name, arguments->prefix); strcat(kmerFreq_name, ".kmerFreq");
    curve = parseKmerFreqFile(kmerFreq_name, arguments);
  
    KMER = arguments->K ; NCPU = arguments->NCPU ;

    DEPTH = arguments->depth ;
    strcpy(kmerFile_name, prefix); strcat(kmerFile_name, ".kmer");
    kmerfp = xopen(kmerFile_name, "r");
    if(NCPU > 7)  hashNumber = 7 ;
    else  hashNumber = NCPU ;
    hashTable = (HashTable*)xcalloc(hashNumber, sizeof(HashTable));
    buf = (HashTable*)calloc(hashNumber, sizeof(HashTable));
    fseek(kmerfp, 0 , SEEK_END); 
    kmerNum = ftell(kmerfp)/ sizeof(KmerInfo) ;
    hashSize = (uint64_t)((kmerNum * 1.7) / hashNumber) ;
    hashSize = find_next_prime_kh(hashSize);
    fseek(kmerfp, 0 , SEEK_SET);
    // initial buffer
    for(int i = 0 ; i < hashNumber ; i++)
    {
        buf[i].size = bufSize ;
        buf[i].table = (kmer_t*)calloc(buf[i].size, sizeof(kmer_t));
        // set as the flag between threads , 0 denote empty, 1 denote full , 2 denote finish work 
        buf[i].loadFactor = 0 ;
    }
 
    time(&timeval); fprintf(arguments->logfp, "Begin construct De Bruijn Graph hash table\ntime:\t%s\n", ctime(&timeval));
    omp_set_nested(2) ;
    omp_set_num_threads(hashNumber);
    int64_t host_delay = 0, co_delay = 0;
    #pragma omp parallel sections 
    {
        #pragma omp section 
        {
            KmerInfo *kmerBuf = (KmerInfo*)calloc(bufSize, sizeof(KmerInfo));
            int kmerLenByWord = (KMER + 32 - 1) / 32 ;
			FILE *tmp = xopen("kmer.out", "w");
            while(kmerNum > 0)
            {
                int  max ;
                kmer_t kmer, rkmer ;
                if(kmerNum >= bufSize) { max = bufSize; kmerNum -= max ; }
                else { max = kmerNum ; kmerNum -= max ; }
                fread(kmerBuf,sizeof(KmerInfo), max , kmerfp); 
                for(int i = 0 ; i < max ; i++)
                {
                    // first kmer 
                    kmer.ku.hi = 0; rkmer.ku.hi = 0;
                    processToKmer(&kmerBuf[i],&kmer, KMER, 1);
                    getRevKmer(kmer.kmer, rkmer.kmer, KMER);
                    pushToBuf(buf, &kmer, &rkmer, kmerLenByWord, kmerBuf[i].freq, hashNumber, &host_delay);
					char *s1 = transformKmerToString(kmer.kmer, 0, KMER);
					char *s2 = transformKmerToString(rkmer.kmer, 0, KMER);
					fprintf(tmp, "%s\t%s\n", s1, s2);
					free(s1); free(s2);
					
                    // second kmer
                    kmer.ku.hi = 0; rkmer.ku.hi = 0;
                    processToKmer(&kmerBuf[i], &kmer, KMER, 2);
                    getRevKmer(kmer.kmer, rkmer.kmer, KMER);
                    pushToBuf(buf, &kmer, &rkmer, kmerLenByWord, kmerBuf[i].freq, hashNumber, &host_delay);
					s1 = transformKmerToString(kmer.kmer, 0, KMER);
					s2 = transformKmerToString(rkmer.kmer, 0, KMER);
					fprintf(tmp, "%s\t%s\n", s1, s2);
					free(s1); free(s2);
                }
            }
			fclose(tmp);
            // finalize the buf 
            for(int i = 0 ; i < hashNumber ; i++)
            {
                if(buf[i].count > 0)
                {
                    ////#pragma omp atomic
                    buf[i].loadFactor = 1 ;
                    while(buf[i].loadFactor != 0) { usleep(100); }
                }
                // denote no any more data for process 
                ////#pragma omp atomic
                buf[i].loadFactor = 2 ;
            }
            free(kmerBuf);
        }

        #pragma omp section 
        {
            #pragma omp parallel for schedule(dynamic)
            for(int i = 0 ; i < hashNumber ; i++)
            {
                hashTable[i].size = hashSize ;
                hashTable[i].loadFactor = 0.75 ;
                hashTable[i].table = (kmer_t*)xcalloc(hashTable[i].size, sizeof(kmer_t)) ;
                assignKmerToHash(&hashTable[i], &buf[i], i, &co_delay);
                setHashInfo(&hashTable[i]);
            }
        }
    }
    time(&timeval); fprintf(arguments->logfp, "Finished construct\ntime:\t%s\n", ctime(&timeval));
    fprintf(arguments->logfp, "host_delay: %lu\t co_delay: %lu\n", host_delay, co_delay);
    fflush(arguments->logfp);
    host_delay = 0; co_delay = 0 ;
    // clean and free work 
    for(int i = 0 ; i < hashNumber ; i++)
    {
        free(buf[i].table);
    }
    free(buf);

    // simplify de bruijn graph  
    // removing tips and merging bubbles
    {
        // parallel process kmer by hashNumber 
        //#pragma omp parallel for schedule(dynamic)
        for(int i = 0; i < hashNumber; i++)
        {
            removingGraphTips(hashTable, hashNumber, i, curve);
        }
        // reset processed flag 
        {
            #pragma omp parallel for schedule(dynamic)
            for(int i = 0; i < hashNumber; i++)
            {
               resetHashTable(&hashTable, hashNumber, i); 
            }
        }
        
	   // simplifyGraph by reverse order
		for(int i = hashNumber -1; i >= 0; i--)
        {
            simplifyGraph(hashTable, hashNumber, i, curve);
        }
        // reset processed flag 
        {
            #pragma omp parallel for schedule(dynamic)
            for(int i = 0; i < hashNumber; i++)
            {
               resetHashTable(hashTable, hashNumber, i); 
            }
        }

		// set kmer_t direction flag
		{
			SetDirectionFlag(hashTable, hashNumber);
		}
        
		// reset processed flag 
        {
            #pragma omp parallel for schedule(dynamic)
            for(int i = 0; i < hashNumber; i++)
            {
               resetHashTable(hashTable, hashNumber, i); 
            }
        }
    }
    time(&timeval); fprintf(arguments->logfp, "Finished simplify\ntime:\t%s\n", ctime(&timeval));
    fflush(arguments->logfp);
    // extract unitig from De Bruijn Graph
	//extractUnitig(arguments, hashTable, hashNumber, curve);
	saveDebruijnGraph(arguments, hashTable, hashNumber, curve);
}

int genGraph_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:     %s genGraph [options]\n\n", PROGRAM_NAME);
	fprintf(stderr, "Options:   -p STR          the Prefix of output files [%s]\n", PROGRAM_NAME);
	fprintf(stderr, "           -c STR          the Configure file that contain library information[%s.ini]\n", PROGRAM_NAME);
	fprintf(stderr, "           -N INT          maximum Number of CPU program used( or the number of threads created) [use all the rest of CPU that system allowed, depend on the system load]\n");
	fprintf(stderr, "           -M INT(Giga)    Maximum RAM memory used(must bigger than X times of the total reads length) [no limitation]\n");
	fprintf(stderr, "           -G INT(Mega)    approximate Genome size by prior known\n");
	fprintf(stderr, "           -o STR          program Output directory by running-time[current directory]\n");
	fprintf(stderr, "           -H INT          Heterozygosity of diploid genome, low heterozygosity set 1, middle set 2, high heterozygosity set 3 [1]\n");
	fprintf(stderr, "           -m INT          the Minimum length of unitig[4 * 'K'(Kmer length)]\n");
	fprintf(stderr, "           -h              print this message to the stderr\n");
	return 0 ;
}

int hwgsa_genGraph(int argc , char *argv[])
{
	int c ;
	char name[PATH_LEN];
	Arguments *arguments ;

	arguments = (Arguments*)xcalloc(1 , sizeof(Arguments));

	if(argc < 2) { genGraph_usage(); return 1 ; }
	while((c = getopt(argc, argv, "c:p:N:M:G:o:q:H:h")) >=0)
	{
		switch(c) {
			case 'c':   strcpy(arguments->conffile, optarg); break ;
			case 'p':   strcpy(arguments->prefix, optarg); break ;
			case 'N':   arguments->NCPU = atoi(optarg); break ;
			case 'M':   arguments->maxMem = atoi(optarg); break ;
			case 'G':   arguments->G = atoi(optarg); break ;
			case 'o':   strcpy(arguments->outputDir, optarg); break ;
			case 'H':   arguments->H = atoi(optarg); break ;
			case 'm':   arguments->min_unitig_len = atoi(optarg); break ;
			case 'h':   genGraph_usage() ; exit(1) ;
			case ':':   fprintf(stderr, "option %c need a value\n", optopt); genGraph_usage(); return 1 ;
			case '?':   fprintf(stderr, "[hwgsa_genGraph] unrecognized option '%c'\n",optopt ); genGraph_usage() ; return 1 ;
			default: genGraph_usage() ; return 1 ;
		}
	}

	// check arguments
	if(arguments->outputDir[0] == '\0')
	{
		fprintf(stderr, "outputDir have not set......\n");
		genGraph_usage(); 
		return 1 ;
	} else {
		char *prefix = (char*)xcalloc(PATH_LEN, sizeof(char));
		if(arguments->prefix[0] == '\0') { strcpy(arguments->prefix, PROGRAM_NAME); }
		strcpy(prefix, arguments->outputDir);
		if(prefix[strlen(prefix)-1] != '/') strcat(prefix, "/");
		strcat(prefix, arguments->prefix);
		strcpy(arguments->prefix, prefix);
		free(prefix);
	}
	// check arguments
	{
		FILE *fp ;
		strcpy(name, arguments->prefix); strcat(name, ".args");
		fp = xopen(name, "r");
		checkGenGraphArgs(arguments, fp);
		fclose(fp);
	}
	KMER = arguments->K ;

	strcpy(name, arguments->prefix); strcat(name, ".log");
	arguments->logfp = xopen(name, "a");

#ifdef DEBUG
	// output arguments to log file
	fprintf(arguments->logfp, "#######################arguments#######################\n");
	fprintf(arguments->logfp, "log file:\t%s\n", name);
	fprintf(arguments->logfp, "c:\t\t%s\n", arguments->conffile);
	fprintf(arguments->logfp, "p:\t\t%s\n", arguments->prefix);
	fprintf(arguments->logfp, "N:\t\t%d\n", arguments->NCPU);
	fprintf(arguments->logfp, "M:\t\t%d\n", arguments->maxMem);
	fprintf(arguments->logfp, "G:\t\t%d\n", arguments->G);
	fprintf(arguments->logfp, "K:\t\t%d\n", arguments->K);
	fprintf(arguments->logfp, "o:\t\t%s\n", arguments->outputDir);
	fprintf(arguments->logfp, "q:\t\t%d\n", arguments->qual);
	fprintf(arguments->logfp, "H:\t\t%d\n", arguments->H);
	fprintf(arguments->logfp, "m:\t\t%d\n", arguments->min_unitig_len);
	fprintf(arguments->logfp, "###########################end#########################\n");

	fflush(arguments->logfp); fflush(stderr); fflush(stdout);
#endif
	// generate De bruijn Graph
	generateKmerGraph(arguments);

	// write arguments to file *.args
	{
		FILE *fp ;
		strcpy(name, arguments->prefix); strcat(name, ".args");
		fp = xopen(name, "w");
		writeArgsToFile(arguments, fp);
		fclose(fp);
	}
	// clean and free work 
	fclose(arguments->logfp);
	free_Arguments(arguments);
	return 0 ;
}

