/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/


#ifndef BWT_H
#define BWT_H

#include <stdint.h>
//#include "bwtindex.h"
#include "bwt_gen64.h"

// requirement: (OCC_INTERVAL%16 == 0)
#define OCC_INTERVAL64 0x80
// occ_interval_major 2^^32
//#define OCC_INTERVAL_MAJOR 0x100000000
#define SA_INTERVAL 64
#define BUF_SIZE 0x100000 // buffer memory size

#ifndef BWA_UBYTE
#define BWA_UBYTE
typedef unsigned char ubyte_t;
#endif
typedef uint64_t bwtint_t;

typedef struct {
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
	bwtint_t L2[5]; // C(), cumulative count
	bwtint_t seq_len; // sequence length
	bwtint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	// occurance array, separated to two parts
	uint32_t cnt_table[256];
    bwtint_t *occMajor ;
    // divide information
    uint8_t isDivide ;
    uint32_t divideNumber ;
    DivideBoundary *divideBoundary ;
    uint64_t *sentinelPosition ;
    uint64_t *sentinelSA ;
	// suffix array
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa; // the leftmost high 8 bit used denote whether SA have been accessed, 1 note has been accessed
} bwt_t;

#define bwt_bwt(b, k) ((b)->bwt[(k)/OCC_INTERVAL64*12 + 4 + (k)%OCC_INTERVAL64/16])

/* retrieve a character from the $-removed BWT string. Note that
 * bwt_t::bwt is not exactly the BWT string and therefore this macro is
 * called bwt_B0 instead of bwt_B */
#define bwt_B0(b, k)  (bwt_bwt(b, k)>>((~(k)&0xf)<<1)&3)                    

#define bwt_occ_intv(b, k) ((b)->bwt + (k)/OCC_INTERVAL64*12)

// inverse Psi function
//#define bwt_invPsi(bwt, k , offset)	not_define							        
#define bwt_invPsi(bwt, k )	13							        


#ifdef __cplusplus
extern "C" {
#endif

	char *transformKmerToString(const uint64_t *buf, const int64_t start, const int64_t end);
    char *transformPacToString(const uint8_t *buf, const int64_t start, const int64_t end);
    void getRevKmer(const uint64_t *kmer, uint64_t *rkmer, const uint32_t len);
    uint64_t calculateSA(const bwt_t *bwt,  uint64_t isa); 
    void dumpAndWriteStorage(const unsigned char *input,const int64_t start, const int64_t end, const char *fileName);
    int64_t pac_seq_len(const char *filename);
    void pac_rev_core(const char *filename);


	void bwt_dump_bwt(const char *fn, const bwt_t *bwt);
	void bwt_dump_sa(const char *fn, const bwt_t *bwt);

	bwt_t *bwt_restore_bwt(const char *fn);
    bwt_t *bwt_restore_bwt_core(const char *fn);
	void bwt_restore_sa(const char *fn, bwt_t *bwt);
    uint8_t *restore_lcp(const char *fn, const int64_t seq_len);

	void bwt_destroy(bwt_t *bwt);

	void bwt_bwtgen(const char *fn_pac, const char *fn_bwt); // from BWT-SW
	void bwt_cal_sa(bwt_t *bwt, int intv);

	void bwt_bwtupdate_core(bwt_t *bwt);
	extern bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c);

	extern void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4]);
	bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k);

	// more efficient version of bwt_occ/bwt_occ4 for retrieving two close Occ values
	void bwt_gen_cnt_table(bwt_t *bwt);
	extern void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol);
	extern void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4]);

	int bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end);
	int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0);

#ifdef __cplusplus
}
#endif

#endif
