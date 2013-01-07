#ifndef KFC_SEQIO_H
#define KFC_SEQIO_H

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>

#ifndef MINBUFSIZE
#define MINBUFSIZE 8192
#endif
#ifndef MAXREADLEN
#define MAXREADLEN 100 
#endif

#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)
#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

extern const unsigned char nst_nt4_table[];
extern const char BIT_NT_CHAR[];
extern const char bit_nt_char[];

typedef struct {
    uint64_t bnt[(MAXREADLEN>>5) + 1] ;
    uint8_t  len ;
    uint8_t highQLen;
} read_bnt ;

typedef struct {
    char *buf ; 
    int begin , end , is_eof ;
    gzFile f;
} kstream_t ;
#ifndef KSTRING_T
#define KSTRING_t kstring_t
typedef struct {
    size_t locPos, max ;
    char *s ;
} kstring_t ;
#endif

typedef struct {
    kstring_t name, comment, seq, qual ;
    int last_char ;
    kstream_t *f ;
} kseq_t ;

typedef struct {
    kstring_t name ;
    kstring_t comment;
    kstring_t seq ;
    kstring_t qual;
    uint32_t strand:1;
} faElem_t ;

typedef struct {
    faElem_t *elem;
    int len;
    int max ;
} faSeq_t ;

void bit64ToChar(char *s, const uint64_t *bt, const int len);
char **readFileList(char *filepath, int *fileNum );
kseq_t *kseq_init(gzFile fd) ;
void kseq_destroy(kseq_t *ks);
int kseq_read(kseq_t *seq);
extern const char BIT_NT_CHAR[] ; 

#endif

