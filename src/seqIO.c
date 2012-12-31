#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "seqIO.h"
#include "utils.h"

#ifndef MAXFILENAMELEN
#define MAXFILENAMELEN 1024
#endif
 
const unsigned char nst_nt4_table[256] = {
    4,4,4,4,  4,4,4,4,  4,4,4,4,  4,4,4,4,
    4,4,4,4,  4,4,4,4,  4,4,4,4,  4,4,4,4,
    4,4,4,4,  4,4,4,4,  4,4,4,4,  4,5/*'-'*/,4,4,
    4,4,4,4,  4,4,4,4,  4,4,4,4,  4,4,4,4,
    4,1,4,0,  4,4,4,3,  4,4,4,4,  4,4,4,4,
    4,4,4,4,  2,4,4,4,  4,4,4,4,  4,4,4,4,
    4,1,4,0,  4,4,4,3,  4,4,4,4,  4,4,4,4,
    4,4,4,4,  2,4,4,4,  4,4,4,4,  4,4,4,4,
    4,4,4,4,  4,4,4,4,  4,4,4,4,  4,4,4,4,
    4,4,4,4,  4,4,4,4,  4,4,4,4,  4,4,4,4,
    4,4,4,4,  4,4,4,4,  4,4,4,4,  4,4,4,4,
    4,4,4,4,  4,4,4,4,  4,4,4,4,  4,4,4,4,
    4,4,4,4,  4,4,4,4,  4,4,4,4,  4,4,4,4,
    4,4,4,4,  4,4,4,4,  4,4,4,4,  4,4,4,4,
    4,4,4,4,  4,4,4,4,  4,4,4,4,  4,4,4,4,
    4,4,4,4,  4,4,4,4,  4,4,4,4,  4,4,4,4,
};

const char bit_nt_char[] = {'c', 'a', 't', 'g' };
const char BIT_NT_CHAR[] = {'C', 'A', 'T', 'G' };

void bit64ToChar(char *s, const uint64_t *bt, const int len)
{
    uint64_t *t = (uint64_t*)xcalloc((len + 32 -1)/32, sizeof(uint64_t));
    memcpy(t, bt, ((len + 32 -1)/32) * sizeof(uint64_t));
    for(int i = len -1; i >= 0; i--)
    {
        uint8_t base = t[i/32] & 0x3;
        s[i] = BIT_NT_CHAR[base];
        t[i/32] >>= 2 ;
    }

    // free and clear
    free(t);
}

static inline kstream_t *ks_init(gzFile f)
{
    kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));
    ks->f = f;
    ks->buf = (char*)calloc(MINBUFSIZE, sizeof(char));
    return ks;
}

static inline void ks_destroy(kstream_t *ks)
{
    if(ks)
    {
        free(ks->buf);
        free(ks);
    }
}

static inline int ks_getc(kstream_t *ks)
{
    if(ks->is_eof && ks->begin >= ks->end) return -1 ; 
    if(ks->begin >= ks->end)
    {
        ks->begin = 0 ;
        ks->end = gzread(ks->f, ks->buf, MINBUFSIZE);
        if(ks->end < MINBUFSIZE) ks->is_eof = 1;
        if(ks->end == 0) return -1 ;
    }
    return (int)ks->buf[ks->begin++];
}

static int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret)
{
    if(dret) *dret = 0;
    str->locPos = 0 ;
    if(ks->begin >= ks->end && ks->is_eof) return -1 ;
    for(;;)
    {
        int i ;
        if(ks->begin >= ks->end)
        {
            if(!ks->is_eof)
            {
                ks->begin = 0;
                ks->end = gzread(ks->f, ks->buf, MINBUFSIZE);
                if(ks->end < MINBUFSIZE) ks->is_eof = 1 ;
                if(ks->end == 0) break ;
            } else break ;
        }
        if(delimiter)
        {
            for(i = ks->begin ; i < ks->end ; ++i)
                if(ks->buf[i] == delimiter) break ;
        } else {
            for(i = ks->begin ; i< ks->end ; ++i)
                if(isspace(ks->buf[i])) break ;
        }
        if(str->max - str->locPos < i - ks->begin +1)
        {
            str->max = str->locPos + (i - ks->begin) + 1 ;
            kroundup32(str->max);
            str->s = (char*)realloc(str->s, str->max);
        }
        memcpy(str->s +str->locPos, ks->buf +ks->begin, i - ks->begin);
        str->locPos = str->locPos + (i - ks->begin);
        ks->begin = i +1 ;
        if(i < ks->end)
        {
            if(dret) *dret = ks->buf[i];
            break ;
        }
    }
    str->s[str->locPos] = '\0';
    return str->locPos ;
}

kseq_t *kseq_init(gzFile fd)
{
    kseq_t *s = (kseq_t*)calloc(1,sizeof(kseq_t));
    s->f = ks_init(fd);
    return s;
}

void kseq_destroy(kseq_t *ks)
{
    if(!ks) return ;
    free(ks->name.s); 
    free(ks->comment.s); 
    free(ks->seq.s); 
    free(ks->qual.s);
    ks_destroy(ks->f);
    free(ks);
}

/* Return value :
    >= 0 length of the sequence (normal)
    -1 end-of-file
    -2 truncated quality string
*/
int kseq_read(kseq_t *seq)
{
    int c;
    kstream_t *ks = seq->f ;
    if(seq->last_char == 0) 
    { /* then jump to the next header line */
        while((c = ks_getc(ks)) != -1 && c != '>' && c != '@');
        if(c == -1) return -1 ; /* end of file*/
        seq->last_char = c;
    } /*  the first header char has been read */
    seq->comment.locPos = seq->seq.locPos = seq->qual.locPos = 0 ;
    if(ks_getuntil(ks, 0, &seq->name, &c) < 0) return -1 ;
    //if(c == '\r') { fprintf(stderr, "The input file maybe the Windows Operating System format, please check, program exit...\n");  exit(1); }
    if(c != '\n') ks_getuntil(ks, '\n', &seq->comment, 0);
    while((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@')
    {
        if(isgraph(c))
        { /* printable non-space character */
            if((seq->seq.locPos +1) >= seq->seq.max) /* double the memory */
            {
                seq->seq.max += (seq->seq.max + 2) ;
                //seq->seq.max = seq->seq.locPos +2 ;
                //kroundup32(seq->seq.max); /* rounded to next closest 2^k*/
                seq->seq.s = (char*)realloc(seq->seq.s,seq->seq.max);
            }
            seq->seq.s[seq->seq.locPos++] = (char)c;
        }
    }
    if(c == '>' || c == '@') seq->last_char = c ; /* the first header char has been read*/
    seq->seq.s[seq->seq.locPos] = '\0' ; /* null terminated string */
    if(c != '+') return seq->seq.locPos ; /* FASTA*/
    if(seq->qual.max < seq->seq.max)
    { /* allocate enough memory */
        seq->qual.max = seq->seq.max ;
        seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.max);
    }
    while((c = ks_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */
    if(c == -1) return -2 ; /* we should not stop here*/
    while((c = ks_getc(ks)) != -1 && seq->qual.locPos < seq->seq.locPos)
        if(c >= 33 && c <= 127) seq->qual.s[seq->qual.locPos++] = (unsigned char)c;
    seq->qual.s[seq->qual.locPos] = 0 ; /* null terminated string */
    seq->last_char = 0 ; /* we have not come to the next header line */
    if(seq->seq.locPos != seq->qual.locPos) return -2 ; /* qual string is shorter than seq string */
    return seq->seq.locPos ;
}
/*
char **readFileList(char *filepath, int *fileNum)
{
    int maxfileNum = 50;
    char **filelist ;
    char filename[MAXFILENAMELEN];
    FILE *fp = xopen(filepath,"r");
    filelist = (char**)calloc(maxfileNum,sizeof(char*));
    while(fgets(filename, sizeof(filename) , fp))
    {
        if(*fileNum > maxfileNum)
        {
            maxfileNum <<= 1 ;
            filelist = (char**)realloc(filelist,maxfileNum);
        }
        filename[strcspn(filename, " \t\n")] = '\0' ;
        filelist[(*fileNum)++] = strdup(filename);
    }
    return filelist;
}
*/
