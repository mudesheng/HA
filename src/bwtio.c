#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "bwt.h"
#include "utils.h"
#include "seqIO.h"

void dumpAndWriteStorage(const unsigned char *input,const int64_t start, const int64_t end, const char *fileName)
{
    int64_t textLength = end - start;
    uint8_t lastByteLen ;
    FILE *fp ;
    uint8_t *textBuffer = (uint8_t*)xcalloc(textLength/CHAR_PER_BYTE + 1, sizeof(uint8_t));

    fp = xopen(fileName, "w");
    for(int64_t i = 0, j = start; i < textLength && j < end ; i++, j++)
    {
        uint8_t c = (input[j>>2] >> ((~j & 0x3)<<1)) & 0x3;
        textBuffer[i>>2] |= (c << ((~i & 0x3)<<1)) ;
    }
#ifdef DEBUG
    /*
    {
        int freq[] = {0, 0, 0, 0 };
        for(int j = 0 ; j < textLength ; j++)
        {
            int c = (textBuffer[j>>2] >> ((~j & 0x3)<<1))  & 0x3 ;
            freq[c]++ ;
        }
        fprintf(stderr, "C:%d A:%d T:%d G:%d\n", freq[0], freq[1], freq[2], freq[3]);
    }
    */
#endif
    fwrite(textBuffer, sizeof(uint8_t) ,(textLength / CHAR_PER_BYTE + 1), fp);
    lastByteLen = textLength % CHAR_PER_BYTE ;
    // left aligned
    fwrite(&lastByteLen, sizeof(uint8_t), 1, fp);
    free(textBuffer);
    fclose(fp);
}

 



int64_t pac_seq_len(const char *filename)
{
    FILE *fp ;
    int64_t pac_len ;
    ubyte_t c ;
    fp = xopen(filename, "rb");
    fseek(fp,-1,SEEK_END);
    pac_len = ftell(fp);
    fread(&c,1,1,fp);
    fclose(fp);
    return (pac_len - 1) * 4 + (int)c ;
}

char *transformPacToString(const uint8_t *buf, const int64_t start, const int64_t end)
{
    char *s = xcalloc(end - start + 1, sizeof(char)) ;
    const char ntc[] = {'C', 'A', 'T', 'G'};
    int64_t i ;
    for(i = 0 ; i < end - start ; i++)
    {
        int c = (buf[(start+i)/4] >> ((~((start+i)%4) & 3)<<1)) & 0x3;
        s[i] = ntc[c] ;
    }
    s[i] = '\0' ;
    return s ;
}

char *transformKmerToString(const uint64_t *buf, const int64_t start, const int64_t end)
{
    char *s = xcalloc(end - start + 1, sizeof(char)) ;
    int64_t i ;
	uint64_t *k = xcalloc(1, sizeof(uint64_t) * ((end - start + 32 -1)/32));
	memcpy(k, buf, sizeof(uint64_t) * ((end - start + 32 -1)/32));
	k[(end - start)/32] <<= ((32 - (end - start)%32)<<1);
    for(i = 0 ; i < end - start ; i++)
    {
        int c = (k[(start+i)/32] >> ((~((start+i)%32) & 31)<<1)) & 0x3;
        s[i] = BIT_NT_CHAR[c] ;
    }
    s[i] = '\0' ;
	free(k);
    return s ;
}

void pac_rev_core(const char *filename)
{
    int64_t seq_len, pac_len, buf_count = 0 , rawbuf_size = 0, last_pos ;
    uint8_t *buf, ct, n, *raw_buf ;
    FILE *fp ;
    int over ;
    seq_len = pac_seq_len(filename);
    pac_len = (seq_len + seq_len%CHAR_PER_BYTE + CHAR_PER_BYTE -1) / CHAR_PER_BYTE ;
    buf = (uint8_t*)xcalloc(pac_len , sizeof(uint8_t));
    raw_buf = (uint8_t*)xcalloc(BUF_SIZE , sizeof(uint8_t));
    fp = xopen(filename, "rb+");
    if(feof(fp)) 
    {
        fprintf(stderr, "[pac_rev_core] The input file : '%s' is NULL file, program exit...\n", filename);
        exit(1);
    }
    // process last byte 
    fseek(fp,-2,SEEK_END);
    fread(&ct , sizeof(uint8_t), 1 , fp);
    fread(&n , sizeof(uint8_t) , 1 , fp);
    buf[buf_count >> 2] = ct ;
    buf_count = n ;
    for(int i = 0 ; i < n ; i++)
    {
        buf[buf_count>>2] <<= BIT_PER_CHAR ;
        buf[buf_count>>2] |= (~ct & 0x3);
        buf_count++ ;
        ct >>= BIT_PER_CHAR ;
    }
    fseek(fp , -2 , SEEK_END);
    last_pos = ftell(fp);
    if(last_pos - BUF_SIZE > 0)
    {
        fseek(fp , -(BUF_SIZE) , SEEK_CUR);
        rawbuf_size = BUF_SIZE ;
        over = 0 ;
    }else {
        fseek(fp, -(last_pos), SEEK_CUR);
        rawbuf_size = last_pos ;
        over = 1 ;
    }
    // process remain bytes
    while((rawbuf_size = fread(raw_buf , sizeof(uint8_t), rawbuf_size , fp)) > 0)
    {
        for(int64_t i = rawbuf_size * CHAR_PER_BYTE -1 ; i >= 0 ; i--)
        {
            buf[buf_count>>2] <<= BIT_PER_CHAR ;
            buf[buf_count>>2] |= ((~(raw_buf[i>>2])) & 0x3);
            buf_count++ ;
            raw_buf[i>>2] >>= BIT_PER_CHAR ;
        }
        if(over == 1) break ;
        else {
            fseek(fp, -(rawbuf_size), SEEK_CUR);
            last_pos = ftell(fp);
            if(last_pos - BUF_SIZE > 0)
            {
                fseek(fp , -(BUF_SIZE) , SEEK_CUR);
                rawbuf_size = BUF_SIZE ;
                over = 0 ;
            }else {
                fseek(fp, -(last_pos), SEEK_CUR);
                rawbuf_size = last_pos ;
                over = 1 ;
            }
        }
    }
    fseek(fp,-2,SEEK_END);
    ct = (seq_len*2) % CHAR_PER_BYTE ;
    if(ct == 0)
    {
        fwrite(buf, sizeof(uint8_t) , pac_len , fp);
        fwrite(&ct , sizeof(uint8_t), 1 , fp);
    } else {
        fwrite(buf , sizeof(uint8_t), pac_len - 1 , fp);
        buf[pac_len -1] <<= ((CHAR_PER_BYTE - ct)<<1);
        fwrite(buf + pac_len - 1 , sizeof(uint8_t), 1 , fp);
    }
    fwrite(&ct, sizeof(uint8_t) , 1 , fp);
    fclose(fp);
    free(buf); free(raw_buf);
}

void bwt_dump_bwt(const char *fn, const bwt_t *bwt)
{
	gzFile fp;
    uint64_t occMajorSize ;
	fp = xzopen(fn, "w");
	gzwrite(fp, &bwt->primary, sizeof(bwtint_t) * 1);
	gzwrite(fp, bwt->L2+1, sizeof(bwtint_t)* 4);
    gzwrite(fp, &bwt->isDivide, sizeof(uint8_t)* 1);
    gzwrite(fp, &bwt->divideNumber, sizeof(uint32_t) * 1);
    if(bwt->isDivide == 1)
    {
        gzwrite(fp, bwt->divideBoundary, sizeof(DivideBoundary) * bwt->divideNumber);
    }
    gzwrite(fp, bwt->sentinelPosition, sizeof(uint64_t) * bwt->divideNumber);
    gzwrite(fp, bwt->sentinelSA, sizeof(uint64_t)* bwt->divideNumber);
    occMajorSize = (bwt->seq_len / OCC_INTERVAL_MAJOR64 + 1) * 4 ;
    gzwrite(fp, &occMajorSize, sizeof(bwtint_t) * 1);
    gzwrite(fp, bwt->occMajor, sizeof(bwtint_t) * occMajorSize);
    gzwrite(fp, &bwt->bwt_size, sizeof(bwtint_t) * 1);
	gzwrite(fp, bwt->bwt, sizeof(uint32_t) * bwt->bwt_size);
	gzclose(fp);
}

void bwt_dump_sa(const char *fn, const bwt_t *bwt)
{
	gzFile fp;
	fp = xzopen(fn, "w");
	gzwrite(fp, &bwt->primary, sizeof(bwtint_t) * 1);
	gzwrite(fp, bwt->L2+1, sizeof(bwtint_t) * 4);
	gzwrite(fp, &bwt->sa_intv, sizeof(bwtint_t) * 1);
	gzwrite(fp, &bwt->seq_len, sizeof(bwtint_t) * 1);
	gzwrite(fp, &bwt->n_sa, sizeof(bwtint_t) * 1);
	gzwrite(fp, bwt->sa, sizeof(uint32_t) * bwt->n_sa);
	gzclose(fp);
}

void bwt_restore_sa(const char *fn, bwt_t *bwt)
{
	char skipped[256];
	FILE *fp;
	bwtint_t primary;

	fp = xopen(fn, "rb");
	fread(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
	fread(skipped, sizeof(bwtint_t), 4, fp); // skip
	fread(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	fread(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

	fread(&bwt->n_sa, sizeof(bwtint_t), 1, fp);
	//bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = (uint32_t*)calloc(bwt->n_sa, sizeof(uint32_t));
	//bwt->sa[0] = -1;

	fread(bwt->sa , sizeof(uint32_t), bwt->n_sa , fp);
	fclose(fp);
}

uint8_t *restore_lcp(const char *fn, const int64_t seq_len)
{
    int64_t len = (seq_len + BITS_IN_BYTE - 1) / BITS_IN_BYTE;
    uint8_t *lcp = (uint8_t*)xcalloc(len, sizeof(uint8_t));
    FILE *fp = xopen(fn, "rb");
    
    fread(lcp, sizeof(uint8_t), len, fp);

    fclose(fp);
    return lcp;
}

/* this version not contain occurence information */
bwt_t *bwt_restore_bwt(const char *fn)
{
	bwt_t *bwt;
	gzFile fp;
    //uint64_t *occMajor ;

	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	fp = xzopen(fn, "r");
	//fseek(fp, 0, SEEK_SET);
	gzread(fp, &bwt->primary, sizeof(bwtint_t) * 1);
	gzread(fp, bwt->L2+1, sizeof(bwtint_t) * 4);
    gzread(fp, &bwt->isDivide, sizeof(uint8_t) * 1);
    gzread(fp, &bwt->divideNumber, sizeof(uint32_t) * 1);
    if(bwt->isDivide == 1)
    {
        bwt->divideBoundary = (DivideBoundary*)calloc(bwt->divideNumber, sizeof(DivideBoundary));
        gzread(fp, bwt->divideBoundary, sizeof(DivideBoundary) * bwt->divideNumber);
    }
    bwt->sentinelPosition = (uint64_t*)calloc(bwt->divideNumber, sizeof(uint64_t));
    bwt->sentinelSA = (uint64_t*)calloc(bwt->divideNumber, sizeof(uint64_t));
    gzread(fp, bwt->sentinelPosition, sizeof(uint64_t) * bwt->divideNumber);
    gzread(fp, bwt->sentinelSA, sizeof(uint64_t) * bwt->divideNumber);
	bwt->bwt_size = (bwt->L2[4] + CHAR_PER_WORD - 1) / CHAR_PER_WORD;
	bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, sizeof(uint32_t));
	gzread(fp, bwt->bwt, sizeof(uint32_t) * bwt->bwt_size);
	bwt->seq_len = bwt->L2[4];
    bwt->occMajor = (uint64_t*)calloc((bwt->seq_len / OCC_INTERVAL_MAJOR64 + 1) * 4 , sizeof(uint64_t));
	gzclose(fp);
	bwt_gen_cnt_table(bwt);

	return bwt;
}

/* this version contain the occurence information */
bwt_t *bwt_restore_bwt_core(const char *fn)
{
    bwt_t *bwt;
    FILE *fp;
    int64_t occMajorSize ;

    bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
    fp = xopen(fn, "rb");
    fseek(fp, 0, SEEK_SET);
    fread(&bwt->primary, sizeof(bwtint_t), 1, fp);
    fread(bwt->L2+1, sizeof(bwtint_t), 4, fp);
    fread(&bwt->isDivide, sizeof(uint8_t), 1, fp);
    fread(&bwt->divideNumber, sizeof(uint32_t), 1, fp);
    if(bwt->isDivide == 1)
    {
        bwt->divideBoundary = (DivideBoundary*)calloc(bwt->divideNumber, sizeof(DivideBoundary));
        fread(bwt->divideBoundary, sizeof(DivideBoundary), bwt->divideNumber, fp);
    }
    bwt->sentinelPosition = (uint64_t*)calloc(bwt->divideNumber, sizeof(uint64_t));
    bwt->sentinelSA = (uint64_t*)calloc(bwt->divideNumber, sizeof(uint64_t));
    fread(bwt->sentinelPosition, sizeof(uint64_t), bwt->divideNumber, fp);
    fread(bwt->sentinelSA, sizeof(uint64_t), bwt->divideNumber, fp);
    fread(&occMajorSize, sizeof(bwtint_t), 1, fp);
    bwt->occMajor = (bwtint_t*)calloc(occMajorSize, sizeof(bwtint_t));
    fread(bwt->occMajor, sizeof(bwtint_t), occMajorSize, fp);
    fread(&bwt->bwt_size, sizeof(bwtint_t), 1, fp);
    bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, sizeof(uint32_t));
    fread(bwt->bwt, sizeof(uint32_t), bwt->bwt_size, fp);

    bwt->seq_len = bwt->L2[4];

    fclose(fp);
    bwt_gen_cnt_table(bwt);

    return bwt;
}

void bwt_destroy(bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->occMajor); free(bwt->sa); free(bwt->bwt); 
    if(bwt->isDivide == 1) free(bwt->divideBoundary);
    free(bwt->sentinelPosition); free(bwt->sentinelSA);
	free(bwt);
}
