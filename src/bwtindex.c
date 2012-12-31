#include "bwtindex.h"
#include "bwt.h"

const uint8_t SET_BIT[8] = { 0X80, 0X40, 0X20, 0X10, 0X08 , 0X04 , 0X02, 0X01 } ;
const uint8_t RESET_BIT[8] = { 0X7F, 0XBF, 0XDF, 0XEF, 0XF7, 0XFB, 0XFD, 0XFE } ;

// parse kmer Frequency characteristic
KmerFreqCurve parseKmerFreqFile(const char *kmerFreq_name, Arguments *arguments)
{
    KmerFreqCurve curve; memset(&curve, 0, sizeof(KmerFreqCurve));
    FILE *kmerfreqF = xopen(kmerFreq_name, "r");
    uint64_t lf = 0,ls = 0,lt = 0,rf = 0,rs = 0,rt = 0 ;
    unsigned long int Frequency;
    uint8_t flag = 0, count = 1 ;
    int peakNum = 0 ;
    int peak[3];
    while((fscanf(kmerfreqF,"%lu",&Frequency)) != EOF)
    {
        lf = ls ;
        ls = lt ;
        lt = rf ;
        rf = rs ;
        rs = rt ;
        rt = Frequency ;
        if(flag == 0)
        {
            if(rf>rs && rs<rt)
            {
                curve.valley = count - 1 ;
                flag = 1 ;
            }
        } else {
            if(lf < ls && ls < lt && lt < rf && rf > rs && rs >rt)
            {
                if(peakNum >= 3 || ( peakNum > 0 && (count > peak[0] * 3))) break ;
                peak[peakNum++] = count - 2 ;
            }
        }
        count++;
    }
    if(peakNum > 0) { curve.peak = peak[peakNum - 1]; }
    
    if(curve.valley == 0 || curve.valley > 20 || curve.peak == 0) 
    {
        fprintf(stderr, "program not got a normal kmer frequency distribution, please check input reads from a metagenomics or not enough WGS coverage, program exit...\n ");
        exit(1);
    }

    if(peakNum >1) {
        if(abs(peak[0]*2 - peak[1]) < peak[1]/4) {
            if(arguments->H < 3) {
                arguments->H = 3 ;
                fprintf(stderr, "genome is a  high heterozygosity of diploid or polyploid, but you not set H with 3, program change set by 3 automatic^_^\n");
            }
            curve.peak = (uint8_t)peak[1] ;
        } else {
            fprintf(stderr, "genome maybe high heterozygosity of diploid or polyploid ,but program can't identified, please check genome is a high heterozygosity diploid and set H is 3, program exit....\n");
            exit(1) ;
        }
    } else curve.peak = (uint8_t)peak[0] ;
    if(arguments->G != 0)
    {
        double cov = (double)arguments->seq_len / arguments->G  * ((150 - arguments->K + 1)/150) ; // default average reads len is 150 
        if(abs(curve.peak - (uint8_t)cov) > (uint8_t)(cov/4)) {
            fprintf(stderr, "the different of curve->peak and coveage bigger than expected, and program exit...\n" );
            exit(1) ;
        }
    }
    return curve ;
}

// find not contained "Nn" region in the read
static inline void findReliableRegion(const kstring_t seq, int *begin, int *end)
{
    int b = 0, e = 0; // begin , end

    for(int i = 0; i < seq.locPos;)
    {
        int len = strcspn(seq.s + i, "Nn");
        if(len > e - b)
        {
            b = i;
            e = b + len; 
        }
        i += (len + 1);
    }
    // check if only contain "ACGTacgt"
    if(strspn(seq.s + b, "ACGTacgt") < e - b)
    {
        fprintf(stderr, "[findReableRegion] the sequence read '%s' contain not 'ACGTacgtNn' symbols, please check, program exit...\n", seq.s);
        exit(1);
    } 
    *begin = b; *end = e;
}

static inline uint8_t adjustQuality(const char quality, const int benchmark)
{
    uint8_t q = (quality - benchmark)/ 16;
    // set [<16] to 0, [16,32] to 1, [32,48] to 2, [>48] to 3
    // recover quality formula: Q = q*16 + 16/2
	if(q < 0) q = 0;
	else if(q >3) q = 3;
    /*if(q > 0) q -= 1 ;
    if(q > 3) 
    { 
        fprintf(stderr, "[bntWritePE] the new adjusted quality score: %u set wrong, must be [0,3], please check, program exit...\n", q);
        exit(1);
    }
	*/
    return q ;
}

// read sequence from raw Pairs End file and transform to bits package
void bntwritepacPE(LIB *libpe, BntWriteArgs *bargs, const  Arguments *arguments , omp_lock_t *lock_log )
{
    const int qualitycutoff = arguments->qual;
    const int K = arguments->K ;
    const int benchmark = libpe->qual_benchmark;
    int64_t bntbuf_count =0, qualbuf_count = 0, max_buf = BUF_SIZE ;
    uint8_t *bnt_buf = (uint8_t*)xcalloc(max_buf,sizeof(uint8_t));
    uint8_t *qual_buf = (uint8_t*)xcalloc(max_buf/2,sizeof(uint8_t));
    int64_t max_length_buf_size = INI_SIZE, processed_rd = 0; 

    for(int i = 0 ; i < bargs->count ; i++)
    {
		int64_t unprocessed_rd = 0;
		int fixed_len = 0;
        gzFile fp_1 = xzopen(bargs->readName[i].f1name , "r");
        gzFile fp_2 = xzopen(bargs->readName[i].f2name , "r");
        kseq_t *seq_1, *seq_2 ;
        xassert(bargs->readName[i].paired == 1 , "not contain the paired files name");
        seq_1 = kseq_init(fp_1); seq_2 = kseq_init(fp_2);
        if(libpe->diverse == 1 && libpe->length == NULL)
        {
            libpe->length = (READ_SIZE*)xcalloc(max_length_buf_size ,sizeof(READ_SIZE));
		}
		// iterater traverse all the reads set
		while(kseq_read(seq_1)>0 && kseq_read(seq_2)>0)
		{
			int noQual = 0;
			int rd1_begin, rd1_end, rd2_begin, rd2_end; 
			// calculate reliable read length
			findReliableRegion(seq_1->seq, &rd1_begin, &rd1_end);
			findReliableRegion(seq_2->seq, &rd2_begin, &rd2_end);
			// set fixed read length 
			if(libpe->diverse == 0 && fixed_len == 0) 
			{
				fixed_len = rd2_end - rd2_begin;
				if(fixed_len < K + MIN_KMER_EXTENSION_LENGTH) fixed_len = K + MIN_KMER_EXTENSION_LENGTH; 
				libpe->read_len = fixed_len ;
			}
			if(( rd1_end - rd1_begin >= K + MIN_KMER_EXTENSION_LENGTH)  && (rd2_end - rd2_begin >= K + MIN_KMER_EXTENSION_LENGTH))
			{
				int rd1_len = 0, rd2_len = 0;
				// set reads length
				if(libpe->diverse == 1)
				{
					rd1_len = rd1_end - rd1_begin;
					rd2_len = rd2_end - rd2_begin;
				} else {
					rd1_len = fixed_len ;  rd2_len = fixed_len ;
				}		

				if((seq_1->qual.locPos == seq_1->seq.locPos) && (seq_2->qual.locPos == seq_2->seq.locPos))
				{
					int totalQual = 0 , num_low = 0, q ;
					for(int i = rd1_begin ; i < rd1_begin + rd1_len ; i++)
					{
						q = seq_1->qual.s[i] - benchmark;
						totalQual += q ;
						if(q <= LOW_QUAL_CUTOFF) num_low++;
					}
					for(int i = rd2_begin ; i < rd2_begin + rd2_len ; i++)
					{
						q = seq_2->qual.s[i] - benchmark;
						totalQual += q ;
						if(q <= LOW_QUAL_CUTOFF) num_low++ ;
					}
					// choose by different quality standard by command line setting
					if(qualitycutoff == 3)
					{
						if(num_low > 0 || totalQual < (35 * (rd1_len + rd2_len))) 
						{ unprocessed_rd += 2 ; continue ;}
					} else if(qualitycutoff == 2 ) {
						if(num_low > 2 || totalQual < (30 * (rd1_len + rd2_len))) 
						{ unprocessed_rd += 2 ; continue ;}
					}else if(qualitycutoff == 1){
						if(num_low > 5 || totalQual < (25 * (rd1_len + rd2_len))) 
						{ unprocessed_rd += 2 ; continue ;}
					}
				} else noQual = 1 ;
				// extending bntbuf and qual_buf at same time
				if(bntbuf_count + (rd1_len + rd2_len)*2 >= (max_buf<<2))
				{
					max_buf <<= 1 ;
					bnt_buf = (uint8_t*)xrecalloc(bnt_buf , max_buf * sizeof(uint8_t));
					qual_buf = (uint8_t*)xrecalloc(qual_buf , (max_buf/2) * sizeof(uint8_t));
				}
				// write to the file and buffer memory
				for(int i = rd1_begin; i < rd1_begin + rd1_len; i++)
				{
					uint8_t c = nst_nt4_table[(int)seq_1->seq.s[i]];
					bnt_buf[bntbuf_count>>2] <<= BIT_PER_CHAR ;
					bnt_buf[bntbuf_count>>2] |=  c ;
					bntbuf_count++ ;
					// change q to [0,3]
					uint8_t q ;
					if(noQual == 0) q = adjustQuality(seq_1->qual.s[i], benchmark);
					else q = 2; // if no quality value to provided, set quality is 2
					qual_buf[qualbuf_count>>2] <<= BIT_PER_CHAR ;
					qual_buf[qualbuf_count>>2] |= q ;
					qualbuf_count++ ;
				}
				for(int i = rd2_begin ; i < rd2_begin + rd2_len; i++)
				{
					uint8_t  c = nst_nt4_table[(int)seq_2->seq.s[i]];
					bnt_buf[bntbuf_count>>2] <<= BIT_PER_CHAR ;
					bnt_buf[bntbuf_count>>2] |= c ;
					bntbuf_count++ ;
					// change q to [0,3]
					uint8_t q ;
					if(noQual == 0) q = adjustQuality(seq_2->qual.s[i], benchmark);
					else q = 2; // if no quality value to provided, set quality is 2
					qual_buf[qualbuf_count>>2] <<= BIT_PER_CHAR ;
					qual_buf[qualbuf_count>>2] |= q ;
					qualbuf_count++ ;
				}
				// write reverse-complement strand to the buffer
				{
					int64_t buf_end = bntbuf_count ;
					uint8_t last_char = bnt_buf[(buf_end -1)>>2] ;
					for(int64_t i = buf_end -1; i >= buf_end - rd1_len - rd2_len ; i--)
					{
						uint8_t c ;
						if(i >= ((buf_end>>2)<<2))
						{
							c = (last_char >> ((buf_end - 1 - i)<<1)) & 0x03 ;
						} else {
							c = (bnt_buf[i>>2] >> ((3 - (i & 0x03))<<1)) & 0x03 ;
						}
						bnt_buf[bntbuf_count>>2] <<= BIT_PER_CHAR ;
						bnt_buf[bntbuf_count>>2] |= (~c & 0x03) ;
						bntbuf_count++;
					}
				}
				//free(rever_com);
				// add read length information
				if(libpe->diverse == 1)
				{
					if(processed_rd + 2 > max_length_buf_size)
					{
						max_length_buf_size <<= 1 ;
						libpe->length = (READ_SIZE*)xrecalloc(libpe->length, max_length_buf_size * sizeof(READ_SIZE));
					}
					libpe->length[processed_rd] = rd1_len;
					libpe->length[processed_rd+1] = rd2_len;
				}
				processed_rd += 2 ;
			}else {unprocessed_rd += 2 ; continue ;}
		}
		// print log to log file
		omp_set_lock(lock_log);
		fprintf(arguments->logfp, "finished process pair files: %s and %s....\tprocessed reads %lu,\tunprocessed reads \t%lu\n", bargs->readName[i].f1name , bargs->readName[i].f2name, processed_rd - libpe->number_rd , unprocessed_rd);
		printf("finished process pair files: %s and %s....\tprocessed reads %lu,\tunprocessed reads \t%lu\n", bargs->readName[i].f1name , bargs->readName[i].f2name, processed_rd - libpe->number_rd , unprocessed_rd);
		omp_unset_lock(lock_log);
		
		// add lib_info content
		libpe->number_rd = processed_rd ;
		
		//free memory
		kseq_destroy(seq_2);kseq_destroy(seq_1);
		gzclose(fp_1); gzclose(fp_2) ;
	}
	
	// realloc memory and efficient memory used
	if((libpe->diverse == 1) && (max_length_buf_size > processed_rd)) libpe->length = (READ_SIZE*)realloc(libpe->length,processed_rd * sizeof(READ_SIZE));

	// finalize the buffer 
	{
		uint64_t bc = (bntbuf_count + CHAR_PER_BYTE - 1) / CHAR_PER_BYTE ;        
		uint64_t qc = (qualbuf_count + CHAR_PER_BYTE - 1) / CHAR_PER_BYTE ;        
		if(max_buf > bc)
		{
			bnt_buf = (uint8_t*)realloc(bnt_buf , bc * sizeof(uint8_t)) ;
			qual_buf = (uint8_t*)realloc(qual_buf , qc * sizeof(uint8_t));
		}
		bargs->bntBuf = bnt_buf ;
		bargs->qualBuf = qual_buf ;
		bargs->seq_len = bntbuf_count ;
		if(bntbuf_count <= 0) 
		{
			fprintf(stderr, "[bntwritepacPE] the length of sequence buffer is NULL....exit\n");
			exit(1);
		}
		xassert(bntbuf_count == 2 * qualbuf_count , "[bntwritepacPE] bntbuf_count != qualbuf_count");
	}
}

/*
// read sequence from raw Pair End file and transform to package
void  bntwritepacSE(LIB *libse , BntWriteArgs *bargs , const  Arguments *arguments , omp_lock_t *lock_log )
{
	const int qualitycutoff = arguments->qual;
	const int K = arguments->K ; 
	const int benchmark = libse->qual_benchmark;
	uint64_t bntbuf_count = 0 , qualbuf_count = 0 , max_buf = BUF_SIZE ;
    uint8_t *bnt_buf = (uint8_t*)xcalloc(max_buf , sizeof(uint8_t));
    uint8_t *qual_buf = (uint8_t*)xcalloc(max_buf , sizeof(uint8_t));
    uint64_t processed_rd = 0 , max_length_buf_size = INI_SIZE  ;
    int fixed_len = 0 ;
    for(int i = 0 ; i < bargs->count ; i++)
    {
        gzFile fp = xzopen(bargs->readName[i].f1name, "r");
        uint64_t unprocessed_rd = 0;
        kseq_t *seq = kseq_init(fp); 
        if(libse->diverse == 1)
        {
            if(libse->length == NULL)   libse->length = (READ_SIZE*)xcalloc(max_length_buf_size , sizeof(READ_SIZE));
            while(kseq_read(seq) > 0)
            {
                int rd_begin, rd_end, noQual = 0 ;
                findReliableRegion(seq->seq, &rd_begin, &rd_end);
                if(rd_end - rd_begin >= K + MIN_KMER_EXTENSION_LENGTH)
                {
                    if(seq->qual.locPos == seq->seq.locPos)
                    {
                        int totalQual = 0 , num_low = 0, q;
                        for(int i = rd_begin; i < rd_end; i++)
                        {
                            q = seq->qual.s[i] - benchmark;
                            totalQual += q ;
                            if(q == LOW_QUAL_CUTOFF) num_low++ ;
                        }
                        // choose by different quality value by command line
                        if(qualitycutoff == 3)
                        {
                            if(num_low > 0 || totalQual > (35 * (rd_end - rd_begin))) {unprocessed_rd++ ; continue ;}
                        } else if(qualitycutoff == 2 ) {
                            if(num_low > 2 || totalQual > (30 * (rd_end - rd_begin))) {unprocessed_rd++ ; continue ;}
                        }else { // qualitycutoff == 1
                            if(num_low > 5 || totalQual > (25 * (rd_end - rd_begin))) {unprocessed_rd++ ; continue ;}
                        }
                    } else noQual = 1 ;
                    // extending buffer
                    if( bntbuf_count + (rd_end - rd_begin) >= (max_buf<<2))
                    {
                        max_buf <<= 1 ;
                        bnt_buf = (uint8_t*)xrecalloc(bnt_buf , max_buf * sizeof(uint8_t));
                        qual_buf = (uint8_t*)xrecalloc(qual_buf , max_buf * sizeof(uint8_t));
                    }
                    // write to buffer memory
                    for(int i = rd_begin; i < rd_end; i++)
                    {
                        uint8_t c = nst_nt4_table[(int)seq->seq.s[i]];
                        bnt_buf[bntbuf_count>>2] <<= BIT_PER_CHAR ;
                        bnt_buf[bntbuf_count>>2] |= c ;
                        bntbuf_count++ ;
                        // change q to [0,3]
                        uint8_t q ;
                        if(noQual == 0) q = adjustQuality(seq->qual.s[i], benchmark);
                        else q = 2; // if no quality value to provided, set quality is 2
                        qual_buf[qualbuf_count>>2] <<= BIT_PER_CHAR ;
                        qual_buf[qualbuf_count>>2] |= q ;
                        qualbuf_count++ ;
                    }
                    if(processed_rd + 1 >= max_length_buf_size)
                    {
                        max_length_buf_size <<= 1 ;
                        libse->length = (READ_SIZE*)xrecalloc(libse->length, max_length_buf_size * sizeof(READ_SIZE));
                    }
                    libse->length[processed_rd] = rd_end - rd_begin;
                    processed_rd++ ;
                }else {unprocessed_rd++ ; continue ;}
            }
        }else {
            while(kseq_read(seq) > 0)
            {
                int rd_begin, rd_end, noQual = 0;
                // locate reliable read region
                findReliableRegion(seq->seq, &rd_begin, &rd_end);
                if(fixed_len == 0)
                {
                    fixed_len = rd_end - rd_begin;
                    if(fixed_len < K + MIN_KMER_EXTENSION_LENGTH) fixed_len = K + MIN_KMER_EXTENSION_LENGTH;
                    libse->read_len = fixed_len;
                }
                if((rd_end - rd_begin) >= fixed_len)
                {
                    if(seq->qual.locPos == seq->seq.locPos)
                    {
                        int totalQual = 0 , num_low = 0, q;
                        for(int i = rd_begin; i < rd_begin + fixed_len; i++ )
                        {
                            q = seq->qual.s[i] - benchmark;
                            totalQual += q;
                            if(q == LOW_QUAL_CUTOFF) num_low++ ;
                        }
                        // choose by different quality value by command line
                        if(qualitycutoff == 3)
                        {
                            if(num_low > 0 || totalQual > (35 * fixed_len) ) {unprocessed_rd++ ; continue ; }
                        } else if(qualitycutoff == 2) {
                            if(num_low > 2 || totalQual > (30 * fixed_len)) { unprocessed_rd++ ; continue ; }
                        } else { // qualitycutoff == 1
                            if(num_low > 5 || totalQual > (25 * fixed_len)) {unprocessed_rd++ ; continue ; }
                        }

                    } else noQual = 1 ;
                    // extending buffer
                    if(bntbuf_count >= (max_buf<<2))
                    {
                        max_buf <<= 1 ;
                        bnt_buf = (uint8_t*)xrecalloc(bnt_buf , max_buf * sizeof(uint8_t));
                        qual_buf = (uint8_t*)xrecalloc(qual_buf , max_buf * sizeof(uint8_t));
                    }
                    for(int i = rd_begin; i < rd_begin + fixed_len; i++)
                    {
                        uint8_t c = nst_nt4_table[(int)seq->seq.s[i]];
                        bnt_buf[bntbuf_count>>2] <<= BIT_PER_CHAR ;
                        bnt_buf[bntbuf_count>>2] |= c ;
                        bntbuf_count++ ;
                        // change q to [0,3]
                        uint8_t q ;
                        if(noQual == 0) q = adjustQuality(seq->qual.s[i], benchmark);
                        else q = 2; // if no quality value to provided, set quality is 2
                        qual_buf[qualbuf_count>>2] <<= BIT_PER_CHAR ;
                        qual_buf[qualbuf_count>>2] |= q ;
                        qualbuf_count++ ;
                    }
                    processed_rd++ ;
                }else { unprocessed_rd++ ; continue ; }
            }
        }
        omp_set_lock(lock_log); 
        fprintf(arguments->logfp, "finished process single file: %s ....\tprocessed reads %lu,\tunprocessed reads \t%lu\n", bargs->readName[i].f1name , processed_rd - libse->number_rd , unprocessed_rd);
        omp_unset_lock(lock_log);
        libse->number_rd = processed_rd ;
        //free memory
        kseq_destroy(seq);
        gzclose(fp); 
    }
    // realloc memory and efficient memory used
    if((libse->diverse == 1) && (max_length_buf_size > processed_rd)) libse->length = (READ_SIZE*)realloc(libse->length,processed_rd * sizeof(READ_SIZE));
    // finalize the buffer 
    {
        uint64_t c = (bntbuf_count + CHAR_PER_BYTE - 1) >>2 ;        
        if(max_buf > c)
        {
            bnt_buf = (uint8_t*)realloc(bnt_buf , c * sizeof(uint8_t)) ;
            qual_buf = (uint8_t*)realloc(qual_buf , c * sizeof(uint8_t));
        }
        bargs->bntBuf = bnt_buf ;
        bargs->qualBuf = qual_buf ;
        bargs->seq_len = bntbuf_count ;
        if(bntbuf_count <= 0) 
        {
            fprintf(stderr, "[bntwritepacSE] the length of sequence buffer is NULL....exit\n");
            exit(1);
        }
        xassert(bntbuf_count == qualbuf_count , "[bntwritepacSE] bntbuf_count != qualbuf_count");
    }
}
*/

// read sequence from raw Pair End file and transform to package
void  bntwritepacSE(LIB *libse , BntWriteArgs *bargs , const  Arguments *arguments , omp_lock_t *lock_log)
{
	const int qualitycutoff = arguments->qual;
	const int K = arguments->K ; 
	const int benchmark = libse->qual_benchmark;
	uint64_t bntbuf_count = 0 , qualbuf_count = 0 , max_buf = BUF_SIZE ;
    uint8_t *bnt_buf = (uint8_t*)xcalloc(max_buf , sizeof(uint8_t));
    uint8_t *qual_buf = (uint8_t*)xcalloc(max_buf/2 , sizeof(uint8_t));
    uint64_t processed_rd = 0 , max_length_buf_size = INI_SIZE  ;
    
	for(int i = 0 ; i < bargs->count ; i++)
    {
        gzFile fp = xzopen(bargs->readName[i].f1name, "r");
        uint64_t unprocessed_rd = 0;
		int fixed_len = 0 ;
        kseq_t *seq = kseq_init(fp); 
        if(libse->diverse == 1 && libse->length == NULL)
        {
            libse->length = (READ_SIZE*)xcalloc(max_length_buf_size , sizeof(READ_SIZE));
		}
		
		while(kseq_read(seq) > 0)
		{
			int rd_begin, rd_end, noQual = 0 ;
			findReliableRegion(seq->seq, &rd_begin, &rd_end);
			// set fixed read length 
			if(libse->diverse == 0 && fixed_len == 0) 
			{
				fixed_len = rd_end - rd_begin;
				if(fixed_len < K + MIN_KMER_EXTENSION_LENGTH) fixed_len = K + MIN_KMER_EXTENSION_LENGTH; 
				libse->read_len = fixed_len ;
			}
			if(rd_end - rd_begin >= K + MIN_KMER_EXTENSION_LENGTH)
			{
				int rd_len = 0 ;
				// set read length
				if(libse->diverse == 1) rd_len = rd_end - rd_begin ;
				else rd_len = fixed_len ;
				if(seq->qual.locPos == seq->seq.locPos)
				{
					int totalQual = 0 , num_low = 0, q;
					for(int i = rd_begin; i < rd_begin + rd_len; i++)
					{
						q = seq->qual.s[i] - benchmark;
						totalQual += q ;
						if(q == LOW_QUAL_CUTOFF) num_low++ ;
					}
					// choose by different quality value by command line
					if(qualitycutoff == 3)
					{
						if(num_low > 0 || totalQual > (35 * rd_len)) {unprocessed_rd++ ; continue ;}
					} else if(qualitycutoff == 2 ) {
						if(num_low > 2 || totalQual > (30 * rd_len)) {unprocessed_rd++ ; continue ;}
					}else if(qualitycutoff == 1){ // qualitycutoff == 1
						if(num_low > 5 || totalQual > (25 * rd_len)) {unprocessed_rd++ ; continue ;}
					}
				} else noQual = 1 ;
				// extending buffer
				if( bntbuf_count + rd_len*2 >= (max_buf<<2))
				{
					max_buf <<= 1 ;
					bnt_buf = (uint8_t*)xrecalloc(bnt_buf , max_buf * sizeof(uint8_t));
					qual_buf = (uint8_t*)xrecalloc(qual_buf , (max_buf/2) * sizeof(uint8_t));
				}
				// write to the file and buffer memory
				//uint8_t *rever_com = xcalloc((rd_len +  CHAR_PER_BYTE - 1)/CHAR_PER_BYTE, sizeof(uint8_t));
				//int base_count = 0;	
				// write to buffer memory
				for(int i = rd_begin; i < rd_begin + rd_len; i++)
				{
					uint8_t c = nst_nt4_table[(int)seq->seq.s[i]];
					bnt_buf[bntbuf_count>>2] <<= BIT_PER_CHAR ;
					bnt_buf[bntbuf_count>>2] |= c ;
					bntbuf_count++ ;
					// add to the temp buffer for reverse complement
					//rever_com[base_count>>2] <<= BIT_PER_CHAR;
					//rever_com[base_count>>2] |= c;
					//base_count++;
					// change q to [0,3]
					uint8_t q ;
					if(noQual == 0) q = adjustQuality(seq->qual.s[i], benchmark);
					else q = 2; // if no quality value to provided, set quality is 2
					qual_buf[qualbuf_count>>2] <<= BIT_PER_CHAR ;
					qual_buf[qualbuf_count>>2] |= q ;
					qualbuf_count++ ;
				}
				// write reverse-complement strand to the buffer
				{
					int64_t buf_end = bntbuf_count ;
					uint8_t last_char = bnt_buf[(buf_end -1)>>2] ;
					for(int64_t i = buf_end - 1; i >= buf_end - rd_len ; i--)
					{
						uint8_t c ;
						if(i >= ((buf_end>>2)<<2))
						{
							c = (last_char >> ((buf_end - 1 - i)<<1)) & 0x03 ;
						} else {
							c = (bnt_buf[i>>2] >> ((3 - (i & 0x03))<<1)) & 0x03 ;
						}
						bnt_buf[bntbuf_count>>2] <<= BIT_PER_CHAR;
						bnt_buf[bntbuf_count>>2] |= (~c & 0x03);
						bntbuf_count++;
					}
				}
				//free(rever_com);

				if(libse->diverse == 1)
				{
					if(processed_rd >= max_length_buf_size)
					{
						max_length_buf_size <<= 1 ;
						libse->length = (READ_SIZE*)xrecalloc(libse->length, max_length_buf_size * sizeof(READ_SIZE));
					}
					libse->length[processed_rd] = rd_len ;
				}
				processed_rd++ ;
			}else {unprocessed_rd++ ; continue ;}
		}
        omp_set_lock(lock_log); 
        fprintf(arguments->logfp, "finished process single file: %s ....\tprocessed reads %lu,\tunprocessed reads \t%lu\n", bargs->readName[i].f1name , processed_rd - libse->number_rd , unprocessed_rd);
        printf("finished process single file: %s ....\tprocessed reads %lu,\tunprocessed reads \t%lu\n", bargs->readName[i].f1name , processed_rd - libse->number_rd , unprocessed_rd);
        omp_unset_lock(lock_log);
        libse->number_rd = processed_rd ;
        //free memory
        kseq_destroy(seq);
        gzclose(fp); 
    }
    // realloc memory and efficient memory used
    if((libse->diverse == 1) && (max_length_buf_size > processed_rd)) libse->length = (READ_SIZE*)realloc(libse->length,processed_rd * sizeof(READ_SIZE));
    // finalize the buffer 
    {
        uint64_t bc = (bntbuf_count + CHAR_PER_BYTE - 1) / CHAR_PER_BYTE ;        
        uint64_t qc = (qualbuf_count + CHAR_PER_BYTE - 1) / CHAR_PER_BYTE ;        
        if(max_buf > bc)
        {
            bnt_buf = (uint8_t*)realloc(bnt_buf , bc * sizeof(uint8_t)) ;
            qual_buf = (uint8_t*)realloc(qual_buf , qc * sizeof(uint8_t));
        }
        bargs->bntBuf = bnt_buf ;
        bargs->qualBuf = qual_buf ;
        bargs->seq_len = bntbuf_count ;
        if(bntbuf_count <= 0) 
        {
            fprintf(stderr, "[bntwritepacSE] the length of sequence buffer is NULL....exit\n");
            exit(1);
        }
        xassert(bntbuf_count == 2 * qualbuf_count , "[bntwritepacSE] bntbuf_count != qualbuf_count");
    }
}

void BntWriteArgs_free(BntWriteArgs *bntWArgs, const int number)
{
    for(int i = 0 ; i < number ; i++)
    {
        free(bntWArgs[i].readName);
        //free(bntWArgs[i].bntBuf);
        //free(bntWArgs[i].qualBuf);
    }
    free(bntWArgs);
}

// transform raw fa/fq file to bit package
void bns2bntseq(lib_info *libIndex_info, BntWriteArgs *bntWriteArgs, const Arguments *arguments)
{
    char pacfn[PATH_LEN], qualfn[PATH_LEN];
    gzFile pacfp , qualfp ;

    strcpy(pacfn, arguments->prefix);   strcat(pacfn, ".pac.gz");
    strcpy(qualfn, arguments->prefix);  strcat(qualfn, ".qual.gz");
    
    pacfp = xzopen(pacfn, "wb+"); qualfp = xzopen(qualfn, "wb+");

    // transform sequence reads to bit representation format file and package 
    {
        // lock log file handle
        omp_lock_t lock_log ;
        omp_init_lock(&lock_log);
        omp_set_num_threads(arguments->NCPU);
        #pragma omp parallel for firstprivate(bntWriteArgs, arguments, libIndex_info) schedule(dynamic)
        for(int i = 0 ; i < (int)libIndex_info->num_lib ; i++)
        {
            if(bntWriteArgs[i].readName[0].paired == 1)
            {
                bntwritepacPE(&(libIndex_info->lib[i]), &(bntWriteArgs[i]), arguments, &lock_log);
            } else {
                bntwritepacSE(&(libIndex_info->lib[i]), &(bntWriteArgs[i]), arguments, &lock_log);
            }
        }
        omp_destroy_lock(&lock_log);
    }
    // write library package  to .pac file 
    {
        // write the first library 
        //uint64_t totalLen = 0 ;
        int64_t bc = bntWriteArgs[0].seq_len / CHAR_PER_BYTE; // remainder
        int64_t qc = (bntWriteArgs[0].seq_len/2) / CHAR_PER_BYTE; // remainder
        uint8_t b = bntWriteArgs[0].bntBuf[bc], q = bntWriteArgs[0].qualBuf[qc] ;
	    uint8_t	bz = bntWriteArgs[0].seq_len % CHAR_PER_BYTE , qz = (bntWriteArgs[0].seq_len/2) % CHAR_PER_BYTE ;
        gzwrite(pacfp, bntWriteArgs[0].bntBuf , sizeof(uint8_t)* bc);
        gzwrite(qualfp, bntWriteArgs[0].qualBuf , sizeof(uint8_t)* qc);
        libIndex_info->lib[0].offset = bntWriteArgs[0].seq_len ;
        for(int i = 1 ; i < (int)libIndex_info->num_lib ; i++)
        {
            uint8_t *bnt , *qual ;
            uint64_t bsize_total = (bntWriteArgs[i].seq_len + bz + CHAR_PER_BYTE - 1 ) / CHAR_PER_BYTE , size_bnt ;
            uint64_t qsize_total = ((bntWriteArgs[i].seq_len/2) + qz + CHAR_PER_BYTE - 1 ) / CHAR_PER_BYTE , size_qual ;
            free(bntWriteArgs[i-1].bntBuf); free(bntWriteArgs[i-1].qualBuf);
            bnt = (uint8_t*)xcalloc(bsize_total , sizeof(uint8_t));
            qual = (uint8_t*)xcalloc(qsize_total , sizeof(uint8_t));
            if(bz > 0)  bnt[0] = b ;
		    if(qz > 0)	qual[0] = q ; 
            size_bnt = (bntWriteArgs[i].seq_len + CHAR_PER_BYTE - 1) / CHAR_PER_BYTE ;
            size_qual = ((bntWriteArgs[i].seq_len/2) + CHAR_PER_BYTE - 1) / CHAR_PER_BYTE ;
            for(uint64_t j = 0 ; j < size_bnt - 1 ; j++)
            {
                bnt[j] <<= ((CHAR_PER_BYTE - bz)<<1) ;
                bnt[j] |= (bntWriteArgs[i].bntBuf[j] >> (bz<<1) );
                bnt[j+1] = bntWriteArgs[i].bntBuf[j] & masklow8[bz] ;
            }
			for(uint64_t j = 0; j < size_qual - 1; j++)
			{
                qual[j] <<= ((CHAR_PER_BYTE -qz)<<1) ;
                qual[j] |= (bntWriteArgs[i].qualBuf[j] >> (qz<<1));
                qual[j+1] = bntWriteArgs[i].qualBuf[j] & masklow8[qz] ;
			}
            if(bsize_total > size_bnt)
            {
                uint8_t t = (bntWriteArgs[i].seq_len + bz)%CHAR_PER_BYTE ;
                bnt[size_bnt - 1] <<= ((CHAR_PER_BYTE - bz)<<1);
                bnt[size_bnt - 1] |= (bntWriteArgs[i].bntBuf[size_bnt-1] >> (t<<1));
                bnt[size_bnt] = bntWriteArgs[i].bntBuf[size_bnt-1] & masklow8[t];
            } else {
                uint8_t t = (bntWriteArgs[i].seq_len + bz)%CHAR_PER_BYTE - bz ;
                bnt[size_bnt -1] <<= (t<<1);
                bnt[size_bnt -1] |= (bntWriteArgs[i].bntBuf[size_bnt -1] & masklow8[t]);
            }
			if(qsize_total > size_qual)
			{
                uint8_t t = ((bntWriteArgs[i].seq_len/2) + qz)%CHAR_PER_BYTE ;
                qual[size_qual -1] <<= ((CHAR_PER_BYTE - qz)<<1);
                qual[size_qual -1] |= (bntWriteArgs[i].qualBuf[size_qual -1] >> (t<<1));
                qual[size_qual] = bntWriteArgs[i].qualBuf[size_qual - 1] & masklow8[t] ;
			} else {
                uint8_t t = ((bntWriteArgs[i].seq_len/2) + qz)%CHAR_PER_BYTE - qz ;
                qual[size_qual -1] <<= (t<<1);
                qual[size_qual -1] |= (bntWriteArgs[i].qualBuf[size_qual-1] & masklow8[t]);
			}
            // write to .pac and .qual file 
            bsize_total = (bntWriteArgs[i].seq_len + bz) / CHAR_PER_BYTE ;
            qsize_total = ((bntWriteArgs[i].seq_len/2) + qz) / CHAR_PER_BYTE ;
            gzwrite(pacfp, bnt , sizeof(uint8_t) * bsize_total);
            gzwrite(qualfp, qual , sizeof(uint8_t) * qsize_total);
            b = bnt[bsize_total] ;
            q = qual[qsize_total] ;
            bz = (bntWriteArgs[i].seq_len + bz) % CHAR_PER_BYTE ;
            qz = ((bntWriteArgs[i].seq_len/2) + qz) % CHAR_PER_BYTE ;
            libIndex_info->lib[i].offset = bntWriteArgs[i].seq_len ;

            // clean and free work
            free(bnt); free(qual);
        }
        free(bntWriteArgs[libIndex_info->num_lib -1].bntBuf); free(bntWriteArgs[libIndex_info->num_lib -1].qualBuf);
        //finalize .pac and .qual file 
        if(bz > 0)  { b <<= ((CHAR_PER_BYTE - bz)<<1); gzwrite(pacfp, &b , sizeof(uint8_t) * 1); }
        else   gzwrite(qualfp, &bz, sizeof(uint8_t)* 1);
        if(qz > 0) { q <<= ((CHAR_PER_BYTE - qz)<<1); gzwrite(qualfp, &q , sizeof(uint8_t) * 1); }
		else gzwrite(qualfp, &qz, sizeof(uint8_t) *1);
        gzwrite(pacfp, &bz , sizeof(uint8_t) * 1);
        gzwrite(qualfp, &qz , sizeof(uint8_t) * 1 );
    }
    gzclose(pacfp) ; gzclose(qualfp);
    // complete libIndex_info structure
    libIndex_info->len_pac = 0 ;
    for(int i = 0 ; i < (int)libIndex_info->num_lib ; i++)
    {
        uint64_t m  = libIndex_info->len_pac ;
        libIndex_info->len_pac += libIndex_info->lib[i].offset ;
        libIndex_info->lib[i].offset = m ;
        libIndex_info->num_seqs += libIndex_info->lib[i].number_rd ;
    }

    // clean and free work
    BntWriteArgs_free(bntWriteArgs, libIndex_info->num_lib);
}

// return isa
uint64_t setBWTSA(bwt_t *bwt, const int64_t start, const int64_t end,  const int offset , const int divideCount )
{
    uint64_t isa ; // S(isa) = sa 
    int64_t i ;
    // calculate SA value
    if(bwt->isDivide == 1) isa = bwt->divideBoundary[divideCount].bwtLocation ;
    else isa = 0 ;
    //sa = start ;
    for(i = start ; i > end ; i-- )
    {
        int c ;
        int64_t bp ;
        if(isa % SA_INTERVAL == 0) bwt->sa[isa/SA_INTERVAL] = i ;
        bp = posSA2PosBWT(bwt, isa);
        c = bwt_B0(bwt, bp);
        isa = bwt->L2[c] + bwt->occMajor[(bp-1)/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp-1, c) + offset;
    }
    if(isa % SA_INTERVAL == 0) bwt->sa[isa/SA_INTERVAL] = i;
    bwt->sa[0] = (bwtint_t)-1 ; // before this line, bwt->sa[0] = bwt->seq_len
#ifdef DEBUG
    {
        int flag = 0 ;
        for(int i = 0 ;  i < bwt->divideNumber; i++)
        {
            if(isa == bwt->sentinelPosition[i]) { flag = 1 ; break; }
        }
        if(flag != 1) fprintf(stderr, "[setBWTSA] isa != bwt->sentinelPosition\n");
    }
#endif
    return isa ;
}

/* // update lcpBound array 
void updateLcpBound(LCPBound *lcpBound, const int lcpIndex,const uint32_t c,const int remainLen, bwt_t *bwt, const uint32_t offset)
{
    int i ;
    
    if(remainLen >= (int)KMER)
    {
        lcpBound[lcpIndex%KMER].lowBound = bwt->L2[c] + offset ; lcpBound[lcpIndex%KMER].highBound =  bwt->L2[c+1] - 1 + offset ;
        if(lcpIndex < (int)KMER)
        {
            for(i = lcpIndex - 1 ; i >= 0 ; i-- )
            {
                int z ;
                for(z = 0 ; z < bwt->divideNumber ; z++)
                {
                    if(lcpBound[i].lowBound < bwt->sentinelPosition[z])   break; 
                }
                lcpBound[i].lowBound -= (z + 1) ;
                lcpBound[i].lowBound = bwt->L2[c] + bwt->occMajor[lcpBound[i].lowBound/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt,lcpBound[i].lowBound , c) + offset;
                for( z = 0 ; z < bwt->divideNumber ; z++)
                {
                    if(lcpBound[i].highBound < bwt->sentinelPosition[z]) break; 
                }
                lcpBound[i].highBound -= (z+1) ;
                lcpBound[i].highBound = bwt->L2[c] + bwt->occMajor[lcpBound[i].highBound/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, lcpBound[i].highBound, c) + offset ;
            }
        } else {
            for(i = lcpIndex - 1 ; i > lcpIndex - (int)KMER ; i-- )
            {
                int z ;
                for(z = 0 ; z < bwt->divideNumber ; z++)
                {
                    if(lcpBound[i%KMER].lowBound < bwt->sentinelPosition[z])  break; 
                }
                lcpBound[i%KMER].lowBound -= (z+1) ;
                lcpBound[i%KMER].lowBound = bwt->L2[c] + bwt->occMajor[lcpBound[i%KMER].lowBound/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt,lcpBound[i%KMER].lowBound , c) + offset ;
                for(z = 0 ; z < bwt->divideNumber ; z++)
                {
                    if(lcpBound[i%KMER].highBound < bwt->sentinelPosition[z]) break; 
                }
                lcpBound[i%KMER].highBound -= (z+1) ;
                lcpBound[i%KMER].highBound = bwt->L2[c] + bwt->occMajor[lcpBound[i%KMER].highBound/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, lcpBound[i%KMER].highBound, c) + offset ;
            }
        }
    }else {
        if(lcpIndex < (int)KMER)
        {
            for(i = lcpIndex - (int)KMER + remainLen; i >= 0 ; i-- )
            {
                int z ;
                for(z = 0 ; z < bwt->divideNumber ; z++)
                {
                    if(lcpBound[i].lowBound < bwt->sentinelPosition[z]) break; 
                }
                lcpBound[i].lowBound -= (z+1) ;
                lcpBound[i].lowBound = bwt->L2[c] + bwt->occMajor[lcpBound[i].lowBound/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt,lcpBound[i].lowBound , c) + offset ;
                for(int z = 0 ; z < bwt->divideNumber ; z++)
                {
                    if(lcpBound[i].highBound < bwt->sentinelPosition[z]) break; 
                }
                lcpBound[i].highBound -= (z+1) ;
                lcpBound[i].highBound = bwt->L2[c] + bwt->occMajor[lcpBound[i].highBound/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, lcpBound[i].highBound, c) + offset ;
            }
        } else {
            for(i = lcpIndex - (int)KMER + remainLen; i > lcpIndex - (int)KMER ; i--)
            {
                int z ;
                for(z = 0 ; z < bwt->divideNumber ; z++)
                {
                    if(lcpBound[i%KMER].lowBound < bwt->sentinelPosition[z])  break; 
                }
                lcpBound[i%KMER].lowBound -= (z+1) ;
                lcpBound[i%KMER].lowBound = bwt->L2[c] + bwt->occMajor[lcpBound[i%KMER].lowBound/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt,lcpBound[i%KMER].lowBound , c) + offset ;
                for(int z = 0 ; z < bwt->divideNumber ; z++)
                {
                    if(lcpBound[i%KMER].highBound < bwt->sentinelPosition[z]) break; 
                }
                lcpBound[i%KMER].highBound -= (z+1) ;
                lcpBound[i%KMER].highBound = bwt->L2[c] + bwt->occMajor[lcpBound[i%KMER].highBound/OCC_INTERVAL_MAJOR64 * 4 + c] +  bwt_occ(bwt, lcpBound[i%KMER].highBound, c) + offset ;
            }
        }
    }
}
*/
/*
static void setLcp(uint8_t *lcp, const int64_t start , const int64_t end)
{
    const uint8_t setArray[8] = { 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 } ;
    int64_t i , bound;
    if(start  == end + 1) 
    {  
        return ; 
    }else if(start > end + 1) {
        fprintf(stderr, "[setLcp] start > end +1\n"); return ;
    }
    // test if has been set  and return 
    if((lcp[start/BITS_IN_BYTE] & setArray[start%BITS_IN_BYTE]) == setArray[start%BITS_IN_BYTE]) return ; 
    if((end - start)/BITS_IN_BYTE > 1) 
    {
        bound = (start + BITS_IN_BYTE - 1)/BITS_IN_BYTE * BITS_IN_BYTE ;
        for(i = start ; i < bound ; i++)
        {    lcp[i/BITS_IN_BYTE] |= setArray[i%BITS_IN_BYTE]; }
        for(i = bound ; i < end / BITS_IN_BYTE * BITS_IN_BYTE ; i += BITS_IN_BYTE)
        {    lcp[i/BITS_IN_BYTE] = 0xff ; }
        for(i = end / BITS_IN_BYTE * BITS_IN_BYTE ; i <= end ; i++)
        {    lcp[i/BITS_IN_BYTE] |= setArray[i%BITS_IN_BYTE] ; }
    } else {
        for(i = start ; i <= end ; i++)
        {    lcp[i/BITS_IN_BYTE] |= setArray[i%BITS_IN_BYTE] ; }
    }
} */

/* // return isa
uint64_t setBWTSAAndLCP(bwt_t *bwt,uint8_t  *lcp, const uint64_t start, uint64_t end, lib_info *libIndex_info, const uint32_t offset , const int divideCount )
{
    uint64_t isa , sa ; // S(isa) = sa 
    int64_t i ;
    uint32_t c ;
    LCPBound *lcpBound = (LCPBound*)calloc(KMER, sizeof(LCPBound));
    int lcpIndex ;
    uint64_t readBound ;
    // calculate SA value
    if(bwt->isDivide == 1) isa = bwt->divideBoundary[divideCount].bwtLocation ;
    else isa = 0 ;
    sa = start ;
    for(i = (int64_t)start ; i > (int64_t)end ; i = readBound )
    {
        readBound = readBoundaryLookUp(libIndex_info, sa - 5 , BACKWARD);
        lcpIndex = 0 ;
        for(int64_t j = i ; j > (int64_t)readBound ; j--)
        {
            int zcount ;
            int64_t bp ;
            if(isa % SA_INTERVAL == 0) bwt->sa[isa/SA_INTERVAL] = sa ;
            --sa ;
            for(zcount = 0 ; zcount < bwt->divideNumber ; zcount++)
            {
                if(isa < bwt->sentinelPosition[zcount]) break ;
                else if(isa == bwt->sentinelPosition[zcount]){
                    fprintf(stderr, "[setBWTSAAndLCP] isa == bwt->sentinelPosition[zcount], %lu\n", bwt->sentinelPosition[zcount]);
                    break ;
                }
            }
            bp = isa - zcount ;
            c = bwt_B0(bwt, bp);
            isa = bwt->L2[c] + bwt->occMajor[bp/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp, c) + offset - 1 ;
            if(lcpIndex >= (int)KMER) setLcp(lcp,lcpBound[lcpIndex%KMER].lowBound + 1, lcpBound[lcpIndex%KMER].highBound - 1);
            updateLcpBound(lcpBound, lcpIndex, c, j - readBound  , bwt,offset);
            lcpIndex++ ;
        }
        if(isa % SA_INTERVAL == 0) bwt->sa[isa / SA_INTERVAL] = sa ;
        bwt->sa[0] = (bwtint_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
        if(lcpIndex >= (int)KMER) setLcp(lcp, lcpBound[lcpIndex%KMER].lowBound + 1, lcpBound[lcpIndex%KMER].highBound - 1); 
    }
#ifdef DEBUG
    {
        int flag = 0 ;
        for(int i = 0 ; i < bwt->divideNumber; i++)
        {
            if(isa == bwt->sentinelPosition[i]) { flag = 1 ; break; }
        }
        if(flag != 1) fprintf(stderr, "[setBWTSAAndLCP] isa != bwt->sentinelPosition\n");
        if(readBound != end)
        {
            fprintf(stderr, "[setBWTSAAndLCP] readBound != end\n");
        }
        assert(readBound == end) ;
    }
#endif
    return isa ;
} */
/*  transform the position of SA to the position in the BWT string ,
 *  if encounter a sentinelPosition, return -1, else return the position of BWT string
 */
inline int64_t posSA2PosBWT(const bwt_t *bwt, const int64_t isa)
{
    int64_t ret_v = -2 ;
    int i ;
    for(i = 0; i < bwt->divideNumber; i++)
    {
        if(isa < bwt->sentinelPosition[i]) break ;
        else if(isa == bwt->sentinelPosition[i]) {
            fprintf(stderr, "[posSA2PosBWT] isa == bwt->sentinelPosition[%d](%lu)\n", i, bwt->sentinelPosition[i]);
            ret_v = -1 ;
			exit(1);
            break;
        }
    }
    if(ret_v == -2) ret_v = isa - i ; 
    return ret_v ;
}

static inline int checkInterLeave(const bwt_t *bwt, const int64_t low, const int64_t high)
{
    int interLeave = -1;
    for(int i = 0; i < bwt->divideNumber; i++)
    {
        if(low < bwt->sentinelPosition[i] && bwt->sentinelPosition[i] <= high) 
        {
            interLeave = bwt->sentinelPosition[i] - low ; 
            break;
        }
    }
    return interLeave;
}

/*
// flag note low or high with the region , low: 0 high: 1
static inline int64_t posSAChangeToPosBWT(const bwt_t *bwt, const int64_t isa, const int flag)
{
    int i ;
    for(i = 0; i < bwt->divideNumber; i++)
    {
        if(isa < bwt->sentinelPosition[i]) break;
        else if(isa == bwt->sentinelPosition[i]) {
            i += flag ; break;
        }
    }
    return isa - i ;
}
*/

void setBWTLCP(const bwt_t *bwt, uint8_t *lcp, const uint64_t word_part, const int depth, const int first_len, const int min_kmerfreq)
{
    // check arguments
    //int mid_depth = depth/2;
    int64_t low = -1, high = -1 ;
    for(int i = 0; i < first_len; i++)
    {
        int64_t bp_l, bp_h ;
        uint32_t c = GET_CHAR_FROM_UINT64(word_part, (32 - first_len + i));
        if(low == -1 && high == -1)
        {
            low = bwt->L2[c] + bwt->divideNumber ;
            high = bwt->L2[c+1] + bwt->divideNumber - 1;
        } else {
            //bp_l = posSAChangeToPosBWT(bwt, low, 0);
            //bp_h = posSAChangeToPosBWT(bwt, high, 1);
			bp_l = posSA2PosBWT(bwt, low);
			bp_h = posSA2PosBWT(bwt, high);
            low = bwt->L2[c] + bwt->occMajor[(bp_l-1)/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_l-1, c) + bwt->divideNumber;
            high = bwt->L2[c] + bwt->occMajor[bp_h/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_h, c) + bwt->divideNumber -1;
        }
        if(low >= high) break;
    }
    if(low >= high) return;
    int s_size = INITIAL_STACK_SIZE, s_count = 0 ;
    LCPBound *lcpBound = (LCPBound*)xcalloc(s_size, sizeof(LCPBound));
    lcpBound[s_count].low = low; lcpBound[s_count].high = high;
    lcpBound[s_count].depth = depth - first_len; s_count++;
    
    while(s_count > 0)
    {
        LCPBound lcp_elem = lcpBound[--s_count];
        if(lcp_elem.depth < 1) { fprintf(stderr, "[setBWTLCP] depth set error, exit...\n");  exit(1); }
        int64_t bp_l, bp_h;
        bp_l = posSAChangeToPosBWT(bwt, lcp_elem.low, 0);
        bp_h = posSAChangeToPosBWT(bwt, lcp_elem.high, 1);
        int64_t region = bp_h - bp_l + 1;
		for(int i = 0; i < DNA_SYMBOLS_NUMBER; i++)
		{
			low = bwt->L2[i] + bwt->occMajor[bp_l/OCC_INTERVAL_MAJOR64 * 4 + i] + bwt_occ(bwt, bp_l, i) + bwt->divideNumber -1;
			high = bwt->L2[i] + bwt->occMajor[bp_h/OCC_INTERVAL_MAJOR64 * 4 + i] + bwt_occ(bwt, bp_h, i) + bwt->divideNumber -1;
			if(high - low >= min_kmerfreq)
			{
				if(lcp_elem.depth > 1)
				{
					if(s_count >= s_size)
					{
						s_size <<= 1;
						lcpBound = (LCPBound*)xrecalloc(lcpBound, s_size * sizeof(LCPBound));
					}
					lcpBound[s_count].low = low; lcpBound[s_count].high = high;
					lcpBound[s_count].depth = lcp_elem.depth - 1; 
					s_count++;
				} else { // lcp_elem.depth == 1
					for(uint64_t j = low + 1; j <= high; j++)
					{
						__sync_or_and_fetch(lcp + (j>>3), SET_BIT[j & 0X7]);
						//#pragma omp atomic
						//SET_LCP_FLAG(lcp, low + j);
					}
				}
			}
			// test region
			if(high > low) region -= (high - low + 1);
			if(region < min_kmerfreq) break;
		}
    }
    free(lcpBound);
}

/*
void setBWTLCP_bak(const bwt_t *bwt, uint8_t *lcp, const uint64_t word_part, const int depth, const int first_len, const int min_kmerfreq)
{
    // check arguments
    int mid_depth = depth/2;
    int64_t low = -1, high = -1 ;
    for(int i = 0; i < first_len; i++)
    {
        int64_t bp_l, bp_h ;
        uint32_t c = GET_CHAR_FROM_UINT64(word_part, (32 - first_len + i));
        if(low == -1 && high == -1)
        {
            low = bwt->L2[c] + bwt->divideNumber ;
            high = bwt->L2[c+1] + bwt->divideNumber - 1;
        } else {
            bp_l = posSAChangeToPosBWT(bwt, low, 0);
            bp_h = posSAChangeToPosBWT(bwt, high, 1);
            low = bwt->L2[c] + bwt->occMajor[bp_l/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_l, c) + bwt->divideNumber -1;
            high = bwt->L2[c] + bwt->occMajor[bp_h/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_h, c) + bwt->divideNumber -1;
        }
        if(low >= high) break;
    }
    if(low >= high) return;
    int s_size = INITIAL_STACK_SIZE, s_count = 0 ;
    LCPBound *lcpBound = (LCPBound*)xcalloc(s_size, sizeof(LCPBound));
    lcpBound[s_count].low = low; lcpBound[s_count].high = high;
    lcpBound[s_count].depth = depth - first_len; s_count++;
    
    while(s_count > 0)
    {
        LCPBound lcp_elem = lcpBound[--s_count];
        if(lcp_elem.depth < 1) { fprintf(stderr, "[setBWTLCP] depth set error, exit...\n");  exit(1); }
        int64_t bp_l, bp_h;
        bp_l = posSAChangeToPosBWT(bwt, lcp_elem.low, 0);
        bp_h = posSAChangeToPosBWT(bwt, lcp_elem.high, 1);
        int64_t region = bp_h - bp_l + 1;
        if(lcp_elem.flag == 1)
        {
            uint32_t c = bwt_B0(bwt, bp_l);
            low = bwt->L2[c] + bwt->occMajor[bp_l/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_l, c) + bwt->divideNumber -1;
            high = bwt->L2[c] + bwt->occMajor[bp_h/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_h, c) + bwt->divideNumber -1;
            if(high - low + 1 == region)
            {
                if(lcp_elem.depth > 1)
                {
                    if(s_count >= s_size)
                    {
                        s_size <<= 1;
                        lcpBound = (LCPBound*)xrecalloc(lcpBound, s_size * sizeof(LCPBound));
                    }
                    lcpBound[s_count].low = low; lcpBound[s_count].high = high;
                    lcpBound[s_count].depth = lcp_elem.depth - 1; lcpBound[s_count].flag = 1; s_count++;
                } else { // lcp_elem.depth == 1
                    for(uint64_t j = low+1; j <= high ; j++)
                    {
                        __sync_or_and_fetch(lcp + (j>>3), SET_BIT[j & 0X7]);
                        //SET_LCP_FLAG(lcp, low + j);
                    }
                }
            } else {
                //if((lcp_elem.depth > mid_depth && high - low > 2 * min_kmerfreq) || 
                //    (lcp_elem.depth <= mid_depth && high - low >= min_kmerfreq))
				if(high - low >= min_kmerfreq)
                {
                    if(lcp_elem.depth > 1)
                    {
                        if(s_count >= s_size)
                        {
                            s_size <<= 1;
                            lcpBound = (LCPBound*)xrecalloc(lcpBound, s_size * sizeof(LCPBound));
                        }
                        lcpBound[s_count].low = low; lcpBound[s_count].high = high;
                        lcpBound[s_count].depth = lcp_elem.depth - 1; lcpBound[s_count].flag = 0; s_count++;
                    } else { // lcp_elem.depth == 1
                        for(uint64_t j = low+1; j <= high ; j++)
                        {
                            __sync_or_and_fetch(lcp + (j>>3), SET_BIT[j & 0X7]);
                            //#pragma omp atomic
                            //SET_LCP_FLAG(lcp, low + j);
                        }
                    }
                }
                // test region
                if(high > low)region -= (high - low + 1);
                if(region <= min_kmerfreq) continue;

                for(int i = 0; i < DNA_SYMBOLS_NUMBER; i++)
                {
                    if(i == c) continue;
                    low = bwt->L2[i] + bwt->occMajor[bp_l/OCC_INTERVAL_MAJOR64 * 4 + i] + bwt_occ(bwt, bp_l, i) + bwt->divideNumber -1;
                    high = bwt->L2[i] + bwt->occMajor[bp_h/OCC_INTERVAL_MAJOR64 * 4 + i] + bwt_occ(bwt, bp_h, i) + bwt->divideNumber -1;
                    if((lcp_elem.depth > mid_depth && high -low > 2 * min_kmerfreq) ||
                        (lcp_elem.depth <= mid_depth && high - low >= min_kmerfreq))
                    {
                        if(lcp_elem.depth > 1)
                        {
                            if(s_count >= s_size)
                            {
                                s_size <<= 1;
                                lcpBound = (LCPBound*)xrecalloc(lcpBound, s_size * sizeof(LCPBound));
                            }
                            lcpBound[s_count].low = low; lcpBound[s_count].high = high;
                            lcpBound[s_count].depth = lcp_elem.depth - 1; lcpBound[s_count].flag = 0; s_count++;
                        } else { // lcp_elem.depth == 1
                            for(uint64_t j = low +1; j <= high ; j++)
                            {
                                __sync_or_and_fetch(lcp + (j>>3), SET_BIT[j & 0X7]);
                                //#pragma omp atomic
                                //SET_LCP_FLAG(lcp, low + j);
                            }
                        }
                    }
                    // test region
                    if(high > low) region -= (high - low + 1);
                    if(region <= min_kmerfreq) break;
                }
            }
        } else {
            int64_t loc_region = region;
            for(int i = 0; i < DNA_SYMBOLS_NUMBER; i++)
            {
                low = bwt->L2[i] + bwt->occMajor[bp_l/OCC_INTERVAL_MAJOR64 * 4 + i] + bwt_occ(bwt, bp_l, i) + bwt->divideNumber -1;
                high = bwt->L2[i] + bwt->occMajor[bp_h/OCC_INTERVAL_MAJOR64 * 4 + i] + bwt_occ(bwt, bp_h, i) + bwt->divideNumber -1;
                if( (lcp_elem.depth > mid_depth && high - low > 2 * min_kmerfreq) || 
                    ( lcp_elem.depth <= mid_depth &&  high - low >= min_kmerfreq))
                {
                    if(lcp_elem.depth > 1)
                    {
                        if(s_count >= s_size)
                        {
                            s_size <<= 1;
                            lcpBound = (LCPBound*)xrecalloc(lcpBound, s_size * sizeof(LCPBound));
                        }
                        lcpBound[s_count].low = low; lcpBound[s_count].high = high;
                        lcpBound[s_count].depth = lcp_elem.depth - 1; 
                        if(high - low + 1 == region) lcpBound[s_count].flag = 1;  
                        else lcpBound[s_count].flag = 0;
                        s_count++;
                    } else { // lcp_elem.depth == 1
                        for(uint64_t j = low + 1; j <= high; j++)
                        {
                            __sync_or_and_fetch(lcp + (j>>3), SET_BIT[j & 0X7]);
                            //#pragma omp atomic
                            //SET_LCP_FLAG(lcp, low + j);
                        }
                    }
                }
                // test region
                if(high > low) loc_region -= (high - low + 1);
                if(loc_region <= min_kmerfreq) break;
            }
        }
    }
    free(lcpBound);
} */

/*
void setBWTLCP(const bwt_t *bwt, uint8_t *lcp, int64_t start, const int64_t len, const int depth)
{
    // check arguments
    if(len < 1000) { fprintf(stderr, "[setBWTLCP] the value of argument 'len' small than 1000, exit...\n"); exit(1); }
    // skip head same to last partition
    if(start > 0)
    {
        uint32_t c1 = bwt_B0(bwt, start-1);
        while(bwt_B0(bwt, start) == c1) start++;
    }
    for(int64_t i = start; i < start+len; i++)
    {
        int64_t bp_s = i, bp_e, pos_s, pos_e;
        uint32_t c1 = bwt_B0(bwt, i);
        while( i+1 < bwt->seq_len && bwt_B0(bwt, i+1) == c1) { i++;}
        bp_e = i;
        if(bp_s >= bp_e) continue;
        pos_s = bwt->L2[c1] + bwt->occMajor[bp_s/OCC_INTERVAL_MAJOR64 * 4 + c1] + bwt_occ(bwt, bp_s, c1) + bwt->divideNumber - 1;
        pos_e = pos_s + (bp_e - bp_s);
        int s_size = INITIAL_STACK_SIZE, s_count = 0 ;
        LCPBound *lcpBound = (LCPBound*)xcalloc(s_size, sizeof(LCPBound)); 
        lcpBound[s_count].low = pos_s ;     lcpBound[s_count].high = pos_e;
        lcpBound[s_count].depth = depth - 1;    s_count++;

        while(s_count > 0)
        {
            LCPBound lcp_elem = lcpBound[--s_count];
            if(lcp_elem.low >= lcp_elem.high) continue;
            if(lcp_elem.depth > 0)
            {
                int64_t cor_low, cor_high;
                if((cor_low = posSA2PosBWT(bwt, lcp_elem.low)) == -1 || (cor_high = posSA2PosBWT(bwt, lcp_elem.high)) == -1
                    || lcp_elem.high - lcp_elem.low != cor_high - cor_low)
                {
                    for(int64_t j = lcp_elem.low; j < lcp_elem.high; j++)
                    {
                        if((cor_low = posSA2PosBWT(bwt, j)) == -1) continue;
                        int64_t low = j ;
                        while(j+1 <= lcp_elem.high && (cor_high = posSA2PosBWT(bwt, j+1)) != -1) { j++; }
                        if(j > low)
                        {
                            if(s_count >= s_size)
                            {
                                s_size <<= 1;
                                lcpBound = (LCPBound*)realloc(lcpBound, s_size * sizeof(LCPBound));
                            }
                            lcpBound[s_count].low = low; lcpBound[s_count].high = j;
                            lcpBound[s_count].depth = lcp_elem.depth;  s_count++;
                        }
                    }
                    continue;
                }
                for(int64_t j = cor_low; j < cor_high; j++)
                {
                    bp_s = j ;
                    c1 = bwt_B0(bwt, j);
                    while(j+1 <= cor_high && bwt_B0(bwt, j+1) == c1) { j++;}
                    bp_e = j;
                    if(bp_s >= bp_e) continue;
                    pos_s = bwt->L2[c1] + bwt->occMajor[bp_s/OCC_INTERVAL_MAJOR64 * 4 + c1] + bwt_occ(bwt, bp_s, c1) + bwt->divideNumber - 1;
                    pos_e = pos_s + (bp_e - bp_s);
                    if(s_count >= s_size)
                    {
                        s_size <<= 1;
                        lcpBound = (LCPBound*)realloc(lcpBound, s_size * sizeof(LCPBound));
                    }
                    lcpBound[s_count].low = pos_s; lcpBound[s_count].high = pos_e;
                    lcpBound[s_count].depth = lcp_elem.depth - 1;  s_count++;
                }
            } else { // lcp_elem.depth == 0
                for(int j = 1; j <= lcp_elem.high - lcp_elem.low; j++)
                {
                    SET_LCP_FLAG(lcp, lcp_elem.low + j);
                }
            }
        }
        free(lcpBound); 
    }
}*/

/*
void setBWTLCP(const bwt_t *bwt, uint8_t *lcp, const int last_width, const int width , const int64_t start, const int64_t len)
{
    int divideNumber = bwt->divideNumber;
    int ext_len = width - last_width ;
    //int64_t bp, location = start;
    while( location < start + len && (bp = posSA2PosBWT(bwt, location)) == -1) { location++; }
    uint32_t c1 = bwt_B0(bwt, location), c2;
     // skip equal to last char of last partition
    if(start >0)
    {
        int64_t loc = start + offset ;
        c2 = bwt_B0(bwt, start - 1); 
        while(c1 == c2 && GET_OLD_LCP_FLAG_FORWARD(lcp, loc) == flag && GET_OLD_LCP_FLAG_BACKWARD(lcp, loc) == flag && 
            ++offset < len) { loc = start + offset; c1 = bwt_B0(bwt, loc); }
    } 
    if(last_width == 0)
    {
        for(int64_t i = start; i < start + len; i++)
        {
            int64_t bp_s, bp_e, pos_s = i, pos_e ;
            if((bp_s = posSA2PosBWT(bwt, i)) == -1) continue ;
            uint32_t c1 = bwt_B0(bwt, bp_s), c2;
            bp_e = bp_s + 1;
            while(bp_e < bwt->seq_len && bwt_B0(bwt,bp_e) == c1) {bp_e++;}
            bp_e--;
            if(bp_e == posSA2PosBWT(bwt, i + (bp_e - bp_s))) 
            {
                i += (bp_e - bp_s);
            } else {
                bp_e = bp_s;
                i++;
                int64_t swap_bp ;
                while(i < bwt_seq_len + divideNumber && (swap_bp = posSA2PosBWT(bwt, i)) != -1) 
                {
                    if(bwt_B0(bwt, swap_bp) != c1) break;
                    else { bp_e = swap_bp ; i++; }
                }
                i--;
            }
            pos_e = i ;
            if(pos_e <= pos_s) continue;
            int s_size = INITIAL_STACK_SIZE, s_count = 0;
            LCPBound *lcpBound = (LCPBound*)xcalloc(s_size, sizeof(LCPBound));
            lcpBound[s_count].low = pos_s ;
            lcpBound[s_count].high = pos_e ;
            lcpBound[s_count].depth = ext_len ;
            lcpBound[s_count].offset = 0; 
            s_count++;
            while(s_count > 0)
            {
                LCPBound lcp_elem = lcpBound[--s_count];
                if(lcp_elem.low >= lcp_elem.high) continue;
                for(int64_t j = lcp_elem.low; j < lcp_elem.high; j++)
                {
                    int64_t ps = j, pe; // position of SA array
                    if((bp_s = posSA2PosBWT(bwt, j)) == -1) continue ;
                    c1 = bwt_B0(bwt, bp_s), c2;
                    bp_e = bp_s + 1;
                    while(bp_e <= bp_s + (lcp_elem.high -j) && bwt_B0(bwt, bp_e) == c1) {bp_e++;}
                    bp_e--;
                    if((bp_e == posSA2PosBWT(bwt, j + (bp_e - bp_s))))
                    {
                        j += (bp_e - bp_s);
                    } else {
                        bp_e = bp_s;
                        j++;
                        int64_t swap_bp ;
                        while(j < bwt_seq_len + divideNumber && (swap_bp = posSA2PosBWT(bwt, j)) != -1)
                        {
                            if(bwt_B0(bwt, swap_bp) != c1) break;
                            else { bp_e = swap_bp; j++ }
                        }
                        j--;
                    }
                    pe = j ;
                    if(pe <= ps) continue;
                    if(lcp_elem.depth > 0)
                    {
                        if(s_size <= s_count)
                        {
                            s_size <<= 1;
                            lcpBound = (LCPBound*)realloc(lcpBound, s_size * sizeof(LCPBound));
                        }
                        lcpBound[s_count].low = ps; lcpBound[s_count].high = pe;
                        lcpBound[s_count].depth = lcp_elem.depth - 1;
                        lcpBound[s_count].offset = lcp_elem.offset + (ps - lcp_elem.low);
                        //lcpBound[s_count].offset_h = lcp_elem.offset_h + (lcp_elem.high - pe);
                        s_count++;
                    } else {  // lcp_elem.depth == 0
                        int64_t t = pos_s + lcp_elem.offset + (ps - lcp_elem.low);
                        for(int z = 1; z <= pe - ps; z++)
                        {
                            SET_NEW_LCP_FLAG_FORWARD(lcp, ps + z);
                            SET_NEW_LCP_FLAG_BACKWARD(lcp, t + z);
                        }
                    }
                }
            }
            free(lcpBound);
        }
    }else { // last_width >0 
        for(int64_t i = start; i < start + len; i++)
        {
            int64_t pos_s = i, pos_e;
            while(i + 1 < bwt_seq_len + divideNumber && GET_OLD_LCP_FLAG_FORWARD(lcp, i+1) == 1 &&
                GET_OLD_LCP_FLAG_BACKWARD(lcp, i+1) == 1) { i++; }
            pos_e = i ;
            if(pos_e <= pos_s) continue;
            int s_size = INITIAL_STACK_SIZE, s_count = 0;
            LCPBound *lcpBound = (LCPBound*)xcalloc(s_size, sizeof(LCPBound));
            lcpBound[s_count].low = pos_s;
            lcpBound[s_count].high = pos_e;
            lcpBound[s_count].depth = ext_len - 1 ;
            lcpBound[s_count].offset = 0 ;
            s_count++;
            while(s_count > 0)
            {
                LCPBound lcp_elem = lcpBound[--s_count];
                if(lcp_elem.low >= lcp_elem.high) continue;
                if(lcp_elem.depth > 0)
                {
                    for(int64_t j = lcp_elem.low; j < lcp_elem.high; j++)
                    {
                        int64_t ps = j, pe, a, b; // position of SA
                        if(lcp_elem.depth >= last_width) 
                        {
                            while( j < lcp_elem.high && GET_OLD_LCP_FLAG_BACKWARD(lcp, j+1) == 1) { j++; }
                        } else {
                            while( j < lcp_elem.high && GET_OLD_LCP_FLAG_FORWARD(lcp, j+1) == 1) { j++; }
                        }
                        pe = j;
                        a = posSA2PosBWT(bwt, ps); b = posSA2PosBWT(bwt, pe);
                        uint32_t c = bwt_B0(bwt, a);
                        ps = bwt->L2[c] + bwt->occMajor[a/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, a, c) + divideNumber - 1;
                        pe = bwt->L2[c] + bwt->occMajor[b/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, b, c) + divideNumber - 1;
                        if(s_size <= s_count)
                        {
                            s_size <<= 1;
                            lcpBound = (LCPBound*)realloc(lcpBound, s_size * sizeof(LCPBound));
                        }
                        lcpBound[s_count].low = ps; lcpBound[s_count].high = pe;
                        lcpBound[s_count].depth
                    }
                }
            }
        }
    }
    for(int64_t i = start ; i < start + len; i++)
    {
        if(i + 1 < bwt->seq_len + bwt->divideNumber && GET_OLD_LCP_FLAG_FORWARD(lcp, i+1) == flag && 
            GET_OLD_LCP_FLAG_BACKWARD(lcp, i+1) == flag)
        {
            int64_t pos_s = i , pos_e, bp ;
            while(i < start+len && (bp = posSA2PosBWT(bwt, i)) == -1) {i++;}
            uint32_t c1 = bwt_B0(bwt, bp), c2;
            i++;
            while(i < start+len && (bp = posSA2PosBWT(bwt,i)) == -1) {i++; }
            while(i < bwt->seq_len + bwt->divideNumber && GET_OLD_LCP_FLAG_FORWARD(lcp, i) == flag && 
                GET_OLD_LCP_FLAG_BACKWARD(lcp, i) == flag && bwt_B0(bwt, bp) == c1) 
            { 
                i++; 
                while(i < start+len && (bp = posSA2PosBWT(bwt,i)) == -1) {i++; }
            }
            pos_e = --i;
            if(pos_e > pos_s)
            {
                int interleave, offset_h = 0; 
                //int64_t low = bwt->L2[c1] + bwt->occMajor[pos_s/OCC_INTERVAL_MAJOR64 * 4 + c1] + bwt_occ(bwt, pos_s, c1) + divideNumber - 1;
                //int64_t high = bwt->L2[c1] + bwt->occMajor[pos_e/OCC_INTERVAL_MAJOR64 * 4 + c1] + bwt_occ(bwt, pos_e, c1) + divideNumber - 1;
                if((interleave = checkInterLeave(bwt, pos_s, pos_e)) != -1)
                {
                    if(interleave > 1) 
                    {
                        offset_h = pos_e - (pos_s + interleave - 1);
                        pos_e = pos_s + interleave -1 ;
                        //high = low + interleave - 1 ;
                        i =  pos_s + interleave - 1 ;
                    } else {
                        i = pos_s + interleave - 1 ;
                        continue;
                    }
                }
                int s_size = INITIAL_STACK_SIZE, s_count = 0 ;
                LCPBound *lcpBound = (LCPBound*)xcalloc(s_size, sizeof(LCPBound)); 
                lcpBound[s_count].low = pos_s ;
                lcpBound[s_count].high = pos_e;
                lcpBound[s_count].depth = ext_len - 1;
                lcpBound[s_count].offset_l = 0; lcpBound[s_count].offset_h = offset_h;
                s_count++;
                while(s_count > 0)
                {
                    LCPBound lcp_elem = lcpBound[--s_count];
                    //if(lcp_elem.depth <= 0) continue;
                    // process encounter sentinelPosition
                    if(lcp_elem.low >= lcp_elem.high) continue;
                    if(lcp_elem.depth > 0)
                    {
                        if(last_width > 0 && lcp_elem.depth >= last_width)
                        {
                            if
                        }
                        int64_t cor_low, cor_high ;
                        while((cor_low = posSA2PosBWT(bwt, lcp_elem.low)) == -1) {lcp_elem.offset_l++; lcp_elem.low++; }
                        while((cor_high = posSA2PosBWT(bwt, lcp_elem.high)) == -1) { lcp_elem.offset_h++; lcp_elem.high--;}
                        for(int64_t j = lcp_elem.low; j < lcp_elem.high; j++)
                        {
                            int64_t begin = j,bp_s, bp_e ;
                            if(last_width > 0 && lcp_elem.depth >= last_width)
                            {
                                while(GET_OLD_LCP_FLAG_BACKWARD(lcp,j+1) == flag) {j++;}
                                if(j > begin)
                                {
                                    while((bp_s = posSA2PosBWT(bwt, begin)))
                                }
                            }
                            c1 = bwt_B0(bwt, j); j++; c2 = bwt_B0(bwt, j); 
                            if(c1 != c2) { j--; continue; } 
                            while(c1 == c2 && j <=cor_high && j < bwt->seq_len) {j++;  c2 = bwt_B0(bwt, j);}
                            j--;
                            if(j > begin)
                            {
                                int64_t a, b ;
                                int inter ;
                                a = bwt->L2[c1] + bwt->occMajor[begin/OCC_INTERVAL_MAJOR64 * 4 + c1] + bwt_occ(bwt, begin, c1) + divideNumber - 1;
                                b = bwt->L2[c1] + bwt->occMajor[j/OCC_INTERVAL_MAJOR64 * 4 + c1] + bwt_occ(bwt, j, c1) + divideNumber - 1;
                                if((inter = checkInterLeave(bwt, a, b)) != -1)
                                {
                                    if(inter > 1)
                                    {
                                        b = a + inter - 1;
                                        j = begin + inter - 1;
                                    } else {
                                        j = begin + inter -1;
                                        continue;
                                    }
                                }
                                if(s_size <= s_count)
                                {
                                    s_size <<= 1 ;
                                    lcpBound = (LCPBound*)realloc(lcpBound, s_size * sizeof(LCPBound));
                                }
                                lcpBound[s_count].low = a ;  lcpBound[s_count].high = b;
                                lcpBound[s_count].depth = lcp_elem.depth - 1;
                                lcpBound[s_count].offset_l = lcp_elem.offset_l + (begin - cor_low);
                                lcpBound[s_count].offset_h = lcp_elem.offset_h + (cor_high - j);
                                s_count++;
                            } 
                        }
                    } else {
                        int64_t t = pos_s + lcp_elem.offset_l;
                        for(int j = 1; j <= lcp_elem.high - lcp_elem.low; j++)
                        {
                            SET_NEW_LCP_FLAG_FORWARD(lcp, lcp_elem.low + j);
                            if(GET_OLD_LCP_FLAG_BACKWARD(lcp, lcp_elem.low + j) == flag) SET_NEW_LCP_FLAG_BACKWARD(lcp, t + j);
                        }
                    }
                }
                // free and clean work
                free(lcpBound);
            }
        } 
    }
}

uint8_t *compressLcp(const bwt_t *bwt, uint8_t *lcp)
{
    uint8_t *p ;
    for(int64_t i = 0; i < bwt->seq_len + bwt->divideNumber; i++)
    {
        GET_OLD_LCP_FLAG_FORWARD(lcp, i) == 1 ? COMPRESS_SET_LCP(lcp, i) : COMPRESS_RESET_LCP(lcp, i) ;
    }
    p = (uint8_t*)realloc(lcp, ((bwt->seq_len + bwt->divideNumber + BITS_IN_BYTE -1)/BITS_IN_BYTE) * sizeof(uint8_t));
    return p ;
} */

void saveKmerSA(const uint8_t *lcp , bwt_t *bwt , uint64_t start , const uint64_t end ,const lib_info *libIndex_info, FILE *kmerSAfp, const int min_kmerfreq, const int lcp_fixed_len)
{
    KmerFreq *kmerSABuf ;
    int64_t totalCount = 0;
    //const uint8_t ArraySet[8] = { 0x80, 0x40, 0x20, 0x10, 0x08 , 0x04 , 0x02, 0x01 } ;
    int kmerSAIndex = 0 , bufSize = 10000000 ;
    int flagCount = 0 , zeroCount = 1 ;
    uint64_t lastRegionISAStart  = 0, lastRegionISAEnd = 0 , i ;
    kmerSABuf = (KmerFreq*)xcalloc(bufSize , sizeof(KmerFreq));
    // skip the head 
    while(start > bwt->divideNumber && GET_LCP_FLAG(lcp, start) > 0) { start++;}
    for( i = start ; i < end ; i++)
    {
        if(GET_LCP_FLAG(lcp, i) > 0)
        {
            flagCount++ ;
            if(zeroCount > 0)
            {
                lastRegionISAStart = i ;
                zeroCount = 0 ;
            }
        } else {
            zeroCount++ ;
            if(flagCount >= min_kmerfreq  - 1)
            {
                int64_t a;
                lastRegionISAEnd = i -1 ;
				while(lastRegionISAStart <= lastRegionISAEnd)
				{
					if(lastRegionISAEnd / SA_INTERVAL * SA_INTERVAL > lastRegionISAStart || lastRegionISAStart%SA_INTERVAL == 0)
					{
						a = bwt->sa[lastRegionISAEnd / SA_INTERVAL] ;
					} else {
						a = calculateSA(bwt , lastRegionISAStart );
					}
					ReadLocation rl = readBoundaryLookup(libIndex_info, a);
					if(checkSA(rl, a, lcp_fixed_len))
					{
						kmerSABuf[kmerSAIndex].sa = a ;
						kmerSABuf[kmerSAIndex++].freq = flagCount + 1 < 255 ? flagCount+1 : 255;
					} else { 
						//printf("[saveKmerSA] rl.bound: %ld, rl.length[0]: %d, a: %ld\n", rl.bound, rl.length[0], a);
					   	lastRegionISAStart++ ;	
						continue;
					}
					totalCount++;
					if(kmerSAIndex >= bufSize)
					{
						fwrite(kmerSABuf, sizeof(KmerFreq) , kmerSAIndex, kmerSAfp);
						kmerSAIndex = 0 ;
					}
					break ;
				}
            }
            flagCount = 0 ;
        }
    }
    // 
    while(i < bwt->seq_len + bwt->divideNumber && flagCount > 0 && GET_LCP_FLAG(lcp, i) > 0)
    {
        i++; flagCount++;
    }
    if(flagCount >= min_kmerfreq -1)
    {
        int64_t a;
        lastRegionISAEnd = i -1 ;
		while(lastRegionISAStart <= lastRegionISAEnd)
		{
			if(lastRegionISAEnd / SA_INTERVAL * SA_INTERVAL > lastRegionISAStart || lastRegionISAStart%SA_INTERVAL == 0)
			{
				a = bwt->sa[lastRegionISAEnd / SA_INTERVAL] ;
			} else {
				a = calculateSA(bwt , lastRegionISAStart );
			}
			ReadLocation rl = readBoundaryLookup(libIndex_info, a);
			if(checkSA(rl, a, lcp_fixed_len))
			{
				kmerSABuf[kmerSAIndex].sa = a ;
				kmerSABuf[kmerSAIndex++].freq = flagCount + 1 < 255 ? flagCount+1 : 255;
			} else { 
				//printf("[saveKmerSA] rl.bound: %ld, rl.length[0]: %d, a: %ld\n", rl.bound, rl.length[0], a);
				lastRegionISAStart++ ;	
				continue ;
			}
			totalCount++;
			break ;
		}
    }
    if(kmerSAIndex > 0)
    {
        fwrite(kmerSABuf, sizeof(KmerFreq) , kmerSAIndex, kmerSAfp);
    }
    fprintf(stderr, "[saveKmerSA] totalCount : %ld\n", totalCount);
	//fprintf(stderr, "no_cyc : %ld\n", no_cyc);
    // clean and free work
    free(kmerSABuf);
}

/*
void saveKmerSAAndLongZeroSA(const uint8_t *lcp , const bwt_t *bwt, const uint64_t start, const uint64_t end, FILE *kmerSAfp, FILE *longZeroSAfp)
{
    KmerFreq *kmerSABuf ;
    LongZeroSA *longZeroSABuf ;
    uint8_t setArray[8] = { 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 } ;
    int kmerSAIndex = 0, longZeroSAIndex = 0 , bufSize = 10000000 ;
    int flagCount = 0 , zeroCount = 0;
    uint64_t lastRegionISAStart, lastRegionISAEnd ;
    uint8_t lastFlag = 0 ;
    kmerSABuf = (KmerFreq*)calloc(bufSize, sizeof(KmerFreq));
    longZeroSABuf = (LongZeroSA*)calloc(bufSize,sizeof(LongZeroSA));
    for(uint64_t i = start ; i < end ; )
    {
        if(lcp[i/BITS_IN_BYTE] & setArray[i%BITS_IN_BYTE] > 0) 
        {
            flagCount++ ; 
            if(zeroCount >= MIN_ZERO_GAP_NUM )
            {
                lastRegionISAEnd = i - 1 ;
                longZeroSABuf[longZeroSAIndex].ISA[0] = lastRegionISAStart + 2;
                longZeroSABuf[longZeroSAIndex].ISA[1] = lastRegionISAEnd - 2 ;
                longZeroSABuf[longZeroSAIndex].SA[0] = calculateSA(bwt, longZeroSABuf[longZeroSAIndex.ISA[0]]);
                longZeroSABuf[longZeroSAIndex].SA[1] = calculateSA(bwt, longZeroSABuf[longZeroSAIndex.ISA[1]]);
                longZeroSAIndex++ ;
                zeroCount = 0 ;
                lastRegionISAStart = i ;
                if(longZeroSAIndex >= bufSize)
                {
                    #pragma omp atomic
                    fwrite(longZeroSABuf, sizeof(LongZeroSA),longZeroSAIndex,longZeroSAfp);
                    longZeroSAIndex = 0 ;
                }
            } else if(zeroCount >0) {
                zeroCount = 0 ;
                lastRegionISAStart = i ;
            }
        } else {
            zeroCount++ ;
            if(flagCount >= MIN_KMERFREQ - 1)
            {
                lastRegionISAEnd = i - 1 ;
                if(lastRegionISAEnd/SA_INTERVAL * SA_INTERVAL >= lastRegionISAStart)
                {
                    kmerSABuf[kmerSAIndex].sa = bwt->sa[lastRegionISAEnd/ SA_INTERVAL] ;
                } else {
                    kmerSABuf[kmerSAIndex].sa = calculateSA(bwt, lastRegionISAStart);
                }
                kmerSABuf[kmerSAIndex].freq = flagCount + 1 ;
                kmerSAIndex++ ;
                flagCount = 0 ;
                lastRegionISAStart = i ;
                if(kmerSAIndex >= bufSize)
                {
                    #pragma omp atomic
                    fwrite(kmerSABuf, sizeof(KmerFreq), kmerSAIndex, kmerSAfp);
                    kmerSAIndex = 0 ;
                }
            } else if(flagCount > 0) {
                flagCount = 0 ;
                lastRegionISAStart = i ;
            }
        }
    }
    if(zeroCount >= MIN_ZERO_GAP_NUM )
    {
        lastRegionISAEnd = i - 1 ;
        longZeroSABuf[longZeroSAIndex].ISA[0] = lastRegionISAStart + 2;
        longZeroSABuf[longZeroSAIndex].ISA[1] = lastRegionISAEnd - 2 ;
        longZeroSABuf[longZeroSAIndex].SA[0] = calculateSA(bwt, longZeroSABuf[longZeroSAIndex.ISA[0]]);
        longZeroSABuf[longZeroSAIndex].SA[1] = calculateSA(bwt, longZeroSABuf[longZeroSAIndex.ISA[1]]);
        longZeroSAIndex++ ;
    } else if(flagCount >= MIN_KMERFREQ){
        lastRegionISAEnd = i - 1 ;
        if(lastRegionISAEnd/SA_INTERVAL * SA_INTERVAL >= lastRegionISAStart)
        {
            kmerSABuf[kmerSAIndex].sa = bwt->sa[lastRegionISAEnd/ SA_INTERVAL] ;
        } else {
            kmerSABuf[kmerSAIndex].sa = calculate(bwt, lastRegionISAStart);
        }
        kmerSABuf[kmerSAIndex].freq = flagCount + 1 ;
        kmerSAIndex++ ;
    }

    if(longZeroSAIndex >= 0)
    {
        #pragma omp atomic
        fwrite(longZeroSABuf, sizeof(LongZeroSA),longZeroSAIndex,longZeroSAfp);
    }
    if(kmerSAIndex >= 0 )
    {
        #pragma omp atomic
        fwrite(kmerSABuf, sizeof(KmerFreq), kmerSAIndex, kmerSAfp);
    }

    free(longZeroSABuf); free(kmerSABuf);
}
*/

void getKmer(uint64_t *kmer, const uint8_t *pac, const uint64_t start, const int len)
{
	for(int64_t i = start, j = 0 ; i < start + len ; i++, j++)
	{
		kmer[j>>5] <<= BIT_PER_CHAR ;      // j/32 ==  j >> 5
		kmer[j>>5] |= (pac[i>>2] >> (((~i)&0x3)<<1));   // i/CHAR_PER_BYTE = i >>2
	}
}

/*
void getKmer(uint64_t *kmer, const uint8_t *pac, const uint64_t start, const int len )
{
    int l = (len + 32 -1) / 32 ;
    int index = 0 ;
    int remain ;
    kmer[l - 1] = 0 ;
    kmer[index/32] = ((pac[start/CHAR_PER_BYTE] << ((start%CHAR_PER_BYTE)<<1)) >> ((start%CHAR_PER_BYTE)<<1));
    for(index = CHAR_PER_BYTE - start % CHAR_PER_BYTE ; index + CHAR_PER_BYTE < len ; index += CHAR_PER_BYTE )
    {
        if((index + CHAR_PER_BYTE - 1)/32 > index/32)
        {
            remain = (index + CHAR_PER_BYTE - 1)/32 * 32 - index ;
            kmer[index/32] <<= (remain<<1) ;
            kmer[index/32] |= (pac[(start + index)/CHAR_PER_BYTE] >> ((CHAR_PER_BYTE - remain)<<1));
            kmer[index/32 + 1] = ((pac[(start + index)/CHAR_PER_BYTE] << (remain<<1)) >> (remain<<1));
        } else {
            kmer[index/32] <<= BITS_IN_BYTE;
            kmer[index/32] |= pac[(start + index)/CHAR_PER_BYTE] ;
        }
    }
    remain = len -index ;
	if((index + remain -1)/32 > index/32)
	{
		int r = (index + remain - 1)/32 * 32 - index ;
		kmer[index/32] <<= (r<<1) ;
		kmer[index/32] |= (pac[(start + index)/CHAR_PER_BYTE] >> ((CHAR_PER_BYTE - remain)<<1));
		kmer[index/32 + 1] = ((pac[(start + index)/CHAR_PER_BYTE] << (remain<<1)) >> (remain<<1));
	}
    kmer[index/32] <<= (remain<<1);
    kmer[index/32] |= (pac[(start + index)/CHAR_PER_BYTE] >> ((CHAR_PER_BYTE - remain)<<1));
} */

/*
void setLcpFromLongZeroSA(LongZeroSA *longZeroSABuf, const uint64_t start,const uint64_t end, uint8_t *pac, uint8_t *lcp)
{
    int kmerLen = (KMER + 1 + 32 - 1 )/32 ;
    uint64_t kmer1[kmerLen], kmer2[kmerLen];
    
    for(int64_t i = start ; i < end ; i++)
    {
        int  breakFlag = 0 ;
        // get KMER + 1 length sequence
        getKmer(kmer1,pac , longZeroSABuf[i].SA[0], KMER + 1);
        getKmer(kmer2,pac, longZeroSABuf[i].SA[1] , KMER + 1);
        for(int j = 0 ; j < kmerLen ; j++)
        {
            if(kmer1[j] != kmer2[j]) {
                breakFlag = 1 ;
                break ;
            }
        }
        if(breakFalg == 0) 
            setLcp(lcp, longZeroSABuf[i].ISA[0] + 1 , longZeroSABuf[i].ISA[1]);
    }
}
*/

int compare64(const void *p, const void *q)
{
    KmerFreq *p1 = (KmerFreq*)p, *q1 = (KmerFreq*)q ;
    if(p1->sa < q1->sa) return -1 ;
    else if(p1->sa > q1->sa) return 1 ;
    else return 0 ;
}

void mergeSort(const KmerFreq *src1, const  uint64_t src1Len, const KmerFreq *src2, const uint64_t src2Len, KmerFreq *dest)
{
    uint64_t destIndex = 0 , src1Index , src2Index = 0;
    for(src1Index = 0 ; src1Index < src1Len ; src1Index++)
    {
        while( src2Index < src2Len && src2[src2Index].sa <= src1[src1Index].sa)
        {
            dest[destIndex++] = src2[src2Index++];
        }
        dest[destIndex++] = src1[src1Index];
    }
    while(src2Index < src2Len)
    {
        dest[destIndex++] = src2[src2Index++];
    }
}

void extractKmerFromPac(KmerFreq *kmerSABuf, uint64_t len ,FILE *kmerfp, uint8_t *pac, omp_lock_t *lock_kmerfp)
{
    int bufSize = 5000000, bufIndex = 0 ;
    KmerInfo *kmerInfoBuf = (KmerInfo*)xcalloc(bufSize, sizeof(KmerInfo));
    for(uint64_t i = 0 ; i < len ; i++)
    {
        getKmer(kmerInfoBuf[bufIndex].kmer, pac , kmerSABuf[i].sa, KMER + 1);
        //getRevKmer(kmerInfoBuf[bufIndex].kmer, kmerInfoBuf[bufIndex].rkmer, KMER + 1);
		char *s1 = transformKmerToString(kmerInfoBuf[bufIndex].kmer, 0, KMER+1);
		//char *s2 = transformKmerToString(kmerInfoBuf[bufIndex].rkmer, 0, KMER+1);
		fprintf(stdout, "%s\n", s1);
		free(s1); 
        kmerInfoBuf[bufIndex].freq = kmerSABuf[i].freq ;
        bufIndex++ ;
        if(bufIndex >= bufSize)
        {
            omp_set_lock(lock_kmerfp);
            fwrite(kmerInfoBuf, sizeof(KmerInfo), bufIndex , kmerfp);
            omp_unset_lock(lock_kmerfp);
            bufIndex = 0 ;
        }
    }
    if(bufIndex > 0)
    {
        omp_set_lock(lock_kmerfp);
        fwrite(kmerInfoBuf, sizeof(KmerInfo), bufIndex, kmerfp);
        omp_unset_lock(lock_kmerfp);
    }

    free(kmerInfoBuf);
}

void calculateKmerFreqAndWriteToFile(uint8_t *lcp, FILE *kmerFreqfp, uint64_t lcpLen)
{
    uint64_t *kmerFreq = (uint64_t*)calloc(256, sizeof(uint64_t));
    //uint8_t setArray[8] = { 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 } ;
    int zeroCount = 0 , flagCount = 0 ;
    // calculate kmerFreq
    for(uint64_t i = 0 ; i < lcpLen ; i++)
    {
        if(lcp[i] == 0xff)
        {
            if(zeroCount >1) { kmerFreq[1] += zeroCount - 1 ; }
            zeroCount = 0;
            flagCount += 8 ;
        } else if(lcp[i] == 0) {
            if(flagCount > 0) flagCount < 255 ? kmerFreq[flagCount + 1]++ : kmerFreq[255]++ ;
            flagCount = 0 ;
            zeroCount += 8 ;
        } else {
            for(int j = 0 ; j < 8 ; j++)
            {
				int64_t p = (i<<3) + j ;
                if(GET_LCP_FLAG(lcp, p) > 0) 
                {
                    if(zeroCount > 1) kmerFreq[1] += zeroCount - 1 ;
                    zeroCount = 0 ; flagCount++ ;
                } else {
                    if(flagCount > 0) flagCount < 255 ? kmerFreq[flagCount + 1]++ : kmerFreq[255]++ ;
                    flagCount = 0 ; zeroCount++ ;
                }
            }
        }
    }
    if(zeroCount >1) kmerFreq[1] += zeroCount - 1 ;
    else if(flagCount >0) flagCount <255 ? kmerFreq[flagCount + 1]++ : kmerFreq[255]++ ;
    
    // writeToFile
    int64_t total = 0;
    for(int i = 1 ; i <256 ; i++)
    {
        fprintf(kmerFreqfp, "%lu\n", kmerFreq[i]/2); // calculate the reversed strand kmer frequency
        total += kmerFreq[i]/2;
    }
    fprintf(stderr, "[calculateKmerFreqAndWriteToFile] total kmerfreq: %ld\n", total);
    free(kmerFreq);
}

static inline void updataLCP(uint8_t *lcp, const int64_t len)
{
    for(int64_t i = 0; i < len; i++)
    {
        lcp[i] <<= 2 ;
        lcp[i] &= 0XCC ; // equal to 11001100
    }
}

static void appendToSAFile(char *dest_fn, char *src_fn)
{
    FILE *dest_fp = xopen(dest_fn, "ab");
    FILE *src_fp = xopen(src_fn, "rb");
    int buf_size = 10000000 ;
    KmerFreq *buf = (KmerFreq*)xcalloc(buf_size, sizeof(KmerFreq));
    while(!feof(src_fp))
    {
        int c = fread(buf, sizeof(KmerFreq), buf_size, src_fp);
        if(c < buf_size && ferror(src_fp))
        {
            fprintf(stderr, "[appendToSAFile] fread occur error, exit..\n" );
            exit(1);
        }
        fwrite(buf, sizeof(KmerFreq), c, dest_fp);
    }

    // free and clean work
    free(buf);
    fclose(src_fp);
    fclose(dest_fp);
}

//Construct BWT for the packed sequence
void generateBWT(const Arguments *arguments)
{
    const char *prefix = arguments->prefix;
    const int NCPU = arguments->NCPU, maxMem = arguments->maxMem ;
    char pac_name[PATH_LEN] ,bwt_name[PATH_LEN], sa_name[PATH_LEN];
    int64_t seqLen ;
    BWTInc64 *bwtInc64 ;
    clock_t t, index_t = clock();
    time_t timeval ;
    lib_info *libIndex_info ;
    char filename[PATH_LEN];
    uint8_t *lcp ;
    char kmerSAF[PATH_LEN] ;
/*  
#ifdef DEBUG
    if(arguments == NULL)
    {
        fprintf(stderr, "arguments is NULL, program exit....\n");
        exit(1);
    } else {
        fprintf(stderr, "prefix: %s\n", arguments.prefix);
    }
#endif
*/
    // restore lib_info 
    strcpy(filename, prefix); strcat(filename, ".ann");
    libIndex_info = lib_info_restore(filename);

    strcpy(pac_name,prefix); strcat(pac_name, ".pac.gz");
    strcpy(bwt_name,prefix); strcat(bwt_name, ".bwt.gz");
    strcpy(sa_name, prefix); strcat(sa_name, ".sa.gz");
	seqLen = libIndex_info->len_pac;
    //seqLen = pac_seq_len(pac_name);
    // check memory available
    if(maxMem < ((seqLen / pow(2 , 30)) * 5 / 8)) // maxMem >= seqLen * 5/8
    {
        fprintf(stderr, "[generateBWT]maxMem: %d < seqLen: %ld, max memory less than the length of package file and exit....\n ",maxMem, seqLen);
        exit(3);
    }
    
    time(&timeval); fprintf(arguments->logfp, "time:\t%s\n", ctime(&timeval));
    
	// Multi-threads BWT construct
	{
		gzFile packedFile ;
		//uint64_t packedFileLen, totalTextLength, processedTextLength ;
		int64_t partLen ;
		// check CPU
		if(NCPU < 2)
		{
			fprintf(stderr, "[generateBWT] not enough CPU for indexing, CPU number must bigger than 2, exit...\n");
			exit(3);
		}

		fprintf(arguments->logfp, "[generateBWT] Begin to construct BWT by multi threads.....\n");
		packedFile = xzopen(pac_name, "r");
		// divide .pac.gz file and parallel construct BWT
		{
			uint8_t *textBuffer ;
			int64_t textBufferLen ;
			int64_t textLocation ;
			int divideCount ;
			// set maximum partition length < 4G(pow(2, 32))
			{
				int i = 2 ;
				while((seqLen + i * NCPU -1) / (i * NCPU) > (pow(2,32) - OCC_INTERVAL_MAJOR - 4 * libIndex_info->max_rd_len)) { i++; }
				
				if(i > 10) 
				{
					fprintf(stderr, "[generateBWT] The number of threads for constructing BWT is too much, Please check! you need use more CPU resource for indexing, program exit...\n ");
					exit(2);
				}
				partLen = ((seqLen + i * NCPU -1) / (i * NCPU)) + libIndex_info->max_rd_len;
			}
				//partLen = pow(2, 32) - CHAR_PER_BYTE ;
			
			textBuffer = (uint8_t*)xcalloc((seqLen + CHAR_PER_BYTE -1)/ CHAR_PER_BYTE, sizeof(uint8_t));
			textBufferLen = seqLen ;
			gzread(packedFile, textBuffer, sizeof(uint8_t) *((seqLen + CHAR_PER_BYTE - 1) / CHAR_PER_BYTE));
			int divideNumber = (seqLen + partLen - 1) / partLen ;
			BWTInc *bwtInc[divideNumber];
			uint8_t finishedFlag[divideNumber];
			memset(finishedFlag, 0 , divideNumber);
			SaIndexRange saIndexRange[divideNumber];
			int i = 0;
			textLocation = textBufferLen ;
			while(textLocation > 0 )
			{
				char fileName[PATH_LEN];
				int64_t textBufferStart;
				if(i >= divideNumber) {fprintf(stderr, "[genarateBWT]divideNumber bigger than seted...exit\n"); exit(2);}
				sprintf(fileName, "%stextBufferDivide_%d.tmp", prefix, i);
				if(textBufferLen - partLen <= 0) textBufferStart = 0 ;
				else {
					// balance the divided package size 
					if((textBufferLen - partLen) < partLen) textBufferStart = textBufferLen/2 ; 
					else textBufferStart = textBufferLen - partLen ;
				}
				ReadLocation rl = readBoundaryLookup(libIndex_info, textBufferStart);
				textLocation = rl.bound ; // pac.gz store "+" and "-" strand 
				dumpAndWriteStorage(textBuffer, textLocation , textBufferLen , fileName);
				saIndexRange[i].startSaIndex = textLocation ; saIndexRange[i].endSaIndex = textBufferLen ;
				textBufferLen = textLocation;
				i++;
			}
			gzclose(packedFile);
			free(textBuffer);
			divideNumber = i;
			// multi-threads construct BWT 
			omp_set_num_threads(NCPU);
			t = clock();
			
			#pragma omp parallel for schedule(dynamic)
			for(int i = 0 ; i < NCPU ; i++)
			{
				char fileName[PATH_LEN];
				sprintf(fileName, "%stextBufferDivide_%d.tmp", prefix, i);
#ifdef DEBUG
				/*
				// check divide file
				{
					FILE *fp, *tfp ;
					char filename[PATH_LEN];
					int64_t loc, lastLoc = libIndex_info->len_pac * 2 - 1 ;
					int64_t seq_len, pac_len ;
					uint8_t *buf ;
					strcpy(filename, arguments->prefix);    strcat(filename, ".pac");
					seq_len = pac_seq_len(filename);
					pac_len = (seq_len + CHAR_PER_BYTE -1)/CHAR_PER_BYTE ;
					fp = xopen(filename, "r");
					tfp = xopen("reverPack.fa", "w");
					buf = (uint8_t*)xcalloc(pac_len, sizeof(uint8_t));
					fread(buf, sizeof(uint8_t), pac_len, fp);
					for(int i = 0 ; i < 5 ; i++)
					{
						loc = readBoundaryLookUp(libIndex_info, lastLoc, BACKWARD);
						transformPacToString(buf, loc, lastLoc+ 1, tfp);
						transformPacToString(buf, libIndex_info->len_pac*2 - lastLoc -1, libIndex_info->len_pac*2 - loc, tfp );
						lastLoc = loc - 1 ;
					}
					fclose(fp); fclose(tfp);
					free(buf);
					
				} */
#endif
				bwtInc[i] = BWTIncConstructFromPacked(fileName, 2.5,10000000,10000000);
				if((bwtInc[i])->bwt->saIndexRange == NULL)
				{
					(bwtInc[i])->bwt->saIndexRange = (SaIndexRange*)xcalloc(1, sizeof(SaIndexRange));
					(bwtInc[i])->bwt->saIndexRange->startSaIndex = saIndexRange[i].startSaIndex ;
					(bwtInc[i])->bwt->saIndexRange->endSaIndex = saIndexRange[i].endSaIndex ;
				}else {
					fprintf(stderr,"bwtInc[%d].bwt->saIndexRange is not NULL...exit\n", i);
					exit(2);
				}
				finishedFlag[i] = 1 ;
				fprintf(arguments->logfp, "finished construct file:%s\n", fileName);
				if( i == 0 && remove(fileName) == -1) fprintf(stderr, "remove file: %s error\n", fileName);
			}
			divideCount = NCPU;
			fprintf(arguments->logfp, "[generateBWT]construct %d files' BWT used %.2f CPU secs\n", NCPU, (float)(clock() - t) / CLOCKS_PER_SEC);
			time(&timeval); fprintf(arguments->logfp, "time:\t%s\n", ctime(&timeval));
			
			omp_set_nested(2) ;
			#pragma omp parallel sections 
			{
				#pragma omp section
				{
					omp_set_num_threads(NCPU/2);
					#pragma omp parallel for schedule(dynamic)
					for(int j = divideCount ; j < divideNumber ; j++)
					{
						char fileName[PATH_LEN];
						sprintf(fileName, "%stextBufferDivide_%d.tmp", prefix, j );
						bwtInc[j] = BWTIncConstructFromPacked(fileName, 2.5,10000000,10000000);
						if((bwtInc[j])->bwt->saIndexRange == NULL)
						{
							(bwtInc[j])->bwt->saIndexRange = (SaIndexRange*)xcalloc(1, sizeof(SaIndexRange));
							(bwtInc[j])->bwt->saIndexRange->startSaIndex = saIndexRange[j].startSaIndex ;
							(bwtInc[j])->bwt->saIndexRange->endSaIndex = saIndexRange[j].endSaIndex ;
						}else {
							fprintf(stderr,"bwtInc[%d].bwt->saIndexRange is not NULL...exit\n", j);
							exit(2);
						}
						finishedFlag[j] = 1 ;
						fprintf(arguments->logfp, "finished construct file:%s\n", fileName);
					}
					time(&timeval); fprintf(arguments->logfp, "finished construct divide files\ntime:\t%s\n", ctime(&timeval)); fflush(arguments->logfp);
				}
				#pragma omp section
				{
					// merge BWT
					bwtInc64 = BWTInc2BWTInc64(bwtInc[0]);
					bwtInc64->bwt->isDivide = 1 ;
					bwtInc64->bwt->divideNumber = 1 ;
					// set sentinel position of the BWT code 
					bwtInc64->bwt->sentinelPosition = (uint64_t*)xcalloc(divideNumber, sizeof(uint64_t));
					bwtInc64->bwt->sentinelSA = (uint64_t*)xcalloc(divideNumber, sizeof(uint64_t));
					bwtInc64->bwt->sentinelPosition[0] = bwtInc64->bwt->inverseSa0 ;
					bwtInc64->bwt->sentinelSA[0] = bwtInc[0]->bwt->saIndexRange->startSaIndex ;
					bwtInc64->bwt->divideBoundary = (DivideBoundary*)xcalloc(divideNumber, sizeof(DivideBoundary));
					bwtInc64->bwt->divideBoundary[0].textLocation = (bwtInc[0])->bwt->saIndexRange->endSaIndex ;
					bwtInc64->bwt->divideBoundary[0].bwtLocation =  0 ; // $ is the small char and the last char in the string is the first symbol in bwt code
					BWTIncFree(bwtInc[0]);
					for(int j = 1 ; j < divideNumber ; j++)
					{
						char fileName[PATH_LEN];
						sprintf(fileName, "%stextBufferDivide_%d.tmp", prefix, j);
						// set number of threads
						int num_threads, finished_number = 0 ;
						for(int z = 0; z < divideNumber; z++)
						{
							if(finishedFlag[z] == 1) finished_number++;
						}
						if(divideNumber - finished_number >= NCPU/2) num_threads = NCPU - NCPU/2;
						else num_threads =  NCPU - (divideNumber - finished_number);
						while(finishedFlag[j] != 1)
						{
							usleep(1000 * 10);
						}
						BWTIncSetMergeSize64(bwtInc64, bwtInc[j]); 
						bwtInc64 = mergeBWT64(bwtInc64,bwtInc[j],j, fileName, num_threads);
						bwtInc64->bwt->divideNumber++ ;
/* #ifdef DEBUG            
						if(j == divideNumber - 1)
						{
							char fn[PATH_LEN];
							int part = 1, len = 1000;
							sprintf(fn, "%stextBufferDivide_%d.tmp", prefix, part);
							FILE *fp = xopen(fn, "r");
							uint8_t str[len];
							fread(str, sizeof(uint8_t), len, fp);
							int64_t rIndex = part ;
							for(int64_t z = bwtInc[part]->bwt->textLength - 1; z >= 0; z--)
							{
								int y;
								for(y = 0; y < bwtInc64->bwt->divideNumber; y++)
								{
									if(rIndex < bwtInc64->bwt->sentinelPosition[y]) break;
									else if(rIndex == bwtInc64->bwt->sentinelPosition[y]) {
										fprintf(stderr, "[generateBWT]rIndex == bwtInc64->bwt->sentinelPosition\n");
										exit(1);
										break;
									}
								}
								int64_t bp = rIndex - y;
								uint32_t c = (bwtInc64->bwt->bwtCode[bp>>4] >> bwtInc64->packedShift[bp & 0xf]) & 0x3 ;
								rIndex = bwtInc64->bwt->cumulativeFreq[c] + BWTOccValue64(bwtInc64->bwt,rIndex, c) + bwtInc64->bwt->divideNumber;
								if(z < 3000 && GET_CHAR_FROM_PACKED(str, z) != c)
								{
									fprintf(stderr, "[generateBWT] position: %ld recover not equal to initial.. %u\t%u\n", z, GET_CHAR_FROM_PACKED(str, z), c);
									exit(1);
								} //else if (z < 3000){ fprintf(stderr, "ok "); }
								const char ntc[] = {'C', 'A', 'T', 'G'};
								if(z < 1000) fprintf(stderr, "%c", ntc[c]);
								if(z < 1000 && z %100 == 0) fprintf(stderr, "\n");
							}
							FILE *tmpfp = xopen("tmpfile", "w");
							fwrite(transformPacToString(str, 0, 1000), sizeof(uint8_t), 1000, tmpfp);
							fclose(tmpfp);
							fclose(fp);
						}
#endif */
						BWTIncFree(bwtInc[j]);
						// remove divide tmp file 
						if(remove(fileName) == -1) fprintf(stderr, "remove file: %s error\n", fileName);
						time(&timeval); fprintf(arguments->logfp, "merged divided BWT: %d\ntime:\t%s\n", j, ctime(&timeval));
						#ifdef DEBUG
						fprintf(stderr, "merged divided BWT: %d\n",j);
						#endif
					}
				}
			}
			time(&timeval); fprintf(arguments->logfp, "merged all the divided BWT\ntime:\t%s\n", ctime(&timeval)); fflush(arguments->logfp);
		}
    }
    // write BWT to file 
    int divideNumber = bwtInc64->bwt->divideNumber;
    BWTSaveBwtCodeAndOcc64(bwtInc64->bwt,bwt_name,0);
    BWTIncFree64(bwtInc64);
    // add Occurrence and SA ,LCP 
    {
        // add Occurrence count and Update BWT...
        bwt_t *bwt ;
        t = clock();
        fprintf(arguments->logfp, "[generateBWT] Update BWT... ");
        bwt = bwt_restore_bwt(bwt_name);
        bwt_bwtupdate_core(bwt);
        fprintf(arguments->logfp, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
        
        // add LCP
        //lcp = (uint8_t*)calloc((bwt->seq_len + bwt->divideNumber + BITS_IN_BYTE - 1) / BITS_IN_BYTE , sizeof(uint8_t));
        /* first use as marker of lcp, every position use 4 bits represent, first bit denote whether that position has
         * certain length lcp before that position of string , second bit denote whether that position has certain length
         * lcp after that position of string. */
        time(&timeval); fprintf(arguments->logfp, "time:\t%s\n", ctime(&timeval));
        fprintf(arguments->logfp, "[generateBWT] Begin set BWT LCP....\n");
        int64_t len_lcp = (bwt->seq_len + bwt->divideNumber + BITS_IN_BYTE -1)/BITS_IN_BYTE;
        //fprintf(stderr, "[generateBWT] bwt->seq_len: %ld, len_lcp: %ld\n",bwt->seq_len, len_lcp);
        lcp = (uint8_t*)xcalloc(len_lcp , sizeof(uint8_t));

/* #ifdef DEBUG            
        {
            char fn[PATH_LEN];
            int part = 2, len = 1000;
            sprintf(fn, "%stextBufferDivide_%d.tmp", prefix, part);
            int64_t part_len = pac_seq_len(fn);
            FILE *fp = xopen(fn, "rb");
            uint8_t str[len];
            fread(str, sizeof(uint8_t), len, fp);
            int64_t rIndex = part ;
            for(int64_t z = part_len - 1; z >= 0; z--)
            {
                int64_t bp = posSA2PosBWT(bwt, rIndex);
                uint32_t c = bwt_B0(bwt, bp) ;
                rIndex = bwt->L2[c] + bwt->occMajor[bp/OCC_INTERVAL_MAJOR64*4 + c] + bwt_occ(bwt, bp, c)+ bwt->divideNumber -1;
                if(z < 3000 && GET_CHAR_FROM_PACKED(str, z) != c)
                {
                    fprintf(stderr, "[generateBWT] recover not equal to initial.. %u\t%u\n", GET_CHAR_FROM_PACKED(str, z), c);
                    exit(1);
                }
                const char ntc[] = {'C', 'A', 'T', 'G'};
                if(z < 1000) fprintf(stderr, "%c", ntc[c]);
                if(z < 1000 && z %100 == 0) fprintf(stderr, "\n");
            }
			FILE *tmpfp = xopen("tmpfile_1", "w");
			fwrite(transformPacToString(str, 0, 1000), sizeof(uint8_t), 1000, tmpfp);
            fclose(tmpfp);
            fclose(fp);
        }
#endif */

        t = clock();
        int64_t width = (bwt->seq_len / DEFAULT_DEPTH) / AVERAGE_CAPACITY;
        int first_len = 6;
		//int fixed_len = 20;
        //while(pow(4, first_len) < width) { first_len++; }
        if(first_len > KMER/2) { fprintf(stderr, "[generateBWT] first length of LCP set wrong, exit...\n"); exit(1); }
        width = pow(4, first_len);
        #pragma omp parallel for schedule(guided)
        for(uint64_t i = 0; i < width; i++)
        {
            setBWTLCP(bwt, lcp, i, KMER+1, first_len, arguments->min_kmerfreq);
        }
        fprintf(arguments->logfp, "[generateBWT] finish setting BWT LCP....\n");
        
        if(bwt->sa != NULL) free(bwt->sa);
        bwt->sa_intv = SA_INTERVAL ;
        bwt->n_sa = (bwt->seq_len + bwt->divideNumber + bwt->sa_intv - 1) / bwt->sa_intv ;
        bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
        time(&timeval); fprintf(arguments->logfp, "time:\t%s\n", ctime(&timeval));
        fprintf(arguments->logfp, "[generateBWT] Begin set BWT SA ....\n");
        if(bwt->isDivide == 0)
        {
            setBWTSA(bwt, bwt->seq_len , 0, bwt->divideNumber , 0);
        }else {
            #pragma omp parallel for schedule(dynamic)
            for(int i = 0 ; i < (int)bwt->divideNumber ; i++)
            {
                int64_t end ;
                if(i == (int)bwt->divideNumber - 1) end = 0 ;
                else end = bwt->divideBoundary[i+1].textLocation;
                setBWTSA(bwt, bwt->divideBoundary[i].textLocation , end, (int)bwt->divideNumber , i);
            }
        }
        time(&timeval); fprintf(arguments->logfp, "finished set\ntime:\t%s\n[generateBWT]Begin save kmer SA....\n", ctime(&timeval));
#ifdef DEBUG
        fprintf(stderr, "finished set BWT SA\n");
        fprintf(arguments->logfp, "Indexing enhanced suffix arrays used %.2f sec\n", (float)(clock() - index_t) / CLOCKS_PER_SEC);
#endif
        // scan lcp array and save kmer SA 
        if(bwt->divideNumber == 0) bwt->divideNumber = 1 ;
        strcpy(kmerSAF , prefix); strcat(kmerSAF , ".kmerSA");
        FILE *kmerSAfp = xopen(kmerSAF, "w");
        if(bwt->isDivide  == 1)
        {
            #pragma omp parallel for schedule(dynamic)
            for(int i = 0 ; i < (int)bwt->divideNumber ; i++)
            {
                if(i == 0) saveKmerSA(lcp , bwt , bwt->divideBoundary[i+1].textLocation , bwt->seq_len + bwt->divideNumber, libIndex_info, kmerSAfp, arguments->min_kmerfreq, KMER+1);
                if( i < (int)bwt->divideNumber - 1)
                    saveKmerSA(lcp , bwt , bwt->divideBoundary[i+1].textLocation , bwt->divideBoundary[i].textLocation, libIndex_info, kmerSAfp, arguments->min_kmerfreq, KMER+1);
                else saveKmerSA(lcp , bwt , bwt->divideNumber , bwt->divideBoundary[i].textLocation , libIndex_info, kmerSAfp, arguments->min_kmerfreq, KMER+1);
            }
        } else {
            saveKmerSA(lcp , bwt , bwt->divideNumber , bwt->seq_len , libIndex_info, kmerSAfp , arguments->min_kmerfreq, KMER+1);
        }
        fclose(kmerSAfp);
        time(&timeval); fprintf(arguments->logfp, "finished saving\ntime:\t%s\n", ctime(&timeval));

        /*
        // scan lcp array , save kmer SA and check long zero range if missing mark 
        if(bwt->divideNumber == 0) bwt->divideNumber = 1 ;
        FILE *kmerSAfp, *longZeroSAfp;
        strcpy(kmerSAF, prefix); strcat(kmerSAF, ".kmerSA");
        strcpy(longZeroSAF, prefix); strcat(longZeroSAF, ".longZeroSA");
        kmerSAfp = xopen(kmerSAF, "wb"); longZeroSAfp = xopen(longZeroSAF, "wb");
        if(bwt->isDivide != 0)
        {
            //#pragma omp parallel for schedule(dynamic)
            for(int i = 0 ; i < bwt->divideNumber ; i++)
            {
            if(i < bwt->divideNumber - 1)
                saveKmerSAAndLongZeroSA(bwt->divideBoundary[i+1].textLocation, bwt->divideBoundary[i].textLocation,kmerSAfp, longZeroSAfp);
            else saveKmerSAAndLongZeroSA(bwt->divideNumber, bwt->divideBoundary[i].textLocation, kmerSAfp,longZeroSAfp); 
            }
        } else {
            saveKmerSAAndLongZeroSA(bwt->divideNumber, bwt->seq_len , kmerSAfp, longZeroSAfp);
        }
        fclose(kmerSAfp); fclose(longZeroSAfp);
        */

        bwt_dump_bwt(bwt_name, bwt);
        bwt_dump_sa(sa_name, bwt);
        bwt_destroy(bwt);
    }
    
    {
        /*
        // calculate  long Zero Gaps lcp  
        uint64_t seq_len = pac_seq_len(pac_name);
        uint8_t *pac = (uint8_t*)calloc(seq_len / CHAR_PER_BYTE + 1 , sizeof(uint8_t));
        FILE *pacfp = xopen(pac_name, "rb");
        FILE *longZeroSAfp = xopen(longZeroSAF, "rb");
        fread(pac, sizeof(uint8_t), seq_len / CHAR_PER_BYTE + 1, pacfp );
        uint64_t longZeroNum = fseek(longZeroSAfp, 0, SEEK_END) / sizeof(LongZeroSA);
        LongZeroSA *longZeroSABuf = (LongZeroSA*)calloc(longZeroNum, sizeof(LongZeroSA));
        fseek(longZeroSAfp, 0, SEEK_SET);
        fread(longZeroSABuf, sizeof(LongZeroSA), longZeroNum, longZeroSAfp);

        //#pragma omp parallel for schedule(dynamic)
        for(int i = 0 ; i < NCPU ; i++)
        {
            if(i == NCPU - 1)
            {
                setLcpFromLongZeroSA(longZeroSABuf, i * (longZeroNum / NCPU),lengZeroNum, pac, lcp);
            } else {
                setLcpFromLongZeroSA(longZeroSABuf, i *(longZeroNum/NCPU), (i+ 1) * (longZeroNum/NCPU), pac,lcp);
            }
        } 
        free(longZeroSABuf); 
        fclose(pacfp); fclose(longZeroSAfp);
        if(remove(longZeroSAF) == -1) fprintf(stderr, "remove file: %s error\n", longZeroSAF);
        */

        // calculate kmerFreq and write to .kmerFreq file 
        char kmerFreq_name[PATH_LEN];
        FILE *kmerFreqfp , *lcpfp ;
        char lcp_name[PATH_LEN];
        strcpy(kmerFreq_name, prefix); strcat(kmerFreq_name, ".kmerFreq");
        kmerFreqfp = xopen(kmerFreq_name, "w");
        calculateKmerFreqAndWriteToFile(lcp, kmerFreqfp, (libIndex_info->len_pac + divideNumber + BITS_IN_BYTE - 1)/BITS_IN_BYTE);
        fclose(kmerFreqfp);
        // write lcp to file
        strcpy(lcp_name, prefix); strcat(lcp_name, ".lcp.gz");
        lcpfp = xzopen(lcp_name, "w");
        gzwrite(lcpfp, lcp, sizeof(uint8_t) * (libIndex_info->len_pac + divideNumber + BITS_IN_BYTE - 1)/BITS_IN_BYTE);
        gzclose(lcpfp);
        free(lcp) ;
    }

    // sort SA ,make less cache missing  and extract kmer and reverse complement kmer 
    {
        FILE *kmerSAfp = xopen(kmerSAF, "r");
        uint64_t SALen ;
        fseek(kmerSAfp, 0 , SEEK_END);
        SALen = ftell(kmerSAfp)/sizeof(KmerFreq);
        fseek(kmerSAfp, 0 , SEEK_SET);
        KmerFreq *kmerSABuf = (KmerFreq*)xcalloc(SALen, sizeof(KmerFreq));
        fread(kmerSABuf, sizeof(KmerFreq), SALen, kmerSAfp);
        t = clock();
        #pragma omp parallel for schedule(dynamic)
        for(int i = 0 ; i < NCPU ; i++)
        {
            if(i == NCPU - 1) qsort(kmerSABuf + i * (SALen/NCPU), SALen - i * (SALen/NCPU), sizeof(KmerFreq), compare64);
            else qsort(kmerSABuf + i * (SALen/NCPU), SALen/NCPU, sizeof(KmerFreq), compare64 );
        }

        fprintf(arguments->logfp, "[generateBWT] qsort() used %.2f CPU secs\n", (float)(clock() - t) / CLOCKS_PER_SEC);
        // merge sort
        {
            uint64_t len = SALen / NCPU ;
            uint64_t sortedNum ;
            KmerFreq *swapBuf = (KmerFreq*)xcalloc(len, sizeof(KmerFreq));
            sortedNum = SALen/NCPU + SALen % NCPU ;
            for(int i = NCPU - 2 ; i >= 0; i-- )
            {
                memcpy(swapBuf, kmerSABuf + i * len , len * sizeof(KmerFreq));
                mergeSort(swapBuf, len, kmerSABuf + (i + 1) * len, sortedNum, kmerSABuf + i * len);
                sortedNum += len ;
            }
            free(swapBuf);
        }
        t = clock();

        // extract kmer 
        {
            uint64_t seq_len = libIndex_info->len_pac ;
            uint8_t *pac = (uint8_t*)xcalloc(seq_len / CHAR_PER_BYTE + 1 , sizeof(uint8_t));
            gzFile pacfp = xzopen(pac_name, "r") ;
			FILE *kmerfp;
            char kmerfile_name[PATH_LEN];
            gzread(pacfp, pac, sizeof(uint8_t) * (seq_len / CHAR_PER_BYTE + 1));
            strcpy(kmerfile_name, prefix); strcat(kmerfile_name, ".kmer");
            kmerfp = xopen(kmerfile_name, "w");
            omp_lock_t lock_kmerfp ;
            omp_init_lock(&lock_kmerfp);

            #pragma omp parallel for firstprivate(kmerSABuf, kmerfp, pac) schedule(dynamic)
            for(int i = 0 ; i < NCPU ; i++)
            {
                if(i == NCPU - 1) extractKmerFromPac(kmerSABuf + i * (SALen/NCPU), SALen - i * (SALen/NCPU), kmerfp, pac , &lock_kmerfp );
                else extractKmerFromPac(kmerSABuf + i * (SALen/NCPU), SALen/NCPU , kmerfp, pac , &lock_kmerfp);
            }
            omp_destroy_lock(&lock_kmerfp);
            free(pac);
            gzclose(pacfp) ; fclose(kmerfp);
            fprintf(arguments->logfp, "[generateBWT] extractKmerFromPac() used %.2f CPU secs\n", (float)(clock() - t) / CLOCKS_PER_SEC);
        }
        free(kmerSABuf);
    }

	// free work
	lib_infoFree(libIndex_info) ;
}

int CheckPacFile(const Arguments *arguments, const lib_info *libIndex_info, const BntWriteArgs *bntWriteArgs)
{
	char fileName[PATH_LEN];
	gzFile pacF ;
	
	strcpy(fileName, arguments->prefix); strcat(fileName, ".pac.gz");
	pacF = xzopen(fileName, "r");
	uint8_t *pac = xcalloc((libIndex_info->len_pac + CHAR_PER_BYTE - 1)/CHAR_PER_BYTE, sizeof(uint8_t)) ;
	gzread(pacF, pac, (libIndex_info->len_pac + CHAR_PER_BYTE -1)/CHAR_PER_BYTE);
	int64_t pos = 0 ;	
	for(int i = 0 ; i < libIndex_info->num_lib ; i++ )	
	{
		if(libIndex_info->lib[i].paired == 1)
		{
			sprintf(fileName, "reverPacFile_%d_1.seq", i);
			FILE *rpf_1 = xopen(fileName, "w");
			sprintf(fileName, "reverPacFile_%d_2.seq", i);
			FILE *rpf_2 = xopen(fileName, "w");
			for(int j = 0; j < libIndex_info->lib[i].number_rd ; j += 2)
			{
				int read_len1 = libIndex_info->lib[i].diverse == 0 ? libIndex_info->lib[i].read_len : libIndex_info->lib[i].length[j] ;
				int read_len2 = libIndex_info->lib[i].diverse == 0 ? libIndex_info->lib[i].read_len : libIndex_info->lib[i].length[j+1] ;
				
            	char *s = transformPacToString(pac, pos, pos + read_len1);
				fprintf(rpf_1, "%s\n", s); free(s);
				pos += read_len1 ;
            	s = transformPacToString(pac, pos, pos + read_len2);
				fprintf(rpf_2, "%s\n", s); free(s);
				pos += read_len2 ;
				pos += (read_len1 + read_len2) ;
			}

			fclose(rpf_1) ;	
			fclose(rpf_2) ;	

		} else {
			sprintf(fileName, "reverPacFile_%d.seq", i);
			FILE *rpf = xopen(fileName, "w");
			for(int j = 0; j < libIndex_info->lib[i].number_rd ; j++)
			{
				int read_len = libIndex_info->lib[i].diverse == 0 ? libIndex_info->lib[i].read_len : libIndex_info->lib[i].length[j] ;
				
            	char *s = transformPacToString(pac, pos, pos + read_len);
				fprintf(rpf, "%s\n", s); free(s);
				pos += read_len * 2;
			}
			fclose(rpf) ;	
		}
	}

	// finalize
	free(pac);
	gzclose(pacF);

#ifdef DEBUG
    fprintf(arguments->logfp, "print out libIndex_info structure:\n");   
    fprintf(arguments->logfp, "#############################################\n");
    for(int i = 0 ; i < (int)(libIndex_info->num_lib); i++)
    {
        fprintf(arguments->logfp,"#####%s#####\n",libIndex_info->lib[i].name);
        if(libIndex_info->lib[i].paired == 1)
            fprintf(arguments->logfp,"insert_size: %u\ninsert_SD: %u\n",libIndex_info->lib[i].pe.insert_size,libIndex_info->lib[i].pe.insert_SD);
        fprintf(arguments->logfp,"offset: %lu\nread_len: %u\n",libIndex_info->lib[i].offset,libIndex_info->lib[i].read_len);
        if(libIndex_info->lib[i].diverse){
            fprintf(arguments->logfp, "length: ");
            for(int j = 0 ; j < 11 ; j++) fprintf(arguments->logfp,"%u ",libIndex_info->lib[i].length[j]);
        }
        fprintf(arguments->logfp, "\nnumber_rd: \%u\ndiverse: %u\n",libIndex_info->lib[i].number_rd,libIndex_info->lib[i].diverse);
    }
    fprintf(arguments->logfp, "#####libIndex_info#####\n");
    fprintf(arguments->logfp, "len_pac: %ld\n",libIndex_info->len_pac);
    fprintf(arguments->logfp, "num_seqs: %ld\n",libIndex_info->num_seqs);
    fprintf(arguments->logfp, "num_lib: %d\n",libIndex_info->num_lib);
    fprintf(arguments->logfp, "#############################################\n");

#endif
	
	return 1 ;
}

int bwtIndex(lib_info *libIndex_info, BntWriteArgs *bntWriteArgs, Arguments *arguments)
{
    clock_t t ;
    char filename[PATH_LEN];
    
    // transform to package
    t = clock();
    bns2bntseq(libIndex_info, bntWriteArgs, arguments);
    fprintf(arguments->logfp, "[bwtindex] Pack reads...");
    fprintf(arguments->logfp, "%.2f seconds elapse\n", (float)(clock() -t)/CLOCKS_PER_SEC);
    
	/*
    // add the complementary reverse package sequence to package file
    t = clock();
    strcpy(filename, arguments->prefix);    strcat(filename, ".pac");
    pac_rev_core(filename);
    fprintf(arguments->logfp, "[bwtindex] Reverse the packed sequence...%.2f seconds elapse\n", (float)(clock() - t)/CLOCKS_PER_SEC);
	*/

    //assert((libIndex_info->len_pac) == pac_seq_len(filename));

    // write libIndex_info content to the .ann file
    strcpy(filename, arguments->prefix);    strcat(filename, ".ann");
    writeInfo2Ann(filename, libIndex_info);
    arguments->seq_len = libIndex_info->len_pac ;

    // DEBUG information
	//CheckPacFile(arguments, libIndex_info, bntWriteArgs);
    
    lib_infoFree(libIndex_info);

    //Construct BWT for the packed sequence
    t = clock();
    generateBWT(arguments);
    fprintf(arguments->logfp, "[bwtindex] Construct BWT for the packed sequence....%.2f seconds elapse.\n", (float)(clock() -t) / CLOCKS_PER_SEC);

    return 0 ;
}

int index_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:     %s index [options]\n\n", PROGRAM_NAME);
    fprintf(stderr, "Options:   -p STR          the prefix of output files [same as program name]\n");
    fprintf(stderr, "           -c STR          the configure file that contain library information[%s.ini]\n", PROGRAM_NAME);
    fprintf(stderr, "           -N INT          Maximum number of CPU that program used( or the number of threads created) [use all the rest of CPU that system allowed, depend on the system load]\n");
    fprintf(stderr, "           -M INT(Giga)    Maximum RAM memory used(must bigger than X times of the total reads length) [no limitation]\n");
    fprintf(stderr, "           -G INT(Mega)    approximate genome size by prior knonw, setting genome size if assemble a high heterozygosity genome\n");
    fprintf(stderr, "           -K INT          The kmer size used for generating graph, must small than minimum read length that used construct contig, program will filter out the reads than small than K+5, the value must be odd and between %d~%d[%d]\n", MIN_KMER, MAX_KMER, DEF_KMER);
    fprintf(stderr, "           -o STR          program output directory by running-time[current directory]\n");
    fprintf(stderr, "           -q INT          the quality of base calling by cutting off, 1 is low quality, 2 is middle value, 3 is high, 4 is unfiltered. [3]\n");
    fprintf(stderr, "           -H INT          heterozygosity of diploid genome, low heterozygosity set 1, middle set 2, high heterozygosity set 3 [1]\n");
    return 0 ;
}

int hwgsa_index(int argc, char *argv[])
{
    int c ;
    char name[PATH_LEN];
    Arguments *arguments ;
    PConfReturn *retV = NULL ;
    lib_info *libIndex_info = NULL ;
    BntWriteArgs *bntWriteArgs = NULL ;

    arguments = (Arguments*)xcalloc(1 , sizeof(Arguments));

    if(argc < 2) { index_usage(); return 1 ; }
    while((c = getopt(argc, argv, "c:p:N:M:G:K:o:q:H:h")) >=0)
    {
        switch(c) {
            case 'c':   strcpy(arguments->conffile, optarg); break ;
            case 'p':   strcpy(arguments->prefix, optarg); break ;
            case 'N':   arguments->NCPU = atoi(optarg); break ;
            case 'M':   arguments->maxMem = atoi(optarg); break ;
            case 'G':   arguments->G = atoi(optarg); break ;
            case 'K':   arguments->K = atoi(optarg); break ;
            case 'o':   strcpy(arguments->outputDir, optarg); break ;
            case 'q':   arguments->qual = atoi(optarg); break;
            case 'H':   arguments->H = atoi(optarg); break ;
            case 'h':   index_usage() ; exit(1) ;
            case ':':   fprintf(stderr, "option %c need a value\n", optopt); index_usage(); return 1 ;
            case '?':   fprintf(stderr, "[hwgsa_all] unrecognized option '%c'\n",optopt ); index_usage(); return 1 ;
            default: index_usage() ; return 1 ;
        }
    }
    preProcessArgs(arguments);
    KMER = arguments->K ;

    strcpy(name, arguments->prefix); strcat(name, ".log");
    arguments->logfp = xopen(name, "w");

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
    fprintf(arguments->logfp, "###########################end#########################\n");

#ifdef DEBUG
    fflush(arguments->logfp); fflush(stderr); fflush(stdout);
#endif
    
    retV = parseConffile(arguments->conffile);
    libIndex_info = retV->libIndex_info ;
    bntWriteArgs = retV->bntWriteArgs ;
    free(retV);
    // construct BWT index 
    bwtIndex(libIndex_info, bntWriteArgs, arguments);

    // write arguments to file *.args
    {
        FILE *fp ;
        strcpy(name, arguments->prefix); strcat(name, ".args");
        fp = xopen(name, "w");
        writeArgsToFile(arguments, fp);
        fclose(fp);
    }
    // clean and free work
    //free_BntWriteArgs(bntWriteArgs);
    //free_lib_info(libIndex_info);
    fclose(arguments->logfp);
    free_Arguments(arguments);

    return 0 ;
}
