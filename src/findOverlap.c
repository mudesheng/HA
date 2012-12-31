#include "findOverlap.h"
#include "bwt.h"


static inline uint8_t *revcompleBnt(const uint8_t *bnt_seq, const int len)
{
	uint8_t *rev_bnt = (uint8_t*)xcalloc(len, sizeof(uint8_t));
	for(int i = 0; i < len; i++)
	{
		rev_bnt[i] = (~(bnt_seq[len - 1 - i])) & 0X3;
	}

	return rev_bnt;
}

static inline uint8_t *revBnt(const uint8_t *bnt_seq, const int len)
{
    uint8_t *rev_bnt = (uint8_t*)xcalloc(len, sizeof(uint8_t));
    for(int i = 0; i <len; i++)
    {
        rev_bnt[len-1-i] = bnt_seq[i];
    }

    return rev_bnt;
}

static inline uint8_t *compleBnt(const uint8_t *bnt_seq, const int len)
{
    uint8_t *comple_bnt = (uint8_t*)xcalloc(len, sizeof(uint8_t));
    for(int i = 0; i < len; i++)
    {
        comple_bnt[i] = (~(bnt_seq[i])) & 0X3;
    }

    return comple_bnt;
}

static inline uint8_t *transformSeq2bnt(const kstring_t seq)
{
    uint8_t *bnt_seq = (uint8_t*)xcalloc(seq.locPos, sizeof(uint8_t));
    for(int i = 0; i < seq.locPos; i++)
    {
        bnt_seq[i] = nst_nt4_table[(int)seq.s[i]];
    }

    return bnt_seq;
}

static inline LCPBound locateLCPBound(const uint8_t *lcp, const int64_t len, const int64_t loc)
{
    LCPBound lcpBound;
    int64_t pos = loc;
    if(GET_LCP_FLAG(lcp, pos) > 0)
    {
        while(pos - 1 >= 0 && GET_LCP_FLAG(lcp, pos-1) > 0) { pos-- ;}
        lcpBound.low = pos - 1;
        while(pos + 1< len && GET_LCP_FLAG(lcp, pos+1) > 0) { pos++; } 
        lcpBound.high = pos ;
    } else {
        lcpBound.low = pos;
        while(pos + 1 < len && GET_LCP_FLAG(lcp, pos+1) > 0) { pos++; }
        lcpBound.high = pos;
    }

    return lcpBound;
}

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

static inline int checkAndAddISABase(ISABase *isaBase, const int64_t isa)
{
	int flag = 0;
	for(int i = isaBase->count - 1; i >= 0; i--)
	{
		if(isaBase->db[i] == isa) { flag = 1; break; }
	}
	if(flag == 0) 
	{
		if(isaBase->count >= isaBase->max)
		{
			isaBase->max <<= 1;
			isaBase->db = (int64_t*)xrecalloc(isaBase->db, isaBase->max * sizeof(int64_t));
		}
		isaBase->db[isaBase->count++] = isa;
	}
	return flag;
}

//flag return if have been processed
static inline int calculateSAAndRecoverBase(const bwt_t *bwt,  ReadInfo *ri, int64_t isa, ISABase *isaBase)
{
    int offset = 0 ;
    uint32_t c = 0;
	// check if have been in the ISABase
	int flag = checkAndAddISABase(isaBase, isa);
	if(flag == 1) return 1;
	if(ri->max == 0)
	{
		ri->max = AVERAGE_READ_LEN;
		ri->count = 0;
		ri->bnt = (uint8_t*)xcalloc(ri->max, sizeof(uint8_t));
	} else {
		fprintf(stderr, "[calculateSAAndRecoverBase] the element bnt of ReadInfo set error, exit...\n");
		exit(1);
	}

    while(isa % SA_INTERVAL != 0)
    {
        int64_t bp ;
        int zcount ;
        for(zcount = 0; zcount < bwt->divideNumber; zcount++)
        {
            if(isa < bwt->sentinelPosition[zcount]) break ;
            else if(isa == bwt->sentinelPosition[zcount]) {
                ri->sa = offset + bwt->sentinelSA[zcount];
				ri->isa = isa;
				return 0;
            }
        }
        bp = isa - zcount ;
        c = bwt_B0(bwt, bp);
		if(ri->count >= ri->max)
		{
			ri->max <<= 1;
			ri->bnt = (uint8_t*)xrecalloc(ri->bnt, ri->max * sizeof(uint8_t));
		}
		ri->bnt[ri->count++] = c ;
        isa = bwt->L2[c] + bwt->occMajor[bp/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp, c) + bwt->divideNumber  - 1;
        offset++ ;
        // check if has been in the ISABase
        if(checkAndAddISABase(isaBase, isa) == 1) { return 1; }
    }

	ri->isa = isa;
	ri->sa = bwt->sa[isa / SA_INTERVAL] + offset ;
	return 0;
}

static inline int recoverBase(const bwt_t *bwt, ReadInfo *ri, const int len, ISABase *isaBase)
{
	int64_t isa = ri->isa;
	for(int i = 0; i < len; i++)
	{
        int zcount ;
		uint32_t c;
		int64_t bp ;
        for(zcount = 0; zcount < bwt->divideNumber; zcount++)
        {
            if(isa < bwt->sentinelPosition[zcount]) break ;
            else if(isa == bwt->sentinelPosition[zcount]) {
                fprintf(stderr, "[recoverBase] isa == bwt->sentinelPosition, exit...\n");
				exit(1);
            }
        }
		bp = isa - zcount;
        c = bwt_B0(bwt, bp);
		if(ri->count >= ri->max)
		{
			ri->max <<= 1;
			ri->bnt = (uint8_t*)xrecalloc(ri->bnt, ri->max * sizeof(uint8_t));
		}
		ri->bnt[ri->count++] = c ;
        isa = bwt->L2[c] + bwt->occMajor[bp/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp, c) + bwt->divideNumber  - 1;
        // check if has been in the ISABase
        if(checkAndAddISABase(isaBase, isa) == 1) 
        { return 1; }
	}
	ri->isa = isa;

    return 0;
}

int checkAndAddToReadBase(const bwt_t *bwt, const lib_info *libIndex_info, ReadBase *rb, ISABase *isaBase, const int64_t isa, const int seed_len)
{

	//rb->ri[rb->count].isa = isa;
	ReadInfo ri; memset(&ri, 0, sizeof(ReadInfo));
	if(calculateSAAndRecoverBase(bwt, &ri, isa, isaBase) == 1)
    {
        if(ri.max >0) free(ri.bnt);
        return 1;
    }
	ri.rl = readBoundaryLookup(libIndex_info, ri.sa);
	// spurious overlap 
    if(ri.sa + seed_len > ri.rl.bound + ri.rl.length) 
    { 
        free(ri.bnt);  
        return 1; 
    }
	if(ri.count < ri.sa - ri.rl.bound)
	{
		int len = ri.sa - ri.rl.bound - ri.count;
		if(recoverBase(bwt, &ri, len, isaBase) == 1)
        {
            free(ri.bnt);
            return 1;
        }
	}else {
        isaBase->count -= (ri.count - (ri.sa - ri.rl.bound));
		ri.count = ri.sa - ri.rl.bound;
	}
    
	if(rb->count >= rb->max)
	{
		rb->max <<= 1;
		rb->ri = (ReadInfo*)xrecalloc(rb->ri, rb->max * sizeof(ReadInfo));
	}
    for(int i = 0; i < rb->count; i++)
    {
        if(abs(rb->ri[i].sa - ri.sa) < 100)
        {
            fprintf(stderr, "[checkAndAddToReadBase]maybe replication recovered\n");
        }
    }
    rb->ri[rb->count] = ri;
    rb->count++;

    // free and clean work

	return 0;
}

static inline int searchReadBaseBySa(const ReadBase *rb, const int64_t sa)
{
    int index = -1 ;
    for(int i = 0; i < rb->count; i++)
    {
        if(rb->ri[i].sa == sa)
        {
            index = i;
            break;
        }
        /*if(abs(rb->ri[i].sa - sa) < 100)
        {
            fprintf(stderr, "%dth sa : %ld\n", i, rb->ri[i].sa);
        }
        if(rb->ri[i].rl.bound == bound)
        {
            fprintf(stderr, "index: %d\n", i);
        }
        fprintf(stderr, "%ld\t", rb->ri[i].sa); */
    }
    /*if(index == -1)
    {
        fprintf(stderr, "[searchReadBaseBySa] index not found, exit...\n");
        exit(1);
    } */
    return index ;
}



OverlapInfo align_seq_core(const uint8_t *bnt1, const uint8_t *bnt2, const int bnt1_len, const int bnt2_len, const int step, const OverlapInfo oi, const int max_diff, const uint8_t last_char)
{
    OverlapInfo oi_new = oi;
    for(int i =  0, j = 0; i < bnt1_len && j < bnt2_len; i++, j++)
    {
        if(bnt1[i] == bnt2[j]) continue;
        else {
            uint8_t c ;
            if(i == 0) c = last_char;
            else c = bnt1[i-1];
            if(i < bnt1_len - step - (max_diff - oi_new.max_diff) && j < bnt2_len - step - (max_diff - oi_new.max_diff))
            {
                int flag1 = 0;
                int z;
                // mismatch
                for(z = 0; z < max_diff - oi_new.max_diff; z++)
                {
                    if(memcmp(bnt1 + i + z + 1, bnt2 + j + z + 1, step) == 0) 
                    {
                        i += (step + z);
                        j += (step + z);
                        flag1 = 1; break; 
                    } 
                }
                if(flag1 == 1) { oi_new.max_diff += (z+1); continue; }

                // indel
                flag1 = 0;
                for(z = 0; z < max_diff - oi_new.max_diff; z++)
                {
                    // deletion
                    int flag2 = 1;
                    for(int y = 0; y <= z; y++)
                    {
                        if(c != bnt1[i+z]) { flag2 = 0; break; }
                    }
                    if(flag2 == 1 && memcmp(bnt1 + i + z + 1, bnt2 + j , step) == 0)
                    {
                        i += (z + step);
                        j += (step -1);
                        flag1 = 1;
                        break;
                    }
                    // insertion
                    flag2 = 1;
                    for(int y = 0; y <= z; y++)
                    {
                        if(c != bnt2[j+z]) { flag2 = 0; break; }
                    }
                    if(flag2 == 1 && memcmp(bnt1 +  i , bnt2 + j + z + 1, step) == 0)
                    {
                        i += (step - 1);
                        j += (z + step);
                        flag1 = 1;
                        break;
                    } 
                }
                if(flag1 == 1) { oi_new.max_diff += (z+1);}
                else { oi_new.flag = 0; break; }
            } else {
                int flag1;
                // mismatch
                int z;
                flag1 = 1;
                for(z = 1; z < bnt1_len - i && z < bnt2_len - j; z++)
                {
                    if(bnt1[i+z] != bnt2[j+z]) { flag1 = 0; break; }
                }
                if(flag1 == 1) { oi_new.max_diff += 1; break;}
                // indel 
                // deletion
                flag1 = 0;
                if(c == bnt1[i])
                {
                    int flag2 = 1;
                    for(z = 0; z < bnt1_len - i - 1 && z < bnt2_len - j; z++)
                    {
                        if(bnt1[i+1+z] != bnt2[j+z]) { flag2 = 0; break;}
                    }
                    if(flag2 == 1) { i += z; j += (z-1); flag1 = 1;  }
                } else if(c == bnt2[j]) { // insertion
                    int flag2 = 1;
                    for(z = 0; z < bnt1_len - i && z < bnt2_len - j - 1; z++)
                    {
                        if(bnt1[i+z] != bnt2[j+1+z]) { flag2 = 0; break; }
                    }
                    if(flag2 == 1) { i += (z-1); j += z; flag1 = 1; }
                } 

                if(flag1 == 1) { oi_new.max_diff += 1;}
                else { oi_new.flag = 0;}
                break;
            }
        }
    }
    oi_new.consis_len += ((bnt1_len <= bnt2_len ? bnt1_len : bnt2_len) - (oi_new.max_diff - oi.max_diff));

    return oi_new;
}

OverlapInfo checkAndConsistent(const bwt_t *bwt, const lib_info *libIndex_info, ReadBase *rb, ISABase *isaBase, const int64_t isa, const uint8_t *rev_bnt, const int seq_len, const int seq_start, const int seed_len, const int max_diff)
{
    OverlapInfo oi;
    int step = 5;
    oi.consis_len = 0; oi.max_diff = 0; oi.flag = 1;
    // check if have been in the ISABase
    //int flag = checkAndAddISABase(isaBase, isa);
    //if(flag == 1) { oi.flag = 0; return oi; }

    ReadInfo *ri = (ReadInfo*)xcalloc(1, sizeof(ReadInfo));
	//rb->ri[rb->count].isa = isa;
	if(calculateSAAndRecoverBase(bwt, ri, isa, isaBase) == 1)
    {
        isaBase->count -= (ri->count + 1);
        if(ri->max > 0) free(ri->bnt); 
        free(ri);
        oi.flag = 0;
        return oi;
    }
    int64_t sa = bwt->seq_len - 1 - (ri->sa + seed_len - 1);
    //int64_t sa = bwt->seq_len - 1 - ri->sa;
    int index ;
    //fprintf(stderr, "sa: %ld\n", sa);
    if((index = searchReadBaseBySa(rb, sa)) == -1)
    {
        isaBase->count -= (ri->count+1);
        free(ri->bnt); free(ri);
        oi.flag = 0;
        return oi;
    }
    ri->rl.bound = bwt->seq_len - rb->ri[index].rl.bound - rb->ri[index].rl.length;
    ri->rl.length = rb->ri[index].rl.length;
	
	if(ri->count < ri->sa - ri->rl.bound)
	{
		int len = ri->sa - ri->rl.bound - ri->count;
		if(recoverBase(bwt, ri, len, isaBase) == 1)
        {
            fprintf(stderr, "[checkAndConsistent] recoverBase have been recover before, exit...\n");
            exit(1);
        }
	}else {
        isaBase->count -= (ri->count - (ri->sa - ri->rl.bound));
		ri->count = ri->sa - ri->rl.bound;
	}
    
    // alignment
    // right part
    uint8_t *bnt1 = revBnt(rev_bnt, seq_start);
    oi = align_seq_core(bnt1, ri->bnt, seq_start, ri->count, step, oi, max_diff, rev_bnt[seq_start]);
    if(oi.flag == 1) // left alignment
    {
        uint8_t *comple_bnt = compleBnt(rb->ri[index].bnt, rb->ri[index].count);
        oi = align_seq_core(rev_bnt + seq_start + seed_len, comple_bnt, seq_len - (seq_start + seed_len), rb->ri[index].count, step, oi, max_diff, rev_bnt[seq_start + seed_len -1]);

        // free and clean work
        free(comple_bnt);
    }

    // free and clean work
    free(bnt1);
    free(ri->bnt);
    free(ri);

    return oi;
}

void free_ReadBase(ReadBase *rb)
{
	for(int i = 0; i < rb->count; i++)
	{
		if(rb->ri[i].max >0) free(rb->ri[i].bnt);
	}
	free(rb->ri);
}

void findOverlap(Arguments *arguments, const char *readfileName, const int seed_len, const int max_diff)
{
     char fn[PATH_LEN];

    strcpy(fn, arguments->prefix); strcat(fn, ".bwt");
    bwt_t *bwt = bwt_restore_bwt_core(fn);
    strcpy(fn, arguments->prefix); strcat(fn, ".lcp");
    uint8_t *lcp = restore_lcp(fn, bwt->seq_len + bwt->divideNumber);
    strcpy(fn, arguments->prefix); strcat(fn, ".sa");
    bwt_restore_sa(fn, bwt);
    strcpy(fn, arguments->prefix); strcat(fn, ".ann");
    lib_info *libIndex_info = lib_info_restore(fn);
    strcpy(fn, arguments->prefix); strcat(fn, ".kmerFreq");
    KmerFreqCurve curve = parseKmerFreqFile(fn, arguments);
    gzFile fp = xzopen(readfileName, "r");
    strcpy(fn, arguments->prefix); strcat(fn, ".overlap");
    gzFile ofp = xzopen(fn, "w");
    int processed_rd = 0;
    kseq_t *seq = kseq_init(fp);
    while(kseq_read(seq) > 0)
    {
        int mapped_count = 0;
		processed_rd++;
        if(strcspn(seq->seq.s, "Nn") != seq->seq.locPos || seq->seq.locPos <= seed_len) continue;
        uint8_t *bnt_seq = transformSeq2bnt(seq->seq);
		uint8_t *rev_bnt = revcompleBnt(bnt_seq, seq->seq.locPos);
        int64_t low = -1, high = -1;
        int flag = 0;
        ISABase isaBase;
        isaBase.count = 0; isaBase.max = 500;
        isaBase.db = (int64_t*)xcalloc(isaBase.max, sizeof(int64_t));
        ReadBase rb;
		rb.count = 0; rb.max = 20;
		rb.ri = (ReadInfo*)xcalloc(rb.max, sizeof(ReadInfo));
		// traverse '+' strand 
        //int count_rb = rb.count;
        for(int i = seq->seq.locPos - 1; i >= 0; i--)
        {
            uint32_t c = bnt_seq[i];
            int64_t bp_l, bp_h;
            if(low == -1 && high == -1)
            {
                low = bwt->L2[c] + bwt->divideNumber ;
                high = bwt->L2[c+1] + bwt->divideNumber - 1 ; 
            } else {
                bp_l = posSAChangeToPosBWT(bwt, low, 0);
                bp_h = posSAChangeToPosBWT(bwt, high, 1);
                if(bp_l > bp_h || (bp_l == bp_h && c != bwt_B0(bwt, bp_l))) { flag = 1; break;}
                low = bwt->L2[c] + bwt->occMajor[bp_l/OCC_INTERVAL_MAJOR64*4 + c] + bwt_occ(bwt,bp_l, c) + bwt->divideNumber - 1;
                high = bwt->L2[c] + bwt->occMajor[bp_h/OCC_INTERVAL_MAJOR64*4 + c] + bwt_occ(bwt, bp_h, c) + bwt->divideNumber - 1;
                
				if(i <= seq->seq.locPos - seed_len)
                {
                    LCPBound lcpBound = locateLCPBound(lcp, bwt->seq_len, low);
                    if(lcpBound.low >= lcpBound.high || lcpBound.high - lcpBound.low >= curve.peak * 2) continue;
                    for(int64_t j = lcpBound.low; j <= lcpBound.high; j++)
                    {
						if(checkAndAddToReadBase(bwt, libIndex_info, &rb, &isaBase, j, seed_len) == 0);
                        /*{
                            for(int z = 0; z < rb.ri[rb.count-1].count; z++)
                            {
                                fprintf(stdout, "%d", rb.ri[rb.count-1].bnt[z]);
                            }
                            fprintf(stdout, "\n");
                            for(int z = i - 1; z >= 0; z--)
                            {
                                fprintf(stdout, "%d", bnt_seq[z]);
                            }
                            fprintf(stdout, "\n");
                        }*/
                    }
                    /*for(int j = count_rb; j < rb.count; j++)
                    {
                        fprintf(stderr, "%ld\t", rb.ri[j].sa);
                    }
                    fprintf(stderr, "\n");
                    count_rb = rb.count; */
                }
            } 
        }
		if(flag != 1)
		{
			// traverse '-' strand
			low = -1; high = -1;
            LCPBound lcpBoundBase[seq->seq.locPos - seed_len + 1];
			for(int i = seq->seq.locPos - 1; i >= 0; i--)
			{
				uint32_t c = rev_bnt[i];
				int64_t bp_l, bp_h;
				if(low == -1 && high == -1)
				{
					low = bwt->L2[c] + bwt->divideNumber;
					high = bwt->L2[c+1] + bwt->divideNumber - 1;
				} else {
					bp_l = posSAChangeToPosBWT(bwt, low, 0);
					bp_h = posSAChangeToPosBWT(bwt, high, 1);
					if(bp_l > bp_h || (bp_l == bp_h && c != bwt_B0(bwt, bp_l)))  
					{ 
						fprintf(stderr, "[findOverlap] the '-' strand not found in Database, exit...\n");
						exit(1);
					}
					low = bwt->L2[c] + bwt->occMajor[bp_l/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_l, c) + bwt->divideNumber - 1;
					high = bwt->L2[c] + bwt->occMajor[bp_h/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_h, c) + bwt->divideNumber - 1;
					
					if(i <= seq->seq.locPos - seed_len)
					{
						lcpBoundBase[i] = locateLCPBound(lcp, bwt->seq_len, low);
                        lcpBoundBase[i].depth = i;
					} 
				}
			}

            // consistent check
            for(int i = 0; i < seq->seq.locPos - seed_len + 1; i++)
            {
                if(lcpBoundBase[i].low >= lcpBoundBase[i].high || lcpBoundBase[i].high - lcpBoundBase[i].low >= curve.peak *2) continue;
                for(int j = lcpBoundBase[i].low; j <= lcpBoundBase[i].high; j++)
                {
                    OverlapInfo oi = checkAndConsistent(bwt, libIndex_info, &rb, &isaBase, j, rev_bnt, seq->seq.locPos, lcpBoundBase[i].depth, seed_len, max_diff);
					if(oi.flag == 1) mapped_count++;
                }
            }
			fprintf(stderr, "mapped_count: %d\n", mapped_count);
		}

        // free and clean work
		free_ReadBase(&rb);
		free(isaBase.db);
		free(rev_bnt);
        free(bnt_seq);
    }
    fprintf(stderr, "processed read number is : %d\n", processed_rd); 

    // free and clean work
	kseq_destroy(seq);
    gzclose(ofp);
    gzclose(fp);
    lib_infoFree(libIndex_info);
    free(lcp);
    bwt_destroy(bwt);
}



int findOverlap_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:     %s findOverlap [options]\n\n", PROGRAM_NAME);
    fprintf(stderr, "Options:   -p STR          the Prefix of output files [same as PROGRAM_NAME]\n");
    fprintf(stderr, "           -f STR          the input File that contain reads sequence(mandatory option)\n");
    fprintf(stderr, "           -N INT          maximum Number of CPU program used( or the number of threads created) [use all the rest of CPU that system allow, depend on the system loading]\n");
    fprintf(stderr, "           -o STR          program Output directory by running-time[current directory]\n");
    fprintf(stderr, "           -s INT          the seed Length of overlap[20]\n");
    fprintf(stderr, "           -d INT          the maximum difference allowed[0]\n");
    fprintf(stderr, "           -h              print this  message to the stderr\n");

    return 0;
}

int hwgsa_findOverlap(int argc, char *argv[])
{
    int c, seed_len = 20, max_diff = 0;
    char *readfileName = (char*)xcalloc(PATH_LEN, sizeof(char));
    Arguments *arguments;
    
    arguments = (Arguments*)xcalloc(1, sizeof(Arguments));
    
    if(argc < 2) { findOverlap_usage(); return 1; }
    while((c = getopt(argc, argv, "p:N:o:s:f:d:h")) >= 0)
    {
        switch(c) {
            case 'p': strcpy(arguments->prefix, optarg); break;
            case 'N': arguments->NCPU = atoi(optarg); break;
            case 'o': strcpy(arguments->outputDir, optarg); break;
            case 's': seed_len = atoi(optarg); break;
            case 'f': strcpy(readfileName, optarg); break;
            case 'd': max_diff = atoi(optarg);
            case 'h': findOverlap_usage(); exit(1);
            case ':': fprintf(stderr, "[hwgsa_findOverlap]option %c need a value\n", optopt); findOverlap_usage(); return 1;
            case '?': fprintf(stderr, "[hwgsa_findOverlap]unrecognized option '%c'\n", optopt); findOverlap_usage(); return 1 ;
            default: findOverlap_usage(); return 1 ;
        }
    }
    
    // check arguments
    if(arguments->outputDir[0] == '\0')
    {
        fprintf(stderr, "outputDir have not set.....\n");
        findOverlap_usage();
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
    if(readfileName[0] == '\0')
    {
        fprintf(stderr, "'f' not setting, 'f' is mandatory option\n");
        findOverlap_usage();
        exit(1);
    }
    
    findOverlap(arguments, readfileName, seed_len, max_diff);

    free(readfileName);
    free_Arguments(arguments);

    return 0;
}

