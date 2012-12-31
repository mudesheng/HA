#include <math.h>
#include "extend.h"

extern int32_t KMER ;
extern int32_t NCPU;


void kseqCopy2faElem(const kseq_t *seq, faElem *fe )
{
    fe->name.locPos = seq->name.locPos ;
    fe->name.max = seq->name.max;
    fe->name.s = strdup(seq->name.s);

    fe->comment.locPos = seq->comment.locPos;
    fe->comment.max = seq->comment.max;
    fe->comment.s = strdup(seq->comment.s);

    fe->seq.locPos = seq->seq.locPos;
    fe->seq.max = seq->seq.max;
    fe->seq.s = strdup(seq->seq.s);

}

/* return value: 0 denote return NULL , others denote return the kmer end position */
static int getBinaryKmer(uint64_t *k, const char *s, int start, const int max, const int count)
{
    int flag = 0;
    for(int i = 0; i < KMER; i++)
    {
        if(start > max) return 0;
        if(s[start] == '{')
        {
            int c = 0;
            while(c < count)
            {
                if(start >max || s[start] == '}') break;
                else if(s[start] == ',') c++;
                start++;
            }
            if(start > max) return 0;
            else if(s[start] == '}') return 0;
            if(s[start] == ',') 
            { 
                while(s[start] != '}' && start <= max) start++;  
                start++;
                if(start > max) return 0;
                else {
                    k[i/32] <<= BITS_PER_CHAR;
                    k[i/32] |= nst_nt4_table[s[start]];
                    start++;
                    continue;
                }
            }
            while(s[start] != ',' && s[start] != '}') 
            {
                if(i < KMER) 
                { 
                    k[i/32] <<= BITS_PER_CHAR; 
                    k[i/32] |= nst_nt4_table[s[start]]; 
                    start++; 
                    i++; 
                }else return 0 ;
            }
        } else {
            k[i/32] <<= BITS_PER_CHAR;
            k[i/32] |= nst_nt4_table[s[start]];
        }
        start++;
    }
    k[i/32] <<= ((32 - i%32)<<1);
    return start ; 
}

static inline int stepLen2Next(const char *s, const int step_len)
{
    int count =0 ;
    for(int i = 0; i < step_len; i++)
    {
        if(s[count] == '{' )
        {
            while(s[count] != '}') count++;
        }
        count++;
    }
    return count;
}

static inline Region calculateRegion(const bwt_t *bwt, const uint64_t *k, const int len)
{
    Region rg ;
    for(int i = 0; i < len; i++)
    {
        int c = (k[i/32] >> ((31 - i%32)<<1)) && 0x3 ;
        if(i = 0) { rg.low = bwt->L2[c] + bwt->divideNumber; rg.high = bwt->L2[c+1] + bwt->divideNumber -1;  }
        else {
            rg.low = posSA2PosBWT(bwt, rg.low); rg.high = posSA2PosBWT(bwt, rg.high);
            rg.low = bwt->L2[c] + bwt->occMajor[rg.low/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, rg.low, c) + bwt->divideNumber - 1;
            rg.high = bwt->L2[c] + bwt->occMajor[rg.high/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, rg.high, c) + bwt->divideNumber - 1;
        }
    }
    
    return rg;
}

void Seq_region_destroy(Seq_region *sr)
{
    for(int i = 0; i < sr->chance_num; i++)
    {
        free(sr->s[i]);
    }
}

Seq_region *getSeqRegion(const char *s, const int start, const int length, const int max, const int direction)
{
    
}

/* return value: 0 denote not in the ISABase, 1 denote in the ISABase */
static inline int searchISABase(const ISABase *isaBase, const int64_t isa)
{
    for(int i = 0; i < isaBase->count; i++)
    {
        if(isa == isaBase->baseIndex[i].isa) return 1 ;
    }

    return 0 ;
}


/* return value : 0 note not necessary add to MapIndex, 1 note need add to MapIndex */
static inline int checkNeedAddToMapIndex(const int64_t bp, const Seq_region *sr, const int len, const int direction, const bwt_t *bwt, const ISABase *isaBase)
{
    int64_t isa = bp; 
    uint8_t read_seq[len];
    for(int i = sr->len - 1 ; i >= len; i--)
    {
        int c = bwt_B0(bwt, isa);
        read_seq[i] = c;
        int flag = 0 ;
        for(int j = 0 ; j < sr->chance_num; j++)
        {
            if(c == sr->s[j][i]) { flag = 1; break; }
        }
        if(flag == 0) return 1 ;
        if(searchISABase(isaBase, isa) == 1) return 0;
        isa = bwt->L2[c] + bwt->occMajor[isa/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, isa, c) + bwt->divideNumber - 1;
        isa = posSA2PosBWT(bwt, isa);
    }
    for(int i = 0; i < sr->chance_num; i++)
    {
        if(memcmp(sr->s[i], read_seq, len * sizeof(uint8_t)) == 0) return 0; 
    }
    return 1 ;
}

static inline void reverseBitBase(uint8_t *s, const int len)
{
    for(int i = 0, j = len -1; i < j; i++, j--)
    {
        uint8_t c = s[i];
        s[i] = s[j]; j[j] = c ;
    }
}

SA_REV calculateSAAndGetBase(const bwt_t *bwt, int64_t isa)
{
    SA_REV sa_rev ;
    sa_rev.len = 0 ; sa_rev.max = 32 ;
    sa_rev.s = (uint8_t*)xcalloc(sa_rev.max, sizeof(uint8_t));
    int offset = 0;
    int c = 0;
    while(isa % SA_INTERVAL != 0)
    {
        int64_t bp ;
        int zcount;
        for(zcount = 0; zcount < bwt->divideNumber; zcount++)
        {
            if(isa < bwt->sentinelPosition[zcount]) break;
            else if(isa == bwt->sentinelPosition[zcount]) {
                sa_rev.sa = bwt->sentinelSA[zcount] + offset;
                sa_rev.isa = isa;
                reverseBitBase(sa_rev.s, sa_rev.len);
                return sa_rev ;
            }
        }
        bp = isa - zcount;
        c = bwt_B0(bwt, bp);
        isa = bwt->L2[c] + bwt->occMajor[bp/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp, c) + bwt->divideNumber - 1;
        if(sa_rev.len >= sa_rev.max)
        {
            sa_rev.max <<= 1;
            sa_rev.s = (uint8_t*)xrecalloc(sa_rev.s, sa_rev.max * sizeof(uint8_t));
        }
        sa_rev.s[sa_rev.len++] = c;
        offset++;
    }
    sa_rev.sa = bwt->sa[isa / SA_INTERVAL] + offset ;
    bwt->sa[isa / SA_INTERVAL]  |= SA_SET_ACCESS_FLAG ;
    sa_rev.isa = isa ;
    reverseBitBase(sa_rev.s, sa_rev.len);
    return sa_rev ;
}

static SA_REV getBaseAndSetSAFlag(const bwt_t *bwt, int64_t isa, const int len)
{
    SA_REV sa_rev ;
    sa_rev.len = 0 ; sa_rev.max = len ;
    sa_rev.s = (uint8_t*)xcalloc(sa_rev.max, sizeof(uint8_t));
    int offset = 0;
    for(int i = 0; i < len; i++)
    {
        int64_t bp ;
        int c = 0;
        bp = posSA2PosBWT(bwt, isa);
        c = bwt_B0(bwt, bp);
        isa = bwt->L2[c] + bwt->occMajor[bp/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp, c) + bwt->divideNumber - 1;
        sa_rev.s[sa_rev.len++] = c;
        offset++;
        if(isa % SA_INTERVAL == 0) 
        {
            bwt->sa[isa / SA_INTERVAL] |= SA_SET_ACCESS_FLAG;
        }
    }
    sa_rev.isa = isa ;
    reverseBitBase(sa_rev.s, sa_rev.len);
    return sa_rev ;
}

static inline void readAlignCore(const uint8_t *ref, const int ref_len, const uint8_t *qur, const int qur_len, AlignLimit alignLimit, read_info *ri)
{
    alignLimit.max_indel_num = ri->read_len * alignLimit.max_indel_num / 100;
    alignLimit.max_mismatch = ri->read_len * alignLimit.max_mismatch /100;
    ri->score = 0XF;
    int cigar_index = MAX_UNMAP_REGION_ALLOW - 2 ;
    int match_len = 0;
    for(int i =qur_len -1, j = ref_len -1; i >= 0, j >= 0; )
    {
        if(qur[i] == ref[j]) { match_len++; i--; j--;}
        else {
            if(ri->read_type == 1 ) // Illumina read
            {
                if(ri->mismatch_num >= alignLimit.max_mismatch)
                {
                    ri->score = 0XF; 
                    return ;
                } else {
                    int flag = 0;
                    int z ;
                    for(z = 1; z < alignLimit.max_mismatch - ri->mismatch_num + 1; z++)
                    {
                        // align to the edge of reference
                        if(j - z < 0) { ri->score = 0XF; return; }
                        else if(i - z < 0 ) {  flag = 1; break;}
                        else if(qur[i -z] == ref[j-z]){ flag = 1; break; }
                    }
                    if(flag == 1)
                    {
                        // write CIGAR
                        if(match_len > 0) { ri->cigar[cigar_index] = match_len ; ri->cigar[cigar_index+1] = 'M'; match_len = 0; cigar_index -= 2; }
                        ri->cigar[cigar_index] = z ; ri->cigar[cigar_index+1] = 'U'; cigar_index -= 2;
                        ri->mismatch_num += z ;
                        i -= z ; j -= z;
                    } else {
                        ri->score = 0XF; return ;
                    }
                }
            } else if(ri->read_type == 2) { // 454 read
                int flag = 0;
                if(ri->indel_num < alignLimit.max_indel_num)
                {
                    for(int z = 1; z < alignLimit.max_indel_len; z++)
                    {
                        int step ;
                        if(i - z < 0) 
                        {
                            if(match_len > 0) { ri->cigar[cigar_index] = match_len; ri->cigar[cigar_index +1] = 'M'; match_len = 0; cigar_index -= 2; }
                            ri->cigar[cigar_index] = z  ; ri->cigar[cigar_index+1] = 'I'; cigar_index -= 2;
                            ri->indel_num++;
                             flag = 1; break; 
                        }
                        if(j - z < 0) break;
                        if(i - z - alignLimit.min_step < 0 || j - z - alignLimit.min_step < 0)
                        {
                            step = i - z - alignLimit.min_step < j - alignLimit.min_step ? i -z : j - z;
                        } else { step = alignLimit.min_step; }
                        if(memcmp(&qur[i - z - step], &ref[j - step], step * sizeof(uint8_t)) == 0)
                        {
                            if(match_len > 0) { ri->cigar[cigar_index] = match_len; ri->cigar[cigar_index + 1] = 'M'; match_len = 0; cigar_index -= 2; }
                            ri->cigar[cigar_index] = z ; ri->cigar[cigar_index+1] = 'D'; cigar_index -= 2;
                            ri->indel_num++;
                            i -= (z + step); j -= step;
                            flag = 1 ; break;
                        } else if(memcmp(&qur[i - step], &ref[j - z - step], step * sizeof(uint8_t)) == 0) {
                            if(match_len > 0) { ri->cigar[cigar_index] = match_len; ri->cigar[cigar_index +1] = 'M'; match_len = 0; cigar_index -= 2; }
                            ri->cigar[cigar_index] = z ; ri->cigar[cigar_index+1] = 'I'; cigar_index -= 2;
                            ri->indel_num++;
                            i -= step; j -= (z + step);
                            flag = 1; break;
                        } else if(ri->mismatch_num + z < alignLimit.max_mismatch && (memcmp(&qur[i -z - step], &ref[j -z - step], step * sizeof(uint8_t)) == 0)) {
                            if(match_len > 0) { ri->cigar[cigar_index] = match_len; ri->cigar[cigar_index+1] = 'M'; match_len = 0; cigar_index -= 2; }
                            ri->cigar[cigar_index] = z ; ri->cigar[cigar_index+1] = 'U'; cigar_index -= 2;
                            ri->mismatch_num++;
                            i -= z ; j -= z;
                            flag = 1; break;
                        }
                    }
                    if(flag == 0)
                    {
                        ri->score = 0XF; return ;
                    }
                }
            } else { // other type read
                fprintf(stderr, "[readAlignCore] unrecognized read type: %u, exit...\n ", ri->read_type);
                exit(1);
            }
        }
    }
    if(match_len > 0) { ri->cigar[cigar_index] = match_len; ri->cigar[cigar_index+1] = 'M'; cigar_index -= 2;}
    ri->cigar_count = MAX_UNMAP_REGION_ALLOW - (cigar_index + 2);
    memmove(ri->cigar, &(ri->cigar[cigar_index + 2]), ri->cigar_count * sizeof(uint16_t)); 
}

static inline int chooseBestMatch(const read_info *ri, const int num)
{
    int min_score = ri[0].score, index = 0;
    for(int i = 1; i < num; i++)
    {
        if(ri[i].score < min_score)
        {
            min_score = ri[i].score;
            index = i;
        }
    }
    return index;
}

static inline void reverseSeq_region(Seq_region *sr)
{
    for(int i =0; i < sr->chance_num; i++)
    {
        reverseBitBase(sr->s[i], sr->len);
        for(int j =0; j < sr->len; j++)
        {
            sr->s[i][j] = (~sr->s[i][j]) & 0x3;
        }
    }
}

int map2Unitig(const Arguments *arguments, const faElem_t *fe, const int start, const int end, omp_lock_t *fpLock, const bwt_t *bwt, const lib_info *libInfo, FILE *readsMap_fp)
{
    int step_len = STEP_LEN ;
    int kmerLenByWord = (KMER + 32 -1)/ 32;
    int max_read_len = MAX_READ_LEN;
    int map_number = 0;
    AlignLimit alignLimit ;
    alignLimit.max_indel_num = MAX_INDEL_NUM; alignLimit.max_indel_len = MAX_INDEL_LEN;
    alignLimit.max_mismacth = MAX_MISMACTH; alignLimit.min_flank = MIN_FLANK;
    for(int i = start; i < end; i++)
    {
        int last_step = step_len;
        ISABase isaBase; 
        isaBase.count = 0; isaBase.max = 1000;
        isaBase.baseIndex = (BaseISA*)xcalloc(isaBase.max, sizeof(BaseISA));
		MapInfo mapInfo;
		mapInfo.max = INITIAL_STACK_SIZE; mapInfo.count = 0;
		mapInfo.mlinfo = (MapLineInfo*)xcalloc(mapInfo.max, sizeof(MapInfo));
        for(int j = step_len; j < fe[i].seq.locPos - KMER; )
        {
            uint64_t k[kmerLenByWord];
            memset(k, 0, kmerLenByWord * sizeof(uint64_t));
            int count = 0;
            int end_position = 0;
            while((end_position = getBinaryKmer(k, fe[i].seq.s, j, fe[i].seq.locPos, count)) > 0)
            {
                Region region = calculateRegion(bwt, k, KMER);
                MapBase mapBase ;
                mapBase.count = 0; mapBase.max = 50;
                mapBase.ri = (read_info*)xcalloc(mapBase.max, sizeof(read_info));
                Seq_region srb = getSeqRegion(fe[i].seq.s, j, max_read_len - KMER, fe[i].seq.locPos, BACKWARD);
                // check if identity the sr.len before position j
                for(int64_t z = region.low; z <= region.high; z++)
                {
                    int64_t bp = posSA2PosBWT(bwt, z) ;
                    int len ;
                    if(sr.len -1 - 2 * last_step > 0) len = sr.len - 1 - 2 * last_step;
                    else len = 0;
                    if(checkNeedAddToMapIndex(bp, &srb, len, BACKWARD, bwt, &isaBase) == 1)
                    {
                        if(mapBase.count >= mapBase.max)
                        {
                            mapBase.max <<= 1 ;
                            mapBase.ri = (read_info*)xrecalloc(mapBase.ri, mapBase.max * sizeof(read_info));
                        }
                        mapBase.ri[mapBase.count++].isa = z;
                    }
                }
                // align the part of fronter KMER 
                for(int z = 0; z < mapBase.count; z++)
                {
                    SA_REV sa_rev = calculateSAAndGetBase(bwt, mapBase.ri[z].isa);
                    getReadInfo(libInfo, sa_rev.sa, &mapBase.ri[z], BACKWARD);
                    //int64_t forwardBoun = readBoundaryLookUp(libInfo, sa_rev.sa, FORWARD);
                    if(sa_rev.sa - mapBase.ri[z].read_head_sa > sa_rev.len)
                    {
                        int len = sa_rev.sa - mapBase.ri[z].read_head_sa - sa_rev.len;
                        SA_REV pad = getBaseAndSetSAFlag(bwt, sa_rev.isa, len);
                        if(sa_rev.len + len > sa_rev.max)
                        {
                            sa_rev.max = sa_rev.len + len;
                            sa_rev.s = realloc(sa_rev.s, sa_rev.max * sizeof(uint8_t));
                        }
                        memmove(&sa_rev.s[len], sa_rev.s, sa_rev.len * sizeof(uint8_t));
                        memcpy(sa_rev.s, pad.s, len * sizeof(uint8_t));
                        sa_rev.len += len;
                        // free and clean work
                        free(pad.s);
                    } else {
                        memmove(sa_rev.s, &sa_rev.s[sa_rev.len - (sa_rev.sa - mapBase.ri[z].read_head_sa)], (sa_rev.sa - mapBase.ri[z].read_head_sa) * sizeof(uint8_t));
                        sa_rev.len = sa_rev.sa - mapBase.ri[z].read_head_sa;
                    }
                    // alignment
                    if(srb.len > sa_rev.len)
                    {
                        read_info ri[srb.chance_num];
                        for(int y = 0; y < srb.chance_num; y++)
                        {
                            ri[y] = mapBase.ri[z];
                            readAlignCore(srb.s[y],srb.len, sa_rev.s,sa_rev.len, alignLimit, &ri[y]);         
                        }
                        // choose best match
                        int bestHit = chooseBestMatch(ri, srb.chance_num);
                        if(ri[bestHit].score == 0XF) { mapBase.ri[z].deleted = 1 ; } 
                        else { mapBase.ri[z] = ri[bestHit]; }
                    } else {
                        mapBase.ri[z].deleted = 1 ; 
                    }
                    // free and clean work
                    free(sa_rev.s);
                }
                // delete that have mark DELETED flag
                int y = 0;
                for(int z = 0; z < mapBase.count; z++)
                {
                    if(mapBase.ri[z].deleted == 0)
                    {
                        if(y < z)  mapBase.ri[y] = mapBase.ri[z];
                        y++;
                    }
                }
                mapBase.count = y;
                // align the part of following KMER seed
                uint64_t rk[kmerLenByWord];
                getRevKmer(k, rk, KMER);
                rk[kmerLenByWord -1] <<= ((32 - KMER%32)<<1);
                region = calculateRegion(bwt, rk, KMER);
                Seq_region srf = getSeqRegion(fe[i].seq.s, j + KMER, max_read_len - KMER, fe[i].seq.locPos, FORWARD);
                reverseSeq_region(&srf); 
                for(int64_t z = region.low; z <= region.high; z++)
                {
                    SA_REV sa_rev = calculateSAAndGetBase(bwt, z);
                    int64_t sa = bwt->seq_len * 2 -1 - sa_rev.sa;
                    int flag = 0, basePos;
					// check if need to map other part of read
                    for(int y = 0; y < mapBase.count; y++)
                    {
                        if(sa > mapBase.ri[y].read_head_sa && sa < mapBase.ri[y].read_head_sa + mapBase.ri[y].read_len)
                        {
                            flag = 1; basePos = y; break;
                        }
                    }
                    if(flag == 1)
                    {
                        if(mapBase.ri[basePos].read_head_sa + mapBase.ri[basePos].read_len - sa > sa_rev.len)
                        {
                            int len = mapBase.ri[basePos].read_head_sa + mapBase.ri[basePos].read_len - sa - sa_rev.len;
                            SA_REV pad = getBaseAndSetSAFlag(bwt, sa_rev.isa, len);
                            if(sa_rev.len + len > sa_rev.max)
                            {
                                sa_rev.max = sa_rev.len + len;
                                sa_rev.s = realloc(sa_rev.s, sa_rev.max * sizeof(uint8_t));
                            }
                            memmove(&sa_rev.s[len], sa_rev.s, sa_rev.len * sizeof(uint8_t));
                            memcpy(sa_rev.s, pad.s, len * sizeof(uint8_t));
                            sa_rev.len += len;
                            // free and clean work
                            free(pad.s);
                        } else {
                            memmove(sa_rev.s, &sa_rev.s[sa_rev.len - (mapBase.ri[basePos].read_head_sa + mapBase.ri[basePos].read_len - sa )], (mapBase.ri[basePos].read_head_sa + mapBase.ri[basePos].read_len -sa) * sizeof(uint8_t));
                            sa_rev.len = mapBase.ri[basePos].read_head_sa + mapBase.ri[basePos].read_len - sa ;
                        }
						// alignment the fronter part of read
						if(srf.len > sa_rev.len)
						{
							read_info ri[srf.chance_num];
							AlignLimit localLimit  = alignLimit;
							localLimit.max_indel_num = alignLimit.max_indel_num - mapBase.ri[basePos].indel_num;
							localLimit.max_mismatch = alignLimit.max_mismatch - mapBase.ri[basePos].mismatch_num;
							for(int y = 0; y < srf.chance_num; y++)
							{
								ri[y] = mapBase.ri[basePos];
								readAlignCore(srf.s[y], srf.len, sa_rev.s, sa_rev.len, localLimit, &ri[y]);
							}
							// choose best match
							int bestHit = chooseBestMatch(ri, srb.chance_num);
							if(ri[bestHit].score == 0XF) { mapBase.ri[basePos].deleted = 1; }
							else { 
								mapBase.ri[basePos].indel_num += ri[bestHit].indel_num;
								mapBase.ri[basePos].mismatch_num += ri[bestHit].mismatch_num;
								mapBase.ri[basePos].score += ri[bestHit].score;
								// padding the CIGAR tags
								mapBase.ri[basePos].cigar[mapBase.ri[basePos].cigar_count++] = KMER ; 
								mapBase.ri[basePos].cigar[mapBase.ri[basePos].cigar_count++] = 'M';
								memcpy(mapBase.ri[basePos].cigar + mapBase.ri[basePos].cigar_count, ri[bestHit].cigar, ri[bestHit].cigar_count * sizeof(uint16_t));
								mapBase.ri[basePos].cigar_count += ri[bestHit].cigar_count;
								char mapLine[256];
								memset(mapline, 0, 256);
								for(int y = 0; y < mapBase.ri[basePos].cigar_count; y++)
								{
									// identify odd and even
									if(y % 2 == 0)
									{
										char s[10];
										memset(s, 0, 10);
										itoa(mapBase.ri[basePos].cigar[y], s, 10);
										strcat(mapLine, s);
									} else {
										char s[2];
										s[0] = mapBase.ri[basePos].cigar[y]; s[1] = '\0' ;
										strcat(mapLine, s);
									}
								}
								if(mapBase.ri[basePos].read_head_sa < bwt->seq_len)
								{
									strcat(mapLine, "\t+");
								} else {
									strcat(mapLine, "\t-");
									mapBase.ri[basePos].read_head_sa = bwt->seq_len * 2 - 1 - mapBase.ri[basePos].read_head_sa - (mapBase.ri[basePos].read_len - 1) ;
									mapBase.ri[basePos].pair_head_sa = bwt->seq_len * 2 - 1 - mapBase.ri[basePos].pair_head_sa - (mapBase.ri[basePos].read_len - 1);
									if(mapBase.ri[basePos].pair_head_sa > mapBase.ri[basePos].read_head_sa)
									{
										mapBase.ri[basePos].pair_head_sa = mapBase.ri[basePos].read_head_sa;
									}
								}


								// write map info to the buffer
								if(mapInfo.count >= mapInfo.max)
								{
									mapInfo.max <<= 1 ;
									mapInfo.mlinfo = (MapLineInfo*)xrecalloc(mapInfo.mlinfo, mapInfo.max * sizeof(MapLineInfo));
								}
								int p = mapBase.ri[basePos].pair_head_sa == mapBase.ri[basePos].read_head_sa ? 1 : 2 ;
								sprintf(mapInfo.mlinfo[mapInfo.count].mapLine, "%ld\t%u\t%s\t%d\t%s", mapBase.ri[basePos].pair_head_sa, p,  fe[i].name.s, mapBase.ri[basePos].ref_position, mapLine) ;
							}
						} else {
							mapBase.ri[basePos].deleted = 1 ;
						}

                    }
                    // free and clean work
                    free(sa_rev.s);
                }
                memset(k, 0, kmerLenByWord * sizeof(uint64_t));
                count++;
                // free and clean work
                free(mapBase.ri); 
                Seq_region_destroy(&sr);
            }
            last_step = stepLen2Next(fe[i].seq.s[j], step_len);
            j += last_step ;
            // update ISABase
            int x = 0;
            for(int z = 0; z < isaBase.count; z++)
            {
                if(isaBase.baseIndex[z].distance + last_step < max_read_len - KMER && x < z)
                {
                    isaBase.baseIndex[x].isa = isaBase.baseIndex[z].isa;
                    isaBase.baseIndex[x].distance = isaBase.baseIndex[z].distance;
                    x++;
                }
            }
            isaBase.count = x ;
            if(isaBase.count + (region.high - region.low + 1) > isaBase.max)
            {
                isaBase.max <<= 1;
                isaBase.baseIndex = (BaseISA*)realloc(isaBase.baseIndex, isaBase.max * sizeof(BaseISA));
            }
            for(int z = 0; z < region.high - region.low + 1; z++)
            {
                isaBase.baseIndex[isaBase.count].isa = region.low + z;
                isaBase.baseIndex[isaBase.count++].distance = last_step;
            }
        }
		// write buffer to the file
		omp_set_lock(fpLock);
		fprintf(readsMap_fp, ">%s\t%s\n", fe[i].name.s, fe[i].comment.s );
		fprintf(readsMap_fp,"%s\n", fe[i].seq.s);
		for(int j = 0; j < mapInfo.count; j++)
		{
			fprintf(readsMap_fp, "%s\n", mapInfo.mlinfo[j].mapLine);
		}
        map_number += mapInfo.count;
		omp_unset_lock(fpLock);
        // free and clean work
		free(mapInfo.mlinfo);
        free(isaBase.baseIndex);
    }

    return map_number;
}

void indexReadsMapFile(const FILE *readsMap_fp, FILE *index_fp)
{
	long int pos = 0;
	fseek(readsMap_fp, SEEK_SET);
	while(readsMap_fp != EOF)
	{
		char line[PATH_LEN];
		fgets(line, PATH_LEN, readsMap_fp);
		if(line[0] == '>')
		{
			char *id;
            id = strtok(line, '\t');
			fprintf(index_fp, "%s\t%ld\n", id + 1, pos);
		} 
		pos = ftell(readsMap_fp);
		// free and clean work
	}
}
        
void rebuildUnitigPInfoHash(const UnitigPInfo *upi, UnitigPInfo *newP, const int64_t new_size)
{
    for(int64_t i = 0; i < upi->max; i++)
    {
        if(upi->unitig_pairs[i].pair_num > 0)
        {
            uint32_t v = upi->unitig_pairs[i].unitigID[0] ^ upi->unitig_pairs[i].unitigID[1];
            uint32_t hashF = 0, hashs = 0;
            uint64_t hashAdd ;
            hashlittle2((void*)&v, sizeof(uint32_t), &hashF, &hashS);
            hashAdd = hashF + (((uint64_t)hashS)<<32);
            hashAdd %= new_size ;
            while(1)
            {
                if(newP[hashAdd].pair_num == 0)
                {
                    newP[hashAdd] = upi->unitig_pairs[i];
                } else {
                    hashAdd++;
                    if(hashAdd >= new_size) hashAdd = 0;
                }
            }
        }
    }
}

void addUnitigPInfoHash(UnitigPInfo *upi, const int unitigID1, const int unitigID2)
{
    uint32_t v = unitigID1 ^ unitigID2 ;
    uint32_t hashF = 0, hashS = 0;
    uint64_t hashAdd ;
    hashlittle2((void*)&v, sizeof(uint32_t), &hashF, &hashS);
    hashAdd = hashF + (((uint64_t)hashS)<<32);
    hashAdd %= upi->max;
    while(1)
    {
        if(upi->unitig_pairs[hashAdd].pair_num == 0)
        {
            upi->unitig_pairs[hashAdd].unitigID[0] = unitigID1;
            upi->unitig_pairs[hashAdd].unitigID[1] = unitigID2;
            upi->unitig_pairs[hashAdd].pair_num++;
            upi->count++;
            break;
        } else if((upi->unitig_pairs[hashAdd].unitigID[0] == unitigID1 && upi->unitig_pairs[hashAdd].unitigID[1] == unitigID2) || (upi->unitig_pairs[hashAdd].unitigID[0] == unitigID2 && upi->unitig_pairs[hashAdd].unitigID[1] == unitigID1)) {
            upi->unitig_pairs[hashAdd].pair_num++;
            break;
        } else {
            hashAdd++;
            if(hashAdd >= upi->max) hashAdd = 0;
        }
    }
    // realloc the hash size
    if((double)upi->count / upi->max > 0.8)
    {
        int64_t new_size = upi->max << 1;
        new_size = find_next_prime_kh(new_size); 
        Unitig_pair *newP =  (Unitig_pair*)xcalloc(new_size, sizeof(Unitig_pair));
        rebuildUnitigPInfoHash(upi, newP);
        free(upi->unitig_pairs); upi->unitig_pairs = newP ;
        upi->max = new_size ;
    }
    
}

static inline int getReadLen(char *cigar)
{
    char *t;
    int len = 0;
    for(int i =0; i < strlen(cigar); i++)
    {
        if(cigar[i] >= 'A' && cigar[i] <= 'Z') cigar[i] = '+';
    }
    t = strtok(cigar, '+');
    len += atoi(t);
    while((t = strtok(NULL, '+')) != NULL)
    {
        len += atoi(t);
    }
    return len ;
}

void constructUnitigPairInfo(ReadPos *readsHash, const int64_t readPair_num, UnitigPInfo *upi, FILE *readsMap_fp, lib_info *libInfo)
{
    //statistic insert size
    Insert_info insertInfo[libInfo->num_lib];
    for(int i =0; i < libInfo->num_lib; i++)
    {
        insertInfo[i].count = 0; insertInfo[i].max = INITIAL_STACK_SIZE;
        insertInfo[i].insert_size = (int*)xcalloc(insertInfo[i].max, sizeof(int));
    }
    while(readsMap_fp != EOF)
    {
        char line[PATH_LEN];
        fgets(line, PATH_LEN, readsMap_fp);
        if(line[0] == '>')
        {
            fgets(line, PATH_LEN, readsMap_fp);
        } else {
            char *id = strtok(line, '\t');
            int64_t pos = atol(id);
            id = strtok(NULL, '\t');
            id = strtok(NULL, '\t');
            int unitigID = atoi(id+6);
            id = strtok(NULL, '\t');
            int ref_pos = atoi(id);
            id = strtok(NULL, '\t');
            int read_len = getReadLen(id);
            uint64_t hashAdd;
            uint32_t hashF = 0, hashS = 0;
            hashlittle2((void*)&pos, sizeof(int64_t), &hashF, &hashS);
            hashAdd = hashF + (((uint64_t)hashS)<<32);
            hashAdd %= readPair_num;
            while(1)
            {
                if(readsHash[hashAdd].flag == 0)
                {
                    readsHash[hashAdd].flag = 0;
                    readsHash[hashAdd].pos = pos;
                    readsHash[hashAdd].unitigID = unitigID;
                    readsHash[hashAdd].ref_pos = ref_pos;
                    readsHash[hashAdd].read_len = read_len;
                    break;
                } else if(readsHash[hashAdd].flag == 1 && readsHash[hashAdd].pos == pos)
                {
                    if(readsHash[hashAdd].unitigID == unitigID)
                    {
                        int insert_s;
                        if(ref_pos > readsHash[hashAdd].ref_pos) {insert_s = ref_pos - readsHash[hashAdd].ref_pos + read_len;}
                        else { insert_s = readsHash[hashAdd].ref_pos - ref_pos + readsHash[hashAdd].read_Len; }
                        int lib_loc = getLibID(libInfo, pos);
                        if(insertInfo[lib_loc].count >= insertInfo[lib_loc].max)
                        {
                            insertInfo[lib_loc].max <<= 1;
                            insertInfo[lib_loc].insert_size = (int*)xrecalloc(insertInfo[lib_loc].insert_size, insertInfo[lib_loc].max * sizeof(int));
                        }
                        insertInfo[lib_loc].insert_size[insertInfo[lib_loc].count++] = insert_s;
                    }
                    addUnitigPInfoHash(upi, unitigID, readsHash[hashAdd].unitigID);
                    memset(&readsHash[hashAdd], 0, sizeof(ReadPos));
                    break;
                }else { 
                    hashAdd++;
                    if(hashAdd >= readPair_num) hashAdd= 0;
                }
            }
            
        }

        // free and clean work
    }

    // free and clean work
    for(int i =0; i < libInfo->num_lib; i++)
    {
        //
        adjust_insert_size();
        free(insertInfo[i].insert_size);
    }
}

int compare_pair_num(const void *p, const void *q)
{
    Unitig_pair *p1 = (Unitig_pair*)p, *q1 = (Unitig_pair*)q;
    if(p1->pair_num < q1->pair_num) return -1 ;
    else if(p1->pair_num > q1->pair_num) return 1 ;
    else return 0;
}

void groupUnitigPInfo(UnitigGroup *unitigGroup, const UnitigPInfo *unitigPInfo)
{
    int64_t i = unitigPInfo->max - 1 ;
    while(unitigPInfo->unitig_pairs[i].pair_num > 0)
    {
        Unitig_pair up = unitigPInfo->unitig_pairs[i];
        int j , z;
        for(j = 0; j < unitigGroup->count; j++)
        {
            for(z = 0; z < unitigGroup->unitigPInfo[j].count; z++)
            {
                if(up.unitigID[0] == unitigGroup->unitigPInfo[j].unitig_pairs[z].unitigID[0] ||
                   up.unitigID[0] == unitigGroup->unitigPInfo[j].unitig_pairs[z].unitigID[1] ||
                   up.unitigID[1] == unitigGroup->unitigPInfo[j].unitig_pairs[z].unitigID[0] ||
                   up.unitigID[1] == unitigGroup->unitigPInfo[j].unitig_pairs[z].unitigID[1]   ) {break; }
            }
        }
        if(j == unitigGroup->count)
        {
            unitigGroup->unitigPInfo[j].max = 10 ;
            unitigGroup->unitigPInfo[j].count = 0 ;
            unitigGroup->unitigPInfo[j].unitig_pairs = (Unitig_pair*)xcalloc(unitigGroup->unitigPInfo[j].max, sizeof(Unitig_pair));
            unitigGroup->unitigPInfo[j].unitig_pairs[unitigGroup->unitigPInfo[j].count++] = up;
            unitigGroup->count++;
        } else {
            if(unitigGroup->unitigPInfo[j].count >= unitigGroup->unitigPInfo[j].max)
            {
                unitigGroup->unitigPInfo[j].max <<= 1 ;
                unitigGroup->unitigPInfo[j].unitig_pairs = (Unitig_pair*)xrecalloc(unitigGroup->unitigPInfo[j].unitig_pairs, unitigGroup->unitigPInfo[j].max * sizeof(Unitig_pair));
            }
            unitigGroup->unitigPInfo[j].unitig_pairs[unitigGroup->unitigPInfo[j].count++] = up;
        }

        i--;
    }
}

void restoreUnitigLocTable(FILE *index_fp, UnitigLocTable *ult)
{
    while(index_fp != EOF)
    {
        UnitigLoc ul;
        char line[PATH_LEN], *t;
        fgets(line, PATH_LEN, index_fp);
        if(ult->count >= ult->max) {fprintf(stderr, "[restoreUnitigLocTable] ult->count >= ult->max, exit....\n"); exit(1); }
        t = strtok(line, '\t');
        ul.unitigID = atoi(t);
        t = strtok(NULL, '\n');
        ul.offset = atol(t);
        uint32_t hashF =0, hashS = 0;
        uint64_t hashAdd ;
        hashlittle2((void*)&ul.unitigID, sizeof(int), &hashF, &hashS);
        hashAdd = hashF + (((uint64_t)hashS)<<32);
        hashAdd %= ult->max;
        while(1)
        {
            if(ult->ul[hashAdd].unitigID == 0)
            {
                ult->ul[hashAdd] = ul;
                break;
            } else {
                hashAdd++;
                if(hashAdd >= ult->max) hashAdd = 0;
            }
        }
        ult->count++;
    }
}

int searchULPos(const UnitigLocTable *ult, const int unitigID)
{
    uint32_t hashF = 0, hashS =0;
    uint64_t hashAdd ;
    int rev ;
    hashlittle2((void*)&uintigID, sizeof(int), &hashF, &hashS);
    hashAdd = hashF + (((uint64_t)hashS)<<32);
    hashAdd %= ult->max;
    while(1)
    {
        if(ult->ul[hashAdd].unitigID == unitigID)
        {
            rev = hashAdd;
            break;
        } else {
            hashAdd++;
            if(hashAdd >= ult->max) hashAdd = 0;
        }
    }
    return rev ;
}

/*  *.unitigGroupReadsMap file format
 *  "@unitigGroupID unitig_pair_info(unitig1--pair_num--unitig2)
 *   >unitigID1 comment
 *   unitigID1 sequence
 *   >unitigID2 comment
 *   unitigID2 sequence
 *   .........
 *   Reads map to this unitig group"
 */
void write2Ugrm(FILE *ugrm_fp, const UnitigLocTable *ult, FILE *readsMap_fp, const UnitigGroup *unitigGroup)
{
    int ugID = 1;
    for(int i = 0; i < unitigGroup->count; i++)
    {
        fprintf(ugrm_fp, "@ug%d", ugID);
        int unitigID[unitigGroup->unitigPInfo[i].count];
        int count = 0;
        for(int j = 0; j < unitigGroup->unitigPInfo[i].count; j++)
        {
            Unitig_pair up = unitigGroup->unitigPInfo[i].unitig_pairs[j];
            fprintf(ugrm_fp, "\t%s-%d-%s", up.unitigID[0], up.pair_num, up.unitigID[1]);
            int flag = 0;
            for(int z = 0; z < count; z++)
            {
                if(unitigID[z] == up.unitig_pairs[0]) { flag = 1 ;break; }
            }
            if(flag == 0) { unitigID[count++] = up.unitig_pairs[0]; }
            flag = 0;
            for(int z = 0; z < count; z++)
            {
                if(unitigID[z] == up.unitig_pairs[1]) { flag = 1; break; }
            }
            if(flag == 0) unitigID[count++] = up.unitig_pairs[1];
        }
        fprintf(ugrm_fp, "\n");
        // write unitig sequence
        for(int j = 0; j < count; j++)
        {
            int ul_pos = searchULPos(ult, unitigID[j]);
            fseek(readsMap_fp, ult->ul[ul_pos].offset, SEEK_SET);
            char *line = NULL;
            getline(line, INT_MAX, readsMap_fp);
            fprintf(ugrm_fp, "%s", line); 
            free(line); line == NULL;
            getline(line, INT_MAX, readsMap_fp);
            fprintf(ugrm_fp, "%s", line);
            // free and clean work
            free(line);
        }
        // write reads map info
        for(int j = 0; j < count; j++)
        {
            int ul_pos = searchULPos(ult, unitigID[j]);
            fseek(readsMap_fp, ult->ul[ul_pos].offset, SEEK_SET);
            char line[PATH_LEN];
            // ignore first two line(unitigID and unitig sequence)
            fgets(line, PATH_LEN, readsMap_fp);
            fgets(line, PATH_LEN, readsMap_fp);
            fgets(line, PATH_LEN, readsMap_fp);
            while(line[0] != '>' && readsMap_fp != EOF)
            {
                fprintf(ugrm_fp, "%s", line);
                fgets(line, PATH_LEN, readsMap_fp);
            }
        }
    }
}

static inline int CIGAR2read_len(const char *cigar)
{
    int len = 0;
    char *s = strdup(cigar), *t;
    for(int i =0; i < strlen(s); i++)
    {
        if(s[i] >= 'A' && s[i] <= 'Z') s[i] = '\t';
    }
    t = strtok(s, '\t');
    len += atoi(t);
    while((t = strtok(NULL, '\t')) != NULL)
    {
        len += atoi(t);
    }
    // free and clean work
    free(s);
    
    return len; 
}

int readUnitigGroupInfo(UGInfo *ugInfo, FILE *ugrm_fp)
{
    char *line = NULL ;
    int c; 
    if((c = fgetc(ugrm_fp) == EOF)return 0;
    if(c != "@") return -1 ;
    int len = getline(&line, INT_MAX, ugrm_fp);
    // set unitig group info
    line[len-1] = '\0';
    char *t = strtok(line, '\t');
    ugInfo->ugID = atoi(t+2);
    ugInfo->upi.count = 0; ugInfo->upi.max = 10;
    ugInfo->unitig_pairs = (Unitig_pair*)xcalloc(ugInfo->upi.max, sizeof(Unitig_pair));
    ugInfo->fs.len = 0; ugInfo->fs.max = 10;
    ugInfo->fs.elem = (faElem_t*)xcalloc(ugInfo->fs.max, sizeof(faElem_t));
    ugInfo->mb.count = 0; ugInfo->mb.max = INITIAL_STACK_SIZE;
    ugInfo->mb.ri = (read_info*)xcalloc(ugInfo->mb.max, sizeof(read_info));
    // ignore a '('
    while(t = strtok( NULL, '\t')) != NULL)
    {
        if(ugInfo->upi.count >= ugInfo->upi.max)
        {
            ugInfo->upi.max <<= 1;
            ugInfo->unitig_pairs = (Unitig_pair*)xrecalloc(ugInfo->unitig_pairs, ugInfo->upi.max * sizeof(Unitig_pair));
        }
        char *ele = strtok(t, '-');
        ugInfo->unitig_pairs[ugInfo->upi.count].unitigID[0] = atoi(ele);
        ele = strtok(NULL, '-');
        ugInfo->unitig_pairs[ugInfo->upi.count].pair_num = atoi(ele);
        ele = strtok(NULL, '-');
        ugInfo->unitig_pairs[ugInfo->upi.count].unitigID[1] = atoi(ele);
        ugInfo->upi.count++;
    }
    free(line); line = NULL
    c = fgetc(ugrm_fp);
    while(c != EOF && c != '@')
    {
        if(ugrm_fp != EOF) fseek(ugrm_fp, -1, SEEK_CUR);
        len = getline(&line, INT_MAX, ugrm_fp);
        if(c == '>')
        {
            if(ugInfo->fs.len >= ugInfo->fs.max)
            {
                ugInfo->fs.max <<= 1;
                ugInfo->fs.elem = (faElem_t*)xrecalloc(ugInfo->fs.elem, ugInfo->fs.max * sizeof(faElem_t));
            }
            char *s = strtok(line, '\t');
            ugInfo->fs.elem[ugInfo->fs.len].name.locPos = strlen(s+1);
            ugInfo->fs.elem[ugInfo->fs.len].name.max = ugInfo->fs.elem[ugInfo->fs.len].name.locPos + 1;
            ugInfo->fs.elem[ugInfo->fs.len].name.s = strdup(s);
            s = strtok(NULL, '\n');
            ugInfo->fs.elem[ugInfo->fs.len].comment.locPos = strlen(s);
            ugInfo->fs.elem[ugInfo->fs.len].comment.max = ugInfo->fs.elem[ugInfo->fs.len].comment.locPos + 1;
            ugInfo->fs.elem[ugInfo->fs.len].comment.s = strdup(s);
            free(line); line = NULL;
            len = getline(&line, INT_MAX, ugrm_fp);
            line[strlen(line)] = '\0';
            ugInfo->fs.elem[ugInfo->fs.len].seq.locPos = strlen(line) ;
            ugInfo->fs.elem[ugInfo->fs.len].seq.max = ugInfo->fs.elem[ugInfo->fs.len].seq.locPos + 1;
            ugInfo->fs.elem[ugInfo->fs.len].seq.s = strdup(line);
            ugInfo->fs.len++;
        } else {
            if(ugInfo->mb.count >= ugInfo->mb.max)
            {
                ugInfo->mb.max <<= 1;
                ugInfo->mb.ri = (read_info*)xrecalloc(ugInfo->mb.ri, ugInfo->mb.max * sizeof(read_info));
            }
            char *s = strtok(line, '\t');
            ugInfo->mb.ri[ugInfo->mb.count].pair_head_sa = atol(s);
            s = strtok(NULL, '\t');
            ugInfo->mb.ri[ugInfo->mb.count].order = atoi(s);
            s = strtok(NULL, '\t');
            ugInfo->mb.ri[ugInfo->mb.count].ref_name = strdup(s);
            s = strtok(NULL, '\t');
            ugInfo->mb.ri[ugInfo->mb.count].ref_position = atoi(s);
            s = strtok(NULL, '\t');
            ugInfo->mb.ri[ugInfo->mb.count].cigar = strdup(s);
            s = strtok(NULL, '\t');
            ugInfo->mb.ri[ugInfo->mb.count].strand = atoi(s);
            // set read_len
            ugInfo->mb.ri[ugInfo->mb.count].read_len = CIGAR2read_len(ugInfo->mb.ri[ugInfo->mb.count].cigar);
            ugInfo->mb.count++;
        }
        c = fgetc(ugrm_fp);
        // free and clean work
        free(line); line = NULL;
    }
    if(ugrm_fp != EOF)  { fseek(ugrm_fp, -1, SEEK_CUR);}
    return 1 ;
}

static inline int searchUnitigLength(const faSeq_t *fs, const int unitigID)
{
	int len = 0;
	for(int i = 0; i < fs->len; i++)
	{
		if(unitigID == atoi(fs->elem[i].name.s+6))
		{
			len = fs->elem[i].seq.locPos;
			break;
		}
	}

	return len;
}

/* about how to estimate gap length, more info please cite "High quality draft assemblies of 
 * mammalian genomes from massively parallel sequence data supproting information"
 */
UnitigFlatten estimateGap(UnitigPairGapInfo *upgi, const faSeq_t *fs)
{
    UnitigFlatten uf;
    uf.unitigID[0] = upgi->unitigID[0];
    uf.unitigID[1] = upgi->unitigID[1];

	// searchUnitigLength
	if((uf.len[0] = searchUnitigLength(fs, uf.unitigID[0])) <= 0) 
	{
		fprintf(stderr, "[UntigFlatten] the length of unitig is error, program exit..,,\n");
		exit(1);
	}
	if((uf.len[1] = searchUnitigLength(fs, uf.unitigID[1])) <= 0)
	{
		fprintf(stderr, "[UnitigFlatten] the length of unitig is error, program exit...\n");
		exit(1);
	}
    
    // estimate SD
    double sd = 0;
    for(int i =0; i < upgi->count; i++)
    {
        sd += 1 / pow(upgi->eg[i].SD, 2);
    }
    sd = sqrt(1 / sd);
    uf.SD = (int)sd;

    // estimate gap_len
    double gap_len = 0;
    for(int i =0; i < upgi->count; i++)
    {
        gap_len += (upgi->eg[i].gap_len / pow(upgi->eg[i].SD), 2);
    }
    gap_len = gap_len / pow(sd, 2);
    uf.gap_len = (int)gap_len;

    // compute score for gap estimated
    double score = 0;
    for(int i = 0; i < upgi->count; i++)
    {
        int k = (upgi->eg[i].gap_len - gap_len) / upgi->eg[i].SD;
        score += pow(k, 2);
    }
    score = sqrt(score/n);
    uf.score = score;
    
    return uf;
}

int compareScore(const void *p, const void *q)
{
    UnitigFlatten *p1 = (UntigFlatten*)p, *q1 = (UnitigFlatten*)q;
    if(p1->score < q1->score) return -1;
    else if(p1->score > q1->score) return 1;
    else return 0;
}

static inline int findNearestUnitig(const OrderUnitigInfo *oui, const int order, const int direction)
{
    int index = -1 ;
    if(direction == FORWARD)
    {
		int new_order = order + 1; 
        for(int i = 0; i < oui->count; i++)
        {
            if(oui->oe[i].order == new_order )
            {
				index = i;
                break;
            }
        }
    } else if(direction == BACKWARD) {
        int new_order = order - 1;
        for(int i = 0; i < us->count; i++)
        {
            if(oui->oe[i].order == new_order)
            {
                index = i;
                break;
            }
        }
    } else {
        fprintf(stderr, "[findNearestUnitig] the direction unrecognized, program exit...\n");
        exit(1);
    }
    return index ;
}

static inline int searchWithInUnitigFlatten(const OrderUnitigInfo *oui, const int unitigID)
{
	int index = -1;
	for(int i = 0; i < oui->count; i++)
	{
		if(oui->oe[i].unitigID == unitigID)
		{
			index = i;
			break;
		}
	}
	return index;
}
                            
static inline void updateOrder(OrderUnitigInfo *oui, const int order, const int direction)
{
    if(direction == FORWARD)
    {
        for(int i = 0; i < oui->count; i++)
        {
            if(oui->oe[i].order >= order) { oui->oe[i].order++; }
        }
    } else if(direction == BACKWARD) {
        for(int i = 0; i < oui->count; i++)
        {
            if(oui->oe[i].order <= order) {oui->oe[i].order--; }
        }
    }
}

int compareOrder(const void *p, const void *q)
{
    OrderElem *p1 = (OrderElem*)p, *q1 = (OrderElem*)q;
    if(p1->order < q1->order) return -1;
    else if(p1->order > q1->order) return 1 ;
    else return 0 ;
}
    
static OrderUnitigInfo orderUnitigFlatten(UnitigScaff *us, const UGInfo *ugInfo)
{
    OrderUnitigInfo oui;
    oui.count = 0; oui.max = us->count;
    oui.oe = (OrderElem*)xcalloc(oui.max, sizeof(OrderElem));
    for(int i =0; i < oui.max; i++) { oui.oe[i].order = INT_MIN;}
     
    oui.oe[0].unitigID = us->uf[0].unitigID[0];
	oui.oe[0].lenth = us->uf[0].len[0];
    //oui.oe[0].gap_len = us->uf[0].gap_len;
    oui.oe[0].order = 1 ;
    oui.count++;
    int loc = 0 ;
	while(1)
	{
		for(int i = 0; i < us->count; i++)
		{
			int uid_1, uid_2;
			if((uid_1 = searchWithInUnitigFlatten(&oui, us->uf[i].unitigID[0])) >= 0 && (uid_2 = searchWithInUnitigFlatten(&oui,  us->uf[i].unitigID[1])) == -1)
			{
				if(oui.oe[uid_1].gap_len != 0)
				{
                    int index, order = oui.oe[uid_1].order, gap_len = us->uf[i].gap_len ;
                    while(1)
                    {
					    index = findNearestUnitig(&oui, order, gap_len, FORWARD);
                        if(index == -1)
                        {
                            oui.oe[oui.count].unitigID = us->uf[i].unitigID[1];
                            oui.oe[oui.count].len = us->uf[i].len[1];
                            oui.oe[oui.count].order = order + 1;
                            oui.oe[uid_1].gap_len = gap_len;
                            oui.count++;
                            break;
                        } else if(oui.oe[uid_1].gap_len >= gap_len) {
                            updateOrder(&oui, oui.oe[index].order, FORWARD);
                            oui.oe[oui.count].unitigID = us->uf[i].unitigID[1];
                            oui.oe[oui.count].len = us->uf[i].len[1];
                            oui.oe[oui.count].order = oui.oe[index].order ;
                            oui.oe[oui.count].gap_len = oui.oe[uid_1].gap_len - gap_len - oui.oe[oui.count].len;
                            oui.oe[uid_1].gap_len = gap_len;
                            oui.count++;
                            break;
                        } else {
                            order = oui.oe[index].order;
                            gap_len = gap_len - oui.oe[uid_1].gap_len - oui.oe[index].len; 
                            uid_1 = index;
                        }
                    }
				} else { 
					oui.oe[oui.count].unitigID = us->uf[i].unitigID[1];
					oui.oe[oui.count].len = us->uf[i].len[1];
					oui.oe[oui.count].order = oui.oe[uid_1].order + 1;
					oui.oe[uid_1].gap_len = us->uf[i].gap_len;
					oui.count++;
				}
                us->uf[i].SD = INT_MIN;
			}else if((uid_1 = searchWithInUnitigFlatten(&oui, us->uf[i].unitigID[0])) == -1 && (uid_2 = searchWithInUnitigFlatten(&oui, us->uf[i].unitigID[1])) >= 0){
		        if(oui.oe[uid_1].gap_len != 0)
                {
                    int index, order = oui.oe[uid_2].order, gap_len = us->uf[i].gap_len;
                    while(1)
                    {
                        index = findNearestUnitig(&oui, order, gap_len, BACKWARD);
                        if(index = -1)
                        {
                            oui.oe[oui.count].unitigID = us->uf[i].unitigID[0];
                            oui.oe[oui.count].len = us->uf[i].len[0];
                            oui.oe[oui.count].order = order - 1;
                            oui.oe[oui.count].gap_len = gap_len;
                            oui.count++; break;
                        } else if(oui.oe[index].gap_len >= gap_len) {
                            updateOrder(&oui, oui.oe[index].order, BACKWARD);
                            oui.oe[oui.count].unitigID = us->uf[i].unitigID[0];
                            oui.oe[oui.count].len = us->uf[i].len[0];
                            oui.oe[oui.count].order = oui.oe[index].order;
                            oui.oe[oui.count].gap_len = gap_len;
                            oui.count++; break;
                        } else {
                            order = oui.oe[index].order;
                            gap_len = gap_len - oui.oe[index].gap_len - oui.oe[index].len;
                            uid_2 = index;
                        }
                    }
                } else {
                    oui.oe[oui.count].unitigID = us->uf[i].unitigID[0];
                    oui.oe[oui.count].len = us->uf[i].len[0];
                    oui.oe[oui.count].order = oui.oe[uid_2].order - 1;
                    oui.oe[oui.count].gap_len = us->uf[i].gap_len;
                    oui.count++;
                }
                us->uf[i].SD = INT_MIN;
			} else {continue; }
		}
        // check have been finished
        int finished_flag = 0;
        for(int i = 0; i <  us->count; i++)
        {
            if(us->uf[i].SD == INT_MIN) finished_flag++;
        }
        if(finished_flag >= us->count) break;
	}
    // sort the order
    qsort(oui.oe, oui.count, sizeof(OrderElem), compareOrder);
    // free and clean work

    return oui;
}

UnitigScaff unitigPairInfoFlatten(const UGInfo *ugInfo, const lib_info *libInfo)
{
    UnitigScaff us;
    us.count = 0; us.max = ugInfo->upi.count;
    us.uf = (UnitigFlatten*)xcalloc(us.max, sizeof(UnitigFlatten));
    UnitigPairGapInfo upgi[ugInfo->upi.count]; memset(upgi, 0, ugInfo->upi.count * sizeof(UnitigPairGapInfo));
    int upgi_count = 0;
    // conform the unitig sequence from same strand of genome
    int flag[ugInfo->fs.len]; memset(flag, 0, ugInfo->fs.len);
    flag[0] = 1;
    for(int i = 1; i < ugInfo->mb.count; i++)
    {
        for(int j = 0; j < i; j++)
        {
            if(ugInfo->mb.ri[i].pair_head_sa == ugInfo->mb.ri[j].pair_head_sa)
            {
                int si_1 = searchSeqIndexByName(ugInfo->mb.ri[i].ref_name,libInfo);
                int si_2 = searchSeqIndexByName(ugInfo->mb.ri[j].ref_name,libInfo);
                if(si_1 == 0)
                {
                    if(ugInfo->mb.ri[i].strand == ugInfo->mb.ri[j].strand && ugInfo->fs[si_2].strand == 0)
                    {
                        ugInfo->fs[si_2].strand = 1;
                        reverseBitBase(ugInfo->fs[si_2].seq.s, ugInfo->fs[si_2].seq.locPos); 
                    }
                    flag[si_2] = 1;
                } else if(si_2 == 0) {
                    if(ugInfo->mb.ri[i].strand == ugInfo->mb.ri[j].strand && ugInfo->fs[si_1].strand == 0)
                    {
                        ugInfo->fs[si_1].strand = 1;
                        reverseBitBase(ugInfo->fs[si_1].seq.s, ugInfo->fs[si_1].seq.locPos);
                    }
                    flag[si_1] = 1;
                } else {
                    if(flag[si_1] == 1 && flag[si_2] == 0)
                    {
                        if(ugInfo->mb.ri[i].strand != ugInfo->mb.ri[j].strand)
                        {
                            ugInfo->fs[si_2].strand = 1;
                            reverseBitBase(ugInfo->fs[si_2].seq.s, ugInfo->fs[si_2].seq.locPos);
                        }
                        flag[si_2] = 1;

                    } else if(flag[si_1] == 0 && flag[si_2] == 1) {
                        if(ugInfo->mb.ri[i].strand != ugInfo->mb.ri[j].strand)
                        {
                            ugInfo-.fs[si_1].strand = 1;
                            reverBitBase(ugInfo->fs[si_1].seq.s, ugInfo->fs[si_1].seq.locPos);
                        }
                        flag[si_1] = 1;
                    } else { continue;}
                }
            }
        }
    }
    // check flag
    for(int i = 0; i < ugInfo->fs.len; i++)
    {
        if(flag[i] != 1) 
        {
            fprintf(stderr, "[unitigPairInfoFlatten] the flag set error, program eixt....\n");
            exit(1);
        }
    }
    // computing every pair reads gap_len and SD
    for(int i = 1; i < ugInfo->mb.count; i++)
    {
        for(j = 0; j < i; j++)
        {
            if(ugInfo->mb.ri[i].pair_head_sa == ugInfo->mb.ri[j].pair_head_sa)
            {
                LIB lib = searchLIBBypair_head_sa(ugInfo->mb.ri[i].pair_head_sa, libInfo);
                faElem_t fe_i = searchSeqByName(ugInfo->mb.ri[i].ref_name, ugInfo->fs);
                faElem_t fe_j = searchSeqByName(ugInfo->mb.ri[j].ref_name, ugInfo->fs);
                if(fe_i.strand == 1)
                {
                    ugInfo->mb.ri[i].strand = ugInfo->mb.ri[i].strand == 1 ? 0 : 1;
                    ugInfo->mb.ri[i].ref_position = fe_i.seq.locPos - ugInfo->mb.ri[i].ref_position - ugInfo->mb.ri[i].read_len;
                } 
                if(fe_j.strand == 1)
                {
                    ugInfo->mb.ri[j].strand = ugInfo->mb.ri[j].strand == 1 ? 0 : 1 ;
                    ugInfo->mb.ri[j].ref_position = fe_j.seq.locPos - ugInfo->mb.ri[j].ref_position - ugInfo->mb.ri[j].read_len;
                }
                if(ugInfo->mb.ri[i].strand == ugInfo->mb.ri[j].strand)
                {
                    fprintf(stderr, "[unitigPairInfoFlatten] the strand of MapBase error, program eixt....\n");
                    exit(1);
                }
                int n = -1;
                for(int z = 0; z < upgi_count; z++)
                {
                    if((upgi[z].unitigID[0] == atoi(ugInfo->mb.ri[i].ref_name + 6) && 
                        upgi[z].unitigID[1] == atoi(ugInfo->mb.ri[j].ref_name + 6)) ||
                       (upgi[z].unitigID[0] == atoi(ugInfo->mb.ri[j].ref_name + 6) &&
                        upgi[z].unitigID[1] == atoi(ugInfo->mb.ri[i].ref_name + 6)))
                    {
                        n = z ;
                        break;
                    }
                }
                if(lib.pe.reverse == 1) // jumping library
                {
                    int first = -1, second = -1;
                    if(ugInfo->mb.ri[i].strand == 1) {first = i; second = j; }
                    else { first = j; second = i; faElem_t t = fe_i; fe_i = fe_j; fe_j = t;}
                    if(n >= 0 && ugInfo->mb.ri[first].strand == 1)
                    {
                        if(upgi[n].unitigID[0] != atoi(ugInfo->mb.ri[first].ref_name +6))
                        {
                            fprintf(stderr, "[unitigPairInfoFlatten] the order of UnitigPairGapInfo error, ignore...\n");
                            continue;
                        }
                    }
                    EstimateGap eg;
                    eg.gap_len = lib.pe.insert_size - ugInfo->mb.ri[first].ref_position - (ugInfo->mb.ri[second].ref_position + ugInfo->mb.ri[second].read_len);
                    eg.SD = lib.pe.insert_SD;
                    if(n >= 0)
                    {
                        if(upgi[n].count >= upgi[n].max)
                        {
                            upgi[n].max <<= 1;
                            upgi[n].eg = (EstimateGap*)recalloc(upgi[n].eg, upgi[n].max * sizeof(EstimateGap));
                        }
                        upgi[n].eg[upgi[n].count] = eg;
                        upgi[n].count++;
                    } else {
                        int d = upgi_count;
                        upgi[d].count = 0; upgi[d].max = 20;
                        upgi[d].unitigID[0] = atoi(ugInfo->mb.ri[first].ref_name + 6);
                        upgi[d].unitigID[1] = atoi(ugInfo->mb.ri[second].ref_name + 6);
                        upgi[d].eg = (EstimateGap*)xcalloc(upgi[d].max, sizeof(EstimateGap));
                        upgi[d].eg[upgi[d].count] = eg;
                        upgi[d].count++;
                        upgi_count++;
                    }
                } else { // fragment library
                    int first = -1, second = -1;
                    if(ugInfo->mb.ri[i].strand == 0) { first = i; second = j; }
                    else { first = j; second = i; faElem t = fe_i; fe_i = fe_j; fe_j = t;}
                    if(n >= 0 && ugInfo->mb.ri[first].strand == 0)
                    {
                        if(upgi[n].unitigID[0] != atoi(ugInfo->mb.ri[first].ref_name + 6))
                        {
                            fprintf(stderr, "[unitigPairInfoFlatten] the order of UnitigPairGapInfo error, program exit...\n");
                            exit(1);
                        } 
                    }
                    EstimateGap eg;
                    eg.gap_len = lib.pe.insert_size - (fe_i.seq.locPos - ugInfo->mb.ri[first].ref_position) - ugInfo->mb.ri[second].ref_position;
                    eg.SD = lib.pe.insert_SD;
                    if(n >= 0)
                    {
                        if(upgi[n].count >= upgi[n].max)
                        {
                            upgi[n].max <<= 1;
                            upgi[n].eg = (EstimateGap*)recalloc(upgi[n].eg, upgi[n].max * sizeof(EstimateGap));
                        }
                        upgi[n].eg[upgi[n].count] = eg;
                        upgi[n].count++;
                    } else {
                        upgi[upgi_count].count = 0; upgi[upgi_count].max = 20;
                        upgi[upgi_count].unitigID[0] = atoi(ugInfo->mb.ri[first].ref_name + 6);
                        upgi[upgi_count].unitigID[1] = atoi(ugInfo->mb.ri[second].ref_name + 6);
                        upgi[upgi_count].eg = (EstimateGap*)xcalloc(upgi[upgi_count].max, sizeof(EstimateGap));
                        upgi[upgi_count].eg[upgi[upgi_count].count] = eg;
                        upgi[upgi_count].count++;
                        upgi_count++;
                    }
                }
                break;
            }
        }
    }
    // estimate every gap between unitigs pair 
    for(int i = 0; i < upgi_count; i++)
    {
        us.uf[i] = estimateGap(&upgi[i], ugInfo->fs);
    }
    qsort(&us, upgi_count, sizeof(UntigFlatten), compareScore);
    OrderUnitigInfo oui = orderUnitigFlatten(&us, ugInfo);

    // free and clean work
    free(us.uf);

    return oui;
}

static inline int findMinLibInsert(const lib_info *libInfo)
{
    int minInsert = INT_MAX;
    for(int i = 0; i < libInfo->num_lib; i++)
    {
        if(libInfo->lib[i].paired == 1 && libInfo->lib[i].pe.insert_size < minInsert)
        {
            minInsert = libInfo->lib[i].pe.insert_size;
        }
    }

    return minInsert;
}

static inline findSeqByID(const faSeq_t *fs, const int unitigID)
{
    int index = -1 ;
    for(int i = 0; i < fs->len; i++)
    {
        if(atoi(fs->elem[i].name.s + 6) == unitigID)
        {
            index = i;
            break;
        }
    }
    return index;
}
static void setUnitigQuality(faElem_t *fe, const int begin, const int end)
{
    for(int i = begin; i < end; i++)
    {
        if(fe->seq.s[i] != '{')
        {
            fe->qual.s[i] = 3 ;
        } else {
            fe->qual.s[i] = '{';
            i++;
            while(i < end && fe->seq.s[i] != '}')
            {
                if(fe->seq.s[i] == ',')
                {
                    fe->qual.s[i] = ',';
                } else {
                    fe->qual.s[i] = 2;
                }
                i++;
            }
            if(fe->seq.s[i] == '}') fe->qual.s[i] = '}';
        }
    }
}

static ReadPool mapFlankContig(faElem_t *fe, const int begin, UGInfo *ugInfo, const bwt_t *bwt, const lib_info *libInfo, const int step_len, ISABase *isaBase)
{
    ReadPool rp;
    rp.count = 0; rp.max = 20;
    int min_lib_insert_size = findMinLibInsert(libInfo);
    int max_read_len = MAX_READ_LEN;
    rp.rs = (ReadSeq*)xcalloc(rp.max, sizeof(ReadSeq));
    int kmerLenByWord = (KMER + 32 -1) / 32;
    uint64_t k[kmerLenByWord];
    memset(k, 0, kmerLenByWord * sizeof(uint64_t));
    int count = 0;
    int end_position = 0;
    while((end_position = getBinaryKmer(k, fe->seq.s, begin, fe->seq.locPos, count)) > 0)
    {
        Region region = calculateRegion(bwt, k, KMER);
        MapBase mapBase;
        mapBase.count = 0; mapBase.max = 50;
        mapBase.ri = (read_info*)xcalloc(mapBase.max, sizeof(read_info));
        Seq_region srb = getSeqRegion(fe->seq.s, begin, max_read_len - KMER, fe->seq.locPos, BACKWARD);
        if(rp.max < rp.count + (region.high - region.low +1))
        {
            rp.max = (rp.count + (region.high - region.low +1)) <<1;
            rp.rs = (ReadSeq*)realloc(rp.rs, rp.max * sizeof(ReadSeq));
            memset(rp.rs + rp.count, 0, (rp.max - rp.count) * sizeof(ReadSeq));
        }
        for(int64_t i = region.low; i <= region.high; i++)
        {
            int64_t bp = posSA2PosBWT(bwt, i);
            int len;
            if(sr.len - 1 - 2 * step_len > 0) len = sr.len - 1 - 2 * step_len;
            else len = 0;
            if(checkNeedAddToMapIndex(bp, &srb, len, BACKWARD, bwt, isaBase) == 1)
            {
                rp.rs[rp.count].pair_head_sa = i;
            } 
        }
        // align the part of fronter KMER
        for(int i = 0; i < rp.count; i++)
        {
            SA_REV sa_rev = calculateSAAndGetBase(bwt, rp.rs[i].pair_head_sa);
            getReadInfo(libInfo, sa_rev.sa, &mapBase.ri[i], BACKWARD);
            rp.rs[i].order = mapBase.ri[i].order;
            if(end_position > 2 * min_lib_insert_size)
            {
                if(checkPairInfo(ugInfo->fs, sa_rev.sa) == 1)
                {
                    rp.rs[i].flag = 1 ;
                } else { rp.rs[i].flag = 2; }
            }
            if(sa_rev.sa - mapBase.ri[i].read_head_sa > sa_rev.len)
            {
                int len = sa_rev.sa - mapBase.ri[i].read_head_sa - sa_rev.len;
                SA_REV pad = getBaseAndSetSAFlag(bwt, sa.rev.isa, len);
                if(sa_rev.len + len > sa_rev.max)
                {
                    sa_rev.max = sa_rev.len + len;
                    sa_rev.s = realloc(sa_rev.s, sa_rev.max * sizeof(uint8_t));
                }
                memmove(&sa_rev.s[len], sa_rev.s, sa_rev.len * sizeof(uint8_t));
                memcpy(sa_rev.s, pad.s, len * sizeof(uint8_t));
                sa_rev.len += len;
                // free and clean work
                free(pad.s);
            } else {
                memmove(sa_rev.s, &sa_rev.s[sa_rev.len - (sa_rev.sa - mapBase.ri[i].read_head_sa)], (sa_rev.sa - mapBase.ri[i].read_head_sa) * sizeof(uint8_t));
                sa_rev.len = sa_rev.sa - mapBase.ri[z].read_head_sa;
            }
            // alignment
            if(srb.len > sa_rev.len)
            {
                read_info ri[srb.chance_num];
                alignLimit al ; 
                al.max_indel_num = alignLimit.max_indel_num * (sa_rev.len + KMER) / mapBase.ri[i].read_len;
                al.max_mismatch = alignLimit.max_mismatch * (sa_rev.eln + KMER) / mapBase.ri[i].read_len;
                for(int y = 0; y < srb.chance_num; y++)
                {
                    ri[y] = mapBase.ri[i];
                    readAlignCore(srb.s[y], srb.len, sa_rev.s, sa_rev.len, al, &ri[y]);
                }
                // choose best match
                int bestHit = chooseBestMatch(ri, srb.chance_num);
                if(ri[bestHit].score == 0XF) { mapBase.ri[i].deleted = 1; }
                else { mapBase.ri[i] = ri[bestHit]; }
            } else {
                // free and clean work
                free(sa_rev.s);
            }
            // delete that have mark DELETED flag
            int y = 0;
            for(int z = 0; z < mapBase.count; z++)
            {
                if(mapBase.ri[z].deleted == 0)
                {
                    if(y < z) mapBase.ri[y] = mapBase.ri[z];
                    y++;
                }
            }
            mapBase.count = y;
            // get the part of following KMER seed
            uint64_t rk[kmerLenByWord];
            getRevKmer(k, rk, KMER);
            SA_REV *sa_rev_array;
            int sa_rev_c = 0, sa_rev_m = 20;
            sa_rev_array = (SA_REV*)xcalloc(sa_rev_m, sizeof(SA_REV));
            rk[kmerLenByWord -1] <<= ((32 - KMER%32)<<1);
            region = calculateRegion(bwt, rk, KMER);
            for(int64_t z = region.low; z <= region.high; z++)
            {
                SA_REV sa_rev = calculateSAAndGetBase(bwt, z);
                int64_t sa = bwt->seq_len * 2 - 1 - sa_rev.sa;
                int flag = 0; basePos;
                // check if need to map other part of read
                for(int y = 0; y < mapBase.count; y++)
                {
                    if(sa > mapBase.ri[y].read_head_sa && sa < mapBase.ri[y].read_head_sa + mapBase.ri[y].read_len)
                    {
                        flag = 1; basePos = y; break;
                    }
                }
                if(flag == 1)
                {
                    if(mapBase.ri[basePos].read_head_sa + mapBase.ri[basePos].read_len - sa > sa_rev.len)
                    {
                        int len = mapBase.ri[basePos].read_head_sa + mapBase.ri[basePos].read_len - sa - sa_rev.len;
                        SA_REV pad = getBaseAndSetSAFlag(bwt, sa_rev.isa, len);
                        if(sa_rev.len + len > sa_rev.max)
                        {
                            sa_rev.max = sa_rev.len + len;
                            sa_rev.s = realloc(sa_rev.s, sa_rev.max * sizeof(uint8_t));
                        }
                        memmove(&sa_rev.s[len], sa_rev.s, sa_rev.len * sizeof(uint8_t));
                        memcpy(sa_rev.s, pad.s, len * sizeof(uint8_t));
                        sa_rev.len += len;
                        // free and clean work
                        free(pad.s);
                    } else {
                        memmove(sa_rev.s, &sa_rev.s[sa_rev.len - (mapBase.ri[basePos].read_head_sa + mapBase.ri[basePos].read_len - sa)], (mapBase.ri[basePos].read_head_sa + mapBase.ri[basePos].read_len - sa) * sizeof(uint8_t));
                        sa_rev.len = mapBase.ri[basePos].read_head_sa + mapBase.ri[basePos].read_len - sa;
                    }
                    int rpIndex = findReadPoolIndex(&rp, mapBase.ri[basePos].pair_head_sa, mapBase.ri[basePos].order);
                    

                }
                // free and clean work
                free(sa_rev.s);
            }
            
        }
        
    }
}

void UGInfo_destroy(UGInfo *ugInfo)
{
    if(ugInfo->upi.max > 0) free(ugInfo->upi.unitig_pairs);
    if(ugInfo->fs.max > 0) free(ugInfo->fs.elem); 
    if(ugInfo->mb.max > 0) free(ugInfo->mb.ri);
}

void extendFirstCloud_core(const OrderUnitigInfo *oui, const UGInfo *ugInfo, const bwt_t *bwt, const lib_info *libInfo)
{
    faElem_t fe; memset(&fe, 0, sizeof(faElem_t));
    fe.seq.max = 5000; fe.qual.max = 5000;
    fe.seq.s = (char*)xcalloc(fe.seq.max, sizeof(char));
    fe.qual.s = (char*)xcalloc(fe.qual.max, sizeof(char));
    int step_len = 5 ;
    AlignLimit alignLimit;
    alignLimit.max_indel_num = MAX_INDEL_NUM; alignLimit.max_indel_len = MAX_INDEL_LEN;
    alignLimit.max_mismatch = MAX_MISMATCH; alignLimit.min_flank = MIN_FLANK;
    if(oui->oe[0].len > fe.seq.max)
    {
        // make sure seq and qual alloc synchronize
        fe.seq.max = oui->oe[0].len << 1; fe.qual.max = fe.seq.max;
        fe.seq.s = (char*)realloc(fe.seq.s, fe.seq.max * sizeof(char));
        fe.qual.s = (char*)realloc(fe.qual.s, fe.qual.max * sizeof(char));
        memset(fe.seq.s + fe.seq.locPos, 0, (fe.seq.max - fe.seq.locPos) * sizeof(char));
        memset(fe.qual.s + fe.qual.locPos, 0, (fe.qual.max - fe.qual.locPos) * sizeof(char));
    }
    int fs_index = findSeqByID(ugInfo->fs, oui->oe[0].unitigID);
    memcpy(fe.seq.s, ugInfo->fs.elem[fs_index].seq.s, ugInfo->fs.elem[fs_index].seq.locPos * sizeof(char));
    fe.seq.locPos = ugInfo->fs.elem[fs_index].seq.locPos;
    fe.qual.locPos = fe.seq.locPos;
    setUnitigQuality(&fe, 0, fe.qual.locPos);
    for(int i = 1 ; i < oui->count; i++)
    {
        int j = fe.seq.locPos - KMER - step_len;
        int end = fe.seq.locPos + oui->oe[i-1].gap_len + oui->oe[i-1].SD;
        int flag = 0 ;
        ISABase isaBase ;
        isaBase.count = 0; isaBase.max = 20;
        isaBase.baseIndex = (BaseISA*)xcalloc(isaBase.max, sizeof(ISABase));
        while(j < end)
        {
            ReadPool rp = mapFlankContig(&fe, j, ugInfo, bwt, libInfo, flag, step_len, &isaBase);
            if(flag == 0) flag = 1; // just set once 

            j += step_len;
        }
    }
}

void extendCloud(const FILE *ugrm_fp, omp_lock_t *fpLock, const bwt_t *bwt, FILE *contig_fp, lib_info *libInfo)
{
    while(ugrm_fp != EOF)
    {
        // read unitig group information from ugrm_fp file
        UGInfo ugInfo;
        omp_set_lock(fpLock);
        if(readUnitigGroupInfo(&ugInfo, ugrm_fp) <= 0) break;
        omp_unset_lock(fpLock);
        
        // link unitigs group to the scaffold
        OrderUnitigInfo oui = unitigPairInfoFlatten(&ugInfo, libInfo);
        // extend first cloud
        extendFirstCloud_core(&oui, &ugInfo, bwt, libInfo);

        // free and clean work
        UGInfo_destroy(&ugInfo);
    }
    

}

int extendUnitig(Arguments *arguments)
{
    time_t timeval;
    char *prefix = arguments->prefix ;
    char kmerFreq_name[PATH_LEN];
    char unitig_name[PATH_LEN], name[PATH_LEN];
    faSeq_t fs ;
    fs.len = 0; fs.max = INITIAL_STACK_SIZE;
    fs.elem = (faElem*)xcalloc(fs.max, sizeof(faElem));
    strcpy(kmerFreq_name, prefix); strcat(kmerFreq_name, ".kmerFreq");
    strcpy(unitig_name, prefix); strcat(unitig_name, ".unitig.gz");
    curve = parseKmerFreqFile(kmerFreq_name, arguments);

    KMER = arguments->K ; NCPU = arguments->NCPU;
    int64_t unitig_len = 0;
    // read untig from *.unitig.fz
    {
        gzFile fp = xzopen(unitig_name, "r");
        kseq_t *seq;
        seq = kseq_init(fp);
        
        while(kseq_read(seq) > 0)
        {
            if(fs.len >= fs.max)
            {
                fs.max <<= 1 ;
                fs.elem = (faElem*)xrecalloc(fs.elem, fs.max * sizeof(faElem));
            }
            kseqCopy2faElem(seq, &fs.elem[fs.len]);
            unitig_len += fs.elem[fs.len].seq.locPos;
            fs.len++;
        }

        // free and clean work
        kseq_destroy(seq);
        gzclose(fp);
    }
    
    // load BWT into the RAM
    strcpy(name, prefix); strcat(name, ".bwt"); 
    bwt_t *bwt = bwt_restore_bwt_core(name);
    strcpy(name, prefix); strcat(name, ".sa"); 
    bwt_restore_sa(name, bwt);
    strcpy(name, prefix); strcat(name, ".ann");
    lib_info *libInfo = lib_info_restore(name);
    
    int total_map_number = 0;
    // mapping reads to the unitig
    {
        int64_t loc = 0 ;
        int64_t thread_len = unitig_len / NCPU;
        int64_t thread_coordinate[NCPU];
        thread_coordinate[0] = 0;
        for(int i = 1; i < NCPU; i++)
        {
            int coor = thread_coordinate[i-1];
            int64_t len = 0;
            while(len <= thread_len || coor < fs.len)
            {
                len += fs.elem[coor++].seq.locPos;
            }
            thread_coordinate[i] = coor;
        }
        
        // initial readsMap2Unitig file
        /* @file format description of *.readsMap2Unitig
         * format:  "read_position  read_pair_id Ref_name Ref_position    CIGAR   strand"
         * read_position: the mapped read(or read pair) 0-based leftmost position of reads concatenated sequence
         * read_pair_id:  if the first read close to read_position denote 1, if the another read of pair denote 2
         * Ref_name: the name of unitig that read mapping
         * Ref_position: 0-based leftmost position of reference 
         * CIGAR: same as SAM format see also <http://samtools.sourceforge.net>
         * strand: the strand(+/-) of the reference that read map
         * It's a TAB-delimited text file format
         */
        char readMapU_name[PATH_LEN];
        strcpy(readMapU_name, prefix); strcat(readMapU_name, ".readsMap2Unitig");
        FILE *readsMap_fp = fopen(readMapU_name, "w");
        omp_lock_t fpLock ;
        omp_init_lock(&fpLock);
        
        int map_number[NCPU]; 
        #pragma omp parallel for schedule(dynamic)
        for(int i = 0; i < NCPU; i++)
        {
            if( i == NCPU - 1) { map_number[i] = map2Unitig(arguments, fs.elem, thread_coordinate[i], fs.len, &fpLock, bwt, libInfo, readsMap_fp);  }
            else map_number[i] = map2Unitig(arguments, fs.elem, thread_coordinate[i], thread_coordinate[i+1], &fpLock, bwt, libInfo, readsMap_fp);
        } 
        for(int i =0; i < NCPU; i++) { total_map_number += map_number[i];}
        // free and clean work
        omp_destroy_lock(&fpLock);
        fclose(readsMap_fp);
    }

    // free memory for hash mapped reads
    bwt_destroy(bwt);
    free(fs.elem);
	// index the reads map file
	{
		FILE *index_fp, *readsMap_fp ;
		char name[PATH_LEN];
		memset(name, 0, PATH_LEN);
		strcpy(name, prefix);  strcat(name, ".unitigIndex");
		index_fp = xopen(name, "w");
		memset(name, 0, PATH_LEN);
		strcpy(name, prefix); strcat(name, ".readsMap2Unitig");
		readsMap_fp = xopen(name, "r");
		indexReadsMapFile(readsMap_fp, index_fp);
		// free and clean work
		fclose(readsMap_fp);
		fclose(index_fp);	
	}
    
	//classify unitig to the unitig group
    {
        FILE *readsMap_fp, *index_fp, *ugrm_fp;
        char name[PATH_LEN];
        memset(name, 0, PATH_LEN);
        strcpy(name, prefix); strcat(name, ".readsMap2Unitig");
        readsMap_fp = xopen(name, "r");
        // reads hash table
        readPair_num = find_next_prime_kh(total_map_number);
        ReadPos *readsHash = (ReadPos*)xcalloc(readPair_num, sizeof(ReadPos));
        UnitigPInfo upi;
        upi.count = 0; upi.max = find_next_prime_kh(fs.count);
        upi.unitig_pairs = (Unitig_pair*)xcalloc(upi.max, sizeof(Unitig_pair));
        
        constructUnitigPairInfo(readsHash, readPair_num, &upi, readsMap_fp, libInfo);
        // free memory
        free(readsHash);
        qsort(upi.unitig_pairs, upi.max, sizeof(Unitig_pair), compare_pair_num);
        // set unitigs group
        UnitigGroup unitigGroup;
        unitigGroup.count = 0; unitigGroup.max = INITIAL_STACK_SIZE;
        unitigGroup.unitigPInfo = (UnitigPInfo*)xcalloc(unitigGroup.max, sizeof(UnitigPInfo));
        groupUntigPInfo(&unitigGroup, &upi);
        UnitigLocTable ult;
        ult.count = 0; ult.max = find_next_prime_kh(fs.len * 1.4);
        ult.ul = (UnitigLoc*)xcalloc(ult.max, sizeof(UnitigLoc));
        strcpy(name, prefix); strcat(name, ".unitigIndex");
        index_fp = xopen(name, "r");
        restoreUnitigLocTable(index_fp, &ult);
        strcpy(name, prefix); strcat(name, ".unitigGroupReadsMap");
        ugrm_fp = xopen(name, "w");
        write2Ugrm(ugrm_fp, &ult, readsMap_fp, &unitigGroup);

        // free and clean work
        fclose(ugrm_fp);
        fclose(index_fp);
        free(unitigGroup.unitigPInfo);
        free(upi.unitig_pairs);
        fclose(readsMap_fp);
    }

    // restore BWT info again
    // load BWT into the RAM
    strcpy(name, prefix); strcat(name, ".bwt"); 
    bwt_t *bwt = bwt_restore_bwt_core(name);
    strcpy(name, prefix); strcat(name, ".sa"); 
    bwt_restore_sa(name, bwt);

    // extending first cloud 
    {
        FILE *ugrm_fp, contig_fp;
        strcpy(name, prefix); strcat(name, ".unitigGroupReadsMap");
        ugrm_fp = xopen(name, "r");
        strcpy(name, prefix); strcat(name, ".contig");
        contig_fp = xopen(name, "w");
        omp_lock_t fpLock;
        omp_init_lock(&fpLock);

        #pragma omp parallel for schedule(dynamic)
        for(int i =0; i < NCPU; i++)
        {
            extendCloud(ugrm_fp, &fpLock, bwt, contig_fp);
        }

        // free and clean work
        omp_destroy_lock(&fpLock);
        fclose(contig_fp);
        fclose(ugrm_fp);
    }

    // free and clean work
    bwt_destroy(bwt);
    free(curve);

    return 0 ;
}

int extend_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:     %s extend [options]\n\n", PROGRAM_NAME);
    fprintf(stderr, "Options:   -p STR          the Prefix of output files [same as PROGRAM_NAME]\n");
    fprintf(stderr, "           -c STR          the Configure file that contain library information[%s.ini]\n", PROGRAM_NAME);
    fprintf(stderr, "           -N INT          maximum Number of CPU program used( or the number of threads created) [use all the rest of CPU that system allow, depend on the system loading]\n");
    fprintf(stderr, "           -M INT(Giga)    Maximum RAM memory used(must bigger than X times of the total reads length) [no limitation]\n");
    fprintf(stderr, "           -o STR          program Output directory by running-time[current directory]\n");
    fprintf(stderr, "           -l INT          the minimum Length of contig output [4 * 'K'(Kmer length)]\n");
    fprintf(stderr, "           -h              print this  message to the stderr\n");

    return 0;
}

int hwgsa_extend(int argc, char *argv[])
{
    int c ;
    char name[PATH_LEN];
    Arguments *arguments ;

    arguments = (Arguments*)xcalloc(1, sizeof(Arguments));

    if(argc < 2) { extend_usage(); return 1 ;  }
    while((c = getopt(argc, argv, "c:p:N:M:o:l:h")) >= 0)
    {
        switch(c) {
            case 'c':   strcpy(arguments->conffile, optarg); break;
            case 'p':   strcpy(arguments->prefix, optarg); break;
            case 'N':   arguments->NCPU = atoi(optarg); break;
            case 'M':   arguments->maxMem = atoi(optarg); break;
            case 'o':   strcpy(arguments->outputDir, optarg); break;
            case 'l':   arguments->min_contig_len = atoi(optarg); break;
            case 'h':   extend_usage(); exit(1);
            case ':':   fprintf(stderr, "option %c need a value\n", optopt); extend_usage(); return 1;
            case '?':   fprintf(stderr, "[hwgsa_extend] unrecognized option '%c'\n", optopt); extend_usage(); return 1;
            default:    extend_usage(); return 1 ;
        }
    }

    // check arguments
    if(arguments->outputDir[0] == '\0')
    {
        fprintf(stderr, "outputDir have not set.....\n");
        extend_usage();
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
        checkExtendArgs(arguments, fp);
        fclose(fp);
    }
    KMER = arguments->K ;

    // extend unitig to construct contig 
    extendUnitig(arguments);


    // clean and free work
    fclose(arguments->logfp);
    free_Arguments(arguments);

    return 0 ;
}

