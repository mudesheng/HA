#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>

#include "bwt.h"
#include "mapToGraph.h"

#define MIN_READS_LEN  0
#define MAX_READS_LEN  1

int MAX_IDENTITY_NUM = 8 ;

int GetReadsExtremumLen(const lib_info *libInfo, int FLAG)
{
	int extremum ;
	if(FLAG == MIN_REDAS_LEN) extremum = INT_LEAST32_MAX ;
	else if(FLAG == MAX_READS_LEN) extremum = 0 ;
	else { fprintf(stderr, "[GetReadsExtremumLen] the FLAG set error, exit....\n"); exit(1);}

	for(int i = 0; i < libInfo->num_lib; i++)
	{
		if(libInfo->lib[i].diverse == 0)
		{
			if(FLAG == MAX_READS_LEN &&  libInfo->lib[i].read_len > extremum) { extremum =  libInfo->lib[i].read_len; }
			else if(FLAG == MIN_READS_LEN && libInfo->lib[i].read_len < extremum) { extremum = libInfo->lib[i].read_len; }
		} else {
			for(int j = 0; j < libInfo->lib[i].number_rd; j++)
			{
				int len = libInfo->lib[i].length[j];
				if(FLAG == MAX_READS_LEN &&  len > extremum) { extremum =  len; }
				else if(FLAG == MIN_READS_LEN && len < extremum) { extremum = len; }
			}
		}
	}

	return extremum ;
}

ISARegion BwtCalWidth(const bwt_t *bwt, const char *seq)
{
	int len = strlen(seq);
	ISAregion isaRegion ;
	uint32_t c = nst_nt4_table[(int)seq[len-1]];
	isaRegion.low = bwt->L2[c] + bwt->divideNumber ;
	isaRegion.high = bwt->L2[c] + bwt->divideNumber - 1 ;
	
	for(int i = len - 2; i >= 0; i--)
	{
		c = nst_nt4_table[(int)seq[i]];
		isaRegion.low = posSA2PosBWT(bwt, isaRegion.low);
		isaRegion.high = posSA2PosBWT(bwt, isaRegion.high);
		isaRegion.low = bwt->L2[c] + bwt->occMajor[(isaRegion.low-1)/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, (isaRegion.low-1), c) + bwt->divideNumber ;
		isaRegion.high = bwt->L2[c] + bwt->occMajor[(isaRegion.high)/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, isaRegion.high, c) + bwt->divideNumber - 1 ;
	}

	return isaRegion ;
}

char *ReverseCompleString(const char *seq)
{
	int len = strlen(seq);
	char *rcSeq = xcalloc(len+ 1, sizeof(char));

	for(int i = len - 1; i >= 0; i++)
	{
		switch(seq[i]) {
			case 'A': rcSeq[len - 1 -i] = 'T' ; break ;
			case 'a': rcSeq[len - 1 -i] = 't' ; break ;
			case 'C': rcSeq[len - 1 -i] = 'G' ; break ;
			case 'c': rcSeq[len - 1 -i] = 'c' ; break ;
			case 'G': rcSeq[len - 1 -i] = 'C' ; break ;
			case 'g': rcSeq[len - 1 -i] = 'c' ; break ;
			case 'T': rcSeq[len - 1 -i] = 'A' ; break ;
			case 't': rcSeq[len - 1 -i] = 'a' ; break ;
			default : fprintf(stderr, "[ReverseCompleString] unknown character %c, exit...\n", seq[i]); exit(1);
		}
	}

	return rcSeq ;
}

int *GetRegionSeq(const bwt_t *bwt, const ISARegion isaRegion)
{
	int64_t bp_l, bp_h ;
	int *seq = xcalloc(isaRegion.high - isaRegion.low + 1, sizeof(int));
	bp_l = posSA2PosBWT(bwt, isaRegion.low); bp_h = posSA2PosBWT(bwt, isaRegion.high);
	if(bp_h - bp_l == isaRegion.high - isaRegion.low)
	{
		for(int i = 0 ; i < bp_h - bp_l + 1 ; i++)
		{
			int64_t pos = bp_l + i ;
			seq[i] = bwt_B0(bwt, pos);	
		}	
	} else {
		for(int i = 0; i < isaRegion.high - isaRegion.low + 1; i++)
		{
			bp_l = posSA2PosBWT(bwt, isaRegion.low + i);
			if(bp_l == -1) // encounter a terminal $ character
			{
				seq[i] = -1 ;
			} else {
				seq[i] = bwt_B0(bwt, bp_l);
			}
		}
	}

	return seq ;
}

/*
void SetReadsStart(const Vertex vtx, const int startPos, const int backExtendLen, const bwt_t *bwt, ISARegion isaRegBran,const lib_info *libInfo, SPArr *spArr)
{
	ISACache *isaCache = (ISACache*)xcalloc(isaRegBran.high - isaRegBran.low +1, sizeof(ISACache));
	for(int i = 0; i < isaRegBran.high - isaRegBran.low + 1, i++)
	{
		isaCache[i].isa = isaRegBran.low + i ;
		isaCache[i].loc = beginPos ;
	}
	int MAX_IDENTITY_NUM = 8 ;
	for(int i = startPos -1, i >= startPos - backExtendLen; i--)
	{
		int *seq  = GetRegionSeq(bwt, isaRegBran);
		uint32_t c = bEdge.seq[i];
		int sinkCount = 0;
		for(int j = 0; j < isaRegBran.high - isaRegBran.low + 1; j++)
		{
			if(c != seq[j] )
			{
				int64_t sa = -1 ;
				ReadLocation rl;
				if(isaCache[j].isa%SA_INTERVAL == 0) 
				{
					sa = bwt->sa[isaCache[j].isa/SA_INTERVAL] ;
					if(isaCache[j].loc > i + MAX_IDENTITY_NUM) { rl = ReadBoundaryLookup(libInfo, sa, BACKWARD); } 
					else { rl = ReadBoundaryLookup(libInfo, sa, BIDIRECTION ); }
					sa -= (isaCache[j].loc - i); 
				} else {
					sa = calculateSA(bwt, isaRegBran.low + j);
					rl = ReadBoundaryLookup(libInfo, sa, FORWARD);
				}
				ReadLocation startLoc = LocateStartSA(sa, rl);
				if(startLoc.bound != -1 )
				{
					spArr->sp[spArr->count].edgeID = edgeID ;
					spArr->sp[spArr->count].rl = startLoc;
					spArr->sp[spArr->count].sa = sa ;
					spArr->sp[spArr->count].loc = startLoc.bound - sa + i;
					spArr->count++;
				} else {
					fprintf(stderr, "[SetReadsStart] startLoc.bound set error, exit...\n", exit(1));
				}
				memmove(isaCache + j, isaCache + j+1, (isaRegBran.high - isaRegBran.low + 1 - (j+1)) * sizeof(ISACache));
				sinkCount++;
			} else {
				if((isaRegBran.low + j) % SA_INTERVAL == 0)
				{
					isaCache[j-sinkCount].isa = isaRegBran.low + j ;
					isaCache[j-sinkCount].loc = i;		
				}
			}
		}

		int64_t bp_l = posSA2PosBWT(bwt, isaRegBran.low), bp_h = posSA2PosBWT(bwt, isaRegBran.high);

		isaRegBran.low = bwt->L2[c] + bwt->occMajor[(bp_l-1)/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_l-1, c) + bwt->divideNumber ;
		isaRegBran.high = bwt->L2[c] + bwt->occMajor[bp_h/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_h, c) + bwt->divideNumber - 1 ;
		if(isaRegBran.high < isaRegBran.low) { break; }

		// final and free work
		free(seq);
	}
} */	

SPArr *BWTCalReadsStartPos(const ISARegion isaRegion, const uint32_t vID, const int beginPos, const int backExtendLen, const lib_info *libInfo, const bwt_t *bwt, const VertexArr *vertexArr)
{
	Vertex vtx = vertexArr->arr[vID];
	SPArr *spArr = (SPArr*)xcalloc(1, sizeof(SPArr));
	spArr->size = isaRegion.high - isaRegion.low + 1 ;
	spArr->sp = (SaPos*)xcalloc(spArr->size, sizeof(SaPos));	
	
	ISACache *isaCache = (ISACache*)xcalloc(isaRegion.high - isaRegion.low +1, sizeof(ISACache));
	for(int i = 0; i < isaRegion.high - isaRegion.low + 1, i++)
	{
		isaCache[i].isa = isaRegion.low + i ;
		isaCache[i].loc = beginPos ;
	}
	int MAX_IDENTITY_NUM = 8 ;
	for(int i = startPos -1; i >= startPos - backExtendLen; i--)
	{
		int *seq  = GetRegionSeq(bwt, isaRegion);
		uint32_t c = vtx.seq[i];
		int sinkCount = 0;
		for(int j = 0; j < isaRegion.high - isaRegion.low + 1; j++)
		{
			if(c != seq[j] )
			{
				int64_t sa = -1 ;
				ReadLocation rl;
				if(isaCache[j].isa%SA_INTERVAL == 0) 
				{
					sa = bwt->sa[isaCache[j].isa/SA_INTERVAL] ;
					if(isaCache[j].loc > i + MAX_IDENTITY_NUM) { rl = ReadBoundaryLookup(libInfo, sa, BACKWARD); } 
					else { rl = ReadBoundaryLookup(libInfo, sa, BIDIRECTION ); }
					sa -= (isaCache[j].loc - i); 
				} else {
					sa = calculateSA(bwt, isaRegion.low + j);
					rl = ReadBoundaryLookup(libInfo, sa, FORWARD);
				}
				ReadLocation startLoc = LocateStartSA(sa, rl);
				if(startLoc.bound != -1 )
				{
					spArr->sp[spArr->count].vID = vtx.ID ;
					spArr->sp[spArr->count].rl = startLoc ;
					spArr->sp[spArr->count].sa = startLoc.bound ;
					spArr->sp[spArr->count].loc = startLoc.bound - sa + i;
					spArr->count++;
				} else {
					fprintf(stderr, "[SetReadsStart] startLoc.bound set error, exit...\n", exit(1));
				}
				memmove(isaCache + j, isaCache + j+1, (isaRegion.high - isaRegion.low + 1 - (j+1)) * sizeof(ISACache));
				sinkCount++;
			} else {
				if((isaRegion.low + j) % SA_INTERVAL == 0)
				{
					isaCache[j-sinkCount].isa = isaRegion.low + j ;
					isaCache[j-sinkCount].loc = i;		
				}
			}
		}

		int64_t bp_l = posSA2PosBWT(bwt, isaRegion.low), bp_h = posSA2PosBWT(bwt, isaRegion.high);

		isaRegion.low = bwt->L2[c] + bwt->occMajor[(bp_l-1)/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_l-1, c) + bwt->divideNumber ;
		isaRegion.high = bwt->L2[c] + bwt->occMajor[bp_h/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_h, c) + bwt->divideNumber - 1 ;
		if(isaRegion.high < isaRegion.low) { break; }

		// final and free work
		free(seq);
	}
	
	// final and free work
	free(isaCache);

	return spArr ;
}	
		
// if have return 1 , else return 0
inline int HasReadEndSPArr(const SPArr *spArr)
{
	int retV = 0;
	for(int i = 0; i < spArr->count; i++)
	{
		if(spArr->sp[i].retainLen <= 0 && spArr->sp[i].retainLen + MAX_IDENTITY_NUM > 0 ) { retV = 1; break; }
	}

	return retV;
}

// if has unfinished element return 1 , else return 0
inline int HasUnifinishedSPArr(const SPArr *spArr)
{
	int retV = 0;
	for(int i = 0; i < spArr->count; i++)
	{
		if(spArr->sp[i].retainLen + MAX_IDENTITY_NUM > 0) { retV = 1; break; }
	}

	return retV ;
}
					

void WriteMapInfo(ReadMapInfo *rmi, SPArr *spArr, ReadLocation startLoc)
{
	// 
}

void FreeRMA(ReadMapArg rma)
{
	free(rma.spArr);
	free(rma.isaCache);
}

// return SPArr *spArr;
ReadMapInfo *BWTCalReadsForwardPos(ISARegion fIsaRegion, const SPArr *spArr, const uint32_t vID, const int beginPos, const lib_info *libInfo, const bwt_t *bwt, const VertexArr *vertexArr)
{
	Vertex vtx = vertexArr->arr[vID];
	int MAX_IDENTITY_NUM = 8 ;
	ReadMapInfo *rmi = (ReadMapInfo*)xcalloc(1, sizeof(ReadMapInfo));
	rmi->size = INITIAL_STACK_SIZE ;
	rmi->cigar = (char*)xcalloc(rmi->size, sizeof(char));
	ISACache *isaCache = (ISACache*)xcalloc(fIsaRegion.high - fIsaRegion.low + 1, sizeof(ISACache));

	for(int i = 0; i < fIsaRegion.high - fIsaRegion.low + 1, i++)
	{
		isaCache[i].isa = fIsaRegion.low + i ;
		isaCache[i].loc = beginPos ;
	}
	
	// adjust the SaPos.retainLen;
	for(int i = 0; i < spArr->count; i++)
	{
		spArr->sp[i].retainLen = spArr->sp[i].rl.size - (beginPos - spArr->sp[i].loc);
		SetSaPosCigar(spArr->sp[i]);
	}
	
	// initial ReadMapStack ;
	ReadMapStack rms ;
	rms.size = 10; rms.count = 0;
	rms.rma = (ReadMapArg*)xcalloc(rms.size, sizeof(ReadMapArg));
	rms.rma[rms.count].vID = vID ;
	rms.rma[rms.count].startLoc = beginPos ;
	rms.rma[rms.count].isaRegion = fIsaRegion ;
	rms.rma[rms.count].spArr = spArr;
	rms.rma[rms.count].isaCache = isaCache ;
	rms.count++;
	
	while(rms.count > 0)
	{
		rms.count-- ;
		ReadMapArg rma = rms.rma[rms.count] ;
		Vertex vtx = vertexArr[rma.vID];
		int i ;
		for(i = rma.startLoc; i < vtx.s_size; i++)
		{
			if(rma.isaRegion.high >= rma.isaRegion.low && HasUnfinishedReads(rma.spArr))
			{
				int *seq = GetRegionSeq(bwt, rma.isaRegion);
				uint32_t c = nst_nt4_table[(int)vtx.seq[i]];
				c = ~c & 0x3 ;
				int sinkCount = 0;
				
				for(int j = 0; j < rma.isaRegion.high - rma.isaRegion.low + 1; j++)
				{
					if(c != seq[j] && HasReadEndSPArr(rma.spArr))
					{
						ReadLocation rl ;
						int64_t sa ;
						if(isaCache[j].isa % SA_INTERVAL == 0)
						{
							sa = bwt->sa[rma.isaCache[j].isa/SA_INTERVAL];
							if(isaCache[j].loc < i - MAX_IDENTITY_NUM) { rl = ReadBoundaryLookup(libInfo, sa, BACKWARD); }
							else { rl = ReadBoundaryLookup(libInfo, sa, BIDIRECTION); }
							sa -= (rma.isaCache[j].loc - i);
						} else 	{
							sa = calculateSA(bwt, rma.isaRegion.low + j);
							rl = ReadBoundaryLookup(libInfo, sa, FORWARD);
			 			}
						ReadLocation startLoc = LocateStartSA(sa, rl);
						WriteMapInfo(rma.spArr, startLoc);
						
						memmove(rma.isaCache + j, rma.isaCache + j+1, (rma.isaRegion.high - rma.isaRegion.low + 1 - (j+1)) * sizeof(ISACache));
						sinkCount++;
					} else {
						if((rma.isaRegion.low + j) % SA_INTERVAL == 0)
						{
							rma.isaCache[j-sinkCount].isa = rma.isaRegion.low + j ;
							rma.isaCache[j-sinkCount].loc = i ;
						}
					}
				}
				
				int64_t bp_l, bp_h ;
				bp_l = posSA2PosBWT(bwt, rma.isaRegion.low); bp_h = posSA2PosBWT(bwt, rma.isaRegion.high);
				rma.isaRegion.low = bwt->L2[c] + bwt->occMajor[(bp_l-1)/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_l-1, c) + bwt->divideNumber ;
				rma.isaRegion.high = bwt->L2[c] + bwt->occMajor[bp_h/OCC_INTERVAL_MAJOR64 * 4 + c] + bwt_occ(bwt, bp_h, c) + bwt->divideNumber - 1 ;

				// reset spArr->sp[i].retainLen
				for(int i = 0; i < rma.spArr->count; i++)
				{
					rma.spArr->sp[i].retainLen-- ;
				}
				// final and free work
				free(seq);
			} else { break; } 
		}

		// add to rms stack
		if(i >= vtx.s_size && rma.isaRegion.high >= rma.isaRegion.low && HasUnfinishedReads(spArr))
		{
			for(int j = 0; j < 4; j++)
			{
				if(vtx.rvertex[j] > 0)
				{
					Vertex v = vertexArr->arr[vtx.rvertex[j]];
					uint32_t c = nst_nt4_table[(int)v.seq[KMER-1]];
					c = ~c & 0x3 ;
					int *seq = GetRegionSeq(bwt, rma.isaRegion);
					
				}
			}
		}	

		// final and free work
		FreeRMA(rma)
	}

	// final and free work
	free(rms.rma);
	free(isaCache);
	
	return rmi ;
}
		
char *GetCigarFromSPArr(SPArr *spArr, const int mapEnd)
{
	int size = INITIAL_STACK_SIZE ;
	char *cigar = (char*)xcalloc(size, sizeof(char));
	
	for(int i = 0; i < spArr->count; i++)
	{
		SaPos sp = spArr->sp[i];
		char *s = (char*)xcalloc(INITIAL_STACK_SIZE, sizeof(char));
		// format specification, "readID\tstrand_orientation\tReadLength\tMapPath(StartPosition:vertexID:mapLength\tSA\n)"
		if(sp.rl.serial == 0) 	{ sprintf(s, "%ld\t+\t%d\t%d:%u:%d\t%ld\n", sp.rl.readID, sp.rl.length[0], sp.loc, sp.vID, mapEnd - sp.loc, sp.sa);	}
		else if(sp.rl.serial == 1) { sprintf(s, "%ld\t+\t%d\t%d:%u:%d\t%ld\n", sp.rl.readID+1, sp.rl.length[1], sp.loc, sp.vID, mapEnd - sp.loc, sp.sa);  }
		else if(sp.rl.serial == 2) { sprintf(s, "%ld\t-\t%d\t%d:%u:%d\t%ld\n", sp.rl.readID+1, sp.rl.length[1], sp.loc, sp.vID, mapEnd - sp.loc, sp.sa); }
		else if(sp.rl.serial == 3) { sprintf(s, "%ld\t-\t%d\t%d:%u:%d\t%ld\n", sp.rl.readID, sp.rl.length[0], sp.loc, sp.vID, mapEnd - sp.loc, sp.sa); }
		else { fprintf(stderr, "[GetCigarFromSPArr] sp.rl.serial set error, exit...\n");  exit(1); }

		if(strlen(cigar) + strlen(s) >= size)
		{
			size <<= 1 ;
			cigar = (char*)xrecalloc(cigar, size * sizeof(char));
		}

		strcat(cigar, s);
		// final and free work
		free(s);
	}
	
	return cigar ;
}

void MapToVertex(const uint32_t vID, const lib_info *libInfo, const bwt_t *bwt, const VertexArr *vertexArr, int maxReadLen, int minReadLen, gzFile gzfp, omp_lock_t *fp_lock)
{
	Vertex v = VertexArr->arr[vID];
	if(v.s_size < KMER) {fprintf(stderr, "[MapToEdge] v.s_size < KMER, exit...\n"); exit(1); }
	else if(v.s_size < KMER  + 5) return ;
	int backExtendLen, i ;
	if(v.s_size <= minReadLen) { i = v.s_size; backExtendLen = v.s_size - KMER; }
	else { i = minReadLen ; backExtend = minReadLen - KMER; }

	for(; i < v.s_size; i += backExtendLen)
	{
		int sSize = KMER + 1 ;
		char *seqFrag = xcalloc(sSize, sizeof(char));
		ISARegion bIsaRegion, fIsaRegion ;
		
		strncpy(seqFrag, v.seq[i - KMER ], KMER);
		bIsaRegion = BWTCalWidth(bwt, seqFrag);
		char *rcstr = ReverseCompleString(seqFrag);
		fIsaRegion = BWTCalWidth(bwt, rcstr);
		free(rcstr);
					
		SPArr *spArr = BWTCalReadsStartPos(bIsaRegion, vID, i - KMER, backExtendLen, libInfo, bwt, vertexArr);
		char *cigar = GetCigarFromSPArr(spArr, i);
		//ReadMapInfo rmi = BWTCalReadsForwardPos(fIsaRegion,spArr, vID, i, libInfo, bwt, vertexArr);

		// write rmi to gzfp
		omp_set_lock(fp_lock);
		gzwrite(gzfp, cigar, strlen(cigar));
		omp_unset_lock(fp_lock);
		if(i >= v.s_size ) break ;
		else if(i + (minReadLen - KMER) > v.s_size) {
			backExtendLen = v.s_size - i ;
		} else {
			backExtendLen =  minReadLen - KMER ;
		}

		// final and free work
		free(cigar);
		free(seqFrag);
		free(spArr.sp); free(spArr);
	}	
}

void MapToGraphCore(const Arguments arguments, const lib_info *libInfo, const bwt_t *bwt, const VertexArr *verArr)
{
	int maxReadLen = arguments->maxReadLen;
	int minReadLen = arguments->minReadLen;
	char name[PATH_LEN];
	strcpy(name, arguments->prefix); strcat(name, ".mapToGraphInfo.gz");
	gzFile gzfp = xzopen(name, "w");

	omp_lock_t fp_lock ;
	omp_init_lock(&fp_lock);
	omp_set_num_threads(arguments->NCPU);

	#pragma omp parallel for schedule(dynamic)
	for(long i = 1 ; i < verArr->count; i++)
	{
		MapToVertex(i, libInfo, bwt, vertexArr, maxReadLen, minReadLen, gzfp, fp_lock);
	}
	
	
	
	// final and free work
	omp_destroy_lock(&fp_lock);
	gzclose(gzfp);
}
	
MapIndex *ReadMapInfoRestore(const Arguments *arguments, const lib_info *libInfo)
{
	char fn[PATH_LEN];
	strcpy(fn, arguments->prefix); strcat(fn, ".mapToGraphInfo.gz");
	gzFile gzfp = xzopen(fn, "r");
	int buf_size = 4 * PATH_LEN;
	char *buf = (char*)xcalloc(buf_size, sizeof(char));
	MapIndex *mi = (MapIndex*)xcalloc(1, sizeof(MapIndex));
	mi->size = libInfo->num_seqs + 1 ;
	mi->arr = (MapInfo*)xcalloc(mi->size, sizeof(MapInfo));	

	while(gzgets(gzfp, buf, buf_size) != EOF)
	{
		if(buf[buf_size -2] != '\0') { fprintf(stderr, "[ReadMapInfoRestore] buf size too smaller, exit...\n");  exit(1); }
		char *s = strtok(buf, " \t\n");
		int64_t readID = atol(s);
		s = strtok(NULL, " \t\n");
		if(*s == "+") mi->arr[readID].strand = 0 ;
		else mi->arr[readID].strand = 1 ;
		s = strtok(NULL, " \t\n");
		mi->arr[readID].len = atoi(s);
		s = strtok(NULL, " \t\n");	
		mi->arr[readID].cigar = strdup(s);
	}

	// final and free work
	free(buf);
	gzclose(gzfp);
	
	return mi ;
}	

void freeMapIndex(MapIndex *mi)
{
	for(int64_t i = 0; i < mi->size; i++)
	{
		if(mi->arr[i].cigar != NULL) { free(mi->arr[i].cigar); }
	}

	free(mi->arr) ;
	free(mi);
}

void freeMapEdgeInfo(MapEdgeInfo *meiArr, const int64_t size)
{
	for(int64_t i = 0; i < size; i++)
	{
		free(meiArr[i].readIDArr);
	}

	free(meiArr);
}
int compareLoc(const void *p, const void *q)
{
	ReadLoc *p1 = (ReadLoc*)p, *q1 = (ReadLoc*)q ;
	if(p1->loc < q1->loc) return -1 ;
	else if(p1->loc > q1->loc) return 1 ;
	else return 0 ;
}

void ConstructMapEdgeInfo(MapEdgeInfo *meiArr, const MapIndex *mi, const uint32_t mei_size)
{
	for(int64_t i = 0; i < mi->count; i++)
	{
		if(mi->arr[i].len > 0)
		{
			int startPos ;
			char *s = strtok(mi->arr[i].cigar, " \t\n:");
			startPos = atoi(s);
			s = strtok(NULL, " \t\n:");
			uint32_t edgeID = atol(s);
			if(meiArr[edgeID].size == 0) 
			{
				meiArr[edgeID].size = 10 ;
				meiArr[edgeID].readLocArr = (ReadLoc*)xcalloc(meiArr[edgeID].size, sizeof(ReadLoc));
			} else if(meiArr[edgeID].count >= meiArr[edgeID].size) {
				meiArr[edgeID].size <<= 1;
				meiArr[edgeID].readLocArr = (ReadLoc*)xrecalloc(meiArr[edgeID].readLocArr, meiArr[edgeID].size * sizeof(ReadLoc));
			}
			meiArr[edgeID].readLocArr[meiArr[edgeID].count].ID = edgeID ;
			meiArr[edgeID].readLocArr[meiArr[edgeID].count].loc = startPos ;
			meiArr[edgeID].count++;
		}
	}

	// sort the ReadID by loc postion of edge
	for(uint32_t i = 1; i < mei_size; i++)
	{
		if(meiArr[i].count > 1)
		{
			qsort(meiArr[i].readLocArr, meiArr[i].count, sizeof(ReadLoc), compareLoc);
		}
	}
}

int IsSelfCycleBubble(const int vID, const VertexArr *verArr)
{
	int l_num = 0, r_num = 0 ;
	for(int i = 0; i < 4; i++)
	{
		if(verArr->arr[vID].ledge[i] > 0) { l_num++; }
		if(verArr->arr[vID].redge[i] > 0) { r_num++; }
	}

	if(l_num == 2 && r_num == 2) return 1 ;
	else return 0 ;
}

void DeconstructSelfCycle(const uint32_t edgeID, const VertexArr *verArr, const EdgeArr *edgeArr, const MapEdgeInfo *meiArr, const MapIndex *mi const Arguments *arguments)
{
	int maxReadLen = arguments->maxReadLen ;
	int min_kmerfreq = arguments->min_kmerfreq ;
	int vID = edgeArr->arr[edgeID].rvertex ;
	int ledgeID = -1, redgeID = -1 ;
	for(int i = 0; i < 4; i++)
	{
		if(verArr->arr[vID].ledge[i] > 0 && verArr->arr[vID].ledge[i] != edgeID) { ledgeID = verArr->arr[vID].ledge[i]; }
		if(verArr->arr[vID].redge[i] > 0 && verArr->arr[vID].redge[i] != edgeID) { redgeID = verArr->arr[vID].redge[i]; }
	}
	if(ledgeID != redgeID)
	{
		int map_num = 0 ;
		int size = 10;
		int64_t *ids = (int64_t*)xcalloc(size, sizeof(int64_t));
		for(int i = meiArr[ledgeID].count - 1; i >= 0; i--)
		{
			ReadLoc rl = meiArr[ledgeID].readLocArr[i];
			if(rl.loc >= edgeArr->arr[ledgeID].s_size - maxReadLen)
			{
				if(IsMapThroughCycle(mi->arr[rl.ID], redgeID, edgeID, ledgeID))
				{
					if(map_num >= size) { size <<= 1; ids = (int64_t*)xrecalloc(ids, size * sizeof(int64_t)); }	
					ids[map_num] = rl.ID ;
					map_num++;
				}
			} else { break; }
		}
		
		if(map_num >= min_kmerfreq)
		{

		}

		// final and free work
		free(ids);
	}
}

inline SaPos GetSP(char *str)
{
	SaPos sp ; memset(&sp, 0, sizeof(SaPos));
	char *s = strtok(str, " \t\n:");		sp.rl.readID = atol(s);	
	s = strtok(NULL, " \t\n:");			sp.strand = *s;
	s = strtok(NULL, " \t\n:");			sp.rl.length[0] = atoi(s);	
	s = strtok(NULL, " \t\n:");			sp.loc = atoi(s);
	s = strtok(NULL, " \t\n:");			sp.vID = atol(s);
	s = strtok(NULL, " \t\n:");			sp.retainLen = sp.rl.length[0] - atoi(s);
	s = strtok(NULL, " \t\n:");			sp.sa = atol(s);

	return sp ;
}

char *GetRetainReadSeq(const uint8_t *pac, const int64_t beginPos, const int len)
{
	char *seq = (char*)xcalloc(len + 1, sizeof(char));
	
	for(int i = 0; i < len; i++)
	{
		uint8_t c = GET_CHAR_FROM_PAC(pac, i + beginPos);
		seq[i] = BIT_NT_CHAR[c]; 
	}

	return seq;
}
			
ReadMapInfo MapEndToVertex(const char *seq, const VertexArr *verArr, const SaPos sp)
{
	ReadMapInfo rmi ;
	rmi.size = PATH_LEN ; 
	rmi.cigar = (char*)xcalloc(rmi.size, sizeof(char));
	sprintf(rmi.cigar, "%ld\t%c\t%d\t%d:", sp.rl.readID, sp.strand, sp.rl.length[0], sp.loc);

	int k = strlen(seq);
	int j = sp.loc;
	int ID = sp.vID ;
	Vertex v = verArr->arr[ID];
	for(int i = 0; i < k; )
	{
		if(k - i <= v.s_size - j)
		{
			if(strncmp(&seq[i], &v.seq[j], k - i) == 0)
			{
				sprintf(rmi.cigar + strlen(rmi.cigar), "%d:%d\n", ID, j + (k - i));		
			}
			i = k;
		} else {
			if(strncmp(&seq[i], &v.seq[j], v.s_size - j) == 0)
			{
				sprintf(rmi.cigar + strlen(rmi.cigar), "%d:", ID);
				i += (v.s_size - j);
				for(int z = 0; z < 4; z++)
				{
					if(v.rvertex[z] > 0)
					{
						int vID = v.rvertex[z];
						if(seq[i] == verArr->arr[vID].seq[KMER -1])
						{
							ID = vID ;
							v = verArr->arr[ID];
							j = KMER - 1;
							break;
						}
					}
				}
			} else { break;}
		}

		// realloc rmi.cigar
		if(strlen(rmi.cigar) + 30 > rmi.size)
		{
			rmi.size <<= 1;
			rmi.cigar = (char*)xrecalloc(rmi.cigar, rmi.size * sizeof(char));
		}
	}


	return rmi ;
}

void MapToGraphEndCore(gzFile gzInfo, gzFile gzMap, omp_lock_t fpInfo_lock, omp_lock_t fpMap_lock, const VertexArr *verArr, const lib_info *libInfo, const uint8_t *pac)
{
	int c = 50;
	char *buf[c];
	
	for(int i = 0; i < c; i++) { buf[i] = (char*)xcalloc(PATH_LEN, sizeof(char)); }
	
	while(1)
	{
		int len = 0;
		// read from gzInfo
		omp_set_lock(fpInfo_lock);
		for(int i = 0; i < c; i++;)
		{
			if(gzgets(gzInfo, buf[i], PATH_LEN) == NULL) break;
			else { len++; }
		}
		omp_unset_lock(fpInfo_lock);
		ReadMapInfo rmi;
		rmi.size = PATH_LEN;
		rmi.cigar = (char*)xcalloc(rmi.size, sizeof(char));
		for(int i = 0; i < len; i++)
		{
			SaPos sp = GetSP(buf[i]);
			char *seq = GetRetainReadSeq(pac, sp.sa + sp.loc, sp.retainLen);
			ReadMapInfo r = MapEndToVertex(seq, verArr, sp);
			if(r.cigar[strlen(r.cigar)-1] == '\n')
			{
				if(strlen(rmi.cigar) + strlen(r.cigar) >= rmi.size)
				{
					rmi.size = (rmi.size<<1) + strlen(r.cigar);
					rmi.cigar = realloc(rmi.cigar, rmi.size * sizeof(char));
					memset(rmi.cigar + strlen(rmi.cigar), 0, (rmi.size - strlen(rmi.cigar)));
				}
				strcat(rmi.cigar, r.cigar);
			}
			// final and free work
			free(r.cigar);
			free(seq);
		}

		// write to gzMap file
		omp_set_lock(fpMap_lock);
		gzputs(gzMap, rmi.cigar);
		omp_unset_lock(fpMap_lock);

		// final and free work
		free(rmi.cigar);

		if(len < c) break;
	}

	// final and free work
	for(int i = 0; i < c; i++) { free(buf[i]); }

}

void MapToGraphEnd(const Arguments *arguments, const lib_info *libInfo, const VertexArr *verArr, const uint8_t *pac)
{
	char fn[PATH_LEN];
	strcpy(fn, arguments->prefix); strcat(fn, ".mapToGraphInfo.gz");
	gzFile gzInfo = xzopen(fn, "r");
	strcpy(fn, arguments->prefix); strcat(fn, ".mapToGraph.gz");
	gzFile gzMap = xzopen(fn, "w");

	omp_lock_t fpInfo_lock, fpMap_lock;
	omp_init_lock(&fpInfo_lock); omp_init_lock(&fpMap_lock);
	omp_set_num_threads(arguments->NCPU);

	#pragma omp parallel for schedule(dynamic)
	for(long i = 0; i < arguments->NCPU; i++)
	{
		MapToGraphEndCore(gzInfo, gzMap, fpInfo_lock, fpMap_lock, verArr, libInfo, pac);
	}
	// final and free work
	omp_destroy_lock(&fpInfo_lock); omp_destroy_lock(&fpMap_lock);
	gzclose(gzMap); gzclose(gzInfo);
}

void SimplifyGraphByMapInfo(const Arguments *arguments, const lib_info *libInfo, const bwt_t *bwt,const VertexArr *verArr, const EdgeArr *edgeArr)
{
	MapIndex *mi =	ReadMapInfoRestore(arguments, libInfo);
	MapEdgeInfo	*meiArr = (MapEdgeInfo*)xcalloc(edgeArr->count, sizeof(MapEdgeInfo));
	
	ConstructMapEdgeInfo(meiArr, mi, edgeArr->count);
	// traverse all edges
	for(uint32_t i = 1 ; i < edgeArr->count; i++)
	{
		if(edgeArr->arr[i].rvertex == edgeArr->arr[i].lvertex) // self cycle branch
		{
			int vID = edgeArr->arr[i].rvertex ;
			if(IsSelfCycleBubble(vID, verArr))
			{
				DeconstructSelfCycle(i, verArr, edgeArr, meiArr, arguments);
			}
		}
	}

	// final and free work
	freeMapEdgeInfo(meiArr, edgeArr->count);
	freeMapIndex(mi);		
}

void MapToGraph(Arguments *arguments)
{
	char name[PATH_LEN];
	strcpy(name, arguments->prefix); strcat(name, ".ann");
	lib_info *libInfo = lib_info_restore(name);
	strcpy(name, arguments->prefix); strcat(name, ".bwt.gz");
	bwt_t *bwt = bwt_restore_bwt_core(name);
	strcpy(name, arguments->prefix); strcat(name, ".sa.gz");
	bwt_restore_sa(name, bwt);
	strcpy(name, arguments->prefx); strcat(name, ".DBGvertices.gz");
	VertexArr *verArr = DBGVerticesRestore(name);
	strcpy(name, arguments->prefix); strcat(name, ".DBGEdges.gz");
	DBGEdgesRestore(verArr, name);

	MapToGraphCore(arguments, libInfo, bwt, verArr); // map head of reads to Graph
	bwt_destroy(bwt);	
	
	strcpy(name, arguments->prefix); strcat(name, ".pac.gz");
	gzFile gzpac = xzopen(name, "r");
	uint8_t *pac = (uint8_t*)xcalloc((libInfo->len_pac + CHAR_PER_BYTE-1)/CHAR_PER_BYTE, sizeof(uint8_t));
	gzread(gzpac, pac, ((libInfo->len_pac + CHAR_PER_BYTE-1)/CHAR_PER_BYTE) * sizeof(uint8_t));
	MapToGraphEnd(arguments, libInfo,verArr, pac); //map end of reads to Graph
	free(pac);

	// construct Map Info by verteices ID
	
	SimplifyGraphByMapInfo(arguments, libInfo, bwt, verArr);

	// final and free work
	gzclose(gzpac);
	VertexArrFree(verArr);

	lib_infoFree(libInfo);	
}

int mapToGraph_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:     %s mapToGraph [options]\n\n", PROGRAM_NAME);
	fprintf(stderr, "Options:   -p STR          the Prefix of output files [%s]\n", PROGRAM_NAME);
	fprintf(stderr, "           -c STR          the Configure file that contain library information[%s.ini]\n", PROGRAM_NAME);
	fprintf(stderr, "           -N INT          maximum Number of CPU program used( or the number of threads created) [use all the rest of CPU that system allowed, depend on the system load]\n");
	fprintf(stderr, "           -M INT(Giga)    Maximum RAM memory used(must bigger than X times of the total reads length) [no limitation]\n");
	fprintf(stderr, "           -o STR          program Output directory by running-time[current directory]\n");
	fprintf(stderr, "           -H INT          Heterozygosity of diploid genome, low heterozygosity set 1, middle set 2, high heterozygosity set 3 [1]\n");
	fprintf(stderr, "           -h              print this message to the stderr\n");
	return 0 ;
}

int hwgsa_mapToGraph(int argc, char *argv[])
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
			case 'h':   mapToGraph_usage() ; exit(1) ;
			case ':':   fprintf(stderr, "option %c need a value\n", optopt); mapToGraphGraph_usage(); return 1 ;
			case '?':   fprintf(stderr, "[hwgsa_genGraph] unrecognized option '%c'\n",optopt ); mapToGraph_usage() ; return 1 ;
			default: mapToGraph_usage() ; return 1 ;
		}
	}

	// check arguments
	if(arguments->outputDir[0] == '\0')
	{
		fprintf(stderr, "outputDir have not set......\n");
		mapToGraph_usage(); 
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
		CheckMapToGraphArgs(arguments, fp);
		fclose(fp);
		arguments->maxReadLen = GetReadsExtremumLen(libInfo, MAX_READS_LEN);
		arguments->minReadLen = GetReadsExtremumLen(libinfo, MIN_READS_LEN);
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
	// map sequence reads to  De bruijn Graph
	MapToGraph(arguments);

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
