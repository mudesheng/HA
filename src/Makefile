CC=     gcc 
CXX=    g++
#Options for development
CFLAGS=     -g -pg -Wall -m64 -fopenmp -DDEBUG -std=gnu99 
#Options for release 
#CFLAGS=     -g -Wall -O2 -m64 -fopenmp -std=gnu99
#CFLAGS=     -Wall -O2 -m64 -fopenmp -std=gnu99
CXXFLAGS=	$(CFLAGS)
DFLAGS=    -D_FILE_OFFSET_BITS=64
#-DHAVE_PTHREAD 
OBJS=		hwgsa_all.o bwtindex.o generateGraph.o hwgsa_para.o seqIO.o bwtio.o bwt.o bwtmisc.o\
			utils.o lookup3.o bwt_merge.o bwt_gen64.o bwt_gen.o QSufSort.o	
PROG=		HWGSA0725
#LIBS=       -lm -lz -Lbwt_gen -lbwtgen -lefence
LIBS=       -lm -lz 

.SUFFIXES:.c .o .cc 

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@


all: $(PROG) 


HWGSA0725: $(OBJS) HWGSA.o
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) HWGSA.o -o $@ $(LIBS)

HWGSA.o: HWGSA.h
hwgsa_all.o: HWGSA.h 
extend.o: extend.h 
generateGraph.o: generateGraph.h bwt.h
#findOverlap.o: findOverlap.h bwt.h
lookup3.o: generateGraph.h bwt.h
bwtindex.o: bwtindex.h bwt.h
hwgsa_para.o: hwgsa_para.h utils.h bwt.h
bwt.o: utils.h bwt.h
bwtio.o: bwt.h utils.h
bwtmisc.o: bntseq.h utils.h bwt.h 
seqIO.o: seqIO.h utils.h
utils.o:  utils.h
bwt_merge.o: bwt_gen.h bwt_gen64.h
bwt_gen64.o: bwt_gen64.h bwt_gen.h QSufSort.h
bwt_gen.o: bwt_gen.h QSufSort.h
QSufSort.o: QSufSort.h bwt_gen.h

clean:
	rm -f gmon.out *.o a.out *.a *~

