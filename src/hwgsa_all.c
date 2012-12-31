#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "HWGSA.h"

int all_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:     %s all [options]\n\n", PROGRAM_NAME);
    fprintf(stderr, "Options:   -p STR          the prefix of output files [same as program name]\n");
    fprintf(stderr, "           -c STR          the configure file that contain library information[%s.ini]\n", PROGRAM_NAME);
    fprintf(stderr, "           -N INT          Maximum number of CPU that program used( or the number of threads created) [use all the rest of CPU that system allowed, depend on the system load]\n");
    fprintf(stderr, "           -M INT(Giga)    Maximum RAM memory used(must bigger than X times of the total reads length) [no limitation]\n");
    fprintf(stderr, "           -G INT(Mega)    approximate genome size by prior knonw, setting genome size if assemble a high heterozygosity genome\n");
    fprintf(stderr, "           -K INT          The kmer size used for generating graph, must small than minimum read length that used construct contig, program will filter out the reads than small than K+5, the value must be odd and between %d~%d[%d]\n", MIN_KMER, MAX_KMER, MAX_KMER);
    fprintf(stderr, "           -o STR          program output directory by running-time[current directory]\n");
    fprintf(stderr, "           -q INT          the quality of base calling by cutting off, 1 is low quality, 2 is middle value, 3 is high [3]\n");
    fprintf(stderr, "           -H INT          heterozygosity of diploid genome, low heterozygosity set 1, middle set 2, high heterozygosity set 3 [1]\n");
    return 0 ;
}

int hwgsa_all(int argc, char *argv[])
{
    int c ;
    char name[PATH_LEN];
    Arguments *arguments ;
    PConfReturn *retV = NULL ;
    lib_info *libIndex_info = NULL ;
    BntWriteArgs *bntWriteArgs = NULL ;

    arguments = (Arguments*)xcalloc(1 , sizeof(Arguments));

    if(argc < 2) { all_usage(); return 1 ; }
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
            case 'h':   all_usage() ; exit(1) ;
            case ':':   fprintf(stderr, "option %c need a value\n", optopt); all_usage(); return 1 ;
            case '?':    fprintf(stderr, "[hwgsa_all] unrecognized option '%c'\n",optopt ); all_usage() ; return 1 ;
            default: all_usage() ; return 1 ;
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
    
    // generate De Bruijn Graph
    //generateKmerGraph(arguments);

    // clean and free work 
    fclose(arguments->logfp);
    free_Arguments(arguments);

    return 0 ;
}
