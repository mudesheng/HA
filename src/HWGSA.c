#include "HWGSA.h"

#ifndef PACKAGE_VERSION 
#define PACKAGE_VERSION "1.0.0-r1"
#endif

int hwgsa_findOverlap(int argc, char *argv[])
{
	// NULL
	return 0;
} 

int hwgsa_genGraph(int argc, char *argv[])
{
	// NULL
	return 0 ;
}

int hwgsa_mapToGraph(int argc, char *argv[])
{
	// NULL
	return 0 ;
}

int hwgsa_extend(int argc, char *argv[])
{
    // NULL
    return 0;
}

int hwgsa_scaffold(int argc, char *argv[])
{
    // NULL
    return 0 ;
}

static int usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: %s(Hybrid Whole_Genome Shotgun Assembler\n", PROGRAM_NAME);
    fprintf(stderr, "         This software is served Whole_Genome Shotgun\n");
    fprintf(stderr, "         sequence assembly,fitting for mixing Next\n");
    fprintf(stderr, "         Generation Sequencing reads and other longer\n");
    fprintf(stderr, "         reads by other sequence technology)\n");
    fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
    fprintf(stderr, "Contact: Desheng Mu <mudesheng@gmail.com>\n\n");
    fprintf(stderr, "Usage:   %s  <command> [options]\n\n", PROGRAM_NAME);
    fprintf(stderr, "Command: index         index the reads\n");
    fprintf(stderr, "         findOverlap   find overlap information bewteen reads\n");
    fprintf(stderr, "         genGraph      generate contig graph\n");
    fprintf(stderr, "         extend        extend the unitig\n");
    fprintf(stderr, "         scaffold      generate scaffold use pair end information\n");
    fprintf(stderr, "         all           run all the steps of program\n");
    fprintf(stderr, "\n");
    return 1 ;
}

int main(int argc, char *argv[])
{
    if(argc < 2) return usage() ;
    if(strcmp(argv[1], "index") == 0) return hwgsa_index(argc-1 , argv+1);
    else if(strcmp(argv[1], "findOverlap") == 0) return hwgsa_findOverlap(argc-1, argv+1);
    else if(strcmp(argv[1], "genGraph") == 0) return hwgsa_genGraph(argc-1, argv+1);
	else if(strcmp(argv[1], "mapToGraph") == 0) return hwgsa_mapToGraph(argc-1, argv+1);
    else if(strcmp(argv[1], "extend") == 0) return hwgsa_extend(argc-1, argv+2);
    else if(strcmp(argv[1], "scaffold") == 0) return hwgsa_scaffold(argc-1, argv+1);
    else if(strcmp(argv[1], "all") == 0) return hwgsa_all(argc-1, argv+1);
    else {
        fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
        return 1 ;
    }
    return 0 ;
}
