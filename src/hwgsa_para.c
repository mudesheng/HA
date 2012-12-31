#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <unistd.h>
#include <sys/sysinfo.h>
#include <time.h>
#include "hwgsa_para.h"
#include "utils.h"
 
// check arguments
void preProcessArgs(Arguments *arguments)
{
    
    if(arguments->prefix[0] == '\0')
    {
        //arguments->prefix = (char*)xcalloc(PATH_LEN, sizeof(char));
        strcpy(arguments->prefix, PROGRAM_NAME);
    }

    if(arguments->conffile[0] == '\0')
    {
        //arguments->conffile = (char*)xcalloc(PATH_LEN, sizeof(char));
        strcpy(arguments->conffile, PROGRAM_NAME); strcat(arguments->conffile, ".ini");
    }

    if(arguments->outputDir[0] == '\0')
    {
        //arguments->outputDir = (char*)xcalloc(PATH_LEN, sizeof(char));
        strcpy(arguments->outputDir, ".");
    } else if(opendir(arguments->outputDir) != NULL)
    {
        fprintf(stderr, "Directory: '%s' has been existed, program exit....\n", arguments->outputDir);
        exit(1) ;
    } else {
        int f ;
        if((f = mkdir(arguments->outputDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) !=0)
        {
            fprintf(stderr, "Directory: '%s' create error,program exit....\n", arguments->outputDir);
            exit(1);
        }
    }
    
    // ckeck 'K'
    if (arguments->K == 0)
    {
        arguments->K = DEF_KMER ;
    } else if((arguments->K < MIN_KMER) || (arguments->K > MAX_KMER)) {
        fprintf(stderr, "argument 'K' must between %d~%d , program exit....\n", MIN_KMER, MAX_KMER);
        exit(1);
    } else if(arguments->K % 2 == 0) {
        (arguments->K)++;
        fprintf(stderr,"argument 'K' must be odd, program automatically add 1....\n"); 
    }
    
    // set system load-balancing
    {
        int totalNCPU = 0 , p = 0 , load = 0 ;
        struct sysinfo info ;
        char tmpname[PATH_LEN], cmd[PATH_LEN];
        FILE *tmpfp ;
        strcpy(tmpname, "hwgsaloadavgtmpXXXXXX");
        strcpy(cmd, "cat /proc/loadavg >"); strcat(cmd, tmpname);
        if(system(cmd) != 0)
        {
            fprintf(stderr, "[preProcessArgs] system call is error, program exit...\n");
            exit(1);
        }
        tmpfp = fopen(tmpname, "r");
        if( tmpfp == NULL || feof(tmpfp))
        {
            fprintf(stderr, "[preProcessArgs] tmp file is error or NULL, program exit...\n");
            exit(1);
        }
        fgets(cmd, PATH_LEN, tmpfp);
        for(int i = 1 ; i < strlen(cmd); i++)
        {
            if((cmd[i-1] != ' ' && cmd[i-1] != '\t') && (cmd[i] == ' ' || cmd[i] == '\t') && p < 2)
            {
                p++ ;
            } else if((p == 2) && (cmd[i] != ' ' && cmd[i] != '\t')) {
                p = i ;
            } else if((p > 2) && (cmd[i] == ' ' || cmd[i] == '\t')) {
                cmd[i] = '\0' ;
            }
        }
        load = (int)atof(cmd+p);
        totalNCPU = sysconf(_SC_NPROCESSORS_ONLN);
        if(sysinfo(&info) != 0) perror("system call 'sysinfo' function error");
        if(arguments->NCPU == 0) { arguments->NCPU = (totalNCPU - load); }
        else if(arguments->NCPU > totalNCPU - load) {
            fprintf(stderr, "[preProcessArgs] system have not enough CPU resource amount to your requirement, please try again a moment later, exit...\n");
            exit(1);
        }
        if(arguments->NCPU < 2) arguments->NCPU = 2 ;

        if(arguments->maxMem == 0)
        {
            arguments->maxMem = (int)(info.totalram/(1024 * 1024 * 1024));
        } else {
            if(arguments->maxMem > (int)(info.totalram/(1024 * 1024 * 1024)))
            {
                fprintf(stderr, "no enough RAM %d (G) be content to program requirement, the system just total have RAM %ld (G) only, program exit....\n", arguments->maxMem, (int64_t)(info.totalram/(1024 * 1024 * 1024)));
                exit(1) ;
            }
        }
        fclose(tmpfp);
        if(remove(tmpname) != 0) 
        {
            fprintf(stderr, "Delete file : %s is error!\n", tmpname);
        }
    }
    
    if(arguments->qual != 0)
    {
        if((arguments->qual < 1) || (arguments->qual > 4))
        {
            fprintf(stderr, "argument 'q' must between 1~4, program exit....\n");
            exit(1);
        }
    } else {
        arguments->qual = 3 ;
    }

    if(arguments->H != 0)
    {
        if((arguments->H < 1) || (arguments->H > 3))
        {
            fprintf(stderr, "argument 'H' must between 1~3, program exit....\n");
            exit(1) ;
        }
    } else {
        arguments->H = 1 ;
    }
    if(arguments->min_kmerfreq == 0)
    {
        arguments->min_kmerfreq = 2 ;
    }


    // check prefix set by command line 
    {
        char *prefix = (char*)xcalloc(PATH_LEN, sizeof(char));
        strcpy(prefix, arguments->outputDir);
        if(prefix[strlen(prefix)-1] != '/') strcat(prefix, "/");
        strcat(prefix, arguments->prefix);
        strcpy(arguments->prefix, prefix);
        free(prefix);
    }
}

PConfReturn *parseConffile(char *conffile)
{
    lib_info *libIndex_info = (lib_info*)xcalloc(1 , sizeof(lib_info));
    FILE *fp = xopen(conffile, "r");
    int count = -1 ;
    PConfReturn *retV = (PConfReturn*)xcalloc(1, sizeof(PConfReturn));
    BntWriteArgs *bntWriteArgs = NULL ;
    char s[PATH_LEN];
    int file_num = 0 ;

    // check if configure file is NULL 
    if(feof(fp) != 0) 
    {
        fprintf(stderr, "[parseConffile]The configure file : '%s' file is NULL, please check, program exit....\n", conffile);
        exit(1);
    }
    while(fgets(s , PATH_LEN, fp))
    {
        // set the begin of line
		char *strip_s = KStrip(s) ;
		char *key = strtok(strip_s, " \t\n=");
	    	
		// check the length of input line
        if(key == NULL) continue; // blank line
        if((key[0] == ';') || (key[0] == '#')) { continue; } // comment line 
        else if(strcmp(key, "[global_setting]") == 0) {
            count = 0 ;
        } else if(strcmp(key , "[LIB]") == 0) {
            count++ ;
            file_num = 0 ;
            if(count >= libIndex_info->num_lib)
            {
                libIndex_info->lib = (LIB*)realloc(libIndex_info->lib, (count<<1) * sizeof(LIB));
                bntWriteArgs = (BntWriteArgs*)realloc(bntWriteArgs, (count<<1)* sizeof(BntWriteArgs));
                memset(libIndex_info->lib + libIndex_info->num_lib, 0, ((count<<1) - libIndex_info->num_lib) * sizeof(LIB));
                memset(bntWriteArgs + libIndex_info->num_lib, 0, ((count<<1) - libIndex_info->num_lib) * sizeof(BntWriteArgs));
                libIndex_info->num_lib = (count<<1);
            }
        } else if(count >= 0){
            char *value = strtok(NULL, " \t\n=");  
            if(value == NULL)
            {
                fprintf(stderr, "[parseConffile]the line '%s' in '%s' file set wrong,program exit....\n", s, conffile);
                exit(1) ;
            }
            if(count == 0) // global setting
            {
                if(strcmp(key, "max_rd_len") == 0)
                {
                    int temp = atoi(value);
                    if(temp >= pow(2,16) || temp <= MIN_KMER)
                    {
                        fprintf(stderr, "[parseConffile] the 'max_rd_len' setting in '%s' must between (%d ~ %d), please check, program exit...\n", s,  MIN_KMER, (int)pow(2, 16));
                        exit(1);
                    } else { libIndex_info->max_rd_len = temp; }
                } else {
                    fprintf(stderr, "[parseConffile]the line: '%s' is unrecognized in '%s' file, please check, program exit...\n", s, conffile);
                    exit(1);
                }
            } else { // sections(library) setting
				if(strcmp(key, "name") == 0) // match name partition
				{
					libIndex_info->lib[count-1].name = strdup(value);
					// check if the name has been existed
					for(int i = 0; i < count - 1; i++)
					{
						if(strcmp(libIndex_info->lib[i].name, libIndex_info->lib[count-1].name) == 0)
						{
							fprintf(stderr, "[parseConffile] the lib name '%s' have been appeared in '%s' file, the lib name must unique, please check, program exit...\n", libIndex_info->lib[count-1].name, conffile);
							exit(1);
						}
					}
				}else if(strcmp(key, "avg_insert_len") == 0) {
                    libIndex_info->lib[count-1].paired = 1 ;
                    int temp = atoi(value);
                    if(temp >= pow(2,31) || temp <= 2 * MIN_KMER)
                    {
                        fprintf(stderr, "[parseConffile] the 'avg_insert_len' setting in '%s' must between (%d ~ %d), please check, program exit...\n", s, 2 * MIN_KMER, (int)pow(2, 31));
                        exit(1);
                    } else { libIndex_info->lib[count-1].pe.insert_size = temp; }
                } else if((strcmp(key, "insert_SD")) == 0) {
                    int temp = atoi(value);
                    if(temp >= pow(2, 31) || temp < 0)
                    {
                        fprintf(stderr, "[parseConffile] the 'insert_SD' setting in '%s' must between [0 ~ %d), please check, program exit...\n", s, (int)pow(2, 31));
                        exit(1);
                    } else { libIndex_info->lib[count-1].pe.insert_SD = temp; }
                } else if((strcmp(key, "diverse_rd_len")) == 0) {
                    int temp = atoi(value);
                    if(temp != 0 && temp != 1)
                    {
                        fprintf(stderr, "[parseConffile] the 'diverse_rd_len' setting in '%s' must be 0 or 1 either, please check, program exit....\n", s);
                        exit(1);
                    } else { libIndex_info->lib[count-1].diverse = temp; } 
                } else if((strcmp(key, "asm_flag")) == 0) {
                    int temp = atoi(value);
                    if(temp != 1 && temp != 2 && temp != 3)
                    {
                        fprintf(stderr, "[parseConffile] the 'asm_flag' setting in '%s' must be [1~3], unrecognize : %d, pleaser check, program exit....\n", s, temp);
                        exit(1);
                    } else { libIndex_info->lib[count-1].asm_flag = temp; }
                } else if((strcmp(key, "seq_profile")) == 0) {
                    int temp = atoi(value);
                    if(temp != 1 && temp != 2)
                    {
                        fprintf(stderr, "[parseConffile] the 'seq_profile' setting in '%s' must be [1~2], the unrecognize symbols : %d, please check, program exit...\n", s, temp);
                        exit(1);
                    } else { libIndex_info->lib[count-1].seq_profile = temp; }
                } else if((strcmp(key, "qual_benchmark")) == 0) {
                    int temp = atoi(value);
                    if(temp != 64 && temp != 33)
                    {
                        fprintf(stderr, "[parseConffile] the 'qual_benchmark' setting in '%s' unrecognized: %d, program just identify phred64/phred33, please check, program exit....\n", s, temp);
                        exit(1);
                    } else { libIndex_info->lib[count-1].qual_benchmark = temp; }
                } else if(strncmp(key, "f", 1) == 0) {
                    file_num++ ;
                    if(bntWriteArgs[count-1].max <= (file_num + 1)/2)
                    {
                        bntWriteArgs[count-1].readName = (ReadName*)realloc(bntWriteArgs[count-1].readName, file_num * sizeof(ReadName));
                        memset(bntWriteArgs[count-1].readName + bntWriteArgs[count-1].max , 0, (file_num - bntWriteArgs[count-1].max) * sizeof(ReadName));
                        bntWriteArgs[count-1].max = file_num;
                    } 
                    if(libIndex_info->lib[count-1].paired == 1)
                    {
                        bntWriteArgs[count-1].readName[bntWriteArgs[count-1].count].paired = 1 ;
                        if(file_num % 2 != 0)
                        {
                            strcpy(bntWriteArgs[count-1].readName[bntWriteArgs[count-1].count].f1name, value);
                        } else  {
                            strcpy(bntWriteArgs[count-1].readName[bntWriteArgs[count-1].count].f2name, value);
                            bntWriteArgs[count-1].count++;
                        }
                    } else {
                        fprintf(stderr, "[parseConffile]the line : '%s' is unrecognized in '%s' file, maybe this is a single end library, please check, program exit....\n", s, conffile);
                        exit(1);
                    }
                } else if((strncmp(key, "sf", 1)) == 0) {
                    file_num++ ;
                    if(libIndex_info->lib[count-1].paired == 0)
                    {
                        if(bntWriteArgs[count-1].max <= file_num)
                        {
                            bntWriteArgs[count-1].readName = (ReadName*)realloc(bntWriteArgs[count-1].readName, (file_num<<1) * sizeof(ReadName));
                            memset(bntWriteArgs[count-1].readName + bntWriteArgs[count-1].max , 0 , ((file_num<<1) - bntWriteArgs[count-1].max) * sizeof(ReadName));
                            bntWriteArgs[count-1].max = (file_num<<1);
                        }
                        strcpy(bntWriteArgs[count-1].readName[bntWriteArgs[count-1].count].f1name, value);
                        bntWriteArgs[count-1].count++ ;
                    } else {
                        fprintf(stderr, "[parseConffile] the line : '%s' is unrecognized in '%s' file, maybe this is a paired end library, please check, program exit...\n", s, conffile);
                        exit(1);
                    }
                }
            }
        } else {
            fprintf(stderr, "[parseConffile]The line: '%s' in the configure file is set wrong, please note that program just recognize the comment line begin with '#' or ';', please check, program exit....\n", s);
            exit(1);
        }
    }
    libIndex_info->num_lib = count;
    if(libIndex_info->max_rd_len == 0) 
    {
        fprintf(stderr, "[parseConffile] the 'max_rd_len' of global_setting section must be set in '%s' file, please check, program exit...\n", conffile);
        exit(1);
    }
    // check LIB and BntWriteArgs information that must be set
    for(int i = 0; i < count; i++)
    {
        if(libIndex_info->lib[i].name == NULL)
        {
            fprintf(stderr, "[parseConffile] the section name set NULL in '%s' file, please check, program exit...\n", conffile);
            exit(1);
        }
        if(libIndex_info->lib[i].asm_flag == 0)
        {
            fprintf(stderr, "[parseConffile] the 'asm_flag' is not set of section '%s' in '%s' file, please check, program exit...\n", libIndex_info->lib[i].name, conffile);
            exit(1);
        }
        if(libIndex_info->lib[i].seq_profile == 0)
        {
            fprintf(stderr, "[parseConffile] the 'seq_profile' is not set of section '%s' in '%s' file, please check, program exit...\n", libIndex_info->lib[i].name, conffile);
            exit(1);
        } 
        if(libIndex_info->lib[i].qual_benchmark == 0)
        {
            fprintf(stderr, "[parseConffile] the 'qual_benchmark' is not set of section '%s' in '%s' file, please check, program exit...\n", libIndex_info->lib[i].name, conffile);
            exit(1);
        }
        if(bntWriteArgs[i].count == 0)
        {
            fprintf(stderr, "[parseConffile] the section '%s' not contain any sequence file in '%s' configure file, please check, program exit...\n", libIndex_info->lib[i].name, conffile);
            exit(1);
        } else {
            if(libIndex_info->lib[i].paired == 1)
            {
                for(int j = 0; j < bntWriteArgs[i].count; j++)
                {
                    if(bntWriteArgs[i].readName[j].f1name[0] == '\0' || bntWriteArgs[i].readName[j].f2name[0] == '\0')
                    {
                        fprintf(stderr, "[parseConffile] the sequence files in section '%s' must be paired in the '%s' configure file, please check, program exit...\n", libIndex_info->lib[i].name, conffile);
                        exit(1);
                    }
                }
            } else {
                for(int j = 0; j < bntWriteArgs[i].count; j++)
                {
                    if(bntWriteArgs[i].readName[j].f1name == '\0')
                    {
                        fprintf(stderr, "[parseConffile] the sequence file name is NULL  in section '%s' in the '%s' configure file, please check, program exit...\n", libIndex_info->lib[i].name, conffile);
                        exit(1);
                    }
                }
            }
        }
    }
    /*
    // check 'asm_flag' and delete that not construct contig library
    {
        int c = 0 ;
        for(int i = 0 ; i < count ; i++)
        {
            if(libIndex_info->lib[i].offset == 2)
            {
                free(libIndex_info->lib[i].name);
                free(bntWriteArgs[i].readName);
            } else if(c < i) {
                memcpy(&(libIndex_info->lib[c]), &(libIndex_info->lib[i]), sizeof(LIB));
                memcpy(&(bntWriteArgs[c]), &(bntWriteArgs[i]), sizeof(BntWriteArgs)); 
                c++ ;
            } else {
                c++ ;
            }
        }
        if(libIndex_info->num_lib > c)
        {
            libIndex_info->lib = (LIB*)realloc(libIndex_info->lib, c * sizeof(LIB));
            bntWriteArgs = (BntWriteArgs*)realloc(bntWriteArgs, c * sizeof(LIB));
        }
        libIndex_info->num_lib = c ;

        //clean 'asm_flag'
        for(int i = 0 ; i < libIndex_info->num_lib ; i++)
        {
            libIndex_info->lib[i].offset = 0 ;
        }
    } */
    retV->libIndex_info = libIndex_info ;
    retV->bntWriteArgs = bntWriteArgs ;
    return retV ;
}

void writeArgsToFile(const Arguments *arguments, FILE *fp)
{
    fprintf(fp, "%s\n", arguments->conffile); 
    fprintf(fp, "%s\n", arguments->prefix);
    fprintf(fp, "%d\n", arguments->NCPU);
    fprintf(fp, "%d\n", arguments->maxMem);
    fprintf(fp, "%u\n", arguments->K);
    fprintf(fp, "%s\n", arguments->outputDir);
    fprintf(fp, "%d\n", arguments->qual);
    fprintf(fp, "%d\n", arguments->G);
    fprintf(fp, "%u\n", arguments->H);
    fprintf(fp, "%lu\n", arguments->seq_len);
    fprintf(fp, "%d\n", arguments->min_kmerfreq);
    fprintf(fp, "%d\n", arguments->depth);
    fprintf(fp, "%d\n", arguments->min_unitig_len);
}

int checkGenGraphArgs(Arguments *arguments, FILE *fp)
{
    char s[PATH_LEN];
    int t ;
    //uint64_t uv ;
    if(arguments->conffile[0] == '\0') {  fscanf(fp, "%s\n", arguments->conffile);  }
    else fscanf(fp, "%s\n", s);
    fscanf(fp, "%s\n", s);
    if(strcmp(s, arguments->prefix) != 0) 
    {
        fprintf(stderr, "[checkGenGraphArgs] arguments->prefix not equal to index step...\n");
        exit(1);
    }
    if(arguments->NCPU == 0) fscanf(fp, "%d\n", &arguments->NCPU);
    else fscanf(fp, "%d\n", &t);
    if(arguments->maxMem == 0) fscanf(fp, "%d\n", &arguments->maxMem);
    else fscanf(fp, "%d\n", &t);
    if(arguments->K == 0) fscanf(fp, "%u\n", &arguments->K);
    else {fprintf(stderr, "[checkGenGraphArgs]the value of argument 'K' have been set first step, exit...\n"); exit(1); }
    fscanf(fp, "%s\n", s);
    if(arguments->qual == 0) fscanf(fp, "%d\n", &arguments->qual);
    else {fprintf(stderr, "[checkGenGraphArgs]the value of argument 'q' have been set first step, exit...\n"); exit(1); };
    if(arguments->G == 0) fscanf(fp, "%d\n", &arguments->G);
    else fscanf(fp, "%d\n", &t);
    fscanf(fp, "%d\n", &t);
    if(arguments->H == 0) arguments->H = (uint8_t)t;
    fscanf(fp, "%lu\n", &arguments->seq_len);
    if(arguments->min_kmerfreq == 0) fscanf(fp, "%d\n", &arguments->min_kmerfreq);
    else fscanf(fp, "%d\n", &t);
    if(arguments->depth == 0) fscanf(fp, "%d\n", &arguments->depth);
    else fscanf(fp, "%d\n", &t);
    if(arguments->min_unitig_len == 0) fscanf(fp, "%d\n", &arguments->min_unitig_len);
    else fscanf(fp, "%d\n", &t);
    
    // check arguments
    xassert(arguments->K >= 33 && arguments->K <= MAX_KMER , "[checkGenGraphArgs]the argument 'K' is beyond the range 33 to 201...\n");
    xassert(arguments->NCPU >= 2 && arguments->NCPU <= 64, "[checkGenGraphArgs]the argument 'N' set wrong^_^\n");
    if(arguments->maxMem == 0) 
    {
        fprintf(stderr, "[checkGenGraphArgs]don't set the argument 'M', exit....\n");
        exit(1);
    }
    if(arguments->depth == 0) { arguments->depth = 4 * arguments->K ; } 
    else {
        if(arguments->depth < 3 * arguments->K || arguments->depth > 10 * arguments->K)
        {
            fprintf(stderr, "[checkGenGraphArgs]Please set 'd' between %d ~ %d, exit....\n ", 3 * arguments->K, 10 * arguments->K);
            exit(1);
        }
    }
    if(arguments->min_unitig_len == 0) { arguments->min_unitig_len = 3 * arguments->K; }
    else {
        if(arguments->min_unitig_len < 2 * arguments->K )
        {
            fprintf(stderr, "[checkGenGraphArgs]Please set 'm' make sure bigger than 2 * 'K'(kmer length),exit...\n");
            exit(1);
        }
    }
    return 0 ;
}

int checkExtendArgs(Arguments *arguments, FILE *fp)
{
    char s[PATH_LEN];
    int t;
    if(arguments->conffile[0] == '\0') { fscanf(fp, "%s\n", arguments->conffile); }
    else fscanf(fp, "%s\n", s);
    fscanf(fp, "%s\n", s);
    if(strcmp(s, arguments->prefix) != 0)
    {
        fprintf(stderr, "[checkExtendArgs] arguments->prefix not equal to index step...\n");
        exit(1);
    }
    if(arguments->NCPU == 0) fscanf(fp, "%d\n", &arguments->NCPU);
    else fscanf(fp, "%d\n", &t);
    if(arguments->maxMem == 0) fscanf(fp, "%d\n", &arguments->maxMem);
    else fscanf(fp, "%d\n", &t);
    if(arguments->K == 0) fscanf(fp, "%u\n", &arguments->K);
    else { fprintf(stderr, "[checkExtendArgs] the value of argument 'K' has been set first step, exit...\n"); exit(1);  }
    if(arguments->qual == 0) fscanf(fp, "%d\n", &arguments->qual);
    else { fprintf(stderr, "[checkExtendArgs] the value of arguments 'q' has been set first step, exit....\n"); exit(1);}
    if(arguments->G == 0) fscanf(fp, "%d\n", &arguments->G);
    else { fprintf(stderr, "[checkExtendArgs] the value of arguments 'G' has been set first step, exit...\n"); exit(1); }
    if(arguments->H == 0) { fscanf(fp, "%u\n", &t) ; arguments->H = (uint8_t)t ; }
    else { fprintf(stderr, "[checkExtendArgs] the value of arguments 'H' has been set first step, exit...\n"); exit(1); }
    
    return 0 ;
}


void free_Arguments(Arguments * arguments)
{
    //free(arguments->conffile);
    //free(arguments->prefix);
    //free(arguments->outputDir);
    free(arguments);
}
// flag note FORWARD , BACKWARD and BIDIRECTION , find the read start position
static inline ReadLocation lookUpBound(const LIB *lib, int64_t position, const int flag)
{
    ReadLocation rl;
    position -= lib->offset ;
	position /= 2; // the offset set by all forward and reverse strand 
    if(lib->diverse == 1)
    {
        int64_t len,lastLen, i;
        int64_t lenMajorSize = (lib->number_rd + LEN_INTERVAL_MAJOR - 1) / LEN_INTERVAL_MAJOR  + 1; 
        int64_t low = 0 ,mid,high  = lenMajorSize - 1;
        while(low <= high)
        {
            mid = (low + high) / 2 ;
            if(position == lib->lenMajor[mid]) break ;
            if(position < lib->lenMajor[mid]) high = mid - 1 ;
            else low = mid + 1 ;
        }
        if(lib->lenMajor[mid] < position) {low = mid; high = mid+1; }
        else {low = mid - 1; high = mid; }
        if((position - lib->lenMajor[low]) <= ( lib->lenMajor[high] - position)) // bi-dirtection search
        {
            len = (lastLen = lib->lenMajor[low]) ;
            i = low * LEN_INTERVAL_MAJOR ;
            while(len <= position)
            {
                lastLen = len ;
                len += lib->length[i];
                i++ ;
            }
			if(flag == FORWARD || (flag == BIDIRECTION && len - position < position - lastLen))
			{ lastLen = len; len += lib->length[i]; i++;  }
			if(lib->paired == 1) 
			{
				rl.paired = 1;
				if((i - 1) % 2 == 0) { rl.bound = lib->offset/2 + lastLen;  rl.length[0] = lib->length[i-1]; rl.length[1] = lib->length[i]; rl.readID = i - 1; }
				else { rl.bound = lib->offset/2 + lastLen - lib->length[i-2]; rl.length[0] = lib->length[i-2]; rl.length[1] = lib->length[i-1]; rl.readID = i - 2;  }
				rl.pe = lib->pe;
			} else {
				rl.bound = lib->offset/2 + lastLen ; 
				rl.length[0] = lib->length[i-1];
				rl.readID = i - 1;
			}
        } else {
            len = lib->lenMajor[high];
            if(high >= lenMajorSize - 1) i = lib->number_rd - 1;
            else  i = high * LEN_INTERVAL_MAJOR - 1;
            while(position < len)
            {
                len -= lib->length[i];
                i-- ;
            }
			if(flag == FORWARD || (flag == BIDIRECTION && len + lib->length[i+1] - position < position - len))
			{ i++; len += lib->length[i];  }
			if(lib->paired == 1) 
			{
				rl.paired = 1;
				if((i + 1) %2 == 0) { rl.bound = lib->offset/2 + len;  rl.length[0] = lib->length[i+1]; rl.length[1] = lib->length[i+2]; rl.readID = i + 1; }
				else { rl.bound = lib->offset/2 + len - lib->length[i]; rl.length[0] = lib->length[i]; rl.length[1] = lib->length[i+1];  rl.readID = i; }
				rl.pe = lib->pe;
			} else {
				rl.bound = lib->offset/2 + len; 
				rl.length[0] = lib->length[i+1];
				rl.readID = i + 1;
			}
        }
    } else {
        int64_t num_offset = position / lib->read_len ;
		if(flag == FORWARD || (flag == BIDIRECTION && (num_offset + 1) * lib->read_len - position < position - num_offset * lib->read_len))
		{ num_offset++ ; }
        rl.bound = lib->offset/2 + num_offset * lib->read_len; rl.length[0] = lib->read_len;
		if(lib->paired == 1) 
		{
			rl.paired = 1;
			if(num_offset % 2 == 0) { rl.bound = lib->offset/2 + num_offset * lib->read_len; rl.readID = num_offset; }
			else { rl.bound = lib->offset/2 + (num_offset - 1) * lib->read_len ; rl.readID = num_offset - 1; }
		   	rl.length[0] = lib->read_len; rl.length[1] = lib->read_len; 
			rl.pe = lib->pe;
		} else {
			rl.bound = lib->offset/2 + num_offset * lib->read_len; 
			rl.length[0] = lib->read_len ;
			rl.readID = num_offset ;
		}
    }
	// reset bound of readLocation
	rl.bound *= 2 ;
    return rl;
}



// flag note FOWARD , BACKWARD and BIDIRECTION , return readID of all sorted library 
ReadLocation ReadBoundaryLookup(const lib_info *libIndex_info, const int64_t position, const int flag)
{
    ReadLocation rl;
    int i ;
   	int64_t num_reads = 0; 
    for(i = 0; i < libIndex_info->num_lib ; i++)
    {
        if(position < libIndex_info->lib[i].offset) break ;
		num_reads += libIndex_info->lib[i].number_rd ;
    }

    rl = lookUpBound(&(libIndex_info->lib[i-1]), position, flag);
	rl.readID += num_reads ;
    
    return rl;
}
 

// if sa < rl region, return rl.bound == -1 , else return read start postion sa, read length set in rl.length[0]
// flag note FORWARD , BACKWARD and BIDIRECTION
ReadLocation LocateStartSA(const int64_t sa, const ReadLocation rl, const int flag)
{
	ReadLocation retRl ; memset(&retRl, 0, sizeof(ReadLocation));
	retRl.readID = rl.readID ;
	if(flag == FORWARD)
	{
		if(rl.paired == 1)
		{
			if(sa >= rl.bound + (rl.length[0] + rl.length[1]) * 2) { retRl.bound = -1; }
			else {
				retRl.paired = 1 ;
				if(sa <= rl.bound) { retRl.bound = rl.bound; retRl.size = rl.length[0]; retRl.serial = 0; }
				else if(sa <= rl.bound + rl.length[0]) { retRl.bound = rl.bound + rl.length[0] ; retRl.size = rl.length[1]; retRl.serial = 1; } 
				else if(sa <= rl.bound + rl.length[0] + rl.length[1]) {
					retRl.bound = rl.bound + rl.length[0] + rl.length[1]; retRl.size = rl.length[1]; retRl.serial = 2 ; 
				} else if(sa <= rl.bound + rl.length[0] + rl.length[1] * 2) {
					retRl.bound = rl.bound + rl.length[0] + rl.length[1] * 2 ; retRl.size = rl.length[0]; retRl.serial = 3;
				}
			}
		} else {
			if(sa >= rl.bound + rl.length[0] * 2) { retRl.bound = -1; }
			else {
				retRl.paired = 0;
				if(sa <= rl.bound) { retRl.bound = rl.bound; retRl.size = rl.length[0]; retRl.serial = 0; }
				else if(sa <= rl.bound + rl.length[0]) { retRl.bound = rl.bound + rl.length[0]; retRl.size= rl.length[0]; retRl.serial = 1; }
			}
		}	
	} else if(flag == BACKWARD) {
		if(rl.paired == 1)
		{
			if(sa < rl.bound || ( sa >= rl.bound + (rl.length[0] + rl.length[1] *2))) {  retRl.bound = -1; }
			else {
				retRl.paired = 1 ;
				if(sa < rl.bound  + rl.length[0]) { retRl.bound = rl.bound; retRl.size = rl.length[0]; retRl.serial = 0; }
				else if(sa < rl.bound + rl.length[0] + rl.length[1]) {retRl.bound = rl.bound + rl.length[0]; retRl.size = rl.length[1]; retRl.serial = 1;  }
				else if(sa < rl.bound + rl.length[0] + rl.length[1] * 2) { retRl.bound = rl.bound + rl.length[0] + rl.length[1]; retRl.size = rl.length[1]; retRl.serial = 2; }
				else { retRl.bound = rl.bound + rl.length[0] + rl.length[1] * 2; retRl.size = rl.length[0]; retRl.serial = 3; }
			}
		} else { // rl.paired == 0
			if(sa < rl.bound || (sa >= rl.bound + rl.length[0] * 2)) { retRl.bound = -1; }
			else {
				retRl.paired = 0 ;
				if(sa < rl.bound + rl.length[0]) { retRl.bound = rl.bound; retRl.size = rl.length[0]; retRl.serial = 0; }
				else { retRl.bound = rl.bound + rl.length[0]; retRl.size = rl.length[0]; retRl.serial = 1; }
			}
		}
	} else { // flag == BIDIRECTION
		if(rl.paired == 1)
		{
			if(sa >= rl.bound + (rl.length[0] + rl.length[1]) * 2) { retRl.bound = -1 ; }
			else {
				retRl.paired = 1 ;
				if(sa <= rl.bound + rl.length[0]/2) { retRl.bound = rl.bound; retRl.size = rl.length[0]; retRl.serial = 0; }
				else if(sa <= rl.bound + rl.length[0] + rl.length[1]/2) { retRl.bound = rl.bound + rl.length[0]; retRl.size = rl.length[1]; retRl.serial = 1; }
				else if(sa <= rl.bound + rl.length[0] + rl.length[1] + rl.length[1]/2) { retRl.bound = rl.bound + rl.length[0] + rl.length[1]; retRl.size = rl.length[1]; retRl.serial = 2; }
				else if(sa <= rl.bound + rl.length[0] + rl.length[1] * 2 + rl.length[0]/2) { retRl.bound = rl.bound + rl.length[0] + rl.length[1] * 2; retRl.size = rl.length[0]; retRl.serial = 3; }
				else { // error
					fprintf(stderr, "[LocateStartSA] sa >rl.bound + rl.length[0] + rl.length[1] * 2 + rl.length[0]/2, exit...\n");
					exit(1);
				}
			}
		} else { // rl.paired == 0
			if(sa >= rl.bound + rl.length[0] *2)  { retRl.bound = -1; }
			else {
				retRl.paired = 0;
				if(sa <= rl.bound + rl.length[0]/2) { retRl.bound = rl.bound; retRl.size = rl.length[0]; retRl.serial = 0; }
				else { retRl.bound = rl.bound + rl.length[0]; retRl.size = rl.length[0]; retRl.serial = 1; }
			}	
		}
	}
	return retRl;
}


int checkSA(const ReadLocation rl, const int64_t a, const int fixed_len)
{
	int64_t retain = a - rl.bound ;
	if(rl.paired == 1)
	{
		if(retain > rl.length[0] - fixed_len && retain < rl.length[0]) { return 0 ; }
		else if(retain > rl.length[0] + rl.length[1] - fixed_len && retain  < rl.length[0] + rl.length[1]) {return 0;}
		else if(retain > rl.length[0] + rl.length[1] * 2 - fixed_len && retain < rl.length[0] + rl.length[1] * 2) { return 0; }
		else if(retain > (rl.length[0] + rl.length[1]) * 2 - fixed_len && retain < (rl.length[0] + rl.length[1]) * 2) { return 0; }
		else return 1 ;
	} else {
		if(retain > rl.length[0] - fixed_len && retain < rl.length[0]) { return 0 ; }
		else if(retain > rl.length[0] * 2 - fixed_len && retain  < rl.length[0] * 2) {return 0;}
		else return 1 ;
	}
}
/*// flag : 0 note find the head of read , 1 note find tail of read 
int64_t lookUpBound(const LIB *lib, int64_t position, const int direction , const int flag)
{
    position -= lib->offset ;
    if(lib->diverse == 1)
    {
        int64_t len,lastLen, i;
        int64_t lenMajorSize = (lib->number_rd + LEN_INTERVAL_MAJOR - 1) / LEN_INTERVAL_MAJOR  + 1; 
        int64_t low = 0 ,mid,high  = lenMajorSize - 1;
        while(low <= high)
        {
            mid = (low + high) / 2 ;
            if(position == lib->lenMajor[mid]) break ;
            if(position < lib->lenMajor[mid]) high = mid - 1 ;
            else low = mid + 1 ;
        }
        if(lib->lenMajor[mid] < position) {low = mid; high = mid+1; }
        else if(lib->lenMajor[mid] > position) {low = mid - 1; high = mid; }
        else {
            if(direction == FORWARD)
            {
                if(flag == 0) return lib->offset + lib->lenMajor[mid] ;
                else { low = mid ; high = mid+1;  }
            } else if(direction == BACKWARD) {
                if(flag == 0 ) return lib->offset + lib->lenMajor[mid] ;
                else return lib->offset + lib->lenMajor[mid] - 1 ;
            } else {
                if(flag == 0) return lib->offset + lib->lenMajor[mid] ;
                else return lib->offset + lib->lenMajor[mid] - 1 ;
            }
        }
        if((position - lib->lenMajor[low]) <= ( lib->lenMajor[high] - position)) // bi-dirtection search
        {
            len = (lastLen = lib->lenMajor[low]) ;
            i = low * LEN_INTERVAL_MAJOR ;
            while(len < position)
            {
                lastLen = len ;
                len += lib->length[i];
                i++ ;
            }
            if(direction == FORWARD) 
            {
                if(flag == 0) return lib->offset +  len ;
                //else return  position == len ? lib->length[i] - 1 : (uint64_t)len - 1 ;
                else return  position == len ? lib->offset + len + lib->length[i] - 1 : lib->offset + len - 1 ;
            }else if(direction == BACKWARD) {
                if(flag == 0) return position == len ? lib->offset + len : lib->offset + lastLen ;
                else return position == len ? lib->offset + len - 1 : lib->offset + lastLen - 1 ;
            } else {
                if((position - lastLen) < (len - position)) 
                {
                    if(flag == 0) return lib->offset + lastLen ;
                    else return lib->offset + lastLen - 1 ;
                }else {
                    if(flag == 0) return lib->offset + len ;
                    else return lib->offset + len - 1 ;
                }
            }
        }else {
            len = (lastLen = lib->lenMajor[high]);
            if(high >= lenMajorSize - 1) i = lib->number_rd - 1 ;
            else  i = high * LEN_INTERVAL_MAJOR - 1 ;
            while(len > position)
            {
                lastLen = len ;
                len -= lib->length[i];
                i-- ;
            }
            if(direction == FORWARD) 
            {
                if(flag == 0) return position == len ? lib->offset + len : lib->offset + lastLen ;
                else return lib->offset + lastLen - 1 ;
            }else if(direction == BACKWARD )
            {
                if(flag == 0) return lib->offset + len ;
                else return lib->offset + len - 1 ;
            }else {
                if((position - len) < (lastLen - position)) 
                {
                    if(flag == 0 ) return lib->offset + len ;
                    else return lib->offset + len - 1 ;
                }else {
                    if(flag == 0) return lib->offset + lastLen ;
                    else return lib->offset + lastLen - 1 ;
                }
            }
        }
    } else {
        int64_t num_offset = position / lib->read_len ;
        if(direction == FORWARD)
        {
            if(flag == 0) return num_offset * lib->read_len == position  ? lib->offset + position : (lib->offset + (num_offset + 1) * lib->read_len) ;
            else return lib->offset + (num_offset + 1) * lib->read_len - 1 ;
        }else if(direction == BACKWARD) {
            if(flag == 0) return lib->offset + num_offset * lib->read_len ; 
            else return position == (num_offset + 1) * lib->read_len -1 ? lib->offset + position : (lib->offset + num_offset * lib->read_len - 1) ;  
        }else {
            if((position - lib->read_len * num_offset) < (lib->read_len * (num_offset + 1) - position))
            {
                if(flag == 0) return  lib->offset + num_offset * lib->read_len ;
                else return lib->offset + num_offset * lib->read_len - 1 ;
            }else {
                if(flag == 0) return lib->offset + (num_offset + 1) * lib->read_len ;
                else return lib->offset + (num_offset + 1) * lib->read_len - 1 ;
            }
        }
    }
}

int64_t readBoundaryLookUp(const lib_info *libIndex_info,const int64_t position,const int direction)
{
    int64_t boundaryPos;
    int i ;
    int flag  ;  // 0 note head of the read , 1 note tail of the read 
    
    if(position < libIndex_info->len_pac) { boundaryPos = position; flag = 0; }
    else { boundaryPos = libIndex_info->len_pac * 2 - 1 - position ; flag = 1 ;  }
    //boundaryPos = position < libIndex_info->len_pac ? position : (libIndex_info->len_pac * 2 - 1 - position)
    for(i = 0; i < libIndex_info->num_lib ; i++)
    {
        if(libIndex_info->lib[i].offset > boundaryPos) break ;
    }

    if((position < libIndex_info->len_pac  && direction == FORWARD) || (position >= libIndex_info->len_pac && direction == BACKWARD))
    {
        boundaryPos = lookUpBound(&(libIndex_info->lib[i-1]), boundaryPos, FORWARD , flag);
    }else if((position < libIndex_info->len_pac && direction == BACKWARD) || (position >= libIndex_info->len_pac && direction == FORWARD)) {
        boundaryPos = lookUpBound(&(libIndex_info->lib[i-1]), boundaryPos, BACKWARD , flag);
    } else {
        boundaryPos = lookUpBound(&(libIndex_info->lib[i-1]), boundaryPos, BIDIRECTION , flag);
    }

    // check if boundaryPos is legal 
    if(boundaryPos < 0 ) boundaryPos = 0 ; 
    
    if(position < libIndex_info->len_pac) { return boundaryPos; }
    else { return libIndex_info->len_pac * 2 - 1- boundaryPos ;  }
    //return position < libIndex_info->len_pac ? boundaryPos : (libIndex_info->len_pac * 2 - 1 - boundaryPos + 1);
} */


// write libIndex_info to .ann file
void writeInfo2Ann(const char *fn, lib_info *libIndex_info)
{
    FILE *fp ;
    int i  ;
    int64_t j ;

    fp = xopen(fn, "w");

    // write libIndex_info entries
    fprintf(fp, "len_pac      = %ld\n", libIndex_info->len_pac);
    fprintf(fp, "num_seqs     = %ld\n", libIndex_info->num_seqs);
    fprintf(fp, "num_lib      = %d\n", libIndex_info->num_lib);

	// search maximum read length of libs and write to .ann file
	{
		int max_rd_len = 0;
		for(int i = 0; i < libIndex_info->num_lib; i++)
		{
			if(libIndex_info->lib[i].diverse == 0) 
			{
				if(libIndex_info->lib[i].read_len > max_rd_len) { max_rd_len = libIndex_info->lib[i].read_len; }
			}else {
				for(int j = 0; j < libIndex_info->lib[i].number_rd; j++)
				{
					if(libIndex_info->lib[i].length[j] > max_rd_len) { max_rd_len = libIndex_info->lib[i].length[j]; }
				}
			}	
		}
		// write to .ann file
		fprintf(fp, "max_rd_len		=	%d\n", max_rd_len);
	}

    // write pair end library info 
    for(i = 0; i < libIndex_info->num_lib; i++)
    {
        fprintf(fp, "#############%s################\n", libIndex_info->lib[i].name);
        fprintf(fp, "lib_name     = %s\n", libIndex_info->lib[i].name);
        fprintf(fp, "paired       = %u\n", libIndex_info->lib[i].paired);
        if(libIndex_info->lib[i].paired == 1)
        {
            fprintf(fp, "insert_size  = %u\n", libIndex_info->lib[i].pe.insert_size);
            fprintf(fp, "insert_SD    = %u\n", libIndex_info->lib[i].pe.insert_SD);
            //fprintf(fp, "reverse      = %u\n", libIndex_info->lib[i].pe.reverse);
        }
        fprintf(fp, "offset       = %lu\n", libIndex_info->lib[i].offset);
        fprintf(fp, "read_len     = %u\n", libIndex_info->lib[i].read_len);
        fprintf(fp, "number_rd    = %u\n", libIndex_info->lib[i].number_rd);
        fprintf(fp, "diverse      = %u\n", libIndex_info->lib[i].diverse);
        if(libIndex_info->lib[i].diverse == 1 )
        {
            // check length buffer
            if(libIndex_info->lib[i].length == NULL)
            {
                fprintf(stderr, "[writeInfo2Ann] the 'libIndex_info->lib[%d].length' is NULL, please check, program exit...\n", i);
                exit(1);
            }

            uint64_t totalLen = 0 ;
            uint64_t lenMajorSize = (libIndex_info->lib[i].number_rd + LEN_INTERVAL_MAJOR -1) / LEN_INTERVAL_MAJOR + 1 ;
            libIndex_info->lib[i].lenMajor = (uint64_t*)xcalloc(lenMajorSize, sizeof(uint64_t));
            // write length array and lenMajor array 
            fprintf(fp, "length       =\t");
            for(j = 0 ; j < libIndex_info->lib[i].number_rd ; j++)
            {
                if(j % LEN_INTERVAL_MAJOR == 0)
                {
                    libIndex_info->lib[i].lenMajor[ j / LEN_INTERVAL_MAJOR ] = totalLen ;
                }
                totalLen += libIndex_info->lib[i].length[j] ;
                fprintf(fp, "%u\t",libIndex_info->lib[i].length[j]);
            }
            fprintf(fp, "\n");
            libIndex_info->lib[i].lenMajor[lenMajorSize - 1] = totalLen ;

            fprintf(fp, "lenMajor     =\t");
            for(j = 0 ; j < lenMajorSize ; j++)
            {
                fprintf(fp, "%lu\t", libIndex_info->lib[i].lenMajor[j]);
            }
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
}

lib_info *lib_info_restore(const char *filename)
{
    lib_info *libIndex_info = (lib_info*)xcalloc(1, sizeof(lib_info));
    FILE *fp ;
    int i ;
    int64_t j ;
    char s[PATH_LEN];

    fp = xopen(filename, "r") ;
    
    // restore lib_info content 
    fgets(s, PATH_LEN, fp); strtok(s, " \t\n"); strtok(NULL, " \t\n"); libIndex_info->len_pac = atol(strtok(NULL, " \t\n")); 
    fgets(s, PATH_LEN, fp); strtok(s, " \t\n"); strtok(NULL, " \t\n"); libIndex_info->num_seqs = atol(strtok(NULL, " \t\n")); 
    fgets(s, PATH_LEN, fp); strtok(s, " \t\n"); strtok(NULL, " \t\n"); libIndex_info->num_lib = atol(strtok(NULL, " \t\n")); 
    fgets(s, PATH_LEN, fp); strtok(s, " \t\n"); strtok(NULL, " \t\n"); libIndex_info->max_rd_len = atol(strtok(NULL, " \t\n")); 
    //fscanf(fp, "%13s%1s %lu\n", tmp, tmp1,  libIndex_info->len_pac);
    //fscanf(fp, "num_seqs=%lu\n", libIndex_info->num_seqs);
    //fscanf(fp, "num_lib       = %u\n", libIndex_info->num_lib);
    libIndex_info->lib = (LIB*)xcalloc(libIndex_info->num_lib, sizeof(LIB));
    // restore pair end library info 
    for(i = 0 ; i < libIndex_info->num_lib ; i++)
    {
        libIndex_info->lib[i].name = (char*)xcalloc(PATH_LEN,1);
        fgets(s, PATH_LEN, fp);
        fgets(s, PATH_LEN, fp); strtok(s, " \t\n"); strtok(NULL, " \t\n"); libIndex_info->lib[i].name = strdup(strtok(NULL, " \t\n")); 
        fgets(s, PATH_LEN, fp); strtok(s, " \t\n"); strtok(NULL, " \t\n"); libIndex_info->lib[i].paired = atoi(strtok(NULL, " \t\n")); 
        //fscanf(fp, "#############%s################\n", libIndex_info->lib[i].name);
        //fscanf(fp, "lib_name     = %s\n", libIndex_info->lib[i].name);
        //fscanf(fp, "paired       = %u\n", libIndex_info->lib[i].paired);
        if(libIndex_info->lib[i].paired == 1)
        {
            fgets(s, PATH_LEN, fp); strtok(s, " \t\n"); strtok(NULL, " \t\n"); libIndex_info->lib[i].pe.insert_size = atoi(strtok(NULL, " \t\n")); 
            fgets(s, PATH_LEN, fp); strtok(s, " \t\n"); strtok(NULL, " \t\n"); libIndex_info->lib[i].pe.insert_SD = atoi(strtok(NULL, " \t\n")); 
            //fgets(s, PATH_LEN, fp); strtok(s, " \t\n"); strtok(NULL, " \t\n"); libIndex_info->lib[i].pe.reverse = atoi(strtok(NULL, " \t\n")); 
            //fscanf(fp, "insert_size  = %u\n", libIndex_info->lib[i].pe.insert_size);
            //fscanf(fp, "insert_SD    = %u\n", libIndex_info->lib[i].pe.insert_SD);
            //fscanf(fp, "reverse      = %u\n", libIndex_info->lib[i].pe.reverse);
        }
        fgets(s, PATH_LEN, fp); strtok(s, " \t\n"); strtok(NULL, " \t\n"); libIndex_info->lib[i].offset = atol(strtok(NULL, " \t\n")); 
        fgets(s, PATH_LEN, fp); strtok(s, " \t\n"); strtok(NULL, " \t\n"); libIndex_info->lib[i].read_len = atoi(strtok(NULL, " \t\n")); 
        fgets(s, PATH_LEN, fp); strtok(s, " \t\n"); strtok(NULL, " \t\n"); libIndex_info->lib[i].number_rd = atoi(strtok(NULL, " \t\n")); 
        fgets(s, PATH_LEN, fp); strtok(s, " \t\n"); strtok(NULL, " \t\n"); libIndex_info->lib[i].diverse = atoi(strtok(NULL, " \t\n")); 
        //fscanf(fp, "offset       = %lu\n", libIndex_info->lib[i].offset);
        //fscanf(fp, "read_len     = %u\n", libIndex_info->lib[i].read_len);
        //fscanf(fp, "number_rd    = %u\n", libIndex_info->lib[i].number_rd);
        //fscanf(fp, "diverse      = %u\n", libIndex_info->lib[i].diverse);
        if(libIndex_info->lib[i].diverse == 1 )
        {
            uint64_t lenMajorSize = (libIndex_info->lib[i].number_rd + LEN_INTERVAL_MAJOR -1) / LEN_INTERVAL_MAJOR + 1 ;
            char *bs = (char*)xcalloc(libIndex_info->lib[i].number_rd * 12 , sizeof(char));
            libIndex_info->lib[i].lenMajor = (uint64_t*)xcalloc(lenMajorSize, sizeof(uint64_t));
            libIndex_info->lib[i].length = (READ_SIZE*)xcalloc(libIndex_info->lib[i].number_rd, sizeof(READ_SIZE));
            
            // restore length array and lenMajor array 
            fgets(bs, libIndex_info->lib[i].number_rd * 10, fp ); strtok(bs, " \t\n"); strtok(NULL, " \t\n");
            for(j = 0 ; j < libIndex_info->lib[i].number_rd ; j++)
            {
                libIndex_info->lib[i].length[j] = atoi(strtok(NULL, " \t\n"));
            }

            fgets(bs, libIndex_info->lib[i].number_rd * 10, fp ); strtok(bs, " \t\n"); strtok(NULL, " \t\n");
            for(j = 0 ; j < lenMajorSize ; j++)
            {
                libIndex_info->lib[i].lenMajor[j] = atol(strtok(NULL, " \t\n"));
            }
            free(bs);
        }
    }
    fclose(fp);

    return libIndex_info ;
}

void lib_infoFree(lib_info *libIndex_info)
{
    for(int i = 0 ; i < libIndex_info->num_lib; i++)
    {
        free(libIndex_info->lib[i].name);
        if(libIndex_info->lib[i].diverse == 1)
        {
            free(libIndex_info->lib[i].length);
            free(libIndex_info->lib[i].lenMajor);
        }
    }
    free(libIndex_info->lib);
    free(libIndex_info);
}
