/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <ctype.h>
#include "utils.h"

char *KStrip(char *s)
{
	char *strip_s = s, *p = s + strlen(s) - 1; // p pointer to end of string s
	if(s == NULL)
	{
		fprintf(stderr, "[KStrip] the input string is NULL, please check, exit.....\n");
		exit(2) ;
	}

	while(isblank(*strip_s)) { strip_s++ ; } // strip start blank character    	
	while(p >= s && (isspace(*p) || *p == ';')) {*p = '\0' ; p--; }

	return strip_s ;
}




const uint64_t masklow64[ ] = {
    0x0000000000000000LLU , 0x0000000000000003LLU , 0x000000000000000FLLU , 0x000000000000003FLLU ,
    0x00000000000000FFLLU , 0x00000000000003FFLLU , 0x0000000000000FFFLLU , 0x0000000000003FFFLLU ,
    0x000000000000FFFFLLU , 0x000000000003FFFFLLU , 0x00000000000FFFFFLLU , 0x00000000003FFFFFLLU ,
    0x0000000000FFFFFFLLU , 0x0000000003FFFFFFLLU , 0x000000000FFFFFFFLLU , 0x000000003FFFFFFFLLU ,
    0x00000000FFFFFFFFLLU , 0x00000003FFFFFFFFLLU , 0x0000000FFFFFFFFFLLU , 0x0000003FFFFFFFFFLLU ,
    0x000000FFFFFFFFFFLLU , 0x000003FFFFFFFFFFLLU , 0x00000FFFFFFFFFFFLLU , 0x00003FFFFFFFFFFFLLU ,
    0x0000FFFFFFFFFFFFLLU , 0x0003FFFFFFFFFFFFLLU , 0x000FFFFFFFFFFFFFLLU , 0x003FFFFFFFFFFFFFLLU ,
    0x00FFFFFFFFFFFFFFLLU , 0x03FFFFFFFFFFFFFFLLU , 0x0FFFFFFFFFFFFFFFLLU , 0x3FFFFFFFFFFFFFFFLLU ,
    0xFFFFFFFFFFFFFFFFLLU
} ;

const uint32_t masklow32[ ] = {
    0X00000000 , 0X00000003 , 0X0000000F , 0X0000003F , 
    0X000000FF , 0X000003FF , 0X00000FFF , 0X00003FFF , 
    0X0000FFFF , 0X0003FFFF , 0X000FFFFF , 0X003FFFFF ,
    0X00FFFFFF , 0X03FFFFFF , 0X0FFFFFFF , 0X3FFFFFFF ,
    0XFFFFFFFF
} ;

const uint8_t masklow8[ ] = {
    0x00 , 0x03 , 0x0F , 0x3F , 0xFF 
} ;

const uint8_t clean_bnt[] = {
	0X3F, 0XCF, 0XF3, 0XFC 
} ;

void *xcalloc(size_t nmemb, size_t size)
{
    void *p ;
    if(nmemb <= 0) return NULL ;
    p = calloc(nmemb , size) ;
    if(p == NULL)
    {
        fprintf(stderr, "calloc memory require failed ,abort.....\n");
        fflush(stdout); fflush(stderr);
        abort();
    }
    return p ;
}

// realloc the memory and clean the last half part of memory
void *xrecalloc(void *ptr, size_t size)
{
    void *p ;
    p = realloc(ptr, size);
    if(p == NULL)
    {
        fprintf(stderr, "recalloc memory require failed, abort.....\n");
        abort();
    }
    memset(p + (size>>1), 0 , (size>>1));
    return p ;
}

FILE *err_xopen_core(const char *func, const char *fn, const char *mode)
{
	FILE *fp = 0;
	if (strcmp(fn, "-") == 0)
		return (strstr(mode, "r"))? stdin : stdout;
	if ((fp = fopen(fn, mode)) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s'. Abort!\n", func, fn);
		abort();
	}
	return fp;
}
FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp)
{
	if (freopen(fn, mode, fp) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s': ", func, fn);
		perror(NULL);
		fprintf(stderr, "Abort!\n");
		abort();
	}
	return fp;
}
gzFile err_xzopen_core(const char *func, const char *fn, const char *mode)
{
	gzFile fp;
	if (strcmp(fn, "-") == 0)
		return gzdopen(fileno((strstr(mode, "r"))? stdin : stdout), mode);
	if ((fp = gzopen(fn, mode)) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s'. Abort!\n", func, fn);
		abort();
	}
	return fp;
}
void err_fatal(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, " Abort!\n");
	va_end(args);
	abort();
}

void err_fatal_simple_core(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s Abort!\n", func, msg);
	abort();
}
