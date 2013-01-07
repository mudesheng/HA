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


#ifndef LH3_UTILS_H
#define LH3_UTILS_H

#include <stdio.h>
#include <zlib.h>
#include <stdint.h>
    
extern const uint64_t masklow64[] ;
extern const uint32_t masklow32[];
extern const uint8_t masklow8[] ;
extern const uint8_t clean_bnt[] ;

#define INITIAL_STACK_SIZE 2000
#define err_fatal_simple(msg) err_fatal_simple_core(__func__, msg)
#define xopen(fn, mode) err_xopen_core(__func__, fn, mode)
#define xreopen(fn, mode, fp) err_xreopen_core(__func__, fn, mode, fp)
#define xzopen(fn, mode) err_xzopen_core(__func__, fn, mode)
#define xassert(cond, msg) if ((cond) == 0) err_fatal_simple_core(__func__, msg)

#define SETBNT(allv, c, freq) (((allv) & clean_bnt[(c)]) | ((freq)<<(((~(c)) & 0x3)<<1)))
#define GETBNT(table, index) (((table)[(index)>>2] >> (((~(index)) & 0x3)<<1)) & 0x3)


#ifndef PATH_LEN
#define PATH_LEN 1024 // file path length
#endif

#ifdef __cplusplus
extern "C" {
#endif
    void *xcalloc(size_t nmemb , size_t size);
    void *xrecalloc(void *ptr , size_t size);
	void err_fatal(const char *header, const char *fmt, ...);
	void err_fatal_simple_core(const char *func, const char *msg);
	FILE *err_xopen_core(const char *func, const char *fn, const char *mode);
	FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp);
	gzFile err_xzopen_core(const char *func, const char *fn, const char *mode);

#ifdef __cplusplus
}
#endif

#endif
