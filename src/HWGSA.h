#ifndef HWGSA_MAIN_H
#define HWGSA_MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwtindex.h"
#include "generateGraph.h"
#include "findOverlap.h"
//#include "extend.h"

#ifdef __cplusplus
extern "C" {
#endif

    int hwgsa_extend(int argc, char *argv[]);
    int hwgsa_scaffold(int argc, char *argv[]);
    int hwgsa_all(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif
