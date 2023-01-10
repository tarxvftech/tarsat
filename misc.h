#ifndef MM_MISC_H
#define MM_MISC_H
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

typedef size_t errorcount;
typedef bool success_or_failure; //where true is success
typedef success_or_failure sof;
typedef int success_or_errorcode; //where 0 is success
typedef success_or_errorcode soe;
#define Pf(t) printf("%s = %f\n", #t, t)

#endif
