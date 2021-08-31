#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include "jr_scatter_gpu.h"

#ifdef __cplusplus
extern "C" {
#endif

__host__ void run_vector_addition(vector *a, vector *b); 

#ifdef __cplusplus
} //extern "C"
#endif

#endif
