#ifndef VECTOR_H
#define VECTOR_H

#include "jr_scatter_gpu.h"

__host__ vector *create_vector(int n, int a, int b);

__host__ void print_vector(vector *v);

#endif
