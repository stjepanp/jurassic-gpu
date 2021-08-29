#ifndef VECTOR_H
#define VECTOR_H

#define MAXN 1000

typedef struct {
  int n;
  int array[MAXN];
  int y;
} vector;

__host__ vector *create_vector(int n, int a, int b);

__host__ void print_vector(vector *v);

#endif
