#include <stdio.h>
#include <stdlib.h>

#include "vector.h"

__host__ vector *create_vector(int n, int a, int b) {
  vector *v;
  v = (vector*) malloc(sizeof(vector));
  v->n = n;
  for(int i = 0; i < n; i++)
    v->array[i] = a * i + b;
  return v;
}

__host__ void print_vector(vector *v) {
  for(int i = 0; i < v->n; i++)
    printf("%d ", v->array[i]); 
  printf("\n");
}
