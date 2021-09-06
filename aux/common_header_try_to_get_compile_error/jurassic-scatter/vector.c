#include "vector.h"

#include <stdio.h>
#include <stdlib.h>

vector *create_vector(int n, int a, int b) {
  vector *v;
  v = (vector*) malloc(sizeof(vector));
  v->n = n;
  for(int i = 0; i < n; i++)
    v->array[i] = a * i + b;
  return v;
}

void print_vector(vector *v) {
  printf("CPU printing vector... ");
  for(int i = 0; i < v->n; i++)
    printf("%d ", v->array[i]); 
  printf("\n");
}

