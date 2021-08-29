#ifndef VECTOR_H
#define VECTOR_H

#define MAXN 1000

typedef struct {
  int n;
  int array[MAXN];
  int x; //because I want structs to differ
} vector;

vector *create_vector(int n, int a, int b);

void print_vector(vector *v);

#endif
