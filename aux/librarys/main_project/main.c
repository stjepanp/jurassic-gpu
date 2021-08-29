#include "stdio.h"
#include "stdlib.h"

#define __host__

#include "interface.h"

#include "vector.h"
#include "mul_vector.h"

__host__ void main_convert_conn_to_normal(vector_connection *v_c, vector *v_n) {
  v_n->n = v_c->n;
  //v_n->array = v_c->array;
  int n = v_n->n;
  for(int i = 0; i < n; i++)
    v_n->array[i] = v_c->array[i];
}

__host__ void main_convert_normal_to_conn(vector *v_n, vector_connection *v_c) {
  v_c->n = v_n->n;
  //v_c->array = v_n->array;
  int n = v_n->n;
  for(int i = 0; i < n; i++)
    v_c->array[i] = v_n->array[i];
}

int main() {
  vector *a = create_vector(7, 1, 1);
  vector *b = create_vector(7, 2, 2);
  print_vector(a);
  print_vector(b);

  printf("CPU multiplication...\n");
  mul_vector(a, b);
  print_vector(a);
  print_vector(b);
  
  printf("GPU addition...\n");
  vector_connection *a_conn = (vector_connection*) malloc(sizeof(vector_connection));
  vector_connection *b_conn = (vector_connection*) malloc(sizeof(vector_connection));
  main_convert_normal_to_conn(a, a_conn);
  main_convert_normal_to_conn(b, b_conn);
  run_vector_addition(a_conn, b_conn);
  main_convert_conn_to_normal(a_conn, a);
  main_convert_conn_to_normal(b_conn, b);
  print_vector(a);
  print_vector(b);
  return 0;
}
