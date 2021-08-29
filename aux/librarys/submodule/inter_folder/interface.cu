#include "interface.h"

#include "../vector.h"
#include "../sum_vector.h"

#include <stdio.h>

__host__ void interface_convert_conn_to_normal(vector_connection *v_c, vector *v_n) {
  v_n->n = v_c->n;
  //v_n->array = v_c->array;
  int n = v_n->n;
  for(int i = 0; i < n; i++)
    v_n->array[i] = v_c->array[i];
}

__host__ void interface_convert_normal_to_conn(vector *v_n, vector_connection *v_c) {
  v_c->n = v_n->n;
  //v_c->array = v_n->array;
  int n = v_n->n;
  for(int i = 0; i < n; i++)
    v_c->array[i] = v_n->array[i];
}

__host__ void run_vector_addition(vector_connection *a_conn, vector_connection *b_conn) {
  vector *a = (vector*) malloc(sizeof(vector));
  vector *b = (vector*) malloc(sizeof(vector));
  interface_convert_conn_to_normal(a_conn, a);
  interface_convert_conn_to_normal(b_conn, b);
  add_vector(a, b);
  interface_convert_normal_to_conn(a, a_conn);
  interface_convert_normal_to_conn(b, b_conn);
}
