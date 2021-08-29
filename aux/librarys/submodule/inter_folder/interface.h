#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#define MAXN_CONN 265

typedef struct {
  int n;
  int array[MAXN_CONN];
} vector_connection;


#ifdef __cplusplus
extern "C" {
#endif

__host__ void run_vector_addition(vector_connection *a_conn, vector_connection *b_conn); 

#ifdef __cplusplus
} //extern "C"
#endif

#endif
