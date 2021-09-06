#include "stdio.h"
#include "stdlib.h"

#define __host__

#include "interface.h"

#include "vector.h"
#include "mul_vector.h"
#include "misc.h"

int main(void) {
  vector *a = create_vector(7, 1, 1);
  vector *b = create_vector(7, 2, 2);
  print_vector(a);
  print_vector(b);

  printf("CPU multiplication...\n");
  mul_vector(a, b);
  print_vector(a);
  print_vector(b);
  
  printf("GPU addition...\n");
  run_vector_addition(a, b);
  print_vector(a);
  print_vector(b);
  return 0;
}
