#include "inter_folder/interface.h"

#include "vector.h"
#include "sum_vector.h"

#include <stdio.h>

__host__ void run_vector_addition(vector *a, vector *b) {
  add_vector(a, b);
  print_vector(a);
  print_vector(b);
}
