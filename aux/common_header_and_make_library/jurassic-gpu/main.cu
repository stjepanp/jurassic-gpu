#include "stdio.h"

#include "jr_scatter_gpu.h"
#include "sum_vector.h"

int main() {
  vector *a = create_vector(10, 2, 3);
  vector *b = create_vector(10, 1, 0);
  print_vector(a); print_vector(b);
  add_vector(a, b);
  print_vector(a); print_vector(b);
  return 0;
}
