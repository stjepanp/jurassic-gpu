#include "mul_vector.h"

#include <assert.h>

void mul_vector(vector *a, const vector *b) {
  assert(a->n == b->n);
  for(int i = 0; i < a->n; i++)
    a->array[i] *= b->array[i];
}
