#include "gpu_summation.h"

#include <cuda.h>
#include <stdio.h>

__global__ void sum(int *a, int *b, int n) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i < n) {
    a[i] += b[i];
  }
}
