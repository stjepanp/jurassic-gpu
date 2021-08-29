#include "sum_vector.h"

#include <assert.h>
#include <cuda.h>
#include <stdio.h>
#include <unistd.h>

#include "gpu_summation.h"

__host__ void add_vector(vector *a, const vector *b) {
  assert(a->n == b->n);
  int *a_gpu, *b_gpu;
  cudaMalloc(&a_gpu, a->n * sizeof(int));
  cudaMalloc(&b_gpu, b->n * sizeof(int));
  cudaMemcpy(a_gpu, a->array, a->n * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(b_gpu, b->array, b->n * sizeof(int), cudaMemcpyHostToDevice); 
  dim3 blockDim(64, 1, 1);
  dim3 gridDim(a->n / 64 + 1);
  sum <<<gridDim, blockDim>>> (a_gpu, b_gpu, a->n);
  cudaError_t err = cudaGetLastError();
  if(err) {
    printf("Error!\n");
    exit(0);
  }
  cudaMemcpy(a->array, a_gpu, a->n * sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(a_gpu);
  cudaFree(b_gpu);
}
