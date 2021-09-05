rm libutil.a
gcc -c jurassic.c -o jurassic.o
nvcc -c interface.cu vector.cu sum_vector.cu gpu_summation.cu -gencode arch=compute_70,code=sm_70 -I ../common_header
ar rc libutil.a *.o
ranlib libutil.a
rm *.o
