nvcc -c inter_folder/interface.cu vector.cu sum_vector.cu gpu_summation.cu
ar rc libutil.a *.o
ranlib libutil.a
rm *.o
