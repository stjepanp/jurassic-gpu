#nvcc main.c -L ../submodule -lutil -I../submodule/inter_folder -o prog
gcc main.c vector.c mul_vector.c -L /p/software/jusuf/stages/2020/modules/all/Compiler/GCC/9.3.0 -lstdc++ -L /p/software/jusuf/stages/2020/software/CUDA/11.0/lib64 -lcudart  -L ../submodule -lutil -I../submodule/inter_folder -o prog
