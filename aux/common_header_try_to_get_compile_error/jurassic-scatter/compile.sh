cd ../jurassic-gpu
./compile.sh
cd -
rm prog
gcc prog.c vector.c misc.c mul_vector.c -L /p/software/jusuf/stages/2020/modules/all/Compiler/GCC/9.3.0 -lstdc++ -L /p/software/jusuf/stages/2020/software/CUDA/11.0/lib64 -lcudart  -L ../jurassic-gpu -lutil -I../jurassic-gpu/inter_folder -I ../common_header -o prog