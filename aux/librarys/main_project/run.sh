#!/bin/bash -x
#SBATCH --account=slmet
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --output=out-gpu.txt
#SBATCH --error=err-gpu.txt
###SBATCH --output=out-gpu.%j
###SBATCH --error=err-gpu.%j
#SBATCH --time=01:59:00
#SBATCH --gres=gpu:1
#SBATCH --partition=gpus

cp out-gpu.txt prev-out.txt
cp err-gpu.txt prev-err.txt

module --force purge
#module load Stages/2020 GCC/9.3.0 CUDA/11.0
module load GCC/9.3.0
module load Python/3.8.5
module load OpenMPI/4.1.0rc1
module load CUDA/11.3
module load Nsight-Systems/2021.1.1
module load Nsight-Compute/2020.3.0
module load ParaStationMPI/5.4.7-1
module load mpi-settings/CUDA
module load GCCcore/.10.3.0
module load GCCcore/.9.3.0
module load GSL/2.6 

#export
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/p/home/jusers/baumeister1/jusuf/tfQMRgpu/build/tfQMRgpu
srun ./prog
