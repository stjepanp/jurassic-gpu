#!/bin/bash -x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --output=out
#SBATCH --error=err
#SBATCH --time=1:00:00
#SBATCH --partition=gpus
#SBATCH --gres=gpu:1
#SBATCH --account=slmet

module load Python/3.8.5
module load GCC/9.3.0
module load OpenMPI/4.1.0rc1
module load CUDA/11.3
module load Nsight-Systems/2021.1.1
module load Nsight-Compute/2020.3.0
module load ParaStationMPI/5.4.7-1
module load mpi-settings/CUDA
module load GCCcore/.10.3.0
module load GCCcore/.9.3.0
module load GSL/2.6 


#! /bin/bash

# Setup...
jurassic=../../src
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../lib/build/lib:../../lib/build/lib64

# Create atmospheric data file...
OMP_NUM_THREADS=1 $jurassic/climatology nadir.ctl atm.tab

# Create observation geometry...
OMP_NUM_THREADS=1 $jurassic/nadir nadir.ctl obs.tab T1 10

rm -f rad.tab
# Call forward model...
# reference=~/Codes/JURASSIC/reference-jurassic/jurassic/src
# OMP_NUM_THREADS=1 $reference/formod nadir.ctl obs.tab atm.tab rad.tab TASK time
OMP_NUM_THREADS=1 $jurassic/formod nadir.ctl obs.tab atm.tab rad.tab TASK time CHECKMODE $1

# Plot results... brightness temperature vs. tangent point latitude
gnuplot <<EOF
set term png enh truecolor font "Helvetica,28" size 1600,1200 crop lw 2
set out "plot.png"
set xla "latitude [deg]"
set yla "brightness temperature [K]"
set mxtics
set mytics
set key spac 1.5
set key outside
plot "rad.org" u 10:11 w lp pt 1 t  "ref (667.8 cm^{-1})", \
     "rad.tab" u 10:11 w lp pt 2 t "test (667.8 cm^{-1})", \
     "rad.org" u 10:12 w lp pt 1 t  "ref (668.5 cm^{-1})", \
     "rad.tab" u 10:12 w lp pt 2 t "test (668.5 cm^{-1})", \
     "rad.org" u 10:13 w lp pt 1 t  "ref (669.8 cm^{-1})", \
     "rad.tab" u 10:13 w lp pt 2 t "test (669.8 cm^{-1})"
EOF

# Plot results... transmittance vs. tangent point latitude
gnuplot <<EOF
set term png enh truecolor font "Helvetica,28" size 1600,1200 crop lw 2
set out "plot2.png"
set xla "latitude [deg]"
set yla "transmittance [%]"
set mxtics
set mytics
set key spac 1.5
set key outside
plot "rad.org" u 10:(\$14*100) w lp pt 1 t  "ref (667.8 cm^{-1})", \
     "rad.tab" u 10:(\$14*100) w lp pt 2 t "test (667.8 cm^{-1})", \
     "rad.org" u 10:(\$15*100) w lp pt 1 t  "ref (668.5 cm^{-1})", \
     "rad.tab" u 10:(\$15*100) w lp pt 2 t "test (668.5 cm^{-1})", \
     "rad.org" u 10:(\$16*100) w lp pt 1 t  "ref (669.8 cm^{-1})", \
     "rad.tab" u 10:(\$16*100) w lp pt 2 t "test (669.8 cm^{-1})"
EOF

# Get differences...
echo -e "\nCheck for differences..."
diff -sq rad.tab rad.org
