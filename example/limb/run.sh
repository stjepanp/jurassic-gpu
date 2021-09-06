#!/bin/bash -x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
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

#!/usr/bin/env bash

# Setup...
jurassic=../../src
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../lib/build/lib:../../lib/build/lib64

# Create atmospheric data file...
OMP_NUM_THREADS=1 $jurassic/climatology limb.ctl atm.tab

# Create observation geometry...
OMP_NUM_THREADS=1 $jurassic/limb limb.ctl obs.tab Z0 3 Z1 68 DZ 1.0

rm -f rad.tab
# Call forward model...
reference=~/Codes/JURASSIC/reference-jurassic/jurassic/src
# OMP_NUM_THREADS=1 $reference/formod limb.ctl obs.tab atm.tab rad.tab TASK time
OMP_NUM_THREADS=1 $jurassic/formod limb.ctl obs.tab atm.tab rad.tab TASK time CHECKMODE $1

# Plot results... tangent height vs. radiance
gnuplot <<EOF
set term png enh truecolor font "Helvetica,28" size 1600,1200 crop lw 2
set out "plot.png"
set xla "radiance [nW/(cm^2 sr cm^{-1})]"
set yla "tangent height [km]"
set mxtics
set mytics
set log x
set key spac 1.5
set key left bot
plot "rad.org" u (\$11*1e5):8 w lp pt 1 t  "ref (792 cm^{-1})", \
     "rad.org" u (\$12*1e5):8 w lp pt 1 t  "ref (832 cm^{-1})", \
     "rad.tab" u (\$11*1e5):8 w lp pt 2 t "test (792 cm^{-1})", \
     "rad.tab" u (\$12*1e5):8 w lp pt 2 t "test (832 cm^{-1})"
EOF

# Plot results... transmittance vs. radiance
gnuplot <<EOF
set term png enh truecolor font "Helvetica,28" size 1600,1200 crop lw 2
set out "plot2.png"
set xla "radiance [nW/(cm^2 sr cm^{-1})]"
set yla "transmittance [%]"
set mxtics
set mytics
set log x
set key spac 1.5
set key left bot
plot "rad.org" u (\$11*1e5):(\$13*100) w lp pt 1 t  "ref (792 cm^{-1})", \
     "rad.org" u (\$12*1e5):(\$14*100) w lp pt 1 t  "ref (832 cm^{-1})", \
     "rad.tab" u (\$11*1e5):(\$13*100) w lp pt 2 t "test (792 cm^{-1})", \
     "rad.tab" u (\$12*1e5):(\$14*100) w lp pt 2 t "test (832 cm^{-1})"
EOF

# Plot results... view point altitude vs. radiance
gnuplot <<EOF
set term png enh truecolor font "Helvetica,28" size 1600,1200 crop lw 2
set out "plot3.png"
set xla "radiance [nW/(cm^2 sr cm^{-1})]"
set yla "view point altitude [km]"
set mxtics
set mytics
set log x
set key spac 1.5
set key left bot
plot "rad.org" u (\$11*1e5):5 w lp pt 1 t  "ref (792 cm^{-1})", \
     "rad.org" u (\$12*1e5):5 w lp pt 1 t  "ref (832 cm^{-1})", \
     "rad.tab" u (\$11*1e5):5 w lp pt 2 t "test (792 cm^{-1})", \
     "rad.tab" u (\$12*1e5):5 w lp pt 2 t "test (832 cm^{-1})"
EOF

# Get differences...
echo -e "\nCheck for differences..."
diff -sq rad.tab rad.org
