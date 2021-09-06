./compile.sh
read -n 1 -s -r -p "Press any key to continue"
sbatch run.sh
watch -n 1 squeue -u pozgaj1
cat out-gpu.txt
