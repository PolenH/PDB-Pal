#!/bin/bash
#SBATCH -A b1039
#SBATCH -p buyin 
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem-per-cpu=650M
#SBATCH -t 01:00:00
#SBATCH --job-name="mm_opt"
#SBATCH --output=mm_opt.txt

module purge
module load orca/5.0.4 

export OMPI_MCA_btl=^vader,tcp,openib

/software/orca/5.0.4/orca_5_0_4_linux_x86-64_openmpi411/orca mm_opt.inp > mm_opt.out
