#!/bin/bash

#SBATCH --account=co2
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4700m
#SBATCH -n 1
#SBATCH --exclusive
#SBATCH --output=/work2/noaa/co2/kaushik/sib4/sib4_corral_v4.2-iso/run_jobs/%x_%a.o%j
#SBATCH --qos=batch
#SBATCH --partition=bigmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aleya.kaushik@noaa.gov

module purge
module load intel-oneapi-compilers/2024.1.0
module load intel-oneapi-mpi/2021.12.0
module load intel-oneapi-mkl
module load hdf5/1.14.3
module load netcdf-c/4.9.0
module load netcdf-fortran/4.6.0

make
