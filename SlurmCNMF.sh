#!/bin/bash

##SBATCH --job-name= PUT_JOB_NAME_HERE (any identifier for your job)
##SBATCH --output=gpu_job.txt
##SBATCH --ntasks=1
##SBATCH --gres=gpu:1
#SBATCH --partition=day
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=60G

# the h5ad file must be in the the directory called adataframes (the same directory that contains runCNMF.py)
# make sure you have the CellBenderEnv environment activated before running it  

python /home/PUT_YOUR_NETID_HERE/adataframes/runCNMF.py adataframes PUT_NAME_OF_h5ad_FILE_HERE
