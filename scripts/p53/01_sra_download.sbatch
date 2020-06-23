#!/bin/bash 

#SBATCH --job-name=01_sra_downloads # Job name
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)

#SBATCH --mail-user=rutendo.sigauke@cuanschutz.edu
#SBATCH --array=1-9 # job array index                                                                                                                       

#SBATCH --nodes=1
#SBATCH --ntasks=1 # Num of CPU

#SBATCH --time=00:10:00 # Time limit
#SBATCH --partition short

#SBATCH --mem=10mb # Memory ##whats the limit in memory required

#SBATCH --output=/scratch/Users/rusi2317/e_and_o/%x_%j.out              
#SBATCH --error=/scratch/Users/rusi2317/e_and_o/%x_%j.err   

module load sra/2.10.5

##load reference srrs
srr_paths=/Users/rusi2317/projects/mumerge_test/data/andrysik2017chip_metadata_20200430.tsv

srr=($(cat ${srr_paths} | cut -f 2))

vdb-config --interactive

prefetch ${srr[${SLURM_ARRAY_TASK_ID}]} 

#prefetch SRR4090089
