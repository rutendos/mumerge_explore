#!/bin/bash 

#SBATCH --job-name=nextflow # Job name        
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)  
#SBATCH --mail-user=rutendo.sigauke@cuanschutz.edu # Where to send mail                                                           
#SBATCH --nodes=1   
#SBATCH --ntasks=1 # Number of CPU (processer cores i.e. tasks)               
#SBATCH --time=04:00:00 # Time limit hrs:min:sec                              
#SBATCH --partition short                                                     
#SBATCH --mem=8gb # Memory limit                 
#SBATCH --output=/scratch/Users/rusi2317/e_and_o/%x_%j.out             
#SBATCH --error=/scratch/Users/rusi2317/e_and_o/%x_%j.err

################### SET VARIABLES ######################################            

export PATH=~:$PATH

MAIN=/Users/rusi2317/bin/ChIP-Flow/main.nf
CONFIG=/Users/rusi2317/bin/ChIP-Flow/conf/fiji.config
PROJECT=/scratch/Users/rusi2317/projects/mumerge_test_p53/analysis


########################################################################                                                   
################### LOAD NECESSARY MODULES #############################                                                           
module purge
module load sra/2.8.0
module load samtools/1.8
module load hisat2/2.1.0
module load bedtools/2.28.0
module load gcc/7.1.0
module load seqkit/0.9.0
module load fastqc/0.11.8
module load bbmap/38.05
module load igvtools/2.3.75
module load preseq/2.0.3
module load openmpi/1.6.4
emodule load python/3.6.3
module load mpich/3.2.1
module load picard/2.6.0
module load python/3.6.3/rseqc/3.0.0

######################################################################## 
################## PRINT JOB INFO ######################################                                                           
printf "Sample ID: $ROOTNAME"
printf "\nDirectory: $PROJECT"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"
##########################################################################                                           
#################### Run across samples in the srr #######################
#slurm_grch38 
nextflow run $MAIN -profile fiji --sras '/scratch/Users/rusi2317/projects/mumerge_test_p53/raw_data/*.sra' \
    --singleEnd \
    --workdir ${PROJECT}/temp/ \
    --outdir ${PROJECT}/processed_data/ \
    --email rutendo.sigauke@cuanschutz.edu \
    --threadfqdump \
    --dedup \
    --singleEnd \
    --savebg \
    --savebw \
    --skipMultiQC \
    -resume


    
