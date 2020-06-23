#!/bin/bash
#SBATCH --job-name=peak_calling_with_macs2  # Job name
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rutendo.sigauke@cuanschutz.edu # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of CPU (processer cores i.e. tasks) In this example I use 32 for bowtie2. I only need one, since none of the commands I run are parallelized.
#SBATCH --time=08:00:00 # Time limit hrs:min:sec
#SBATCH -p short
#SBATCH --mem=10gb # Memory limit
#SBATCH --output=/scratch/Users/rusi2317/e_and_o/%x.%j.out
#SBATCH --error=/scratch/Users/rusi2317/e_and_o/%x.%j.err

####################################################################
################### LOAD MODULES ###################################

module load python/2.7.14/MACS/2.1.1
module load bedtools/2.25.0

####################################################################
################## SET VARIABLES ###################################
BAM_DIR='/scratch/Users/rusi2317/projects/mumerge_test_p53/analysis/processed_data/mapped/bams'
ROOTNAME=$1

##1->HCT116, 2->MCF7, 3->SJSA
##input samples
CTRL1=SRR4090089
CTRL2=SRR4090092
CTRL3=SRR4090095

##DMSO samples
DMSO1=SRR4090090
DMSO2=SRR4090093
DMSO3=SRR4090096

##Nutlin treated samples
NTLN1=SRR4090091
NTLN2=SRR4090094
NTLN3=SRR4090097

OUT='/scratch/Users/rusi2317/projects/mumerge_test_p53/peak_calls'
BLACKLIST='/scratch/Shares/dowell/genomes/hg38/wgEncodeDacMapabilityConsensusExcludable_grch38.bed'

####################################################################
################## JOB INFO ########################################

printf "\nControl ID: $CTRL1, $CTRL2, $CNTL3"
printf "\nDMSO ID: $DMSO1, $DMSO2, DMSO3"
printf "\nNutlin ID: $NTLN1, $NTLN2, $NTLN3"
printf "\nBAM Directory: $BAM_DIR"
printf "\nOutput Directory: $OUT"
printf "\nOutput Rootname: $ROOTNAME"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

#######################################################################
#######################Call Peaks with MACS2###########################

echo macs2
date
date

##name ids
HCT116_1=HCT116_DMSO_TP53
HCT116_2=HCT116_Nutlin_TP53

MCF7_1=MCF7_DMSO_TP53
MCF7_2=MCF7_Nutlin_TP53

SJSA_1=SJSA_DMSO_TP53
SJSA_2=SJSA_Nutlin_TP53

##peak calling
##HCT116_1
macs2 callpeak \
 -c $BAM_DIR/${CTRL1}.sorted.bam \
 -t $BAM_DIR/${DMSO1}.sorted.bam \
 --outdir $OUT \
 -n ${HCT116_1} \
 -g hs \
 -f BAM \
 -q 0.001 
 
echo removing blacklist regions
date
date
 
bedtools intersect \
 -a $OUT/${HCT116_1}_peaks.broadPeak \
 -b ${BLACKLIST} \
 -v \
 > $OUT/${HCT116_1}_peaks_clean.broadPeak

##HCT116_2
macs2 callpeak \
 -c $BAM_DIR/${CTRL1}.sorted.bam \
 -t $BAM_DIR/${NTLN1}.sorted.bam \
 --outdir $OUT \
 -n ${HCT116_2} \
 -g hs \
 -f BAM \
 -q 0.001 
 
echo removing blacklist regions
date
date
 
bedtools intersect \
 -a $OUT/${HCT116_2}_peaks.broadPeak \
 -b ${BLACKLIST} \
 -v \
 > $OUT/${HCT116_2}_peaks_clean.broadPeak

##MCF7
##MCF7_1
macs2 callpeak \
 -c $BAM_DIR/${CTRL2}.sorted.bam \
 -t $BAM_DIR/${DMSO2}.sorted.bam \
 --outdir $OUT \
 -n ${MCF7_1} \
 -g hs \
 -f BAM \
 -q 0.001 
 
echo removing blacklist regions
date
date
 
bedtools intersect \
 -a $OUT/${MCF7_1}_peaks.broadPeak \
 -b ${BLACKLIST} \
 -v \
 > $OUT/${MCF7_1}_peaks_clean.broadPeak

##MCF7_2
macs2 callpeak \
 -c $BAM_DIR/${CTRL2}.sorted.bam \
 -t $BAM_DIR/${NTLN2}.sorted.bam \
 --outdir $OUT \
 -n ${MCF7_2} \
 -g hs \
 -f BAM \
 -q 0.001 
 
echo removing blacklist regions
date
date
 
bedtools intersect \
 -a $OUT/${MCF7_2}_peaks.broadPeak \
 -b ${BLACKLIST} \
 -v \
 > $OUT/${MCF7_2}_peaks_clean.broadPeak


##SJSA
##SJSA_1
macs2 callpeak \
 -c $BAM_DIR/${CTRL3}.sorted.bam \
 -t $BAM_DIR/${DMSO3}.sorted.bam \
 --outdir $OUT \
 -n ${SJSA_1} \
 -g hs \
 -f BAM \
 -q 0.001 
 
echo removing blacklist regions
date
date
 
bedtools intersect \
 -a $OUT/${SJSA_1}_peaks.broadPeak \
 -b ${BLACKLIST} \
 -v \
 > $OUT/${SJSA_1}_peaks_clean.broadPeak

##HCT116_2
macs2 callpeak \
 -c $BAM_DIR/${CTRL3}.sorted.bam \
 -t $BAM_DIR/${NTLN3}.sorted.bam \
 --outdir $OUT \
 -n ${SJSA_2} \
 -g hs \
 -f BAM \
 -q 0.001 
 
echo removing blacklist regions
date
date
 
bedtools intersect \
 -a $OUT/${SJSA_2}_peaks.broadPeak \
 -b ${BLACKLIST} \
 -v \
 > $OUT/${SJSA_2}_peaks_clean.broadPeak


echo macs2 done
date
date



