#!/bin/sh                                                        
#SBATCH --job-name=merging_significant_regions
#SBATCH -p short                    
#SBATCH --nodes=1                                                                                                                              
#SBATCH --ntasks=1                                                                                                                            
#SBATCH --mem=150mb                                                                                                                            
#SBATCH --time=00:05:00                                                                                                                        
#SBATCH --output /scratch/Users/rusi2317/e_and_o/%x.out                                                               
#SBATCH --error /scratch/Users/rusi2317/e_and_o/%x.err                                                                
#SBATCH --mail-type=ALL                                                                    
#SBATCH --mail-user=rutendo.sigauke@cuanschutz.edu

module purge
module load python/3.6.3
module load python/3.6.3/numpy
module load python/3.6.3/scipy
module load bedtools/2.28.0


###################################################################
################### Merge all samples #############################

echo '--------------Combining bed files with merge----------------'
bed=/scratch/Users/rusi2317/projects/mumerge_test_rest/processed_data/significant_regions_whole_seq
outdir=/scratch/Users/rusi2317/projects/mumerge_test_rest/processed_data/merged_regions_full_len

cat ${bed}/ENCFF001UN* | bedtools sort -i stdin | bedtools merge -i stdin > ${outdir}/MCF7_merged.bed

##################################################################
######################Intersect then Beds  #######################

echo '---------Combining bed files with bedtools intersect--------'
bedtools intersect -a ${bed}/ENCFF001UNH.bed_sigFimo.bed -b ${bed}/ENCFF001UNI.bed_sigFimo.bed > ${outdir}/MCF7_intersect.bed

##################################################################
######################## MuMerge  ################################
mumerge_path=/scratch/Users/rusi2317/bin/TFEA/TFEA
data_dir=/scratch/Users/rusi2317/projects/mumerge_test_rest/processed_data

echo '-----------Combining bed files with muMerge---------------'
echo 'MCF7 samples'
python3 ${mumerge_path}/mumerge.py -i ${data_dir}/01_mcf7_samples.txt -o ${outdir}/mcf7_mumerged.bed -v
 
