#!/bin/sh                                                        
#SBATCH --job-name=merging_significant_regions
#SBATCH -p short                    
#SBATCH --nodes=1                                                                                                                              
#SBATCH --ntasks=1                                                                                                                            
#SBATCH --mem=1gb                                                                                                                            
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
bed=/scratch/Users/rusi2317/projects/mumerge_test_p53/processed_data/significant_regions_whole_seq
outdir=/scratch/Users/rusi2317/projects/mumerge_test_p53/processed_data/merged_regions_full_len

cat ${bed}/*.bed | bedtools sort -i stdin | bedtools merge -i stdin > ${outdir}/Andrysik2017_merged.bed


##################################################################
######################Intersect then Merge #######################

echo '----------Combining bed files with intersect merge---------'
bedtools intersect -a ${bed}/HCT116_DMSO_TP53_sigFimo.bed -b ${bed}/MCF7_DMSO_TP53_sigFimo.bed ${bed}/SJSA_DMSO_TP53_sigFimo.bed > ${outdir}/dmso_intersect.bed
bedtools intersect -a ${bed}/HCT116_Nutlin_TP53_sigFimo.bed -b ${bed}/MCF7_Nutlin_TP53_sigFimo.bed ${bed}/SJSA_Nutlin_TP53_sigFimo.bed > ${outdir}/nutlin_intersect.bed
                                                                                                                                                                                           
cat ${outdir}/*intersect.bed | bedtools sort -i stdin | bedtools merge -i stdin > ${outdir}/Andrysik2017_intersect_merged.bed

rm ${outdir}/*_intersect.bed

##################################################################
######################Intersect then Merge #######################

echo '----------Combining bed files with intersect all-----------'
bedtools intersect -a ${bed}/HCT116_DMSO_TP53_sigFimo.bed -b ${bed}/MCF7_DMSO_TP53_sigFimo.bed ${bed}/SJSA_DMSO_TP53_sigFimo.bed ${bed}/MCF7_Nutlin_TP53_sigFimo.bed ${bed}/SJSA_Nutlin_TP53_sigFimo.bed ${bed}/HCT116_Nutlin_TP53_sigFimo.bed > ${outdir}/Andrysik2017_intersect_all.bed

echo '----------Combining bed files with intersect merge---------'
bedtools intersect -a ${bed}/HCT116_DMSO_TP53_sigFimo.bed -b ${bed}/MCF7_DMSO_TP53_sigFimo.bed > ${outdir}/ab_intersect.bed
bedtools intersect -a ${outdir}/ab_intersect.bed -b ${bed}/SJSA_DMSO_TP53_sigFimo.bed > ${outdir}/bc_intersect.bed
bedtools intersect -a ${outdir}/bc_intersect.bed -b ${bed}/HCT116_Nutlin_TP53_sigFimo.bed > ${outdir}/cd_intersect.bed
bedtools intersect -a ${outdir}/cd_intersect.bed -b ${bed}/MCF7_Nutlin_TP53_sigFimo.bed > ${outdir}/ef_intersect.bed
bedtools intersect -a ${outdir}/ef_intersect.bed -b ${bed}/SJSA_Nutlin_TP53_sigFimo.bed > ${outdir}/Andrysik2017_intersect_recursive.bed

rm ${outdir}/*_intersect.bed

##################################################################
######################## MuMerge  ################################
mumerge_path=/scratch/Users/rusi2317/bin/TFEA/TFEA
data_dir=/scratch/Users/rusi2317/projects/mumerge_test_p53/processed_data

echo '-----------Combining bed files with muMerge---------------'
##run the mumerge                                                                                                                                                                           
python3 ${mumerge_path}/mumerge.py -i ${data_dir}/andrysik2017_samples.txt -o ${outdir}/Andrysik2017_mumerged.bed -v

