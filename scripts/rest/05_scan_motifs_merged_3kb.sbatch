#!/bin/sh                                                        
#SBATCH --job-name=fimo_scan_rest_merged3kb
#SBATCH -p short                    
#SBATCH --nodes=1                                                                                                                              
#SBATCH --ntasks=1                                                                                                                            
#SBATCH --mem=2gb                                                                                                                            
#SBATCH --time=00:05:00                                                                                                                        
#SBATCH --output /scratch/Users/rusi2317/e_and_o/%x.out                                                               
#SBATCH --error /scratch/Users/rusi2317/e_and_o/%x.err                                                                
#SBATCH --mail-type=ALL                                                                    
#SBATCH --mail-user=rutendo.sigauke@cuanschutz.edu

module load meme/5.1.1
module load bedtools/2.28.0

############################################################
###################get fasta sequences###################### 
hg19=/scratch/Shares/public/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
chip_fa=/scratch/Users/rusi2317/projects/mumerge_test_rest/processed_data/chip_fasta
merged=/scratch/Users/rusi2317/projects/mumerge_test_rest/processed_data/merged_regions_full_len
merged3kb=/scratch/Users/rusi2317/projects/mumerge_test_rest/processed_data/merged_regions3kb

for bed in `ls ${merged} | grep -v "log"`
do
    echo 'Getting fasta sequencies for  --------------------------->' ${bed}
    
    ##get new coordinates and remove regions with -                                                                                                       
    cat ${merged}/${bed} | awk '{printf ("%s %d %d\n",$1,(($3-$2)/2)+$2-1500,(($3-$2)/2)+1+$2+1500)}' | grep -v '-' | tr ' ' '\t' > ${merged3kb}/${bed}

    cat ${merged3kb}/${bed} | bedtools getfasta -fi ${hg19} -bed stdin -fo ${chip_fa}/${bed}_3kb.fa

done


#############################################################
######################scan fimo motifs#######################
outdir=/scratch/Users/rusi2317/projects/mumerge_test_rest/processed_data/fimo_scans_merged3kb
hocomoco=/scratch/Users/rusi2317/projects/mumerge_test_p53/hocomoco 
motiffile=REST_HUMAN.H11MO.0.A

for bed in `ls ${merged} | grep -v "log"`
do
    echo 'Scanning the REST motif for  ------------------------------>' ${bed}

    fimo -max-stored-scores 1000000 --thresh 1e-4 -oc ${outdir}/${bed} -motif $motiffile ${hocomoco}/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme ${chip_fa}/${bed}_3kb.fa

done

