#!/bin/sh                                                        
#SBATCH --job-name=fimo_scan
#SBATCH -p short                    
#SBATCH --nodes=1                                                                                                                              
#SBATCH --ntasks=1                                                                                                                            
#SBATCH --mem=2gb                                                                                                                            
#SBATCH --time=00:07:00                                                                                                                        
#SBATCH --output /scratch/Users/rusi2317/e_and_o/%x.out                                                               
#SBATCH --error /scratch/Users/rusi2317/e_and_o/%x.err                                                                
#SBATCH --mail-type=ALL                                                                    
#SBATCH --mail-user=rutendo.sigauke@cuanschutz.edu

module load meme/5.1.1 #4.12.0
module load bedtools/2.28.0

############################################################
###################get fasta sequences###################### 

hg38=/scratch/Shares/public/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
p53_chip=/scratch/Users/rusi2317/projects/mumerge_test_p53/peak_calls
chip_fa=/scratch/Users/rusi2317/projects/mumerge_test_p53/processed_data/chip_fasta


for bed in `ls ${p53_chip} | grep "_peaks.narrowPeak"`
do
    echo 'Getting fasta sequencies for  --------------------------->' ${bed}

    cat ${p53_chip}/${bed} | bedtools getfasta -fi ${hg38} -bed stdin -fo ${chip_fa}/${bed}_whole_seq_macs.fa

done


#############################################################
######################scan fimo motifs#######################
outdir=/scratch/Users/rusi2317/projects/mumerge_test_p53/processed_data/fimo_scans_whole_seq
hocomoco=/scratch/Users/rusi2317/projects/mumerge_test_p53/hocomoco 
motiffile=P53_HUMAN.H11MO.0.A

for bed in `ls ${p53_chip} | grep "_peaks.narrowPeak"`
do
    echo 'Scanning the p53 motif for  ------------------------------>' ${bed}

    fimo -max-stored-scores 1000000 --thresh 1e-4 -oc ${outdir}/${bed} -motif $motiffile ${hocomoco}/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme ${chip_fa}/${bed}_whole_seq_macs.fa

done

