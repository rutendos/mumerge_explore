#!/bin/sh                                                        
#SBATCH --job-name=chip_rest_summit
#SBATCH -p short                    
#SBATCH --nodes=1                                                                                                                              
#SBATCH --ntasks=1                                                                                                                            
#SBATCH --mem=50mb                                                                                                                            
#SBATCH --time=00:02:00                                                                                                                        
#SBATCH --output /scratch/Users/rusi2317/e_and_o/%x.out                                                               
#SBATCH --error /scratch/Users/rusi2317/e_and_o/%x.err                                                                
#SBATCH --mail-type=ALL                                                                    
#SBATCH --mail-user=rutendo.sigauke@cuanschutz.edu

############################################################
###################get fasta sequences###################### 
rest_chip=/scratch/Users/rusi2317/projects/mumerge_test_rest/raw_data

for bed in `ls ${rest_chip} | grep ".bed"`
do
    echo 'Getting theoretical summits --------------------------->' ${bed}
    cat ${rest_chip}/${bed} | awk 'BEGIN{FS=OFS="\t"} {print $1, $2+int(($3 - $2)/2), $2+int(($3 - $2)/2)+1, $4, int(($3 - $2))}' > ${rest_chip}/${bed}_summits.bed 

done

