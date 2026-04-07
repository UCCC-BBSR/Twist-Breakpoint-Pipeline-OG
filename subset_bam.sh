#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --output=job_outs/BAMsubset-%j.out
#SBATCH --error=job_outs/BAMsubset-%j.err
#SBATCH --qos normal

# Goals
# Different than my March analysis, I will first subset to IGV regions and then to discordant reads

# Load samtools module
module load samtools


### Source config to get paths
. config_loose.sh # Using $project, $pipeoutdir    
path_main_directory=${pipeoutdir}"/preprocessing/recalibrated"
output=${project}"/bam_subset_igv_discordant"
bed_igv_ranges=${project}"/igv_ranges.bed"

# Create output directories if they don't exist
mkdir -p "$output"

# List subfolders within path_main_directory and write to variable samples
samples=$(ls -d "$path_main_directory"/*/)

# Loop through each sample                                                                                                                                                                                                     
for sample_dir in $samples; do

    # Get basename                                                                                                                                                                                                                                                          
    sample=$(basename "$sample_dir")
    echo "Processing sample: $sample"
    
    # For testing only                                                                                                                                                                                                         
    #samtools view -b -h -L "${bed_igv_ranges}" "${sample_dir}/${sample}.recal.bam" > "${output}/postIGV_${sample}.bam"
    #samtools index "${output}/postIGV_${sample}.bam"
    
    # Pipe V2                                                                                                                                                                                                                  
    samtools view -b -h -L "${bed_igv_ranges}" "${sample_dir}/${sample}.recal.bam" | \
        samtools view -h | awk '$7 != "=" || $1 ~ /^@/' | samtools view -bS > "${output}/${sample}.bam"

    samtools index "${output}/${sample}.bam"
    
done
