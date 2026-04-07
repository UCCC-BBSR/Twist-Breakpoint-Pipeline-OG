# This file needs to be soured, so the variables are available to the calling script
set -eu							# die on any error

# AG added sample csv:
samplefile=/pl/active/bbsrinformatics/analysts/mkaufman/projects/Hayashi_Twist_Jan2026/pipeline/scripts/samples.csv

# directories:
projbn=pipeline 
scratch=/scratch/alpine/$USER
plhome=/pl/active/bbsrinformatics/analysts/mkaufman/projects/Hayashi_Twist_Jan2026
scriptdir=~tdanhorn@xsede.org/Programs
project="$plhome/$projbn"
scrproj="$scratch/$projbn"
dataroot=/pl/active/bbsrinformaticsdata/bbsr
#fqroot="$dataroot/Hayashi"
#fastq_raw=("$fqroot"/WGS_sarcoma_compiled_Mar2024) #Path to raw fastq files
guidir="$plhome/proj/_gui-files_/${project##*/proj/}"
pipeoutdir="$project/nf-core_sarek_loose"

# reference files:
refpath=/pl/active/bbsrinformatics/references
fasta="$refpath/nf-core/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
fastQscreen_config="$refpath/genomes/FastQ_Screen_Genomes/alpine_default-with-Mycoplasma_fastq_screen.conf" #Location of config file
annovardb="$refpath/ANNOVAR/hg38"
## auxiliary files for making results tables:
formatconfig="$project/SNP-indel_TN_comparison_format-config.yaml"
# from ~danhornt/proj/_gui-files_/_templates_/SNP-indel_TN_comparison_format-config.yaml
keyfile="$project/SNP-indel_TN_comparison_column-key.txt"
# from ~danhornt/proj/_gui-files_/_templates_/SNP-indel_TN_comparison_column-key.txt

# commands for tools:
swpath=/projects/tdanhorn@xsede.org/bbsr/software
contdir="$swpath/singularity_containers"
container="$contdir/VCtools_2023-06.sif"
# bind mount of various directories into container, so we can see reference files, etc.
containercmd="singularity exec -B $scratch,$plhome,$refpath $container"
rcontainercmd="${containercmd%/*}/R-4.2.2_bioc-3.16.sif"
fqscreencmd="${containercmd%/*}/RNAseq_2021-08.sif fastq_screen" # FastQscreen tool
multiqccmd="${containercmd%/*}/multiqc_v1.11.sif multiqc"
bcftcmd="$containercmd bcftools"
annovarcmd="$containercmd table_annovar.pl"
pipelinedir="$swpath/nf-core/nf-core-sarek-3.3.2/3_3_2"

# Filtering parameters:
minreads=10
minvafdiff=0.05 # minimum difference in normal & tumor VAF for FreeBayes results
minvafratio=3

# Setting for FastQ Screen:
maxfastqscreenallreads=1000000 # scan at most this many reads from full sample
