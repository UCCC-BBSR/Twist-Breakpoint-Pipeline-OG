#!/bin/bash
#SBATCH -J sarek
#SBATCH --partition amilan
#SBATCH --account=amc-general
#SBATCH -o log/sarek_%j.out  # Output file with the job ID
#SBATCH -e log/sarek_%j.err  # Error file with the job ID
#SBATCH -t 00-10:00:00   # Set the wall time: D-HH:MM:SS
#SBATCH --qos=normal
#SBATCH -n 1 -c 2  # ask for number of nodes/cores
#SBATCH --mem=18GB  # Specify memory allocation

set -euo pipefail

# make software accessible:
module load nextflow/23.04
module load singularity/3.7.4

echo ____________________________________
nextflow -version
echo ____________________________________
singularity --version
echo ____________________________________

## Source configuration
. config_loose.sh

## Override defaults set by nextflow module:
export NXF_WORK="$scrproj/work"
export NXF_TEMP="$scrproj/tmp"
export NXF_HOME="/projects/$USER/software/nextflow_config"
export _JAVA_OPTIONS='-Xmx16G' # increase heap to avoid java.lang.OutOfMemoryError

## Make sample file:
#samplefile="$project/all_samples.csv"
#[ -f "$samplefile" ] ||
#./make_samplefile.sh "${fastq_raw[@]}" >"$samplefile" || exit 1

## nextflow puts everything in $PWD, so make an output dir and go there:
mkdir -p "$pipeoutdir"
cd "$pipeoutdir"

nextflow run "$pipelinedir" -profile curc_alpine -ansi-log false \
	--genome GATK.GRCh38 --input "$samplefile" --trim_fastq \
	--save_trimmed --save_mapped --save_output_as_bam --outdir "$pipeoutdir" \
	--tools tiddit,snpeff -resume \
	-c "$project/scripts/fastp.config" \
	-c "$project/scripts/tiddit_loose_v3.config" \
	-c /pl/active/bbsrinformatics/analysts/andrew/scripts/nextflow_configs/lower_cpu.config
# manta doesn't work on alpine
# ascat,controlfreec need sex info
# msisensorpro needs T/N
