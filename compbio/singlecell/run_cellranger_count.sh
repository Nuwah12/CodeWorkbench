#!/bin/bash

####################
# Script to run cellranger count on multiple (10X) fastq files
####################

USAGE="Usage: ./run_cellranger_count -f <fastq dir> -a <cellranger assembly dir> -c <path to cellranger binary, not binary itself>"

function print_usage_exit()
{
	if [[ ! -z $1 ]]; then
		echo -e "Unknown parameter '$1' \n${USAGE}"
		exit 1
	fi
	echo ${USAGE}
	exit 1
}

if [[ $# -eq 0 ]]; then
	echo $USAGE
	exit 1
fi

while [[ $# -gt 0 ]]; do
	case "$1" in
		-f | --fastq-directory ) 
			dir="$2" 
			shift 2 # Shift past both argument and value 
			;; 
		-a | --assembly ) 
			assembly="$2"
			shift 2 
			;;
		-c | --cellranger ) 
			cranger="$2"
			shift 2 
			;;
		-* | --* ) 
			print_usage_exit $1
			;;
	esac
done

declare -a sample_names=() # declare empty array for sample names

echo "Gather sample names"
for file in $dir/*.fastq.gz; do
	name=$(basename $file) # Strip path
	sample=${name%%_*} # Strip suffix after first _
	sample_names+=("$sample") # Append to array
done

eval uniq_samples=($(printf "%q\n" "${sample_names[@]}" | sort -u))

echo "${uniq_samples[@]}"

echo "Running cellranger count on each sample separately..."
echo "Cellranger location: $cranger"
for s in "${uniq_samples[@]}"; do
	echo "Aligning sample $s"
	${cranger}/cellranger count --id ${s} \
		--transcriptome ${assembly}\
		--sample ${s} \
		--fastqs ${dir} \
		--localcores 32 \
		--localmem 128 \
		--create-bam false
done

