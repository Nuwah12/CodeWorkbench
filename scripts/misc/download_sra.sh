#!/bin/bash

####################
# Script to dpwnload a bunch of files from SRA (just a simple loop through a text file of accession numbers)  
####################

USAGE="Usage: ./download_sra.sh --samplesheet/-s <samplesheet> --sratools/-r <SRAtools bin location>"

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
		-s | --samplesheet )
			s="$2"
			shift 2 # Shift past both argument and value
			;;
		-r | --sratools )
			r="$2"
			shift 2
			;;
		-* | --* )
			print_usage_exit $1
			;;
	esac
done

if [[ ! -d $r ]]; then
        echo "Cannot find provided path to SRAtools bin: $r"
        echo $USAGE
        exit 1
fi

while read srr; do
	echo "$(date): Downloading $srr"
	$r/fasterq-dump $srr --progress
	echo "$(date): Finished downloading $srr"
done < $s


