#!/bin/bash

####################
# Script to make a ChromSizes file for a supplied fasta file
###################

USAGE="Usage: $0 -f/--fasta <fasta> [--3col] [-o/--outfile <ChromSizes.txt>]"

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

use_2col=1
use_3col=0
out_name="ChromSizes.txt"

while [[ $# -gt 0 ]]; do
	case "$1" in
		-f | --fasta )
			fasta="$2"
			shift 2 # Shift past both argument and value
			;;
		--3col )
			use_2col=0
			use_3col=1
			shift
			;;
		-o | --outfile )
			out_name="$2"
			shift 2
			;;
		-* | --* )
			print_usage_exit $1
			;;
	esac
done


awk -v use_2col="$use_2col" -v use_3col="$use_3col" '
    BEGIN { seqname=""; seqlen=0; }
    /^>/ {
        if (seqname != "") {
            if (use_2col == 1)
                print seqname "\t" seqlen;
            else if (use_3col == 1)
                print seqname "\t" 1 "\t" seqlen;
        }
        split($0, header, " ");
        sub(/^>/, "", header[1]);
        seqname = header[1];
        seqlen = 0;
        next;
    }
    {
        gsub(/\s+/, "", $0);
        seqlen += length($0);
    }
    END {
        if (seqname != "") {
            if (use_2col == 1)
                print seqname "\t" seqlen;
            else if (use_3col == 1)
                print seqname "\t" 1 "\t" seqlen;
        }
    }
' "$fasta"
