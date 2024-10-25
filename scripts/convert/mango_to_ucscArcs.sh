#!/bin/bash

if [[ $# -lt 4 ]]; then
	echo "Usage: ./mango_to_ucscArcs.sh <.mango file> <output file name> <interact.as> <chromSizes>"
	exit 1
fi

mango=$1
out=$2
ucscTools="/mnt/data0/apps/ucscTools/exec"
formatFile=$3
chromSizes=$4

#header='track type=interact name="Arc Interactions" description="Description useScore=on maxHeightPixels=200:100:50 visibility=full'
#region='chr3:64,562,440-64,642,288'

### Text conversion
awk -v region="${region}" -v header="${header}" -F"\t" 'BEGIN {OFS="\t"} NR>1 {print $1,$2,$6,$7,$10,$17,".",0,$1,$2,$3,$1,".",$4,$5,$6,$4,"."}' ${mango} > ${out}

### Make binary "bigInteract"
${ucscTools}/bedToBigBed -as=${formatFile} -type=bed5+13 ${out} ${chromSizes} ${out}.bb

