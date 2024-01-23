#!/bin/bash

###############
# Script to convert FitHiChIP (significant) interactions to a .cool matrix
# 1/23/24
# Noah Burget
###############

if [[ $# -lt 4 ]]; then
	echo "Usage: ./FitHiChIP_to_cool.sh <chromSizes> <resolution> <FitHiChIP interactions file> <output cool file>"
	exit 1
fi

chromSizes=$1
res=$2
fithichip_out=$3
out_cool=$4

cut -f1-7 ${fithichip_out} > ${fithichip_out}.tmp.bg2

awk 'BEGIN {OFS="\t"} NR == 1 {print ""} NR != 1' ${fithichip_out}.tmp.bg2 > ${fithichip_out}.bg2
rm ${fithichip_out}.tmp.bg2

cooler load -f bg2 ${chromSizes}:${res} ${fithichip_out}.bg2 ${out_cool}

echo "$(date): done!"
