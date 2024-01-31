#!/bin/bash

###############
# Script to convert allValidPairs (HiC-Pro) output to HiCSummary format (HOMER makeTagDirectory input)
# Noah Burget
# 1/21/24
###############

if [[ $# -lt 3 ]]; then
	echo "Usage: ./allValidPairs_to_HiCsummary.sh <allValidPairs> <outDir> <outPrefix>"
	exit 1
fi

allValidPairs=$1
outDir=$2
outPrefix=$3

cut -f1-7  ${allValidPairs} > ${outDir}/${outPrefix}.HiCsummary

echo "$(date): done!"
