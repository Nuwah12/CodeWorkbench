#!/bin/bash

JUICER="/mnt/data0/noah/software/juicer/CPU/juicer.sh"

Top_Dir="/mnt/data0/noah/processing/032524_hic/large_Pre_B_MboI/"
Restriction_Enzyme="MboI"
ChromSizes="/mnt/data0/noah/references/Mus_musculus/GRCm38/Annotation/ChromSizes_ordered.tsv"
BWA_Genome="/mnt/data0/noah/references/Mus_musculus/GRCm38/Sequence/BWAIndex/genome.fa"
RestrictionSites="/mnt/data0/noah/references/Mus_musculus/GRCm38/GRCm38_MboI.txt"
Juicer_Dir="/mnt/data0/noah/software/juicer"

CoolFile_Prefix="myCoolFile"

source /mnt/data0/apps/anaconda/anaconda2/bin/activate noah
echo "Running command bash ${JUICER} -d ${Top_Dir} -s ${Restriction_Enzyme} -p ${CromSizes} -z ${BWA_Genome} -y ${RestrictionSites} -D ${Juicer_Dir} -t 12"
bash ${JUICER} -d ${Top_Dir} -s ${Restriction_Enzyme} -p ${ChromSizes} -z ${BWA_Genome} -y ${RestrictionSites} -D ${Juicer_Dir} -t 12

hic2cool convert -r 0 "${Top_Dir}/aligned/inter_30.hic" "${CoolFile_Prefix}.mcool"

