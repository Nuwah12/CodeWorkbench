#!/bin/bash                                                                                                                                                                      
# 2016-12-13
# Yeqiao
#### simplify and integrate SE scripts to pipeline
#### take narrow peaks; assume island filtered BED/BAM is done at peak calling
#### quantify RPM
#### sort and write to bed6
#### call SE using Babak's R code
 
#### to run the script, suppliy 1) path to peaks 2) peak file 3) island filtered bam/bed 4) SHIFT 5) output directory

#### SET UP for Plutus #---------------------------------------------------------------------------------------------//
#module() { eval `/usr/bin/modulecmd bash $*`; }    
#module use /mnt/data0/modulefiles
#module load bwa-0.7.13
#module load bamtools-2.5.1
#module load samtools-1.9
#module load bedtools-2.27.1
#module load macs-2.0.9
#module load ucsc369
#module load homer-4.9

GENOME="/mnt/data0/fdb/igenome/Homo_sapiens/Ensembl/GRCh37.75/Annotation/Genes/ChromInfo_chr.txt"
config_PICARD_BASE="/mnt/data0/apps/picard"
BLACK="/mnt/data0/fdb/igenome/Homo_sapiens/Ensembl/GRCh37.75/Annotation/blackList/newBlackList/hg19_blackList_ATAC_170206.bed"
PIPELINE="/mnt/data0/yeqiao/analysis/pipeline/SuperEnhancer" 
 
#### Check and pass arguments #-------------------------------------------------------------------------------------//
 
if [ $# -eq 5 ]; then
	EXPDIR=${1}
	PEAK=${EXPDIR}/${2}
	IFFILE=${EXPDIR}/${3}
	SHIFT=${4} 
	OUTDIR=${5}

	#### output prefix
	PREFIX=${EXPDIR##/*/}
else 
	echo "please specify 1)path to peaks, 2)narrow peak file 3)island filtered BAM/BED, 4) SHIFT 5)output directory"
fi

#### FUNCTIONS #----------------------------------------------------------------------------------------------------//

#### RPM #---------------------------------
#### use coverage -counts to get raw read counts of slopped combined islands in each island filtered bed and normalize on total tag count only  
echo " "
get_normalized_count_rpm(){
        IFBED="$1"        #reads filtered on islands
        ISLD="$2"         #combined island file
        SLOP="$3"         #extend each read in island filtered bed
        TEMP="$4"         #output file
        SLOPISLAND="temp.slopped.islands.bed" #slopped island

	#### slop the islands by shift size and make sure there is no black list
	bedtools sort -i ${ISLD} | bedtools slop -i stdin -b ${SLOP} -g ${GENOME} | bedtools intersect -a stdin -b ${BLACK} -v | bedtools sort -i stdin > ${SLOPISLAND}        

	#### detect island filtered file: bed or bam
	if [[ "$IFBED" =~ ".bam" ]]; then
		TAGCOUNT=$(samtools view -c ${IFBED}  | awk '{print $1}')
		samtools sort -m 15000000000 ${IFBED} | bedtools coverage -a ${SLOPISLAND} -b stdin -g ${GENOME} -counts | awk -v S="$TAGCOUNT" -F "\t" 'BEGIN{OFS=FS};{print $0 "\t" $4*1000000/S}'  > ${TEMP}
	else 
		if [[ "$IFBED" =~ ".bed" ]]; then
			TAGCOUNT=$(wc -l ${IFBED}  | awk '{print $1}')
			sort -k 1,1 -k2,2n ${IFBED} | bedtools coverage -a ${SLOPISLAND} -b stdin -g ${GENOME} -counts | awk -v S="$TAGCOUNT" -F "\t" 'BEGIN{OFS=FS};{print $0 "\t" $4*1000000/S}'  > ${TEMP}
		
		else 
			echo "wrong island filtered file type!"
			exit
		fi
	fi
        echo "processed ${IFBED} ..."
        rm ${SLOPISLAND} 
}

#### calling SE ----------------------------------------------------------------------------------------------------//

#### make output dirctory
mkdir -p ${OUTDIR}/${PREFIX}
cd ${OUTDIR}/${PREFIX}

#### quantify RPM
get_normalized_count_rpm ${IFFILE} ${PEAK} ${SHIFT} ${PREFIX}_RPM.txt

##### sort RPM files -- rank -- write bed6
sort -nrk 5 ${PREFIX}_RPM.txt | awk -F"\t" 'BEGIN {OFS=FS} {print $1"\t"$2"\t"$3"\t""RPM_"FNR"\t"$5"\t""+"}' > ${PREFIX}_RPM_rank.bed 

#### call SE
sh ${PIPELINE}/Quantitate_MapSEs_YoungMethod.sh ${OUTDIR}/${PREFIX} ${PREFIX}_RPM_rank.bed
