#!/bin/bash
PIPELINE="/mnt/data0/yeqiao/analysis/pipeline/SuperEnhancer"
BLACK="/mnt/data0/fdb/igenome/Homo_sapiens/Ensembl/GRCh37.75/Annotation/blackList/newBlackList/hg19_blackList_170206.bed"
#windowsize=50
mergesize=12500

#### pass directory 
MYDIR=""
if [ "$#" -eq 2 ]
then
    MYDIR=$1
    INFILE=$2
else
    echo "Please specify working directory and island bed file"
    exit 1
fi


mydir=$MYDIR
inputFile=$INFILE
cd ${mydir}


processed_DIR=${mydir}
Tempdir=${mydir}/TempFiles
SEdir=${mydir}/SEs
TEdir=${mydir}/TEs

mkdir ${Tempdir}
mkdir ${SEdir}
mkdir ${TEdir}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1) ------------------- merge peaks -------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++

cd $processed_DIR
#loop for each library you need to find SE.
#for i in `ls *.bedgraph`;
#do
# I used .bed in the loop to find the filename properly
#  filename="${i%.*}"
  filename="${inputFile%%.*}"
  echo $filename
  #output=${filename}.text
  ##note: ${filename}_Th0WCE_peaks.bedgraph is the SICER/MACs output, you can change the name to reflect that
  #awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $4}' ${processed_DIR}/${filename}.bedgraph > ${mydir}/${filename}_5col.bedgraph
  #mergeBed -scores sum -i ${mydir}/${filename}_5col.bedgraph -d $mergesize > ${mydir}/${filename}_merged.bedgraph
  #rm ${mydir}/${filename}_5col.bedgraph
  #mv ${mydir}/${filename}_merged.bedgraph $Tempdir
  awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' ${processed_DIR}/${inputFile} | bedtools sort -i stdin | bedtools merge -c 5 -o sum -d ${mergesize} -i stdin > ${mydir}/${filename}_merged.bedgraph
  mv ${mydir}/${filename}_merged.bedgraph ${Tempdir}
#done

# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # 2) ------------------- find SEs in R -------------------
# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++
cd ${mydir}

Rscript --vanilla  ${PIPELINE}/MapSEs_YoungMethod.r $mydir

# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # 3) ---------- do some processing ---------------------
# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++
cd ${SEdir}


# to remove the sex chromosomes
#sort -rgk 4 ${filename}_SE.bedgraph | grep -v -E "^chrX|^chrY" | awk '{print $1 "\t" $2 "\t" $3 "\t" "SE_"FNR "\t" $4}' | bedtools sort -i stdin > ${filename}_SE_sort_5col.bed


sort -rgk 4 ${filename}_SE.bedgraph | awk '{print $1 "\t" $2 "\t" $3 "\t" "SE_"FNR "\t" $4}' | bedtools sort -i stdin | bedtools intersect -a stdin -b $BLACK -v > ${filename}_SE_sort_5col.bed

cut -f 1-3 ${filename}_SE_sort_5col.bed > ${filename}_SE_sort_3col.bed

# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # 4) ---------- do some processing ---------------------
# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++
cd ${TEdir}

# to remove the sex chromosomes
# cat ${filename}_TE.bed | grep -v -E "^chrX|^chrY" | bedtools sort -i stdin > ${filename}_TE_sort_3col.bed

cat ${filename}_TE.bed | bedtools sort -i stdin | bedtools intersect -a stdin -b $BLACK -v > ${filename}_TE_sort_3col.bed

