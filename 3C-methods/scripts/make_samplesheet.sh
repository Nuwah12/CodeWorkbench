#!/bin/bash
if [[ -f "samplesheet.tsv" ]]; then
	rm samplesheet.tsv
fi
MODE="single"
while true; do
	case "$1" in
		-p | --paired ) MODE="paired"; shift ;;
		-d | --directory ) dir="$2"; shift ;;
		* ) break ;;
	esac
done

if [[ ! -f "samplesheet.tsv" ]]; then
	touch samplesheet.tsv
fi
#echo "samples" >> samplesheet.csv
if [[ $MODE == "single" ]]; then
	echo "Making single-end samplesheet"
	for i in $(ls ${dir}); do
		i="${i##*/}"
		i="${i%%.*}"
		echo $i >> samplesheet.tsv
	done
elif [[ $MODE == "paired" ]]; then
	echo "Making paired-end samplesheet"
	for i in $(ls ${dir}); do
		i="${i##*/}"
		echo $i
		if [[ $i =~ "_R2" ]]; then
			continue
		fi
		i=${i%%_R1*}
		echo $i >> samplesheet.tsv
	done
fi
