#!/bin/bash

USAGE="Usage: $0 --mapping/-m <name mapping file with form [oldname, newname]> --directory/-d <directory with files to rename>"

function print_usage_exit()
{
	if [[ ! -z $1 ]]; then
		echo -e "Unknown parameter '$1' \n${USAGE}"
		exit 1
	fi
	echo ${USAGE}
	exit 1
}

if [[ $# -lt 2 ]]; then
	echo $USAGE
	exit 1
fi

while [[ $# -gt 0 ]]; do
	case "$1" in
		-m | --mapping )
			map="$2"
			shift 2 # Shift past both argument and value
			;;
		-d | --directory )
			dir="$2"
			shift 2
			;;
		-* | --* )
			print_usage_exit $1
			;;
	esac
done



while read old new; do
    for file in "$dir/${old}"*; do
        [[ -e "$file" ]] || continue
        base="${file##*/}"                    # Strip directory
        suffix="${base#$old}"                 # Remove old name from filename
        mv "$file" "$dir/${new}${suffix}"     # Rename in same directory
    done
done < $map
