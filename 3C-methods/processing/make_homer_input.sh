########## Script for converting Juicer merged_nodups.txt to homer's HiCsummary format
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <Juicer merged_nodups.txt> <output tag directory> ..."
    exit 1
fi

echo "Constructing HiCsummary format..."
# Extract only relevant columns from the merged_nodups file
awk -v OFS='\t' '{ print " ", $2, $3, $1, $6, $7, $5}' merged_nodups > $(basename $1)_noSeq.tmp

# Replace 0s/1s in columns 4 and 7 with +/-, respectively
awk -v OFS='\t' '{gsub(/0/, "+", $3); gsub(/1/, "-", $3); gsub(/0/, "+", $6); gsub(/1/, "-", $6); print}' $(basename $1)_noSeq.tmp > $(basename $1)_hicsummary.tmp

# Add an empty column the file
awk -v OFS='\t' '{print "", $0}' $(basename $1)_hicSummary.txt > $(basename $1)_noSeq.tmp > $(basename $1)_hicsummary.txt

# Make a Homer tag directory
makeTagDirectory ${2} -format HiCsummary $(basename $1)_hicsummary.txt 
