#!/bin/bash
if [[ -z $1 ]]; then
	echo "Pass in path to config file!"
	exit 1
fi
source $1

#conda
if [[ ! -z ${conda} ]]; then
	source ${conda}
fi

mkdir -p dotfinder
mkdir -p mustache
mkdir -p chromosight

echo -e "Calling dots with all 3 methods, using matrix ${coolfile}"

#Run cooltools dotfinder
echo "$(date): Running cooltools dotfinder..."
nohup python3 ${dotfinder} ${coolfile} ${genome} "${dots_out}" --res ${resolution} --expected_val_col ${expected_value_col} --clr_weight_name ${clr_weight_name} --max_sep ${max_distance} --cluster_rad ${cluster_radius} --max_nans ${max_nans} --lambda_bin ${lambda_bin} --tile_size ${tile_size} --fdr ${fdr_threshold} --nproc ${threads} &> "./dotfinder.log" &
cool_pid=$!

#Run chromosight
echo "$(date): Running chromosight..."
nohup chromosight detect --threads ${threads} --min-dist ${min_distance} --max-dist ${max_distance} ${coolfile}::/resolutions/${resolution} "${chromosight_out}" &> "./chromosight.log" &
chromo_pid=$!

#Run mustache
echo "$(date): Running mustache..."
nohup ${mustache} -f ${coolfile} -r ${resolution} -o "${mustache_out}" -p ${threads} -pt ${pval_threshold} -st ${sparsity_threshold} &> "./mustache.log" &
mustache_pid=$!

#Wait for processes to finish before continuing
echo -e "$(date): Waiting for all methods to finish . . .\nCooltools dotfinder PID = ${cool_pid}\nChromosight PID = ${chromo_pid}\nMustache PID = ${mustache_pid}"
wait -f ${cool_pid} ${chromo_pid} ${mustache_pid}

#Generate overlap stats and consensus loop set
echo "$(date): Comparing loopsets, generating consensus set. The consensus loopset with be the ${overlap} of datasets ${datasets}"
Rscript --vanilla bin/compare_called.R "${chromosight_out}.tsv" "${mustache_out}" "${dots_out}" ${overlap} ${datasets}

#Extract raw and balanced scores for each called loop
echo "$(date): Extracting scores of consensus loops from cooler matrix..."
python3 bin/extract_scores.py --dots "./tmp_consensus.txt" --cool ${coolfile} --res ${resolution}

#Make final df
echo "$(date): Making final data structure..."
Rscript --vanilla bin/make_final_output.R "./tmp_consensus.txt" "./tmp_dotscores.txt" ${final_out} ${filter_na}

#Cluster consensus loopset
echo "$(date): Clustering aggregated loops, using same clustering radius as dotfinder..."
python3 bin/dotfinder_clustering.py --dots ${final_out} --radius ${cluster_radius}

rm "./tmp_consensus.txt" "./tmp_dotscores.txt"

echo "$(date): DONE!"

