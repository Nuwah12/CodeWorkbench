#!/bin/bash

# Define what counts as an analysis job
KEYWORDS="python|Rscript|samtools|bwa|bowtie|bedtools|snakemake|nextflow|perl|java|awk|gzip|sed|dorado|porechop"

# Optional: only show your user’s jobs
USER_FILTER=$(whoami)

# Header
printf "%-8s %-10s %-8s %-10s %s\n" "PID" "USER" "CPU%" "MEM%"  "STARTED" "COMMAND"

# Use ps to list all user processes, filter by keyword
ps -eo pid,user,pcpu,pmem,lstart,args --sort=-pcpu | \
  grep -E "$KEYWORDS" | grep "$USER_FILTER" | grep -v grep | \
  awk '{pid=$1; user=$2; cpu=$3; start=$4" "$5" "$6" "$7; cmd=""; for(i=8;i<=NF;++i)cmd=cmd $i" "; printf "%-8s %-10s %-8s %-10s %s\n", pid, user, cpu, start, cmd}'

