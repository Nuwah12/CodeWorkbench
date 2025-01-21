##############################
# Script to rename a directory of SRA files whose names are just accession numbers (i.e. SRR123456) given a table of metadata
# Noah Burget
# 1/21/25
##############################

import pandas as pd
import os
import shutil
import argparse

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--sra", type=str, required=True, help="Path to .fastq or .bam files with SRA accession numbers in names.")
    parser.add_argument("--metadata", type=str, required=True, help="Path to and name of tab-delimited file containing metadata for each SRA accession.")
    parser.add_argument("--filetype", type=str, required=False, default='fastq', help="Filetype to be renamed. Use fastq for .fastq files [default], bam for .bam files, and bw for bigWig files")

    args = parser.parse_args()

    if not os.path.exists(args.sra):
        raise Exception(f"SRA .fastq directory at {args.sra} not found!")
    if not os.path.exists(args.metadata):
        raise Exception(f"Metadata table at {args.metadata} not found!")

    return args

def main():
    args = parse_args()

    if args.filetype == "fastq":
        extension = ".fastq.gz"
    elif args.filetype == "bam":
        extension = ".bam"
    elif args.filetype == "bw":
        extension = ".bw"

    # How many .fastq files are in specified SRA directory?
    c = 0
    for i in os.listdir(args.sra):
        if not i.endswith(extension):
            continue
        elif "_1" in i or "_2" in i and extension == ".fastq.gz":
            c+=0.5
        else:
            c+=1

    # open metadata table
    meta = pd.read_csv(args.metadata, sep=',')

    num_samples = len(meta)

    # If we have mismatching number of samples in data dir and metadata table, thats a problem
    if c != num_samples:
        raise Exception(f"Found {num_samples} rows in metadata table, but {c} files in specified directory {args.sra}.")

    # Check that each of the SRA accessions in the metadata table have a .fastq.gz file in the specified directory
    for i in os.listdir(args.sra):
        if "_1" not in i and "_2" not in i:
            acc = i.split(".")[0]
        else:
            acc = i.split("_")[0]
        if acc not in meta.Run.to_list():
            raise Exception(f"Found file {i} in SRA directory {args.sra}, but it is not present in the metadata table.")
        else:
            continue
    
    fields = meta.columns.to_list()
    
    # Rename the files
    os.mkdir(f"{args.sra}/renamed")
    for fq in os.listdir(args.sra):
        if not fq.endswith(extension):
            continue
        rename = ""
        read = None
        if "_1" not in fq and "_2" not in fq:
            acc = fq.split(".")[0]
        elif "_1" in fq:
            acc = fq.split("_")[0]
            read = "1"
        elif "_2" in fq:
            acc = fq.split("_")[0]
            read = "2"
        data_row = meta[meta.Run==acc]
    
        for i, f in enumerate(fields):
            value = data_row[f].iloc[0]
            # Replace all spaces with underscores
            value = value.replace(" ","_")
            if i == len(fields)-1:
                if "_1" in f:
                    rename += f"{value}_R1"
                elif "_2" in f:
                    rename += f"{value}_R2"
                else:
                    rename += f"{value}"
            else:
                rename += f"{value}_"

        if read is None:
            rename = rename + extension
        else:
            rename = rename + f"_{read}" + extension
        
        print(f"{fq} will be renamed to {rename}")

        shutil.copy(f"{args.sra}/{fq}", f"{args.sra}/renamed/{rename}")

if __name__ == "__main__":
    main()





