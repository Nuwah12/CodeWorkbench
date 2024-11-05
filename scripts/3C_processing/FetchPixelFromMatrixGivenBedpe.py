import cooler
import cooltools
import pandas as pd
import numpy as np
import argparse
import re

### Get regions from cl
parser = argparse.ArgumentParser()

parser.add_argument("--cool", type=str, nargs='+', required=True, help="Path to cool file(s) to extract values from.")
parser.add_argument("--res", type=int, required=True, help="Resolution (bin size) of your regions and the cool file.")
parser.add_argument("--bedpe", type=str, required=True, help="Bedpe file of regions to extract (must have column names : chrom1, start1, end1, chrom2, start2, end2)")
parser.add_argument("--int_chromnames", action='store_true', required=False, help="When set, use integer chromosome names (1,2,3,...) and not chr1,chr2,etc.")
args = parser.parse_args()
# Make output filename based on input bedpe name
name = args.bedpe.split("/")
name = name[len(name)-1]

### Parse input and create cooler object
bedpe = pd.read_csv(args.bedpe, sep = "\t", header=0)

# Read in coolers and make names for columns 
cools = []
cooler_names = []
for i in args.cool:
    fname = i.split("/")
    cooler_names.append(fname[len(fname)-1].split('.')[0])
    if ".mcool" in i:
        cools.append(cooler.Cooler("{}::resolutions/{}".format(i, args.res)))
    elif ".cool" in i:
        cools.append(cooler.Cooler(i))

# Make column names for the final output
columns = ['chrom1','start1','end1','chrom2','start2','end2']
for i in range(len(cooler_names)):
    columns.append('{}_count'.format(i))
    columns.append('{}_balanced'.format(i))

print("Fetching regions...")
res_df = pd.DataFrame(columns = columns)
vals = []

def get_sum_of_contacts(x, cools, cooler_names, cols, intnames):
    (chrom1) = x.loc["chrom1"]
    (start1) = int(x.loc["start1"])
    (end1) = int(x.loc["end1"])
    (chrom2) = x.loc["chrom2"]
    (start2) = int(x.loc["start2"])
    (end2) = int(x.loc["end2"])
    if intnames:
        chrom1 = re.findall(r'\d+', chrom1)[0]
        chrom2 = re.findall(r'\d+', chrom2)[0]
    # Loop through coolers
    for i in range(len(cools)):
        pixel = cools[i].matrix(balance='VC', as_pixels=True, join=True).fetch("{}:{}-{}".format(chrom1, start1, end1), "{}:{}-{}".format(chrom2, start2, end2))
        balanced = sum(pixel['balanced'])
        count = sum(pixel['count'])
        if pixel.shape[0] == 0:
            x["{}_count".format(cooler_names[i])] = 0
            x["{}_balanced".format(cooler_names[i])] = 0
        else:
            x["{}_count".format(cooler_names[i])] = count
            x["{}_balanced".format(cooler_names[i])] = balanced
    #print(regions)
    return(x)

res_df = bedpe.apply(get_sum_of_contacts, axis=1, cools=cools, cooler_names=cooler_names, cols=columns, intnames=args.int_chromnames)
res_df.to_csv('{}_out.tsv'.format(name),sep='\t')







