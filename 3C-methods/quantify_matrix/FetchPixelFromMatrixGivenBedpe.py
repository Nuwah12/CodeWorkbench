import cooler
import cooltools
import pandas as pd
import numpy as np
import argparse
import re
from functools import reduce

### Get regions from cl
parser = argparse.ArgumentParser()

parser.add_argument("--cool", type=str, nargs='+', required=True, help="Path to cool file(s) to extract values from.")
parser.add_argument("--res", type=int, required=True, help="Resolution (bin size) of your regions and the cool file.")
parser.add_argument("--bedpe", type=str, required=True, help="Bedpe file of regions to extract (must have column names : chrom1, start1, end1, chrom2, start2, end2)")
parser.add_argument("--norm", type=str, required=False, default='VC', help="Name of weight column to apply to contact counts. (Default = VC). VC and VC_SQRT weights are treated as divisive.")
parser.add_argument("--int_chromnames", action='store_true', required=False, help="When set, use integer chromosome names (1,2,3,...) and not chr1,chr2,etc.")
parser.add_argument("--sum", action='store_true', required=False, help="Sum intensities over all pixels in the current region, else return all bins for each region")
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

def get_bins(x,i,clr):
    if i==0:
        return clr.bins()._fetch('{}:{}-{}'.format(x['chrom1'],x['start1'],x['end1']))[0]
    else:
        return clr.bins()._fetch('{}:{}-{}'.format(x['chrom2'],x['start2'],x['end2']))[0]

def get_sum_of_contacts(x, cools, cooler_names, cols, intnames):
    # Define coordinates
    (chrom1) = x.loc["chrom1"]
    (start1) = int(x.loc["start1"])
    (end1) = int(x.loc["end1"])
    (chrom2) = x.loc["chrom2"]
    (start2) = int(x.loc["start2"])
    (end2) = int(x.loc["end2"])
    #print(x)
    # If intnames option selected, make the chromosome names 1,2,3... and not chr1,chr2.chr3... (assuming the contacts are intra-chromosomal)
    if intnames:
        if 'X' in chrom1:
            chrom1 = 'X'
            chrom2 = 'X'
        elif 'Y' in chrom1:
            chrom1 = 'Y'
            chrom2 = 'Y'
        elif 'MT' in chrom1:
            chrom1 = 'MT'
            chrom2 = 'MT'
        else:
            chrom1 = re.findall(r'\d+', chrom1)[0]
            chrom2 = re.findall(r'\d+', chrom2)[0]
    # Loop through coolers and extract raw and balanced counts
    c = 0
    allPixels = pd.DataFrame()
    for i in range(len(cools)):
        # Fetch pixels for current coordinate region via cooler's matrix selector
        pixel = cools[i].matrix(balance=args.norm, as_pixels=True, join=True).fetch("{}:{}-{}".format(chrom1, start1, end1), "{}:{}-{}".format(chrom2, start2, end2))
        if pixel.shape[0] == 0: # if no pixels returned, just set counts to 0
            x["{}_count".format(cooler_names[i])] = 0
            x["{}_balanced".format(cooler_names[i])] = 0
        else:
            count = sum(pixel['count'])
            balanced = sum(pixel['balanced'])
            x["{}_count".format(cooler_names[i])] = count
            x["{}_balanced".format(cooler_names[i])] = balanced
    return(x)

def quant_by_bin(x, cools, cooler_names):
    f = pd.DataFrame(columns=['chrom1','start1','end1','chrom2','start2','end2'])
    for i, v in enumerate(x.values):
        (chrom1) = v[0]
        (start1) = int(v[1])
        (end1) = int(v[2])
        (chrom2) = v[3]
        (start2) = int(v[4])
        (end2) = int(v[5])
        if args.int_chromnames:
            if 'X' in chrom1:
                chrom1 = 'X'
                chrom2 = 'X'
            elif 'Y' in chrom1:
                chrom1 = 'Y'
                chrom2 = 'Y'
            elif 'MT' in chrom1:
                chrom1 = 'MT'
                chrom2 = 'MT'
            else:
                chrom1 = re.findall(r'\d+', chrom1)[0]
                chrom2 = re.findall(r'\d+', chrom2)[0]
        dfs = []

        for i in range(len(cools)):    
            bins = cools[i].matrix(balance=args.norm, as_pixels=True, join=True).fetch("{}:{}-{}".format(chrom1, start1, end1), "{}:{}-{}".format(chrom2, start2, end2))
            bins.rename({'balanced':'{}_balanced'.format(cooler_names[i]), 'count':'{}_count'.format(cooler_names[i])}, inplace=True, axis=1)
            bins['{}_bin1'.format(cooler_names[i])] = bins.apply(get_bins,i=0,clr=cools[i],axis=1)
            bins['{}_bin2'.format(cooler_names[i])] = bins.apply(get_bins,i=1,clr=cools[i],axis=1)
            
            dfs.append(bins)
        # Merge cooler dfs
        bins = reduce(lambda x, y: x.merge(y, on='Date'), dfs)
        f = pd.concat([bins,f])
    return f


if args.sum:
    res_df = bedpe.apply(get_sum_of_contacts, axis=1, cools=cools, cooler_names=cooler_names, cols=columns, intnames=args.int_chromnames)
    res_df.to_csv('{}_out.tsv'.format(name),sep='\t')
else:
   res_df =  quant_by_bin(bedpe, cools, cooler_names).reset_index(drop=True).sort_values(by='chrom1')
   res_df.to_csv('{}_out.tsv'.format(name),sep='\t',index=False)
