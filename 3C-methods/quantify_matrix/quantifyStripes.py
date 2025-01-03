import argparse
import cooler
import cooltools
import pandas as pd
import numpy as np
from progress.bar import *

"""
Script for quantifying 'Stripe' regions defined for a contact matrix from a 3C experiment. Tested with 'Stripenn' tool.
In general, this will quantify some given area for a .cool file's matrix, a generalization of this is available in QuantifyMatrixGivenBedpe.py
"""

parser = argparse.ArgumentParser()
parser.add_argument('--stripes', type=str, help="Path to bed-style file with dot coordinates", required = True)
parser.add_argument('--cool', type=str, help="Cool file containing interaction matrix with base resolution of 1000bp", required = True)
parser.add_argument('--prefix', type=str, help="Prefix for output file: <prefix>_dotscores.tsv", required = True)
parser.add_argument('--norm', type=str, help="Type of normlization to use on raw counts (depth, balance)", choices = ['depth', 'balance'], required = False)
parser.add_argument('--matnorm', type=str, help="Type of normlization to use on raw counts that is stored in the matrix itself (in cooler/HDF5 file)", choices = ['VC', 'VC_SQRT', 'KR'], required = False, default = True)
#parser.add_argument('--skinny', type=bool, help="Use only center bin in width range as width of stripe", required = False)

args = parser.parse_args()

raw_counts = []
balanced_counts = []
peak_sizes = []

stripes = pd.read_csv(args.stripes, header=None)

if ".mcool" in args.cool:
    clr = cooler.Cooler("{}::resolutions/{}".format(args.cool, resolution))
elif ".cool" in args.cool:
    clr = cooler.Cooler(args.cool)

bar = ChargingBar('Getting stripe scores', max = stripes.shape[0])
for i in range(0, stripes.shape[0]):
    s = stripes.iloc[i].values[0].split("\t")
    if int(s[1]) < int(s[4]):
        s1 = [s[0], s[1], s[2]]
        s2 = [s[3], s[4], s[5]]
    else:
        s2 = [s[0], s[1], s[2]]
        s1 = [s[3], s[4], s[5]]
    #print("s1: {} s2: {}".format(s1,s2))
    b = clr.matrix(balance = args.matnorm, as_pixels = True, join = True).fetch("{}:{}-{}".format(s1[0], s1[1], s1[2]), "{}:{}-{}".format(s2[0], s2[1], s2[2]))
    #print(b) 
    if len(b) == 0:
        raw_counts.append(-1)
        balanced_counts.append(-1)

    else:
        #### Get sum of bins in stripe
        count_sum = b["count"].sum()
        raw_counts.append(count_sum)
        balanced = b["balanced"].sum()
        balanced_counts.append(balanced)
        #print(stripes)
    peak_sizes.append(int(s[2]) - int(s[1])) 
    bar.next()

counts_df = pd.DataFrame({"Raw_counts":raw_counts})

if args.norm == "depth":
    rpm = []
    rpkm = []
    count_sum = counts_df["Raw_counts"].sum()
    ### RPM
    for i in range(0, counts_df.shape[0]):
       rpm.append((counts_df.iloc[i,0] * 10**6)/count_sum)
    counts_df["rpm"] = rpm
    ### RPKM - 'length' parameter is length of SMC1 peak
    for i in range(0, counts_df.shape[0]):
        size = peak_sizes[i]
        #print("Peak size is {}".format(size))
        rpkm.append((counts_df.iloc[i,0] * 10**9)/(count_sum * size))
    counts_df["rpkm"] = rpkm

counts_df.to_csv("{}_stripeScores.tsv".format(args.prefix), index = False, header = False)
bar.finish()

#           _   
#       .__(.)< (MEOW)
#        \___)
