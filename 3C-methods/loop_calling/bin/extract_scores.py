import argparse
import cooler
import cooltools
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('--dots', type=str, help="Path to temp file with string-encoded dot coordinates")
parser.add_argument('--cool', type=str, help="Cool file containing interaction matrix")
parser.add_argument('--res', type=int, help="Resolution of matrix")

args = parser.parse_args()

#print("Loading matrix and called loops...")
if ".mcool" in args.cool:
    clr = cooler.Cooler("{}::resolutions/{}".format(args.cool, args.res))
elif ".cool" in args.cool:
    clr = cooler.Cooler(args.cool)

raw_counts = []
balanced_counts = []

dots = pd.read_csv(args.dots, header=None)

#print("Extracting scores...")
for i in np.arange(0, dots.shape[0]):
    j = dots.iloc[i,0].split("_")
    b = clr.matrix(balance = True, as_pixels = True, join = True).fetch("{}:{}-{}".format(j[0], j[1], j[2]), "{}:{}-{}".format(j[3], j[4], j[5]))

    if len(b) == 0:
        #print("Could not extract region {}".format(j))
        #print("Processing loop {}".format(i))
        raw_counts.append(-1)
        balanced_counts.append(-1)

    else:
        #print("Processing loop {}".format(i))
        raw_counts.append(b.iloc[0,6])
        balanced_counts.append(b.iloc[0,7])
        #print("Loop {}, raw value {}, balanced value {}".format(j, b.iloc[0,6], b.iloc[0,7]))

print("Constructing output...")
counts_df = pd.DataFrame({"Raw_counts":raw_counts, "Balanced_counts":balanced_counts})
counts_df.to_csv("tmp_dotscores.txt", index = False, header = False)


