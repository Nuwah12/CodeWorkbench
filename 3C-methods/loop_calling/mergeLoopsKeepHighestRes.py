from pybedtools import BedTool
import pandas as pd
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(description="Merge loops from different resolutions. The representative dot is from the highest resolution.")
    parser.add_argument("dots", type=str, help="Space-separated list of paths to dot files", nargs="+")
    parser.add_argument("resolution", type=int, help="Highest Resolution")
    parser.add_argument("prefix", type=str, help="Output: <prefix>_mergedResolutions.dots")
    args = parser.parse_args()

    dfs = {}
    for i in args.dots:
        print("Reading dot file {}".format(i))
        f = pd.read_csv(i, sep = "\t")
        res = f.iloc[0,2] - f.iloc[0,1]
        dfs[res] = BedTool(i)
    
    # repeated inverse intersection
    for i in sorted(dfs):
        if i == args.resolution:
            continue
        for j in sorted(dfs):
            if int(j) > int(i) or i == j:
                continue
            #print("Intersecting {} and {} resolutions".format(i, j))
            #print(dfs[i])
            dfs[i] = dfs[i].pair_to_pair(dfs[j], type='neither')
        print("{} bin size had {} unique dots".format(i, dfs[i].to_dataframe().shape[0]))
    print(dfs)

    # concat bedtool objs
    dots = dfs[args.resolution]
    #print(dots)
    for i in dfs:
        dfs[i] = dfs[i].to_dataframe()

    dots = pd.concat(dfs)
    dots.to_csv("{}_mergedResolutions.dots".format(args.prefix), sep = "\t")

if __name__ == "__main__":
    main()

    



