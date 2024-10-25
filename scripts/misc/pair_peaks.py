import argparse
import numpy as np
import pandas as pd
import bioframe
import itertools

###################################
# Script for pairing all entries of a .bed file within a distance threshold
# 6/10/24
# Noah Burget
##################################

def _parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", type=str, required=True, help="Path to bed file")
    parser.add_argument("--lower", type=int, required=False, help="Lower distance threshold (default = 20,000)", default=20000)
    parser.add_argument("--upper", type=int, required=False, help="Upper distance threshold (default = 2,000,000)", default=2000000)
    parser.add_argument("--out", type=str, required=False, help="Output file name (default = out_paired.bedpe)", default="out_paired.bedpe")
    return parser.parse_args()

def main():
    CHROMS = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
    args = _parseArgs()
    bed = pd.read_csv(args.bed, sep="\t", names=['chr','start','end'])
    
    groupByChrom = bed.groupby('chr')
    master = []
    for c in CHROMS:
        try:
            group = groupByChrom.get_group(c)
            print(c)
            groupStr = group.astype(str).apply(lambda x: '_'.join(x), axis=1).tolist()
            keep = []
            for i in itertools.combinations(groupStr, 2):
                #print(i)
                dist = np.abs(int(i[0].split('_')[1]) - int(i[1].split('_')[2]))
                if dist < args.upper and dist > args.lower:
                    keep.append(i)
            master.append(keep)
        except:
            continue
    # Make final output
    master_flat = [i for chromList in master for i in chromList]
    df = pd.DataFrame(master_flat)
    df = df.astype(str)
    a1 = df[0].str.split('_', expand=True)
    a2 = df[1].str.split('_', expand=True)
    pd.concat([a1, a2], axis=1).to_csv(args.out, sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()
