#!/home/noah/.conda/envs/noah/bin/python3
#########################
# Script using polychrom's contact map generation implementation
# Noah Burgt
# 3/26/24
#########################
import sys
import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

def main(p):
    dist_mats = []

    filenames = os.listdir(p)
    print("Starting distance calculations..")
    for i,x in enumerate(filenames):
        f = "{}/{}".format(p,x)
        f_df = np.loadtxt(f, delimiter=' ', skiprows=1)
        dists = pdist(np.array(f_df), metric="euclidean")
        dist_mat = squareform(dists)
        dist_mats.append(dist_mat)
    
    print("Averaging distances..")
    avg_dist = np.mean(dist_mats, axis=0)
    print("Aggregating distance matrix...")
    agg_dist = avg_dist.reshape(30,30,30,30).mean(axis=(2,3))
    print("Writing to file..")
    np.savetxt('distMatrix.txt', agg_dist, delimiter='\t')

    return agg_dist

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print('Usage: ./make_distanceMap <dir of confs>')
        sys.exit()
    main(sys.argv[1])
