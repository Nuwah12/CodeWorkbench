#!/home/noah/.conda/envs/noah/bin/python3
###############
# Using polychrom polymerutils module to transform h5 --> text x,y,z
# Noah Burget
# 3/26/24
###############
import os
import sys
import polychrom.polymerutils as polu
from pathlib import Path

def main(d, n=10000):

    p = Path(f"{d}/confs_txt")
    if not p.exists():
        p.mkdir()
    else:
        print('Delete ./confs_txt and re-run.')
        os._exit(1)
    for i in range(n):
        t = polu.fetch_block(f"{d}/sim_outs", i)
        # Write as text file
        polu.save(t, filename="{}/conf{}.txt".format(p, i), mode="txt")
