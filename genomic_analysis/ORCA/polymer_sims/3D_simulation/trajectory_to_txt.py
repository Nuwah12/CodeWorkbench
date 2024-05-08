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
if len(sys.argv) != 3:
    print('Usage: python3 trajectory_to_txt.py <total number of confs> <path to h5>')
    os._exit(1)

d = "confs_txt"
p = Path(d)
if not p.exists():
    p.mkdir()
else:
    print('Delete ./confs_txt and re-run.')

n=int(sys.argv[1])
for i in range(n):
    t = polu.fetch_block(sys.argv[2], i)
    # Write as text file
    polu.save(t, filename="{}/conf{}.txt".format(d, i), mode="txt")
