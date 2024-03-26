###############
# Using polychrom polymerutils module to transform h5 --> text x,y,z
# Noah Burget
# 3/26/24
###############
import os
import sys
import polychrom.polymerutils as polu

if len(sys.argv) != 2:
    print('Usage: python3 trajectory_to_txt.py <path to h5>')
    os._exit(1)
n = 1000
for i in range(n):
    t = polu.fetch_block(sys.argv[1], i)
    # Write as pdb
    polu.save(t, filename="blocks_0-99_conf{}.txt".format(i), mode="txt")
