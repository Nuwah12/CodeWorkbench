#!/home/noah/.conda/envs/noah/bin/python3
#########################
# Script using polychrom's contact map generation implementation
# Noah Burgt
# 3/26/24
#########################
import sys
import os
import numpy as np
from polychrom import contactmaps as cm
if len(sys.argv) != 2:
    print('Usage: python3 make_contactMaps <dir of confs>')
    os._exit(1)

filenames = os.listdir(sys.argv[1])
for i,x in enumerate(filenames):
    filenames[i] = "{}/{}".format(sys.argv[1],x)

out = cm.monomerResolutionContactMap(filenames, cutoff=10)
np.savetxt('matrix.txt', out, delimiter='\t')
