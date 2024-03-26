##########
# Script for transforming simulation outs to PDB files
# First fetches coordinates from ONE h5 trajectory file, and then converts to .pdb format
# Directory in `data` must contain one `.h5` file with `n` blocks 
##########
import sys
sys.path.append('/mnt/data0/noah/software/polychrom/polychrom')
import polymerutils

data = "/mnt/data0/noah/analysis/misc-analysis/ORCA/polymer_sims/3D_simulation/sim_outs/"
n = 100

# Fetch blocks
for i in range(n):
    t = polymerutils.fetch_block(data, i)
    # Write as pdb
    polymerutils.save(t, filename="blocks_0-99_conf{}.pdb".format(i), mode="pdb")