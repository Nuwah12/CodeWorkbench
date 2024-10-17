#!/usr/bin/env python
# coding: utf-8

# ## 1-D loop extrusion simulation
# ### !!! Make sure you run the cells IN ORDER !!!

# In[6]:


from extruder import Extruder
import numpy as np
from pathlib import Path
import seaborn as sns
import h5py
import sys


# ### Translocating and writing trajectories using custom extruder class

# In[7]:


RUN_NAME = "MYC_SemiPerturb_highBlock_lowLEF" #### Define a name for this simulation run, sowe can save parameters
N1 = 900 # Size of 1 system
M = 1 # No. of systems
N = N1*M # Total size of system, in momomers
occupied = np.zeros(N) # List to tell if current monomer is occupied by an extruder
occupied[0] = 1 
occupied[-1] = 1
steps = 50000 # Timesteps for 1D sim.
LEFNum = 8 # No. of RANDOMLY LOADED extruders
num_chunks = 50 # No. of chunks to write trajectories in

### Blockers (i.e. CTCF) - {pos. : prob.}
left_blockers_capture = {}
right_blockers_capture = {}
left_blockers_release = {}
right_blockers_release = {}
## Here we define strong, mild, and weak [blocking, release]
# WT blocking strengths
#STRONG_BLOCK = [0.999, 0.0001]
#MEDIUM_BLOCK = [0.75, 0.001]
#WEAK_BLOCK = [0.25, 0.025]
# Perturbation blocking strengths
#STRONG_BLOCK_p = [0.75, 0.01]
#MEDIUM_BLOCK_p = [0.5, 0.02]
#WEAK_BLOCK_p = [0.2, 0.03]
TIER1_BLOCK = [0.99, 0.0001]
TIER2_BLOCK = [0.75, 0.01]
TIER3_BLOCK = [0.5, 0.015]
#TIER3_BLOCK_V2 = [0.4, 0.02]
TIER4_BLOCK = [0.25, 0.03]
# Manually assigning blockers in dict
# Boundary is ~half a probe away from E1.1, so we will put it at the END of Probe 6 (30*6=180)
left_blockers_capture[165] = TIER1_BLOCK[0] # Boundary -- CENTER OF PROBE 6 -- SHOULD NOT CHANGE **
left_blockers_release[165] = TIER1_BLOCK[1]
left_blockers_capture[180] = TIER1_BLOCK[0] # E1.1 -- BEGINNING OF OF PROBE 7 (BiD) -- SHOULD NOT CHANGE **
left_blockers_release[180] = TIER1_BLOCK[1]
left_blockers_capture[210] = TIER1_BLOCK[0] # E1.2 -- EBF1 -- BEGINNING OF PROBE 8 (BiD) -- DECREASES
left_blockers_release[210] = TIER1_BLOCK[1]
left_blockers_capture[315] = TIER3_BLOCK[0] # E2 -- EBF1 -- CENTER OF PROBE 11 (BiD) - DRAMATICALLY DECREASES
left_blockers_release[315] = TIER3_BLOCK[1]
left_blockers_capture[405] = TIER4_BLOCK[0] # B1 -- CENTER OF PROBE 14 -- SHOULD NOT CHANGE **
left_blockers_release[405] = TIER4_BLOCK[1]
left_blockers_capture[555] = TIER2_BLOCK[0] # B3 -- EBF1 -- CENTER OF PROBE 19 (BiD) -- DECREASES
left_blockers_release[555] = TIER2_BLOCK[1]

right_blockers_capture[210] = TIER1_BLOCK[0] # E1.2 -- EBF1 -- BEGINNING OF PROBE 8 (BiD) -- DECREASES
right_blockers_release[210] = TIER1_BLOCK[1]
right_blockers_capture[315] = TIER3_BLOCK[0] # E2 -- EBF1 -- CENTER OF PROBE 11 (BiD) -- DRAMATICALLY DECREASES
right_blockers_release[315] = TIER3_BLOCK[1]
right_blockers_capture[495] = TIER4_BLOCK[0] # B2 -- CENTER OF PROBE 17 -- SHOULD NOT CHANGE **
right_blockers_release[495] = TIER4_BLOCK[1]
right_blockers_capture[555] = TIER2_BLOCK[0] # B3 -- EBF1 -- CENTER OF PROBE 19 (BiD) -- DECREASES
right_blockers_release[555] = TIER2_BLOCK[1]
right_blockers_capture[735] = TIER1_BLOCK[0] # MYC promoter -- CENTER OF PROBE 25 -- SHOULD NOT CHANGE **
right_blockers_release[735] = TIER1_BLOCK[1]

#for i in range(M):
#    for locs in left_blocker_locs:
#        pos = i * N1 + locs
#        left_blockers_capture[pos] = 0.90
#        left_blockers_release[pos] = 0.01
#    for locs in right_blocker_locs:
#        pos = i * N1 + locs
#        right_blockers_capture[pos] = 0.90
#        right_blockers_release[pos] = 0.01


# In[8]:


EXTRUDERS = []
targeted = []
# Generate random spots
LOADING_SPOTS = []
for i in range(LEFNum):
    LOADING_SPOTS.append(np.random.randint(N))
    targeted.append(False)
# Targeted loading
TARGETED_LOADING_SPOTS = [] # 5/16/24 - No targeted loading!!!
LEFNum = LEFNum + len(TARGETED_LOADING_SPOTS)
for i in range(len(TARGETED_LOADING_SPOTS)):
    LOADING_SPOTS.append(TARGETED_LOADING_SPOTS[i])
    targeted.append(True)
LOADING_PROBS = [0.9,0.1]
print(targeted)
#if len(LOADING_SPOTS) != LEFNum:
#    print('ERROR - There must be as many loading spots as there are Extruders (LEFNum).')
#    sys.exit(1)
LIFETIME = 800
LIFETIME_STALLED = LIFETIME // 10
targeted_idx = 0
for i in range(len(LOADING_SPOTS)):
    #leg = np.random.randint(N)
    #print('Loading extruder at {}'.format(LOADING_SPOTS[i]))
    leg = LOADING_SPOTS[i]
    print('Loading extruder at {}'.format(leg))
    if not targeted[i]:
        EXTRUDERS.append(Extruder(
                extruder_index = i,
                leg1 = leg,
                leg2 = leg+1,
                left_blockers_capture = left_blockers_capture,
                right_blockers_capture = right_blockers_capture,
                left_blockers_release = left_blockers_release,
                right_blockers_release = right_blockers_release,
                extrusion_occupancy = occupied,
                lifetime = LIFETIME,
                lifetime_stalled = LIFETIME_STALLED,
                )
            )
    else:
            EXTRUDERS.append(Extruder(
                extruder_index = targeted_idx,
                leg1 = leg,
                leg2 = leg+1,
                left_blockers_capture = left_blockers_capture,
                right_blockers_capture = right_blockers_capture,
                left_blockers_release = left_blockers_release,
                right_blockers_release = right_blockers_release,
                extrusion_occupancy = occupied,
                targeted_loading = True,
                loading_spots = TARGETED_LOADING_SPOTS,
                loading_dist = LOADING_PROBS,
                lifetime = LIFETIME,
                lifetime_stalled = LIFETIME_STALLED,
                )
            )
            targeted_idx+=1


# In[9]:


### Write parameters to text file
with open('/mnt/data0/noah/analysis/misc-analysis-local/ORCA/polymer_sims/MYC/3D_simulation/{}_params.txt'.format(RUN_NAME),'w') as pf:
    pf.write("N: {}\n1D Steps: {}\nTotal LEF: {}\nTargeted: {} Probs: {}\nLifetime: {} Lifetime stalled: {}\nLeft cap: {}, Left rel: {}\nRight cap: {}, Right rel: {}".format(
                                                                                            N1,steps,LEFNum,
                                                                                            TARGETED_LOADING_SPOTS,LOADING_PROBS,
                                                                                            LIFETIME,LIFETIME_STALLED,left_blockers_capture,left_blockers_release,
                                                                                            right_blockers_capture,right_blockers_release))


# In[10]:


outf = "trajectory/LEFPositions.h5"
p = Path(outf)
if p.exists():
    p.unlink()
with h5py.File(outf, mode='w') as f:
    dset = f.create_dataset("positions", 
            shape=(steps, LEFNum, 2), 
            dtype=np.int32, 
            compression="gzip")
    bins = np.linspace(0, steps, num_chunks, dtype=int)
    for st,end in zip(bins[:-1], bins[1:]): # Loop through bins
        cur = []
        for i in range(st,end): # For bin in bins 
            positions = [(extruder.leg1.pos, extruder.leg2.pos) for extruder in EXTRUDERS] # Get both leg positions for all extruders
            cur.append(positions)
            for extruder in EXTRUDERS:
                occupied = extruder.translocate(occupied) # Translocate extruder
        cur = np.array(cur)
        dset[st:end] = np.array(cur)
        #print(dset[st:end])
    f.attrs["N"] = N
    f.attrs["LEFNum"] = LEFNum
del EXTRUDERS

