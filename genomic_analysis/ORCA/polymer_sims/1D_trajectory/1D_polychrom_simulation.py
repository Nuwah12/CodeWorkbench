#!/usr/bin/env python
# coding: utf-8

# ## 1-D loop extrusion simulation
# ### !!! Make sure you run the cells IN ORDER !!!

# In[76]:


from extruder import Extruder
import numpy as np
from pathlib import Path
import seaborn as sns
import h5py
import sys


# ### Translocating and writing trajectories using custom extruder class

# In[77]:


RUN_NAME = "scenario19_v8" #### Define a name for this simulation run, sowe can save parameters
N1 = 900 # Size of 1 system
M = 1 # No. of systems
N = N1*M # Total size of system, in momomers
occupied = np.zeros(N) # List to tell if current monomer is occupied by an extruder
occupied[0] = 1 
occupied[-1] = 1
steps = 50000 # Timesteps for 1D sim.
LEFNum = 15 # No. of RANDOMLY LOADED extruders
num_chunks = 50 # No. of chunks to write trajectories in

### Blockers (i.e. CTCF) - {pos. : prob.}
left_blockers_capture = {}
right_blockers_capture = {}
left_blockers_release = {}
right_blockers_release = {}
## Here we define strong, mild, and weak [blocking, release]
STRONG_BLOCK = [0.999, 10**-10]
MEDIUM_BLOCK = [0.5, 0.02]
WEAK_BLOCK = [0.2, 0.025]
# Manually assigning blockers in dict
#left_blockers_capture[10] = STRONG_BLOCK[0] # 5' CTCF
#left_blockers_release[10] = STRONG_BLOCK[1]
left_blockers_capture[195] = STRONG_BLOCK[0] # E1a
left_blockers_release[195] = STRONG_BLOCK[1]
left_blockers_capture[225] = STRONG_BLOCK[0] # E1b
left_blockers_release[225] = STRONG_BLOCK[1]
left_blockers_capture[315] = MEDIUM_BLOCK[0] # E2
left_blockers_release[315] = MEDIUM_BLOCK[1]
left_blockers_capture[405] = WEAK_BLOCK[0] # B1
left_blockers_release[405] = WEAK_BLOCK[1]
left_blockers_capture[555] = MEDIUM_BLOCK[0] # B3
left_blockers_release[555] = MEDIUM_BLOCK[1]
#left_blockers_capture[735] = STRONG_BLOCK[0] # MYC promoter
#left_blockers_release[735] = STRONG_BLOCK[1]

#right_blockers_capture[10] = STRONG_BLOCK[0] # 5' CTCF
#right_blockers_release[10] = STRONG_BLOCK[1]
#right_blockers_capture[195] = STRONG_BLOCK[0] # E1
#right_blockers_release[195] = STRONG_BLOCK[1]
right_blockers_capture[315] = WEAK_BLOCK[0] # E2
right_blockers_release[315] = WEAK_BLOCK[1]
right_blockers_capture[495] = WEAK_BLOCK[0] # B2
right_blockers_release[495] = WEAK_BLOCK[1]
#right_blockers_capture[555] = MILD_BLOCK[0] # B3
#right_blockers_release[555] = MILD_BLOCK[1]
right_blockers_capture[735] = STRONG_BLOCK[0] # MYC promoter
right_blockers_release[735] = STRONG_BLOCK[1]

#for i in range(M):
#    for locs in left_blocker_locs:
#        pos = i * N1 + locs
#        left_blockers_capture[pos] = 0.90
#        left_blockers_release[pos] = 0.01
#    for locs in right_blocker_locs:
#        pos = i * N1 + locs
#        right_blockers_capture[pos] = 0.90
#        right_blockers_release[pos] = 0.01


# In[78]:


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


# In[79]:


### Write parameters to text file
with open('/mnt/data0/noah/analysis/misc-analysis-local/ORCA/polymer_sims/3D_simulation/{}_params.txt'.format(RUN_NAME),'w') as pf:
    pf.write("N: {}\n1D Steps: {}\nTotal LEF: {}\nTargeted: {} Probs: {}\nLifetime: {} Lifetime stalled: {}\nLeft cap: {}, Left rel: {}\nRight cap: {}, Right rel: {}".format(
                                                                                            N1,steps,LEFNum,
                                                                                            TARGETED_LOADING_SPOTS,LOADING_PROBS,
                                                                                            LIFETIME,LIFETIME_STALLED,left_blockers_capture,left_blockers_release,
                                                                                            right_blockers_capture,right_blockers_release))


# In[80]:


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

