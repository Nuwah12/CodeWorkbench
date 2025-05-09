###################
# Noah Burget
# Last updated: 12/7/24
# Driver script for the 1-Dimensional Portion of the loop extrusion simulations
# Originally written as a jupyter notebook, but converted to .py script for simplicity
###################

### 1-D loop extrusion simulation
from extruder import Extruder
import numpy as np
from pathlib import Path
import seaborn as sns
import h5py
import sys


def le_oneD(name, B1_cap, B1_rel, prom_cap, prom_rel, B2_cap, B2_rel):
    #### This is a function for performing 1-D loop extrusion simulation and writing output to h5 file
    # Parameters
    # ----------
    # name: str
    #   a name for the simulation run, output parameters file will be named like this
    # params: dict
    #   a dictionary of parameter {name:value}. AS OF 12/20/24, the only parameters we are working to optimize are capture/release probabilities
    #   it is organized as {'blockerName':[capture_prob, release_prob]}
    ###############################################################
    ### Translocating and writing trajectories using custom extruder class
    # Basic parameters for simulation run
    ###############################################################
    RUN_NAME = f"{name}_simulation" #### Define a name for this simulation run, sowe can save parameters
    N1 = 900 # Size of 1 system
    M = 1 # No. of systems
    N = N1*M # Total size of system, in momomers
    occupied = np.zeros(N) # List to tell if current monomer is occupied by an extruder
    occupied[0] = 1 
    occupied[-1] = 1
    steps = 50000 # Timesteps for 1D sim.
    num_chunks = 50 # No. of chunks to write trajectories in
    LIFETIME = 800 # Cohesin lifetime
    LIFETIME_STALLED = LIFETIME // 10 # Cohesin lifetime when stalled

    ### Blockers (i.e. CTCF) - {pos. : prob.}
    left_blockers_capture = {}
    right_blockers_capture = {}
    left_blockers_release = {}
    right_blockers_release = {}

    #####################################################################
    ### Defining 'blocking regions' - consecutive monomers that can block
    # Each call below defines the blocking region as (low, high), EXCLUSIVE on the upper bound
    #####################################################################
    B1_blockingRegion = np.arange(220, 222) # B1, EBF1 involved
    promoter_blockingRegion = np.arange(405, 411) # EBF1 involved, ZEB2 promoter
    B2_blockingRegion = np.arange(592, 594) # B2, EBF1 involved

    ###########################################################################
    ### Any EBF1-associated blockers below should have '_EBF1' as the suffix in the keys of the dictionary below
    ###########################################################################
    blockingRegions = {'B1_EBF1':B1_blockingRegion, 'promoter_EBF1':promoter_blockingRegion, 'B2_EBF1': B2_blockingRegion}

    ### Assign each blocking region its capture and release probabilities (0<=x<=1)
    for br in blockingRegions:
        if 'B1' in br:
            cap = B1_cap
            rel = B1_rel
        elif 'promoter' in br:
            cap = prom_cap
            rel = prom_rel
        elif 'B2' in br:
            cap = B2_cap
            rel = B2_rel
        # Now, also according to the dictionary key, we assign the directional blockers with the previously defined probabilities
        if 'EBF1' in br: # IF EBF1 blocker, bidirectional (in the EBF1 loading system, and the perturbed EBF1 blocking system, this will never happen)
            for loc in blockingRegions[br]:
                left_blockers_capture[loc] = cap
                left_blockers_release[loc] = rel
                right_blockers_capture[loc] = cap
                right_blockers_release[loc] = rel
        # The below two conditions will only be evaluated if the blocker is NOT EBF1-associated
        ### Other conditionals here for blocker orientation based on dict keyi
        if 'EBF1' not in br:
            for loc in blockingRegions[br]:
                left_blockers_capture[loc] = cap
                left_blockers_release[loc] = rel
        # ...
   
    ########################################
    ### Now we define the cohesion LOADING regions
    ########################################
    wholePolymer_loading = [0, N1-1] # Entire polymer
    ZEB2_5prime_loading = [411, 592]
    ### List of loading regions for EBF1 loading
    #loading_regions = [E1_loading_1, E1_loading_2, E2_loading, B3_loading_1, MYC_loading, wholePolymer_loading]

    ### List of loading regions for EBF1 blocking (just the entire polymer)
    loading_regions = [wholePolymer_loading, ZEB2_5prime_loading]

    ### Number of cohesins per loading region (each index is associated with the corresponding index of the loading_regions list)
    loading_region_freqs = [11, 4]
    LEFNum = np.sum(loading_region_freqs)

    ### Beyond this point, you need not change any variable values.
    ###############################################################

    # Generate (initial) loading region spots for the loading region LEFs and list of loading site for each LEF
    LOADING_SPOTS = []
    REGIONS_INDEX = []
    for i, region in enumerate(loading_regions):
        print(i, region)
        for j in range(loading_region_freqs[i]):
            while True:
                spot = np.random.randint(low=region[0], high=region[1])
                if spot not in LOADING_SPOTS and spot+1 not in LOADING_SPOTS and spot-1 not in LOADING_SPOTS:
                    print(spot)
                    break
            LOADING_SPOTS.append(spot)
            REGIONS_INDEX.append(loading_regions[i])
    print(REGIONS_INDEX)

    EXTRUDERS = []
    targeted_idx = 0
    for i in range(len(LOADING_SPOTS)):
        #leg = np.random.randint(N)
        #print('Loading extruder at {}'.format(LOADING_SPOTS[i]))
        leg = LOADING_SPOTS[i]
        print('Loading extruder at {}, its loading region is {}'.format(leg, REGIONS_INDEX[i]))
        EXTRUDERS.append(Extruder(
                extruder_index = i,
                leg1 = leg,
                leg2 = leg+1,
                left_blockers_capture = left_blockers_capture,
                right_blockers_capture = right_blockers_capture,
                left_blockers_release = left_blockers_release,
                right_blockers_release = right_blockers_release,
                extrusion_occupancy = occupied,
                loading_region = REGIONS_INDEX[i],
                lifetime = LIFETIME,
                lifetime_stalled = LIFETIME_STALLED)
        )

    ### Write parameters to text file
    with open('/mnt/data0/noah/analysis/misc-analysis-local/ORCA/polymer_sims/optimizeSim/simulation_run_{}_dir/{}_params.txt'.format(name, RUN_NAME),'w') as pf:
        pf.write("N: {}\n1D Steps: {}\nTotal LEF: {}\nLoading Regions: {} LEFs per loading region: {}\n Lifetime: {} Lifetime stalled: {}\nBlocking Regions: {}\nLeft cap: {}, Left rel: {}\nRight cap: {}, Right rel: {}".format(
                                                                                                N1,steps,LEFNum,
                                                                                                REGIONS_INDEX,loading_region_freqs,
                                                                                                LIFETIME,LIFETIME_STALLED,blockingRegions,left_blockers_capture,left_blockers_release,
                                                                                                right_blockers_capture,right_blockers_release))
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
        f.attrs["N"] = N
        f.attrs["LEFNum"] = LEFNum
    del EXTRUDERS

if __name__ == '__main__':
    main()


#          _
#      .__(.)< (MEOW)
#       \___)
