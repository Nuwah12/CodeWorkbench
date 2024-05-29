#!/usr/bin/env python
# coding: utf-8

# ### Simple MD (3D) simulation

# In[41]:


import time
import numpy as np
import h5py
from bondUpdater import bondUpdater
from polychrom.starting_conformations import grow_cubic
from polychrom.simulation import Simulation
from polychrom.hdf5_format import HDF5Reporter
import polychrom.forcekits as forcekits
import polychrom.forces as forces
import matplotlib.pyplot as plt


# ### Define parameters for 3D portion

# In[42]:


### Gather parameters from the 1D portion
trajectories = h5py.File("../1D_trajectory/trajectory/LEFPositions.h5") # Saved trajectories from 1D siumulation
N = trajectories.attrs["N"] # Length of polymer
LEFNum = trajectories.attrs["LEFNum"] # Number of extruders
LEFpositions = trajectories["positions"] # Positions of extruders at each 1D step
Nframes = LEFpositions.shape[0] # Number of 1D steps (= number of extruder steps)

print("""
Polymer is {} monomers long. There are {} Extruders loaded. 
It was processed in {} steps.
""".format(N,LEFNum,Nframes))

### Set molecular dynamics parameters
steps = 500 # MD steps PER STEP OF EXTRUDER
box = (N / 0.1) ** 0.35 # Dimensions of bounding box with Periodic Boundary Conditions (PBC)
data = grow_cubic(N, int(box)) # Initialize random-walk chains for our polymers
# SMC (Extruder) parameters
smcBondWiddleDist = 0.2
smcBondDist = 0.5

### Simulation saving parameters
saveEveryBlocks = 5 # Write coordinates every this many blocks
restartSimulationEveryBlocks = 100 # 
# Checks
assert Nframes % restartSimulationEveryBlocks == 0 # So we don't have leftover steps that won't get saved
assert (restartSimulationEveryBlocks % saveEveryBlocks) == 0

savesPerSim = restartSimulationEveryBlocks // saveEveryBlocks
simInitsTotal = Nframes // restartSimulationEveryBlocks # Number of simulation initializations

print("""
There will be {} MD steps done for every step of the extruder (aka there will be {} steps PER ""BLOCK"")
Simulation restarts every {} blocks, for a total of {} initializations.
Each simulation run will produce {} conformations, for a total of {} conformations.
      """.format(steps,steps,restartSimulationEveryBlocks,simInitsTotal,int(restartSimulationEveryBlocks/saveEveryBlocks),int((restartSimulationEveryBlocks/saveEveryBlocks)*simInitsTotal)))


# ### The Simulation Loop

# In[43]:


milker = bondUpdater(LEFpositions)

reporter = HDF5Reporter(folder="sim_outs", # Save data location
                        max_data_length=100, # Write data in chunks of this size - THIS CONTROLS HOW MANY CONFIGS ARE IN EACH BLOCK!
                        overwrite=True, # overwrite existing file in out location
                        blocks_only=True) # only save simulation blocks

for iter in range(simInitsTotal):
    # Create the simulation object
    a = Simulation(
            platform="cuda", # platform to do computations on
            integrator="variableLangevin", # Integrator from OpenMM
            error_tol=0.01, # error rate parameter for variableLangevin integrator
            GPU="0", # GPU index
            collision_rate=0.03, # collision rate of particles in inverse picoseconds
            N=len(data), # no. of particles
            reporters=[reporter], # list of reporter objects to use
            PBCbox=[box,box,box], # Periodic Boundary Conditions (PBC) box dimensions (x,y,z)
            precision="mixed" # GPU calculation precision, mixed is slow on 3080 and newer GPUs
    )
    # Loads the polymer we created, and puts center of mass at (0,0,0)
    a.set_data(data) 
    # Add a force to the simulation object - since we are doing polymer simulation, we add a 'forcekit' that describes all the forces in a polymer chain and the interactions between them
    a.add_force(
        forcekits.polymer_chains(
            a, # Simulation object
            chains=[(0, None, 0)], # List of tuples desctibing 1 chain each - this is the default value, i.e. one chain of length N that is not a ring (i.e. a chain)
            bond_force_func=forces.harmonic_bonds, # Define the bonded force as harmonic bonds
            bond_force_kwargs={'bondLength':1.0, 'bondWiggleDistance':0.05}, # Parameters for harmonic bonds
            angle_force_func=forces.angle_force, # Angle force 
            angle_force_kwargs={'k':1.5}, # Angle force parameters. k = stiffness bond (8=very stiff, k=1.5 is "realistically flexible")
            nonbonded_force_func=forces.grosberg_repulsive_force, # Nonbonded force
            nonbonded_force_kwargs={'trunc':1.5, # Allows chains to cross, the energy value at dist=0
                                    'radiusMult':1},
            except_bonds=True # Nonbonded forces do not affect bonded pieces
        )
    )
    #a.add_force(forces.spherical_confinement(a,density=0.2)) # Confine polymer in a sphere
    # Calculate bond parameters for extruder contact
    kbond = a.kbondScalingFactor / (smcBondWiddleDist**2)
    bondDist = smcBondDist * a.length_scale
    activeParams = {"length":bondDist, "k":kbond}
    inactiveParams = {"length":bondDist, "k":0}
    # Set up bond manager object ("milker")
    milker.setParams(activeParams, inactiveParams)
    milker.setup(bondForce=a.force_dict["harmonic_bonds"], blocks=restartSimulationEveryBlocks)

    # During the first simulation initiation, minimize energy of conformations
    if iter == 0:
        a.local_energy_minimization()
    else:
        a._apply_forces()
    ########## Start of the actual physics/MD calculations ##########
    for i in range(restartSimulationEveryBlocks): # Loop for our simulation length
        if i % saveEveryBlocks == (saveEveryBlocks-1): ### THIS IS WHERE WE SAVE A BLOCK!!! At the last step of the simulation before we restart
            a.do_block(steps=steps) # do steps AND GET new monomer positions consisting of <steps> steps
        else:
            a.integrator.step(steps) # do steps WITHOUT getting new monomer positions (faster)
        if i < restartSimulationEveryBlocks - 1: # if this is not the final block...
            curBonds, pastBonds = milker.step(a.context) # Update bonds with the milker
    data = a.get_data() # Fetch new polymer positions 
    del a 

    reporter.blocks_only = True # Write only blocks, not individual steps in block
    time.sleep(0.2) # wait so garbage collector can clean up

reporter.dump_data() # Output
