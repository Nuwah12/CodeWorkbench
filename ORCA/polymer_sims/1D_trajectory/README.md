### 1D Loop-extusion simulation
##### The functions defined in this notebook describe the behavior or cohesins and their interactions with eachother and with CTCFs in terms of the polymer's trajectory (how the polymer's elements move relative to one another)
##### Once the simulation finishes, the necessary calculated information is stored in a `trajectories/` folder, which is then used as input into the 3D (molecular dynamics) part of the simulation.

#### Workflow:
1. Load cohesin 
    * Choose random point [1,N) to load a cohesin onto
    * If `point` and `point+1` are unoccupied, mark them as occupied, else, choose a different point
    * Append a `cohesin(point, point+1)` object to the cohesin list
2. Translocate cohesin \
    a. **Unload cohesin** - With probability `q`, mark occupancy at legs as 0 \
    b. **Load cohesin** - With probability `p`, mark occupancy at legs as 1 \
    c. **Translocate cohesin** - As long as the cohesin can move (i.e. it is not located next to another cohesin and is not at the end of the polymer), mark legs +/- 1 as occupied, and current legs as unoccupied, 'extruding' the loop through the cohesin ring.

Files in this directory:
`extrusion_1D_trajectory_example.ipynb` - An example setting up necessary functions for cohesin behavior and demonstrating 1D extrusion with a small polymer, printing progress messages so that an inutition of the process can be formed. No trajectory is saved.
    
