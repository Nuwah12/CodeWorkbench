## 3D (Mol. Dynamics) simulation
This is the 3D portion of the loop extrusion simulation. Recorded LEF positions at each timestep are used to track progression and movement of the LEF with respect to its polymer positions.

#### Bond Identification
As described in the 1D simulation portion, extrusion is measured by the movement of an extruder along a polymer. At each timestep, this extruder is contacting two monomers within the polymer. These two monomers are considered **bonded**, and thus have the "active" bonded parameters set. Monomers that are bonded may or may not also be affected by nonbonding forces (`except_bonds` parameter). All other monomers are given the "inactive" parameters.
#### Adding forces
The forces added to our simulation are packaged into a `forcekit` (polychrom) describing a polymer chain. The forcekit contains 3 forces - the bonded force (Harmonic Bonds), angle force, and nonbonded force. 
* **Bonded force** (`forces.harmonic_bonds`, `openmm.HarmonicBondForce`) 
  * **k** (double) - harmonic force constant of the bond in $kJ/mol/nm^2$ 
  * **length** (double) - equilibrium length of the bond in nm
  \
These are the parameter values that are varied for 'active' and 'inactive' bonds. For *active* bonds: \
  `k = (2 * self.kT / self.conlen**2) / (simtk.unit.kilojoule_per_mole / simtk.unit.nanometer**2) * 1/bondWiggleDistance**2` \
  `length = bondDistance * length_scale` \
And for *inactive* bonds: \
  `k = 0` \
  `length = bondDistance * length_scale`

* **Angle force** describes interactions between triplets of monomers
  * **k** (float) - stiffness of the bond
  * **$\theta_0$** (float) - equilibrium angle of the bond \
  The energy of the bond at angle $\theta$ is described by a custom function, defined here as `kT * k * (theta - theta_0) * (theta - theta_0) * 0.5`

* **Nonbonded force** describes the interactions of monomers that are not regarded to be participating in a bond. It is defined by the module `forces.polynomial_repulsive` (polychrom wrapper) which wraps `openmm.CustomNonbondedForce` \
  * `trunc` (float) - energy value at dist=0
  * The repulsion energy is described by the algerbraic expression:
    ```
        rsc12 * (rsc2 - 1.0) * (trunc * kT) / emin12 + (trunc * kT);
        rsc12 = rsc4 * rsc4 * rsc4;
        rsc4 = rsc2 * rsc2;
        rsc2 = rsc * rsc;
        rsc = r / radius * rmin12;
    ```
### Visualizing conformations
This is the workflow I use to go from simulation outputs (conformations) to an NxN aggregate matrix showing contact frequency.
1. Run a simulation. This will output as many 'blocks' as the number of simulation initializations you did in `.h5` format.
2. Use `trajectory_to_txt.py <total number of confs> <path to h5> <matrix our DIR>` to dump `.h5` files to text.
3. Run `python3 make_contactMaps <dir of confs> <matrix out file>` to obtain the aggregate matrix, where the specified input directory holds the previously created text files. Plot with your favorite plotting library.
