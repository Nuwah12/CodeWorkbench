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
#### Integration
A [**Variable Langevin Integrator**](http://docs.openmm.org/latest/userguide/theory/04_integrators.html?highlight=variablelangevin) is used to calculate future states for the 3D simulation. A Langevin Integrator similates a system in contact with a heat bat by integrating the Langevin equation of motion:
```math
m_i\frac{d\textbf{v}_i}{dt} = f_i - \gamma m_i v_i + R_i
```
Where $v_i$ is the velocity of particle $i$, $f_i$ is the force acting on it, $m_i$ its mass, and $\gamma$ the friction cpefficient, $R_i$ an uncorrelated random forse chosen from a normal distribution with 0 mean and unit variance. \
The *variable* part comes into play when we consider the timestep $\Delta t$. Instead of using a fixed timestep, it continuously updates the step size to keep the integration error below a user-defined threshold. To estimate the integration error, it compares positions calculated by verlet integration with those that would be given by an explicit [Eruler integration](https://physics.umd.edu/hep/drew/numerical_integration/):
```math
\text{error}=(\Delta t)^2 \sum_i \frac{|\bf{f}_i|}{m_i}
```
It selects the value of $\Delta t$ that makes the error exactly equal to the specified error tolerance, i.e. it solves for $\Delta t$ in the above equation.
**Why use a variable time step integrator**? These integrators are usually superior to fixed time step integrators in both stabvility and efficiency. Step sizes are automatically reduced to preserve accuracy and avoid instability when large forces occur. Read more on the benefits [here](http://docs.openmm.org/latest/userguide/theory/04_integrators.html?highlight=variablelangevin#variableverletintegrator).
