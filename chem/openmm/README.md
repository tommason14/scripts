# OpenMM scripts

Note that scripts may include the `ommhelper` module - this can be found [here](https://github.com/z-gong/ms-tools/tree/master/mstools/ommhelper). Download and add the directory to your `PYTHONPATH` shell environmental variable.

## Creating an OPLS-based system

[create_openmm.py](create_openmm.py) uses the [mstools](https://github.com/z-gong/ms-tools) package to create simulation files based on
the OPLS forcefield. Polarisable forcefields such as the [CL&Pol FF](https://github.com/paduagroup/clandpol) are
supported.

## Running simulations with Gromacs files (forcefield agnostic)

[omm_run_gmx.py](omm_run_gmx.py) loads in Gromacs simulation files using Parmed
and then runs an OpenMM simulation. Options include:

- charge scaling
- position restraints based on residue names
- cosine acceleration based on the periodic-perturbation method of [Hess](https://aip.scitation.org/doi/10.1063/1.1421362), used for viscosity predictions. (Requires the velocity-verlet plugin found [here](https://github.com/z-gong/openmm-velocityVerlet)

Example usage of a run at 300 K and 1 bar, with a 1 fs timestep for 10 ns,
restraining residues named 'pol', with minimisation before the run:
`omm_run_gmx.py --gro conf.gro --top topol.top -t 300 -p 1 -dt 1 -n 10_000_000 -r pol -min`
