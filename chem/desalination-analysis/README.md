# Desalination analysis

Analysing trajectories where a polymer chain is placed in between two
water reservoirs - one salinated on the left, one freshwater on the right.

Simulations are assumed to have been perfomed with Gromacs, as tpr files are
needed when unwrapping. 
If a trajectory has already been processed, this module can be used with any 
files that are readable by MDAnalysis.

Analysis is performed purely in Python through MDAnalysis, however Gromacs is used to
unwrap large trajectories, so try to have access to a gmx / gmx_mpi executable.

# Example workflows

To compute ion counts:
```python
from polfuncs import PolSim, plot_ion_counts
nvt = PolSim('nvt.tpr', 'nvt.xtc')
nvt.unwrap()
counts = nvt.compute_ion_counts()  # counts is a pandas DataFrame
plot_ion_counts(counts, fname='ion_counts.png')

# or using the pandas DataFrame pipe:
nvt.compute_ion_counts().pipe(plot_ion_counts, fname='ion_counts.png')
```

To compute partial densities:
```python
from polfuncs import PolSim, plot_partial_densities
nvt = PolSim('nvt.tpr', 'nvt.xtc')
nvt.compute_partial_densities().pipe(plot_partial_densities, fname='partial_densities.png')
```

> When computing partial densities, do not unwrap the trajectory - for a 4 ns
> test run, the computation took 4.7 seconds on a wrapped trajectory, and 4 min
> 57 seconds when unwrapped. The reason is that MDAnalysis has to then wrap the
> coordinates of each atom back inside the box before computing densities of
> each slice.
