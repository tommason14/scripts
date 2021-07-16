# Transferable MD Analysis

Analysis tools written using the python MDAnalysis library.
As such, the tools accept a multitude of file formats and should work
with Gromacs/Amber/OpenMM simulation output with no issues.

For LAMMPS simulations, ensure the trajectory is written in the 'atom' style (`dump ... atom`).
This is the only LAMMPS trajectory format accepted by MDAnalysis. Alternatively, write trajectories
to a DCD format (`dump ... dcd`), which are handled completely.
and then edit the `MDAnalysis.Universe` constructor in each of the scripts if necessary.

## Converting trajectories (traj_convert.py)

- Pass in coordinates and topology files, with the output format detected from
  the filetype of the output file

## Diffusion coefficients (diff_coeffs_from_md.py)

- Computed via the Einstein relation (1/2N * gradient of MSD against time, N = dimensionality)
- The middle 80 % of trajectory is used by default to compute diffusion coefficients, but easily changed
- Diffusion coefficients printed to the screen, MSD vs time plot saved as pdf

## Dipole moments (dipole_moments_from_md.py)

- Calculated relative to the centre of mass of each molecule in the simulation, so one moment per molecule
per frame
- Moments in each direction along with the absolute value saved for every molecule in each frame
- For simpler plotting, total dipole moments for each molecule are also grouped into intervals and data is saved
in order for a histogram to be plotted
