# Transferable MD Analysis

Analysis tools written using the python MDAnalysis library.
As such, the tools accept a multitude of file formats and should work
with Gromacs/Amber/OpenMM simulation output with no issues.

For LAMMPS simulations, ensure the trajectory is written in the 'atom' style (`dump ... atom`).
This is the only LAMMPS trajectory format accepted by MDAnalysis. Alternatively, write trajectories
to a DCD format (`dump ... dcd`), which are handled completely.
and then edit the `MDAnalysis.Universe` constructor in each of the scripts if necessary.

All scripts have help messages available - pass in the `-h` flag to see the available options.


## Converting trajectories (traj_convert.py)

- Pass in coordinates and topology files, with the output format detected from
  the filetype of the output file

## Extracting a selection of atoms (extract_selection.py)

- Alternative to `traj_convert.py` that writes a `gro` coordinate file and
  `xtc` binary trajectory of a selection provided by the user

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

## Radial distribution functions (rdf_from_md.py)

- Given reference and selection atom groups, an rdf is calculated within 15 angstroms of the reference

## RMSD (rmsd_from_md.py)

- Root mean squared deviation of a given selection
- Optional mass-weighting
- CSV saved, optional plot created

## RMSF (rmsf_from_md.py)

- Root mean squared fluctuations for each residue of a given selection
- CSV saved, optional plot created

## Radius of gyration (rgyr_from_md.py)

- Radius of gyration found over time for a given selection, and optionally
  for each residue of that selection
- CSV saved, optional plot created
