#!/usr/bin/env python3
import MDAnalysis as mda
from MDAnalysis import transformations
import sys


def unwrap_traj(universe):
    """
    Unwrap trajectories of all frames in Universe.
    Apply to universe defined in the global scope. 
    i.e.
    unwrap_traj(u)
    not
    u = unwrap_traj(u)
    """
    a = universe.atoms
    transform = mda.transformations.unwrap(a)
    universe.trajectory.add_transformations(transform)


if len(sys.argv) != 3 or '-h' in sys.argv:
    print('Syntax: mdanalysis_unwrap_traj.py data.lmps traj.dcd')
    sys.exit(1)

u = mda.Universe(sys.argv[1], sys.argv[2], topology_format='DATA')
unwrap_traj(u)
with mda.Writer('.'.join(sys.argv[2].split('.')[:-1] + ['unwrapped', 'dcd']),
                len(u.atoms)) as w:
    for ts in u.trajectory:
        w.write(u)
