#!/usr/bin/env python3
import MDAnalysis as mda
import sys

if len(sys.argv) != 3 or any("-h" in x for x in sys.argv[1:]):
    print("Syntax: last_frame.py coords traj")
    sys.exit(1)

u = mda.Universe(*sys.argv[1:])
u.trajectory[-1]
u.atoms.write("last_frame.gro")
