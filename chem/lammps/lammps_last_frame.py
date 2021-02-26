#!/usr/bin/env python3
import sys
from ovito.io import import_file, export_file

if len(sys.argv) != 2 or sys.argv[1] == '-h':
    print('Syntax: lammps_last_frame.py lammps_dump')
    sys.exit(1)

data = import_file(sys.argv[1])
last = data.compute(data.source.num_frames)
export_file(last, 'last_frame.xyz', 'xyz',
columns=['Particle Type', 'Position.X', 'Position.Y', 'Position.Z'])
