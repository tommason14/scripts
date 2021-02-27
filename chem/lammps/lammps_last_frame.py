#!/usr/bin/env python3
import sys
from ovito.io import import_file, export_file

if len(sys.argv) not in [2, 3] or sys.argv[1] == "-h":
    print("Syntax: lammps_last_frame.py lammps_dump [output]")
    sys.exit(1)

try:
    output = sys.argv[2]
except IndexError:
    output = "last_frame.xyz"

data = import_file(sys.argv[1])
last = data.compute(data.source.num_frames)
export_file(
    last,
    output,
    "xyz",
    columns=["Particle Type", "Position.X", "Position.Y", "Position.Z"],
)
