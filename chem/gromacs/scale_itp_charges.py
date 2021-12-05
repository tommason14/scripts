#!/usr/bin/env python3
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("files", nargs="+", help="itp file(s) to alter")
parser.add_argument(
    "-c", "--charge", help="Scaling factor, default=0.8", default=0.8, type=float
)

args = parser.parse_args()


def scale(itp, scaling_factor):
    newname = itp.replace(".itp", "-scaled.itp")
    with open(itp) as f:
        lines = f.readlines()

    newfile = []
    found = False
    for line in lines:
        if "[ atoms ]" in line:
            found = True
            newfile.append(line)
            continue
        if re.search("^[\[#]", line):
            found = False
        if found and not line.startswith(";") and not re.search("^\s*$", line):
            line = line.split()
            line[6] = f"{float(line[6]) * scaling_factor:.8f}"
            line.append("\n")
            newfile.append(" ".join(line))
        else:
            newfile.append(line)

    with open(newname, "w") as new:
        for line in newfile:
            new.write(line)


if __name__ == "__main__":
    for f in args.files:
        scale(f, args.charge)
