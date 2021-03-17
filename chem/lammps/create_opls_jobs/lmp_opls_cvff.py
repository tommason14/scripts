#!/usr/bin/env python3
from get_coeffs import getAtomData, getBond, getImproperCVFF
from get_coeffs import getAngle, getDihedral
from get_coeffs import getAtomPartialCharge, getMass
import subprocess as sp
import glob, re, sys
""" Substitute in different labels to *-l.data files and create
topology file with VMD topo. After VMD original names are
substituted back in and coefficients, box size and atom data are
added. This includes Pair Coeffs, Bond Coeffs, Angle Coeffs,
Dihedral Coeffs, Improper Coeffs and Atoms (partial charges)
sections and box dimensions - xlo, xhi, ylo, yhi, zlo, zhi. """

# PATTERN MATCHING
pattern_labelled_xyz = "^\s*[A-Z]*[a-z]?[0-9]*\s+-?[0-9]\."
pattern_coeffs = "^# [0-9]+\s+[A-Z]{1,2}[a-z]?[0-9]*"

if len(sys.argv) < 3:
    sys.exit("Syntax: lmp_opls.py xyzfile ff_file [output]")

File = sys.argv[1]
ff = sys.argv[2]
try:
    output = sys.argv[3]
except IndexError:
    output = File.replace(".xyz", ".data")

# NUMBER OF NEW ELEMENTS
count = 0

# NEW ELEMENT NAMES
elemDict = {}

# NEW LINES
newLines = []
with open(File, "r") as f:
    for line in f.readlines()[2:]:
        sym, coords = line.strip().split(" ", 1)
        if sym not in elemDict.keys():
            # unique atoms as C001, H002 etc...
            count += 1
            elemDict[sym] = sym[0] + f"{count:03d}"

        newID = elemDict[sym]
        newLines.append(newID + "   " + coords + "\n")
with open("topo-in.xyz", "w") as new:
    new.write(f"{len(newLines)}\n\n")
    for line in newLines:
        new.write(line)

##### RUN VMD TOPO ----------------------------------------

# COMMANDS FOR VMD
lines = ("package require topotools\n" + "topo retypebonds\n" +
         "topo guessangles\n" + "topo guessdihedrals\n" +
         "topo guessimpropers\n" + "set sel [atomselect top all]\n" +
         "$sel set resid [$sel get fragment]\n" +
         "topo writelammpsdata topo.out\n" + "exit")

open("tempfile", "w+").write(lines)

# OPEN VMD WITH EDITTED XYZ AND COMMAND FILE
user = sp.check_output("echo $USER", shell=True).decode("utf8").strip()
if user == "tommason":
    vmd = "/Applications/VMD\ 1.9.4a38.app/Contents/vmd/vmd_MACOSXX86_64"
elif user == "tmas0023":
    vmd = "/Applications/VMD\ 1.9.3.app/Contents/vmd/vmd_MACOSXX86"
else:
    vmd = "vmd"  # module loaded
cmd = f"{vmd} -dispdev none -m topo-in.xyz -e tempfile"
sp.check_output(cmd, shell=True)

### SUBSTITUTE IN ORIG NAMES ------------------------------

newLines = []

# SWITCH VALUES AND KEYS IN DICT FOR EASY ACCESS
elemDict = {y: x for x, y in elemDict.items()}

# FOR EACH LINE IN TOPO FILE
lines = open("topo.out", "r+").readlines()
for line in lines:

    # FOR EACH NEW SYMBOL TO REPLACE
    for newID in elemDict.keys():

        # IF SYMBOL AFTER A HASH
        if re.search("#.*" + newID, line):
            # GET ORIG NAME
            sym = elemDict[newID]

            # REPLACE NEW WITH OLD STRING
            line = line.replace(newID, sym)

    newLines.append(line)

### EDIT DATA FILE ----------------------------------------

# READ LINES
lines = newLines[:]

# SWITCHES - TURN ON WHEN FIND IN FILE
pcoef = False  # Pair Coeffs
bcoef = False  # Bond Coeffs
acoef = False  # Angle Coeffs
dcoef = False  # Dihedral Coeffs
icoef = False  # Improper Coeffs
atoms = False  # Atoms

# SAVE X, Y, Z VALUES
xvals = []
yvals = []
zvals = []

# SAVE PARTIAL CHARGES
pcharges = []

# START EDITTING LINES
newLines = []

for line in lines:

    # masses???

    # PAIR COEFFS SWITCH
    if "Pair Coeffs" in line:
        pcoef = True
        newLines.append("Pair Coeffs\n")

    # BOND COEFFS SWITCH
    elif "Bond Coeffs" in line:
        pcoef = False
        bcoef = True
        newLines.append("Bond Coeffs\n")

    # ANGLE COEFFS SWITCH
    elif "Angle Coeffs" in line:
        bcoef = False
        acoef = True
        newLines.append("Angle Coeffs\n")

    # DIHEDRAL COEFFS SWITCH
    elif "Dihedral Coeffs" in line:
        acoef = False
        dcoef = True
        newLines.append("Dihedral Coeffs\n")

    # IMPROPER COEFFS SWITCH
    elif "Improper Coeffs" in line:
        dcoef = False
        icoef = True
        newLines.append("Improper Coeffs\n")

    # ATOMS PARTIAL CHARGES SWITCH
    elif "Atoms" in line:
        dcoef = False
        icoef = False
        atoms = True
        newLines.append(line)

    # ADD PAIR COEFFS
    elif pcoef and re.search(pattern_coeffs, line):
        h, num, name = line.split()
        pars, pot = getAtomData(name, ff)
        newline = "{:4} {:>9.5f} {:>10.3f}    # {}\n".format(
            num, pars, pot, name)
        newLines.append(newline)

    # ADD BOND COEFFS
    elif bcoef and re.search(pattern_coeffs, line):
        h, num, atms = line.split()
        atm1, atm2 = atms.split("-")
        kr, Re = getBond(atm1, atm2, ff)
        newline = "{:4} {:>9.5f} {:>10.3f}    # {}\n".format(num, kr, Re, atms)
        newLines.append(newline)

    # ADD ANGLE COEFFS
    elif acoef and re.search(pattern_coeffs, line):
        h, num, atms = line.split()
        atm1, atm2, atm3 = atms.split("-")
        ka, th = getAngle(atm1, atm2, atm3, ff)
        newline = "{:4} {:>9.5f} {:>10.3f}    # {}\n".format(num, ka, th, atms)
        newLines.append(newline)

    # ADD DIHEDRAL COEFFS
    elif dcoef and re.search(pattern_coeffs, line):
        h, num, atms = line.split()
        atm1, atm2, atm3, atm4 = atms.split("-")
        v1, v2, v3, v4 = getDihedral(atm1, atm2, atm3, atm4, ff)
        newline = "{:4} {:>9.5f} {:>10.5f} {:>10.5f} {:>10.5f}   # {}\n".format(
            num, v1, v2, v3, v4, atms)
        newLines.append(newline)

    # ADD IMPROPER COEFFS
    elif icoef and re.search(pattern_coeffs, line):
        h, num, atms = line.split()
        atm1, atm2, atm3, atm4 = atms.split("-")
        v1, v2, v3 = getImproperCVFF(atm1, atm2, atm3, atm4, ff)
        newline = "{:4} {:>9.5f} {:>4d} {:>4d}  # {}\n".format(
            num, v1, v2, v3, atms)
        newLines.append(newline)

    # ADD PARTIAL CHARGES AND FIND BOX MIN & MAX
    elif atoms and re.search("^[0-9]{0,5} [0-9]{0,5} [0-9]{0,5}", line):
        ord, num, ID, pc, x, y, z, h, name = line.split()

        # SAVE VALUES OF X, Y, Z
        xvals.append(float(x))
        yvals.append(float(y))
        zvals.append(float(z))

        # GET NEW PARTIAL CHARGE
        pc = getAtomPartialCharge(name, ff)
        newline = "{:4} {:3} {:3} {:>6.3f} {:>11.6f} {:>11.6f} {:>11.6f}   # {}\n".format(
            ord, num, ID, pc, float(x), float(y), float(z), name)
        newLines.append(newline)
        pcharges.append(pc)

    # FINISHED EDITTING LINES
    elif re.search("^\s*Bonds\s*$", line):
        atoms = False
        newLines.append(line)

    # SAVE UNCHANGED LINES
    else:
        newLines.append(line)

# SORT X, Y, Z, VALS FOR MIN MAX VALUES
xvals.sort()
yvals.sort()
zvals.sort()

# SUB IN BOX MAX AND MIN AND REMOVE HASH ON BLANK LINES
for i in range(len(newLines)):

    if "xlo" in newLines[i]:
        newLines[i] = "{:6.3f} {:6.3f}  xlo xhi\n".format(
            xvals[0] - 1, xvals[-1] + 1)

    elif "ylo" in newLines[i]:
        newLines[i] = "{:6.3f} {:6.3f}  ylo yhi\n".format(
            yvals[0] - 1, yvals[-1] + 1)

    elif "zlo" in newLines[i]:
        newLines[i] = "{:6.3f} {:6.3f}  zlo zhi\n".format(
            zvals[0] - 1, zvals[-1] + 1)

    elif re.search("^\s*#\s*$", newLines[i]):
        newLines[i] = "\n"

# Manually add masses so that odd atom names don't result in the wrong mass
found_masses = False
rm_start = 0
rm_end = 0
for ind, line in enumerate(newLines):
    if 'Masses' in line:
        found_masses = True
        rm_start = ind
        continue
    if found_masses:
        if re.search('^\s*[A-Z]', line):
            rm_end = ind
            break

newLines = newLines[:rm_start] + newLines[rm_end:]

found_pairs = False
masses = []
for line in newLines:
    if "Pair Coeffs" in line:
        found_pairs = True
        continue
    if "Bond Coeffs" in line:
        break
    if found_pairs and not re.search("^\s*$", line):
        line = line.split()
        newline = [line[0], getMass(line[-1], ff), line[-1]]
        newline = "{:4} {:>9.4f}    # {}\n".format(*newline)
        masses.append(newline)
masses = ["\n", "Masses\n", "\n"] + masses

# add after box length
zlo = 0
for ind, line in enumerate(newLines, 1):
    if "zlo" in line:
        zlo = ind
        break
newLines = newLines[:zlo] + masses + newLines[zlo:]

# REMOVE EXCESS FILES
sp.check_output("rm topo-in.xyz topo.out tempfile", shell=True)

# add new line below first line- polymatic's pack.pl needs it
# lammps doesn't mind it
newLines.insert(1, "\n")

# WRITE .data FILE
open(output, "w+").writelines(newLines)

# PRINT SUM OF PARTIAL CHARGES
print(f"{File}   Charge: {sum(pcharges):5.5}")
