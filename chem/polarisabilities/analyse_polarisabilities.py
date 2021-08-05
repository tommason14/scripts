#!/usr/bin/env python3
import numpy as np
import sys
import re
from glob import glob

# From http://www.mdy.univie.ac.at/python-stuff/atomic_polarizabilities_charge.py


print("---------------------------------------------------------")
print("---------------------Welcome-----------------------------")
print("---------------------------------------------------------")

# ==============================================================================================
# This part reads the input files for the connectivity
# ==============================================================================================
read_input = open("connected.in", "r")
tmp = "".join([n for n in read_input.readlines()][:])  # Read input file into variable
read_input.close()
bond = re.findall(r"bond\s*\d+\s+\d+\s*", tmp)  # Search for all bonds
bonds = len(bond)
ring = re.findall(r"ring[\s\d]+", tmp)  # Search for all rings
rings = len(ring)

if bonds == 0:
    sys.exit("No bonds found, check input format.")
if len(re.findall(r"[\S ]+\n", tmp)) > bonds + rings:
    answer = input(
        "Not all lines in the input file contained readable information. Continue nevertheless?  yes/no   "
    )
    if answer == "yes":
        print("---------------------------------------------------------")
        print(
            "Continue... if program fails at a later point, check connectivity input file."
        )
    else:
        if answer == "no":
            print("---------------------------------------------------------")
            sys.exit("Program stopped. Check connectivity input file.")
        else:
            print("---------------------------------------------------------")
            sys.exit("Unknown option")
print("Molecule contains " + str(bonds) + " bonds and " + str(rings) + " rings.")

# =============================================================================================
# This part reads the GDMA output for six files (field in positive and
# negative x, y and z direction) and creates arrays for atoms names ("name"),
# coordinates ("c"), charges ("q") and dipoles ("d")
# =============================================================================================
print("Reading files...")

read_input_x = open(glob("plusX/**/gdma.out", recursive=True)[0], "r")
read_input_y = open(glob("plusY/**/gdma.out", recursive=True)[0], "r")
read_input_z = open(glob("plusZ/**/gdma.out", recursive=True)[0], "r")
read_input_mx = open(glob("minusX/**/gdma.out", recursive=True)[0], "r")
read_input_my = open(glob("minusY/**/gdma.out", recursive=True)[0], "r")
read_input_mz = open(glob("minusZ/**/gdma.out", recursive=True)[0], "r")

tmp = [
    "".join([n for n in read_input_x.readlines()][:]),  # Read input files into variable
    "".join([n for n in read_input_y.readlines()][:]),
    "".join([n for n in read_input_z.readlines()][:]),
    "".join([n for n in read_input_mx.readlines()][:]),
    "".join([n for n in read_input_my.readlines()][:]),
    "".join([n for n in read_input_mz.readlines()][:]),
]

read_input_x.close()
read_input_y.close()
read_input_z.close()
read_input_mx.close()
read_input_my.close()
read_input_mz.close()

coordinates = re.findall(
    r"[A-Z]\S* *x =\s*\S*\s*y =\s*\S*\s*z =\s*\S*", tmp[0], flags=0
)
atoms = len(coordinates)  # Number of atoms
q = np.zeros(
    (atoms + rings, 6)
)  # Create empty charge array (+zero rows for ring conditions)
c = np.zeros((atoms, 3))  # Create empty coordinates array
d = np.zeros((atoms, 3, 6))  # Create empty atomic dipole array
name = [0] * atoms  # Create empty array for atom names
for i in range(atoms):
    name[i] = coordinates[i][
        : int(re.search(r"\S*", coordinates[i]).end())
    ]  # Read atom names from file
    c[i, 0] = float(
        re.findall(r"[ -]\d*[.]\d*", coordinates[i])[0]
    )  # Read x coordinates from file
    c[i, 1] = float(
        re.findall(r"[ -]\d*[.]\d*", coordinates[i])[1]
    )  # Read y coordinates from file
    c[i, 2] = float(
        re.findall(r"[ -]\d*[.]\d*", coordinates[i])[2]
    )  # Read z coordinatas from file
for j in range(6):
    charges = re.findall(r"Q00\s*\S*\s*\S*", tmp[j], flags=0)[:atoms]
    dipoles = re.findall(r"\|Q1\| =.*", tmp[j], flags=0)[:atoms]
    for i in range(atoms):
        q[i, j] = float(re.findall(r"[ -]\d*[.]\d*", charges[i])[0])  # Read charges
        if re.search(
            r"Q11c =\s*[ -]\d*[.]\d*", dipoles[i]
        ):  # Look for x, y and z output (if zero, there is not respective output in GMDA!)
            d[i, 0, j] = float(
                re.findall(
                    r"[ -]\d*[.]\d*",
                    re.findall(r"Q11c =\s*[ -]\d*[.]\d*", dipoles[i])[0],
                )[0]
            )  # atomic dipole x
        if re.search(r"Q11s =\s*[ -]\d*[.]\d*", dipoles[i]):
            d[i, 1, j] = float(
                re.findall(
                    r"[ -]\d*[.]\d*",
                    re.findall(r"Q11s =\s*[ -]\d*[.]\d*", dipoles[i])[0],
                )[0]
            )  # atomic dipole y
        if re.search(r"Q10  =\s*[ -]\d*[.]\d*", dipoles[i]):
            d[i, 2, j] = float(
                re.findall(
                    r"[ -]\d*[.]\d*",
                    re.findall(r"Q10  =\s*[ -]\d*[.]\d*", dipoles[i])[0],
                )[0]
            )  # atomic dipole z
# print("Substract overall charge (" + str(q[:, 0].sum()) + ")?")
# answer=input("............ yes/no   ")
# yes, want to subtract any charge
for j in range(6):
    ch = q[:, j].sum() / atoms
    for i in range(atoms):
        q[i, j] -= ch

print(
    "%4s %7s %7s %7s %7s %7s %7s  %7s %7s %7s "
    % ("Name", "x", "y", "z", "q(Fx)", "q(Fy)", "q(Fz)", "q(-Fx)", "q(-Fy)", "q(-Fz)")
)
for i in range(atoms):
    print(
        "%4s %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f "
        % (
            name[i],
            c[i, 0],
            c[i, 1],
            c[i, 2],
            q[i, 0],
            q[i, 1],
            q[i, 2],
            q[i, 3],
            q[i, 4],
            q[i, 5],
        )
    )
print("")
print(
    "%4s %7s %7s %7s %7s %7s %7s  %7s %7s %7s "
    % (
        "Name",
        "dx(Fx)",
        "dy(Fx)",
        "dz(Fx)",
        "dx(Fy)",
        "dy(Fy)",
        "dz(Fy)",
        "dx(Fz)",
        "dy(Fz)",
        "dz(Fz)",
    )
)
for i in range(atoms):
    print(
        "%4s %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f "
        % (
            name[i],
            d[i, 0, 0],
            d[i, 1, 0],
            d[i, 2, 0],
            d[i, 0, 1],
            d[i, 1, 1],
            d[i, 2, 1],
            d[i, 0, 2],
            d[i, 1, 2],
            d[i, 2, 2],
        )
    )
print("")
print(
    "%4s %7s %7s %7s %7s %7s %7s  %7s %7s %7s "
    % (
        "Name",
        "dx(-Fx)",
        "dy(-Fx)",
        "dz(-Fx)",
        "dx(-Fy)",
        "dy(-Fy)",
        "dz(-Fy)",
        "dx(-Fz)",
        "dy(-Fz)",
        "dz(-Fz)",
    )
)
for i in range(atoms):
    print(
        "%4s %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f "
        % (
            name[i],
            d[i, 0, 3],
            d[i, 1, 3],
            d[i, 2, 3],
            d[i, 0, 4],
            d[i, 1, 4],
            d[i, 2, 4],
            d[i, 0, 5],
            d[i, 1, 5],
            d[i, 2, 5],
        )
    )
print("")


# =============================================================================================
#   CHARGE TRANSFER
# =============================================================================================
# This part calculates the connectivity matrix "n", the matrix with indices
# of each bond charge "index" the bond charges "b" and then the atomic
# dipoles "mu" and the polarizability "a" for the charge transfer
# =============================================================================================
print("")
print("---------------------------------------------------------")
print("------------------Starting calculation-------------------")
print("---------------------------------------------------------")


n = np.zeros((atoms, atoms))  # Create empty connectivity matrix
for i in range(bonds):
    first = int(re.findall(r"\d+", bond[i])[0]) - 1  # Read atom numbers of each bond
    second = int(re.findall(r"\d+", bond[i])[1]) - 1
    if first < second:  # Write entry in antisymmetric connectivity matrix
        n[first, second] = 1
        n[second, first] = -1
    else:
        n[first, second] = -1
        n[second, first] = 1
index = np.zeros(
    (bonds, 2)
)  # Create index matrix (which bond charge contains which atoms)
ctr = 0
for i in range(atoms):
    for j in range(i, atoms):
        if n[i, j] == 1:
            index[ctr, 0] = i
            index[ctr, 1] = j
            ctr += 1
a = np.zeros((atoms + rings, bonds))  # Create matrix for linear equations ab=q
for i in range(atoms):
    for j in range(bonds):
        if i == index[j, 0]:
            a[i, j] = 1
        if i == index[j, 1]:
            a[i, j] = -1
for i in range(rings):  # Write ring conditions into matrix a
    structure = re.findall(r"\d+", ring[i])
    for j in range(len(structure)):
        if j == len(structure) - 1:
            k = 0
        else:
            k = j + 1
        first = int(structure[j]) - 1
        second = int(structure[k]) - 1

        factor = 1
        if second < first:
            save = first
            first = second
            second = save
            factor = -1
        # Find the correct bond charge which is involved in the ring
        element = np.where(np.all(index == np.array([[first, second]]), axis=1))[0][0]
        a[atoms + i, element] = factor  # Write entry in a

# Calculate bond charges and mu
mu = np.zeros((atoms, 3, 6))
for k in range(6):
    b = np.linalg.lstsq(a, q[:, k])[0]  # Solve linear equations ab=q for b

    for i in range(
        atoms
    ):  # Calculate atomic dipole from charge transfer as sum of bond charges times vector of the
        for j in range(bonds):  # atom in direction of each bond
            if i == index[j, 0]:
                mu[i, :, k] += (c[i] - (c[i] + c[int(index[j, 1])]) / 2) * b[j]
            if i == index[j, 1]:
                mu[i, :, k] += (c[i] - (c[i] + c[int(index[j, 0])]) / 2) * b[j] * (-1)

# Calculate polarizability
field = float(input("Enter field strength (in au):   "))
# answer = input(
#     "Do you want to change the default units (au.) of the output to angstrom^3?  yes/no   "
# )
# if answer == "yes":
unit = 0.529177249 ** 3  # transform au^3 into A^3
# else:
#     if answer == "no":
#         unit = 1
#     else:
#         print("---------------------------------------------------------")
#         sys.exit("unknown option")

a_xx = (
    (mu[:, 0, 0] - mu[:, 0, 3]) / (2 * field) * 1.889725989 * (-1) * unit
)  # numeric differentiation, yields polarizability
a_yy = (mu[:, 1, 1] - mu[:, 1, 4]) / (2 * field) * 1.889725989 * (-1) * unit
a_zz = (mu[:, 2, 2] - mu[:, 2, 5]) / (2 * field) * 1.889725989 * (-1) * unit
a_tot_c = (a_xx + a_yy + a_zz) / 3
print("")
print("Results:")
print("")
print("Charge transfer contribution:")
print("%4s %7s %7s %7s %7s" % ("Name", "a_xx", "a_yy", "a_zz", "a_tot"))
for i in range(atoms):
    print(
        "%4s %7.2f %7.2f %7.2f %7.2f "
        % (name[i], a_xx[i], a_yy[i], a_zz[i], a_tot_c[i])
    )
print("Summed up contributions:")
print(a_tot_c[:].sum())

# ========================================================================================
#   Atomic polarizability
# ========================================================================================
# This part calculates the polarizability arising from non-uniform
# distribution of electrons around the core
# ========================================================================================

a_xx = (
    (d[:, 0, 0] - d[:, 0, 3]) / (2 * field) * (-1) * unit
)  # numeric differentiation, yields polarizability
a_yy = (d[:, 1, 1] - d[:, 1, 4]) / (2 * field) * (-1) * unit
a_zz = (d[:, 2, 2] - d[:, 2, 5]) / (2 * field) * (-1) * unit
a_tot_p = (a_xx + a_yy + a_zz) / 3
print("Polarization contribution:")
print("%4s %7s %7s %7s %7s" % ("Name", "a_xx", "a_yy", "a_zz", "a_tot"))
for i in range(atoms):
    print(
        "%4s %7.2f %7.2f %7.2f %7.2f "
        % (name[i], a_xx[i], a_yy[i], a_zz[i], a_tot_p[i])
    )
print("Summed up contributions:")
print(a_tot_p[:].sum())

print("Total polarizability:")
print("%4s %7s" % ("Name", "a"))
for i in range(atoms):
    print("%4s %7.2f " % (name[i], a_tot_p[i] + a_tot_c[i]))
print("Summed up contributions:")
print(a_tot_p[:].sum() + a_tot_c[:].sum())

np.savetxt("atomic_pol.dat", np.c_[a_tot_p, a_tot_c, a_tot_p + a_tot_c])

print("---------------------------------------------------------")
print("-------------------------Goodbye-------------------------")
print("---------------------------------------------------------")
