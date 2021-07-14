#!/usr/bin/env python3
import numpy as np
import math
import re
import sys
import argparse
import csv


## Autochem imports ##
def eof(file, percFile):
    # OPEN IN BYTES
    with open(file, "rb") as f:
        f.seek(0, 2)  # Seek @ EOF
        fsize = f.tell()  # Get size
        Dsize = int(percFile * fsize)
        f.seek(max(fsize - Dsize, 0), 0)  # Set pos @ last n chars lines
        lines = f.readlines()  # Read to end

    # RETURN DECODED LINES
    for i in range(len(lines)):
        try:
            lines[i] = lines[i].decode("utf-8")
        except:
            lines[i] = "CORRUPTLINE"
            print("eof function passed a corrupt line in file ", file)
        # FOR LETTER IN SYMBOL
    return lines


def check_user_input(user_input, condition, if_error):
    """
    Uses a try/except statement to create a scenario where the 
    end user cannot give unexpected input. 
    Give the condition referring to an item in a lambda expression 
    i.e. lambda item: item.endswith('.csv'), or lambda item: item in range(...)

    Usage:
        >>> check_user_input('Filename of output', lambda item: item.endswith('.csv'), "Please print a name ending in '.csv'")
        # Produces:
        while not correct:
            try:
                item = input('Filename of output: ')
            except ValueError:
                print("Please enter a filename ending in '.csv'")
            if not filename.endswith('.csv'):
                print("Please enter a filename ending in '.csv'")
            else:
                correct = True
        return item
    """
    f = condition
    correct = False
    while not correct:
        try:
            item = input(user_input + ": ")
        except ValueError:
            print(if_error)
        if not f(item):
            print(if_error)
        else:
            correct = True
    return item


def sort_elements(lst):
    """
    Sort a list of |Atom| objects by atomic number. 
    Returns a list of tuples- [(symbol, atomic number), (symbol, atomic number), ...]

    TODO: Extend to giving back the objects- more useful than just for formatting of symbols
    """
    els = []
    elements = set([atom.symbol for atom in lst])
    for i in elements:
        atom = Atom(i)
        els.append((i, float(PT.get_atnum(atom))))
    sorted_els = sorted(els, key=lambda val: val[1])
    return sorted_els


def write_csv_from_dict(data, filename=None, autosave=False):
    """Write to file from dictionary"""
    write = True if autosave else False
    if not autosave:
        done = False
        while not done:
            to_file = input("Print to csv? [Y/N] ")
            if to_file.lower() in ("y", "n"):
                done = True
                if to_file.lower() == "y":
                    write = True
                    if filename is None:
                        filename = check_user_input(
                            "Filename",
                            lambda item: item.endswith(".csv"),
                            "Please give a filename ending in '.csv'",
                        )
            else:
                print("Please select 'Y' or 'N'")
    if write:
        with open(filename, "w", encoding="utf-8-sig") as f:
            writer = csv.writer(f)
            writer.writerow(data.keys())
            content = zip(*[data[key] for key in data.keys()])
            writer.writerows(content)


def responsive_table(data, strings, min_width=13, decimal_places=5):
    """
    Returns a table that is responsive in size to every column.
    Requires a dictionary to be passed in, with the keys referring to
    the headers of the table.
    Also pass in the number of each column that should be a string, starting
    from 1.

    Usage:
        >>> d = {'col1': [1,2,3,4],
                 'col2': ['One', 'Two', 'Three', 'Four']}
        >>> responsive_table(d, strings = [2])

    Can also give a minimum width, defaults to 13 spaces. A `decimal_places` parameter
    can be passed in to define the number of decimal places of floats.
    """
    num_cols = len(data.keys())
    content = zip(*[data[key] for key in data.keys()])  # dict values into list of lists
    # unknown number of arguments
    max_sizes = {}
    try:
        for k, v in data.items():
            max_sizes[k] = len(max([str(val) for val in v], key=len))
    except ValueError:
        sys.exit("Error: No data is passed into autochem.core.utils.responsive_table")

    # create the thing to pass into .format()- can't have brackets like zip gives
    formatting = []
    index = 0
    all_sizes = []
    for val in zip(data.keys(), max_sizes.values()):
        entry, size = val
        if size < min_width or index + 1 not in strings:
            size = min_width
        # also check dict key length
        if len(entry) > size:
            size = len(entry)
        formatting += [entry, size]
        all_sizes.append(size)
        index += 1
    line_length = sum(all_sizes) + num_cols * 3 - 1  # spaces in header
    print("+" + "-" * line_length + "+")
    output_string = "|" + " {:^{}} |" * len(data.keys())
    print(output_string.format(*formatting))
    print("+" + "-" * line_length + "+")
    for line in content:
        formatting = []
        for val in zip(line, all_sizes):
            entry, size = val
            # if not isinstance(entry, str):
            if isinstance(entry, float):
                size = f"{size}.{decimal_places}f"
            formatting.append(entry)
            formatting.append(size)
        print(output_string.format(*formatting))
    print("+" + "-" * line_length + "+")


class PT:
    """Helper class to allow for lookup of atomic properties. Can convert between symbol and atomic number"""

    ptable = {}
    # [symbol, mass, radius, connectors, vdw radii]
    # atomic weights from: http://www.ciaaw.org/atomic-weights.htm
    # need to add more vdw radii from molden oglmol file (vdwr array)
    ptable[0] = ["Xx", 0.00000, 0.00, 0, 0.000]
    ptable[1] = ["H", 1.00798, 0.30, 1, 0.430]
    ptable[2] = ["He", 4.00260, 0.99, 0, 0.741]
    ptable[3] = ["Li", 6.96750, 1.52, 8, 0.880]
    ptable[4] = ["Be", 9.01218, 1.12, 8, 0.550]
    ptable[5] = ["B", 10.81350, 0.88, 6, 1.030]
    ptable[6] = ["C", 12.01060, 0.77, 4, 0.900]
    ptable[7] = ["N", 14.00685, 0.70, 3, 0.880]
    ptable[8] = ["O", 15.99940, 0.66, 2, 0.880]
    ptable[9] = ["F", 18.99840, 0.64, 1, 0.840]
    ptable[10] = ["Ne", 20.17970, 1.60, 0, 0.815]
    ptable[11] = ["Na", 22.98977, 1.86, 8, 1.170]
    ptable[12] = ["Mg", 24.30550, 1.60, 8, 1.300]
    ptable[13] = ["Al", 26.98154, 1.43, 8, 1.550]
    ptable[14] = ["Si", 28.08500, 1.17, 8, 1.400]
    ptable[15] = ["P", 30.97376, 1.10, 8, 1.250]
    ptable[16] = ["S", 32.06750, 1.04, 2, 1.220]
    ptable[17] = ["Cl", 35.45150, 0.99, 1, 1.190]
    ptable[18] = ["Ar", 39.94800, 1.92, 0, 0.995]
    ptable[19] = ["K", 39.09830, 2.31, 8, 1.530]
    ptable[20] = ["Ca", 40.07800, 1.97, 8, 1.190]
    ptable[21] = ["Sc", 44.95591, 1.60, 8, 1.640]
    ptable[22] = ["Ti", 47.86700, 1.46, 8, 1.670]
    ptable[23] = ["V", 50.94150, 1.31, 8, 1.530]
    ptable[24] = ["Cr", 51.99610, 1.25, 8, 1.550]
    ptable[25] = ["Mn", 54.93804, 1.29, 8, 1.555]
    ptable[26] = ["Fe", 55.84500, 1.26, 8, 1.540]
    ptable[27] = ["Co", 58.93319, 1.25, 8, 1.530]
    ptable[28] = ["Ni", 58.69340, 1.24, 8, 1.700]
    ptable[29] = ["Cu", 63.54600, 1.28, 8, 1.720]
    ptable[30] = ["Zn", 65.38000, 1.33, 8, 1.650]
    ptable[31] = ["Ga", 69.72300, 1.41, 8, 1.420]
    ptable[32] = ["Ge", 72.63000, 1.22, 8, 1.370]
    ptable[33] = ["As", 74.92159, 1.21, 8, 1.410]
    ptable[34] = ["Se", 78.97100, 1.17, 8, 1.420]
    ptable[35] = ["Br", 79.90400, 1.14, 1, 1.410]
    ptable[36] = ["Kr", 83.79800, 1.97, 0, 1.069]
    ptable[37] = ["Rb", 85.46780, 2.44, 8, 1.670]
    ptable[38] = ["Sr", 87.62000, 2.15, 8, 1.320]
    ptable[39] = ["Y", 88.90584, 1.80, 8, 1.980]
    ptable[40] = ["Zr", 91.22400, 1.57, 8, 1.760]
    ptable[41] = ["Nb", 92.90637, 1.41, 8, 1.680]
    ptable[42] = ["Mo", 95.95000, 1.36, 8, 1.670]
    ptable[43] = ["Tc", 98.00000, 1.35, 8, 1.550]
    ptable[44] = ["Ru", 101.07000, 1.33, 8, 1.600]
    ptable[45] = ["Rh", 102.90550, 1.34, 8, 1.650]
    ptable[46] = ["Pd", 106.42000, 1.38, 8, 1.700]
    ptable[47] = ["Ag", 107.86820, 1.44, 8, 1.790]
    ptable[48] = ["Cd", 112.41400, 1.49, 8, 1.890]
    ptable[49] = ["In", 114.81800, 1.66, 8, 1.830]
    ptable[50] = ["Sn", 118.71000, 1.62, 8, 1.660]
    ptable[51] = ["Sb", 121.76000, 1.41, 8, 0.000]
    ptable[52] = ["Te", 127.60000, 1.37, 8, 0.000]
    ptable[53] = ["I", 126.90447, 1.33, 1, 0.000]
    ptable[54] = ["Xe", 131.29300, 2.17, 0, 0.000]
    ptable[55] = ["Cs", 132.90545, 2.62, 8, 0.000]
    ptable[56] = ["Ba", 137.32700, 2.17, 8, 0.000]
    ptable[57] = ["La", 138.90547, 1.88, 8, 0.000]
    ptable[58] = ["Ce", 140.11600, 1.818, 8, 0.000]
    ptable[59] = ["Pr", 140.90766, 1.824, 8, 0.000]
    ptable[60] = ["Nd", 144.24200, 1.814, 8, 0.000]
    ptable[61] = ["Pm", 145.00000, 1.834, 8, 0.000]
    ptable[62] = ["Sm", 150.36000, 1.804, 8, 0.000]
    ptable[63] = ["Eu", 151.96400, 2.084, 8, 0.000]
    ptable[64] = ["Gd", 157.25000, 1.804, 8, 0.000]
    ptable[65] = ["Tb", 158.92535, 1.773, 8, 0.000]
    ptable[66] = ["Dy", 162.50000, 1.781, 8, 0.000]
    ptable[67] = ["Ho", 164.93033, 1.762, 8, 0.000]
    ptable[68] = ["Er", 167.25900, 1.761, 8, 0.000]
    ptable[69] = ["Tm", 168.93422, 1.759, 8, 0.000]
    ptable[70] = ["Yb", 173.04500, 1.922, 8, 0.000]
    ptable[71] = ["Lu", 174.96680, 1.738, 8, 0.000]
    ptable[72] = ["Hf", 178.49000, 1.57, 8, 0.000]
    ptable[73] = ["Ta", 180.94788, 1.43, 8, 0.000]
    ptable[74] = ["W", 183.84000, 1.37, 8, 0.000]
    ptable[75] = ["Re", 186.20700, 1.37, 8, 0.000]
    ptable[76] = ["Os", 190.23000, 1.34, 8, 0.000]
    ptable[77] = ["Ir", 192.21700, 1.35, 8, 0.000]
    ptable[78] = ["Pt", 195.08400, 1.38, 8, 0.000]
    ptable[79] = ["Au", 196.96657, 1.44, 8, 0.000]
    ptable[80] = ["Hg", 200.59200, 1.52, 8, 0.000]
    ptable[81] = ["Tl", 204.38350, 1.71, 8, 0.000]
    ptable[82] = ["Pb", 207.20000, 1.75, 8, 0.000]
    ptable[83] = ["Bi", 208.98040, 1.70, 8, 0.000]
    ptable[84] = ["Po", 209.00000, 1.40, 8, 0.000]
    ptable[85] = ["At", 210.00000, 1.40, 1, 0.000]
    ptable[86] = ["Rn", 222.00000, 2.40, 0, 0.000]
    ptable[87] = ["Fr", 223.00000, 2.70, 8, 0.000]
    ptable[88] = ["Ra", 226.00000, 2.20, 8, 0.000]
    ptable[89] = ["Ac", 227.00000, 2.00, 8, 0.000]
    ptable[90] = ["Th", 232.03770, 1.79, 8, 0.000]
    ptable[91] = ["Pa", 231.03588, 1.63, 8, 0.000]
    ptable[92] = ["U", 238.02891, 1.56, 8, 0.000]
    ptable[93] = ["Np", 237.00000, 1.55, 8, 0.000]
    ptable[94] = ["Pu", 244.00000, 1.59, 8, 0.000]
    ptable[95] = ["Am", 243.00000, 1.73, 8, 0.000]
    ptable[96] = ["Cm", 247.00000, 1.74, 8, 0.000]
    ptable[97] = ["Bk", 247.00000, 1.70, 8, 0.000]
    ptable[98] = ["Cf", 251.00000, 1.86, 8, 0.000]
    ptable[99] = ["Es", 252.00000, 1.86, 8, 0.000]
    ptable[100] = ["Fm", 257.00000, 2.00, 8, 0.000]
    ptable[101] = ["Md", 258.00000, 2.00, 8, 0.000]
    ptable[102] = ["No", 259.00000, 2.00, 8, 0.000]
    ptable[103] = ["Lr", 266.00000, 2.00, 8, 0.000]
    ptable[104] = ["Rf", 267.00000, 2.00, 8, 0.000]
    ptable[105] = ["Db", 268.00000, 2.00, 8, 0.000]
    ptable[106] = ["Sg", 269.00000, 2.00, 8, 0.000]
    ptable[107] = ["Bh", 270.00000, 2.00, 8, 0.000]
    ptable[108] = ["Hs", 277.00000, 2.00, 8, 0.000]
    ptable[109] = ["Mt", 278.00000, 2.00, 8, 0.000]
    ptable[110] = ["Ds", 281.00000, 2.00, 8, 0.000]
    ptable[111] = ["Rg", 282.00000, 2.00, 8, 0.000]
    ptable[112] = ["Cn", 285.00000, 2.00, 8, 0.000]
    ptable[113] = ["Nh", 286.00000, 2.00, 8, 0.000]
    ptable[114] = ["Fl", 289.00000, 2.00, 8, 0.000]
    ptable[115] = ["Mc", 290.00000, 2.00, 8, 0.000]
    ptable[116] = ["Lv", 293.00000, 2.00, 8, 0.000]
    ptable[117] = ["Ts", 294.00000, 2.00, 8, 0.000]
    ptable[118] = ["Og", 294.00000, 2.00, 8, 0.000]

    def __init__(self):
        raise AttributeError("The PeriodicTable class cannot be instantiated.")

    @classmethod
    def get_atnum(cls, atom):
        """Converts symbol to atomic number"""
        for key in cls.ptable.keys():
            if cls.ptable[key][0] == atom.symbol.capitalize():
                return key

    @classmethod
    def get_symbol(cls, atom):
        """Converts atomic number to symbol"""
        return cls.ptable[atom.atnum][0]

    @classmethod
    def get_radius(cls, atom):
        """Returns atomic radius for a given element"""
        return float(cls.ptable[atom.atnum][2])

    @classmethod
    def get_mass(cls, atom):
        """Returns atomic mass for a given element"""
        return float(cls.ptable[atom.atnum][1])

    @classmethod
    def get_connectors(cls, atom):
        """Returns number of possible attachments to a given element"""
        return int(cls.ptable[atom.atnum][-2])

    @classmethod
    def get_vdw(cls, atom):
        """Returns van der waals radius of a given element"""
        return float(cls.ptable[atom.atnum][-1])


class Atom:
    """A class representing a atom in 3 dimensional euclidean space.
    An instance has the following attributes:

    * ``atnum`` -- atomic symbol, equal to zero for a dummy atom
    * ``coords`` -- tuple of x,y,z coordinates
    * ``bonds`` -- list of bonds that this atom is a part of
    * ``mol`` -- molecule this atom is a part of. Assigned programmatically when a molecule is separated using the *mol.separate* method, or can be assigned manually if building up a molecule from scratch

    Access these properties directly:
    * ``x``, ``y``, ``z`` -- for atom coordinates
    * ``symbol`` -- read or write the atom symbol directly

    >>> a = Atom('H', coords = (1,2,3))

    """

    def __init__(self, symbol=None, atnum=0, coords=None, mol=None, bonds=None):
        if symbol is not None:
            self.symbol = symbol
            self.atnum = PT.get_atnum(self)
        else:
            self.atnum = atnum
            self.symbol = PT.get_symbol(self)

        self.mass = PT.get_mass(self)
        self.mol = mol
        self.bonds = bonds or []
        self.connected_atoms = []
        self.h_bonded_to = []
        self.fragment = None

        if coords is None:
            self.coords = [0, 0, 0]
        elif len(coords) == 3:
            self.coords = [float(i) for i in coords]
        else:
            raise TypeError("Atom: Invalid coordinates given")

    @property
    def x(self):
        return self.coords[0]

    @property
    def y(self):
        return self.coords[1]

    @property
    def z(self):
        return self.coords[2]

    def __repr__(self):
        """Unambiguous representation of an |Atom| instance"""
        if self.fragment is not None:
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Index: {self.index} Mol: {self.fragment}"

        if hasattr(self, "index") and not hasattr(self, "mol"):
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Index: {self.index} "

        if hasattr(self, "index") and hasattr(self, "mol"):
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Index: {self.index} Mol: {self.mol}"

        elif hasattr(self, "number"):
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Mol: {self.mol} Atom: {self.number}"

        elif hasattr(self, "index") and len(self.h_bonded_to) > 0:
            h_bonded = [
                (atom.symbol, {"mol": atom.mol, "atom": atom.index})
                for atom in self.h_bonded_to
            ]
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Mol: {self.mol} Index: {self.index} Number: {self.number} H-Bonds: {h_bonded}"

        return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}"

    def __iter__(self):
        """Iterates through coordinates when called"""
        return iter(self.coords)

    def translate(self, vector):
        """Move atom in space by passing a vector in angstroms"""
        self.coords = [i + j for i, j in zip(self, vector)]

    def move_to(self, vector):
        """Move atom in space to the values, in angstroms, given in this vector. The vector passed represents a point in euclidean space"""
        self.coords = [i for i in vector]

    def distance_to(self, vector):
        """Measure the distance between the atom and a point in space, given as a vector in angstroms"""
        # pythagoras in 3D
        return sum((i - j) ** 2 for i, j in zip(self, vector)) ** 0.5

    def vector_to(self, point):
        """Returns a vector from the atom to a given point, in angstroms"""
        return [(i - j) for i, j in zip(point, self)]

    def angle_between(self, pos1, pos2):
        """Returns an angle between positions 1 and 2 in degrees, with this atom lying at the centre"""
        # dot product, angle = cos^-1([vec(a).vec(b)] / [dist(a) * dist(b)])
        num = np.dot(self.vector_to(pos1), self.vector_to(pos2))
        denom = self.distance_to(pos1) * self.distance_to(pos2)
        return math.acos(num / denom) * (180 / math.pi)

    def as_xyz(self, end_of_line=""):
        """
        Return atom in xyz format: symbol x y z. Can also give an optional 
        end of line character such as a newline
        """
        return f"{self.symbol:5s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}{end_of_line}"

    __str__ = __repr__


## Start of actual script ##
def get_dispersion_and_induction(log):
    """
    Returns SAPT dispersion and induction values from a given log file
    """
    disp = 0
    ind = 0
    found = False
    with open(log) as f:
        for line in f.readlines():
            if "SAPT Results" in line:
                found = True
                continue
            if found:
                if "Dispersion" in line and not "sSAPT0" in line:
                    disp = float(line.split()[-2])
                if "Induction" in line and not "sSAPT0" in line:
                    ind = float(line.split()[-2])
    return disp, ind


def kij(disp, ind):
    return 1 / (1 + (ind / disp))


class Molecule:
    def __init__(self, atomlist):
        self.coords = atomlist
        self.split_into_molecules()

    def distance_matrix(self):
        """
        Creates an N x N matrix of interatomic distances
        between every atom in the system. N = number of 
        atoms in system.
        """
        num_atoms = len(self.coords)
        matrix = np.zeros((num_atoms, num_atoms))

        for i, atom_i in enumerate(self.coords):
            for j, atom_j in enumerate(self.coords):
                matrix[i, j] = atom_i.distance_to(atom_j)
        return matrix

    def split_into_molecules(self):
        """
        Split a system into fragments using van der waals radii. Modifies
        attributes of atoms in self.coords directly in loop, instead of creating
        a dictionary and appending to the dictionary as we go. This
        significantly speeds up the fragmentation.
        """
        mol_count = 0
        dists = self.distance_matrix()
        for i, atom_i in enumerate(self.coords):
            connected = False
            for j, atom_j in enumerate(self.coords):
                if i != j:
                    vdw_dist = PT.get_vdw(atom_i) + PT.get_vdw(atom_j)
                    if dists[i, j] < vdw_dist:
                        # connected
                        connected = True
                        if atom_i not in atom_j.connected_atoms:
                            atom_j.connected_atoms.append(atom_i)
                        if atom_j not in atom_i.connected_atoms:
                            atom_i.connected_atoms.append(atom_j)
                        if atom_i.mol is None and atom_j.mol is None:
                            atom_i.mol = mol_count
                            atom_j.mol = mol_count
                            mol_count += 1
                        elif atom_i.mol is None and atom_j.mol is not None:
                            atom_i.mol = atom_j.mol
                        elif atom_j.mol is None and atom_i.mol is not None:
                            atom_j.mol = atom_i.mol
                        # if different assignments, remove original assignment
                        # combine the two fragments together, as they are connected
                        elif atom_i.mol is not None and atom_j.mol is not None:
                            if atom_i.mol != atom_j.mol:
                                orig = atom_j.mol
                                for atom in self.coords:
                                    if atom.mol is orig:
                                        atom.mol = atom_i.mol
            if not connected:
                atom_i.mol = mol_count
                mol_count += 1

        nums = set([atom.mol for atom in self.coords])
        self.molecules = [
            [atom for atom in self.coords if atom.mol == val] for val in nums
        ]


class com_distance:
    """
    Assumes that the xyz only contains two fragments
    """

    def __init__(self, molecule):
        self.mol = molecule.molecules

    def centre_of_mass(self, atom_list):
        """
        1/M * sum(m_i * r_i)
        """
        total_mass = sum(PT.get_mass(atom) for atom in atom_list)
        com_x = 1 / total_mass * sum(PT.get_mass(atom) * atom.x for atom in atom_list)
        com_y = 1 / total_mass * sum(PT.get_mass(atom) * atom.y for atom in atom_list)
        com_z = 1 / total_mass * sum(PT.get_mass(atom) * atom.z for atom in atom_list)
        return com_x, com_y, com_z

    @property
    def distance_between_centre_of_masses(self):
        coms = [self.centre_of_mass(frag) for frag in self.mol]
        return (
            (coms[1][0] - coms[0][0]) ** 2
            + (coms[1][1] - coms[0][1]) ** 2
            + (coms[1][2] - coms[0][2]) ** 2
        ) ** 0.5


class calc:
    def __init__(self, log):
        self.log = log
        self.create_molecule()
        self.check_if_completed()

    def check_if_completed(self):
        lastline = eof(self.log, 0.02)[-1]
        if "Psi4 exiting successfully. Buy a developer a beer!" not in lastline:
            self.completed = False
        self.completed = True

    def create_molecule(self):
        atomlist = []
        with open(self.log) as f:
            found = False
            for line in f:
                if (
                    "Center              X                  Y                   Z"
                    in line
                ):
                    found = True
                    continue
                if found and re.search("^\s+$", line):
                    break
                if found and "---" not in line:
                    line = line.split()
                    sym = line[0]
                    coords = [float(i) for i in line[1:4]]
                    atomlist.append(Atom(sym, coords=coords))
        self.mol = Molecule(atomlist)

    @property
    def centre_of_mass_separation(self):
        com = com_distance(self.mol)
        return com.distance_between_centre_of_masses

    @property
    def k_ij(self):
        disp, ind = get_dispersion_and_induction(self.log)
        return kij(disp, ind)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", help="SAPT2+ log files to read.", nargs="+")
    parser.add_argument(
        "-o",
        "--output",
        help="filename of csv to write to.",
        default="sapt_analysis.csv",
    )
    args = parser.parse_args()

    datadict = {"Name": [], "k_ij": [], "r_COM": []}
    for log in args.files:
        results = calc(log)
        datadict["Name"].append(results.log)
        datadict["k_ij"].append(results.k_ij)
        datadict["r_COM"].append(results.centre_of_mass_separation)
    responsive_table(datadict, strings=[1])
    write_csv_from_dict(datadict, filename=args.output)


if __name__ == "__main__":
    main()
