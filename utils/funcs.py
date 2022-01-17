import numpy as np
import seaborn as sns
import re


def datafile_elements(datafile):
    """
    Extracting elements from
    1      12.0100    # C
    2      12.0100    # C3
    3       1.0080    # HC
    4       1.0080    # H1
    5       1.0080    # HO
    6       1.0080    # HN
    """
    masses = []
    found = False
    with open(datafile) as f:
        for line in f:
            if "Masses" in line:
                found = True
                continue
            if len(masses) > 0 and re.search("^\s*[A-Z]", line):
                break
            if found and len(line.split()) > 1:  # element or #element
                masses.append(float(line.split()[1]))
    ptable = {}
    # atomic weights from: http://www.ciaaw.org/atomic-weights.htm
    ptable["Xx"] = 0.00000
    ptable["H"] = 1.00798
    ptable["He"] = 4.00260
    ptable["Li"] = 6.96750
    ptable["Be"] = 9.01218
    ptable["B"] = 10.81350
    ptable["C"] = 12.01060
    ptable["N"] = 14.00685
    ptable["O"] = 15.99940
    ptable["F"] = 18.99840
    ptable["Ne"] = 20.17970
    ptable["Na"] = 22.98977
    ptable["Mg"] = 24.30550
    ptable["Al"] = 26.98154
    ptable["Si"] = 28.08500
    ptable["P"] = 30.97376
    ptable["S"] = 32.06750
    ptable["Cl"] = 35.45150
    ptable["Ar"] = 39.94800
    ptable["K"] = 39.09830
    ptable["Ca"] = 40.07800
    ptable["Sc"] = 44.95591
    ptable["Ti"] = 47.86700
    ptable["V"] = 50.94150
    ptable["Cr"] = 51.99610
    ptable["Mn"] = 54.93804
    ptable["Fe"] = 55.84500
    ptable["Co"] = 58.93319
    ptable["Ni"] = 58.69340
    ptable["Cu"] = 63.54600
    ptable["Zn"] = 65.38000
    ptable["Ga"] = 69.72300
    ptable["Ge"] = 72.63000
    ptable["As"] = 74.92159
    ptable["Se"] = 78.97100
    ptable["Br"] = 79.90400
    ptable["Kr"] = 83.79800
    ptable["Rb"] = 85.46780
    ptable["Sr"] = 87.62000
    ptable["Y"] = 88.90584
    ptable["Zr"] = 91.22400
    ptable["Nb"] = 92.90637
    ptable["Mo"] = 95.95000
    ptable["Tc"] = 98.00000
    ptable["Ru"] = 101.07000
    ptable["Rh"] = 102.90550
    ptable["Pd"] = 106.42000
    ptable["Ag"] = 107.86820
    ptable["Cd"] = 112.41400
    ptable["In"] = 114.81800
    ptable["Sn"] = 118.71000
    ptable["Sb"] = 121.76000
    ptable["Te"] = 127.60000
    ptable["I"] = 126.90447
    ptable["Xe"] = 131.29300
    ptable["Cs"] = 132.90545
    ptable["Ba"] = 137.32700
    ptable["La"] = 138.90547
    ptable["Ce"] = 140.11600
    ptable["Pr"] = 140.90766
    ptable["Nd"] = 144.24200
    ptable["Pm"] = 145.00000
    ptable["Sm"] = 150.36000
    ptable["Eu"] = 151.96400
    ptable["Gd"] = 157.25000
    ptable["Tb"] = 158.92535
    ptable["Dy"] = 162.50000
    ptable["Ho"] = 164.93033
    ptable["Er"] = 167.25900
    ptable["Tm"] = 168.93422
    ptable["Yb"] = 173.04500
    ptable["Lu"] = 174.96680
    ptable["Hf"] = 178.49000
    ptable["Ta"] = 180.94788
    ptable["W"] = 183.84000
    ptable["Re"] = 186.20700
    ptable["Os"] = 190.23000
    ptable["Ir"] = 192.21700
    ptable["Pt"] = 195.08400
    ptable["Au"] = 196.96657
    ptable["Hg"] = 200.59200
    ptable["Tl"] = 204.38350
    ptable["Pb"] = 207.20000
    ptable["Bi"] = 208.98040
    ptable["Po"] = 209.00000
    ptable["At"] = 210.00000
    ptable["Rn"] = 222.00000
    ptable["Fr"] = 223.00000
    ptable["Ra"] = 226.00000
    ptable["Ac"] = 227.00000
    ptable["Th"] = 232.03770
    ptable["Pa"] = 231.03588
    ptable["U"] = 238.02891
    ptable["Np"] = 237.00000
    ptable["Pu"] = 244.00000
    ptable["Am"] = 243.00000
    ptable["Cm"] = 247.00000
    ptable["Bk"] = 247.00000
    ptable["Cf"] = 251.00000
    ptable["Es"] = 252.00000
    ptable["Fm"] = 257.00000
    ptable["Md"] = 258.00000
    ptable["No"] = 259.00000
    ptable["Lr"] = 266.00000
    ptable["Rf"] = 267.00000
    ptable["Db"] = 268.00000
    ptable["Sg"] = 269.00000
    ptable["Bh"] = 270.00000
    ptable["Hs"] = 277.00000
    ptable["Mt"] = 278.00000
    ptable["Ds"] = 281.00000
    ptable["Rg"] = 282.00000
    ptable["Cn"] = 285.00000
    ptable["Nh"] = 286.00000
    ptable["Fl"] = 289.00000
    ptable["Mc"] = 290.00000
    ptable["Lv"] = 293.00000
    ptable["Ts"] = 294.00000
    ptable["Og"] = 294.00000

    elements = []
    for mass in masses:
        for el, reference in ptable.items():
            # int rounds down, so int(15.999) = 15 not 16, fix by adding 0.5
            if int(mass + 0.5) == int(reference + 0.5):
                elements.append(el)
    return elements


def lammps_mdanalysis(data, traj):
    """
    Import LAMMPS datafile and trajectory,
    adding atom names by looking up the mass from the datafile
    """
    import MDAnalysis as mda

    # deal with incorrect filenames
    format_needed = any(traj.endswith(x) for x in ("lmp", "lammpstrj"))
    top_needed = not data.endswith("data")
    if format_needed and top_needed:
        u = mda.Universe(data, traj, topology_format="DATA", format="LAMMPSDUMP")
    elif format_needed:
        u = mda.Universe(data, traj, format="LAMMPSDUMP")
    elif top_needed:
        u = mda.Universe(data, traj, topology_format="DATA", format="LAMMPSDUMP")
    else:
        u = mda.Universe(data, traj)

    u.add_TopologyAttr("names")
    atoms = datafile_elements(data)

    types = sorted([int(i) for i in np.unique(u.atoms.types).tolist()])

    for a, t in zip(atoms, types):
        u.select_atoms(f"type {t}").names = a
    # convert to ints
    u.atoms.types = u.atoms.types.astype(int)

    # for hydrogen bonding analysis to work, need to assign residues (meaningless for
    # these systems)
    u.add_TopologyAttr("resnames")
    u.atoms.residues.resnames = "UNK"  # set to unknown
    return u


def add_elements_to_universe(universe, datafile="pack.lmps"):
    universe.add_TopologyAttr("names")
    atoms = datafile_elements(datafile)

    types = sorted([int(i) for i in np.unique(universe.atoms.types).tolist()])

    for a, t in zip(atoms, types):
        universe.select_atoms(f"type {t}").names = a
    # convert to ints
    universe.atoms.types = universe.atoms.types.astype(int)
    # for hydrogen bonding analysis to work, need to assign residues (meaningless for
    # these systems)
    universe.add_TopologyAttr("resnames")
    universe.atoms.residues.resnames = "UNK"  # set to unknown
    return universe


def unwrap_traj(universe):
    """
    Unwrap trajectories of all frames in Universe.
    Apply to universe defined in the global scope.
    i.e.
    unwrap_traj(u)
    not
    u = unwrap_traj(u)
    """
    import MDAnalysis as mda
    from MDAnalysis import transformations  # used as mda.transformations

    a = universe.atoms
    transform = mda.transformations.unwrap(a)
    universe.trajectory.add_transformations(transform)


def remove_legend_title(plot):
    """
    Removes legend title from any seaborn plot
    """
    if isinstance(plot, sns.FacetGrid):
        plot.legend.set_title(None)
    else:
        plot.legend().set_title(None)


def boltz(diffs):
    """
    Pass in energy differences in kJ/mol
    """
    R = 8.3145
    T = 298.15
    kJ_to_J = 1000
    numerator = np.exp((-1 * kJ_to_J * diffs) / (R * T))
    return numerator / numerator.sum()


def completion(iterable):
    """
    Write percentage completion of a loop to the screen
    """
    import sys

    total = len(iterable)
    for num, val in enumerate(iterable):
        yield val
        sys.stdout.write(f"\r{num/total*100:.2f} % complete")
    sys.stdout.write("\n")


def pipe(original):
    """
    Use as a decorator on any function.
    For example:
    @pipe
    def func(*args, **kwargs):
        ...

    obj >> func()
    """

    class PipeInto(object):
        data = {"function": original}

        def __init__(self, *args, **kwargs):
            self.data["args"] = args
            self.data["kwargs"] = kwargs

        def __rrshift__(self, other):
            return self.data["function"](
                other, *self.data["args"], **self.data["kwargs"]
            )

    return PipeInto
