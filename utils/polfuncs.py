import sys
import MDAnalysis as mda
import MDAnalysis.transformations as trans
import nglview as nv
import numpy as np
import matplotlib.pyplot as plt


def completion(iterable):
    """
    Write percentage completion of a loop to the screen
    """
    total = len(iterable)
    for num, val in enumerate(iterable):
        yield val
        sys.stdout.write(f"\r{num/total*100:.2f} % complete")
    sys.stdout.write("\n")


def plot_ion_counts(counts, ion="CL", save=False, fname=None):
    """
    Plot the ion counts in the all regions as a function of time (assumes ns).
    Pass in residue name of the ion of interest.

    If a filename is given, the plot is saved to that file without having to set save=True.
    If save=True, the plot is saved to a file named {ion.lower()}_counts.png
    """
    ION_NAMES = {"cl": r"Cl$^-$", "na": r"Na$^+$"}
    _ion = ION_NAMES.get(ion.lower(), ion) # if ion is not in ION_NAMES, use resname given
    fig, ax = plt.subplots()
    ax.plot(counts[:, 0], counts[:, 1], label="Saltwater")
    ax.plot(counts[:, 0], counts[:, 2], label="Membrane")
    ax.plot(counts[:, 0], counts[:, 3], label="Freshwater")
    ax.legend()
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(f"{_ion} ion count")
    max_count = counts[:, 1:].max()
    ax.set_yticks(np.arange(0, max_count + 1, 2))
    if save or fname is not None:
        if fname is None:
            fname = f"{ion.lower()}_counts.png"
        fig.savefig(fname, dpi=300)
    else:
        plt.show()


class PolSim:
    """
    Class to read in a desalination-type simulation with a polymer in the centre of the box,
    surrounded by two solvent reservoirs - saltwater on the left of the membrane, freshwater
    on the right.
    """

    def __init__(self, coords, traj):
        self.coordname = coords
        self.trajname = traj
        self.load()

    def load(self):
        self.universe = mda.Universe(self.coordname, self.trajname)
        self.polymer = self.universe.select_atoms("resname pol")
        self.solvent = self.universe.select_atoms("not resname pol")
        self.traj = self.universe.trajectory  # for convenience

    def unwrap(self):
        """
        Unwrap the trajectory so that the polymer chain is centred.
        Be aware that for large trajectories this can take a long time, and
        using gmx trjconv like so will be significantly faster:
        $ echo -e "2\n0" | gmx trjconv -f traj.xtc -s topol.tpr -o traj_unwrapped.xtc -center
        -pbc mol (the polymer chain is the group labelled 2 in this case)
        """
        transforms = [
            trans.unwrap(self.universe.atoms),
            trans.center_in_box(self.polymer),
            trans.wrap(self.solvent, compound="residues"),
        ]
        self.universe.trajectory.add_transformations(*transforms)

    def show_traj(self, hide_virtual_sites=True):
        """
        View the trajectory using NGLView
        """
        _selection = self.universe.atoms
        if hide_virtual_sites:
            _selection = _selection.select_atoms("not name MW")
        view = nv.show_mdanalysis(_selection)
        view.add_unitcell()
        # solvent not loaded in later versions, so do this explicitly
        view.add_representation("ball+stick", "not resname pol")
        return view

    def compute_ion_counts(self, ion="CL", save=True):
        """
        Tracking ion numbers in each portion of the simulation box.
        The polymer is defined by the minimum and maximum z-coordinate of the polymer
        in each frame. Anything to the left is designated as saltwater and any space
        to the right is designated as freshwater.
        Pass in the residue name of the ion of interest.

        Ion counts are saved per frame (given as time in ns) as a numpy ndarray
        to a filename of {ion.lower()}_counts.dat
        """
        counts = np.zeros((self.traj.n_frames, 4))
        for i, ts in enumerate(completion(self.traj)):
            # these may change, so can't really pre-compute bounds...
            # of course, it would be a lot faster if bounds were only computed once
            min_z = self.polymer.positions[:, 2].min()
            max_z = self.polymer.positions[:, 2].max()
            ions_in_saltwater = self.universe.select_atoms(f"resname {ion} and prop z < {min_z}")
            ions_in_membrane = self.universe.select_atoms(
                f"resname {ion} and prop z > {min_z} and prop z < {max_z}"
            )
            ions_in_freshwater = self.universe.select_atoms(f"resname {ion} and prop z > {max_z}")
            counts[i] = [
                ts.time,
                ions_in_saltwater.n_atoms,
                ions_in_membrane.n_atoms,
                ions_in_freshwater.n_atoms,
            ]
        # convert time from ps to ns
        counts[:, 0] /= 1000
        if save:
            np.savetxt(f"{ion.lower()}_counts.dat", counts)
        return counts
