import MDAnalysis as mda
import MDAnalysis.transformations as trans
from MDAnalysis.analysis.lineardensity import LinearDensity
import matplotlib.pyplot as plt

try:
    import nglview as nv
except ImportError:
    nv = None
import numpy as np
import pandas as pd
import seaborn as sns
import subprocess as sp
import sys

PS_TO_NS = 1e-3


def completion(iterable):
    """
    Write percentage completion of a loop to the screen
    """
    total = len(iterable)
    for num, val in enumerate(iterable):
        yield val
        sys.stdout.write(f"\r{num/total*100:.2f} % complete")
    sys.stdout.write("\n")


def plot_ion_counts(counts, fname=None):
    """
    Plot the ion counts in the all regions as a function of time (assumes ns).

    If a filename is given, the plot is saved.
    """
    ION_NAMES = {"cl": r"Cl$^-$", "na": r"Na$^+$"}
    REGIONS = {"salt": "Saltwater", "fresh": "Freshwater", "mem": "Membrane"}
    tidy = counts.melt(id_vars="time_ns", var_name="ion_region", value_name="count")
    tidy[["ion", "region"]] = tidy["ion_region"].str.split("_", expand=True)
    tidy["ion"] = tidy["ion"].map(ION_NAMES)
    tidy["region"] = tidy["region"].map(REGIONS)
    sns.set(
        style="ticks",
        font_scale=1,
        font="DejaVu Sans",
        rc={"mathtext.default": "regular", "figure.figsize": (8, 6)},
    )
    p = sns.relplot(
        data=tidy,
        x="time_ns",
        y="count",
        hue="region",
        col="ion",
        col_order=[r"Na$^+$", r"Cl$^-$"],
        kind="line",
        height=3,
        aspect=1,
        palette="Set2",
        facet_kws={"sharey": False},
    )
    p.set(xlabel="Time (ns)", ylabel="Ion Count")
    sns.move_legend(p, loc="lower center", ncol=3, title=None, borderaxespad=0)
    # no direction param, so control orientation with number of columns -
    # reducing the border-axes padding gives more space around the legend
    p.set_titles(col_template="{col_name}")
    plt.tight_layout()

    if fname is not None:
        plt.savefig(fname, dpi=300)
    else:
        plt.show()


def plot_partial_densities(df, fname=None):
    """
    Plot the partial densities of molecules across the simulation box.
    If a filename is given, the plot is saved.
    Using PolSim.compute_partial_densities() returns a dataframe with the
    first column as slices in the desired dimension.
    This function uses that and obtains the dimension from the name of that column.
    """
    sns.set(
        style="ticks",
        font_scale=1,
        font="DejaVu Sans",
        rc={"mathtext.default": "regular"},
        palette="Set2",
    )
    NAMES = {"pol": "Polymer", "SOL": r"H$_2$O", "NA": r"Na$^+$", "CL": r"Cl$^-$"}
    dim = df.columns[0].split("_")[0]
    df.set_index(f"{dim}_coord").rename(columns=NAMES).plot().set(
        xlabel=f"{dim.upper()}-coordinate (Å)", ylabel="Density (g/cm$^3$)"
    )
    plt.tight_layout()
    if fname is not None:
        plt.savefig(fname, dpi=300)
    else:
        plt.show()


class PolSim:
    """
    Class to read in a desalination-type simulation with a polymer in the centre of the box,
    surrounded by two solvent reservoirs - saltwater on the left of the membrane, freshwater
    on the right.

    Expects Gromacs files - tpr files as coordinates, xtc files as trajectories.
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
        self.dt = self.traj[1].time - self.traj[0].time  # picoseconds

    def unwrap(self, gmx=False):
        """
        Note: only works with tpr files, not gro files.
        Unwrap the trajectory, centre the polymer and then wrap molecules so that the centre of mass
        lies in the periodic box.
        For large trajectories, MDAnalysis' frame-by-frame unwrapping and analysis takes a long time
        so gmx trjconv can be used, and the universe is re-created using the unwrapped trajectory.
        Note that gmx trjconv is used for any trajectory with over 10000 frames.
        """
        if gmx or self.traj.n_frames > 10000:
            # Assumes simulation was created with the polymer first, then solvated afterwards -
            # this means that the polymer is the group labelled 2
            _unwrapped = self.trajname.replace(".xtc", "_unwrapped.xtc")
            if sp.getstatusoutput("gmx")[0] == 0:
                gmxexe = "gmx"
            elif sp.getstatusoutput("gmx_mpi")[0] == 0:
                gmxexe = "gmx_mpi"
            else:
                raise RuntimeError(
                    "gmx or gmx_mpi not found. Please load GROMACS and try again."
                )
            print(f"Unwrapping trajectory using {gmxexe} trjconv...")
            sp.call(
                f'printf "2\n0" | {gmxexe} trjconv -f {self.trajname} -s {self.coordname} -o {_unwrapped} -center -pbc mol',
                shell=True,
            )
            self.trajname = _unwrapped
            self.load()
        else:
            transforms = [
                trans.unwrap(self.universe.atoms),
                trans.center_in_box(self.polymer),
                trans.wrap(self.solvent, compound="residues"),
            ]
            self.universe.trajectory.add_transformations(*transforms)

    def show_traj(self, hide_virtual_sites=True, hide_graphene=False):
        """
        View the trajectory using NGLView
        """
        if nv is None:
            raise ImportError("NGLView not installed")
        _selection = self.universe.atoms
        if hide_virtual_sites:
            _selection = _selection.select_atoms("not name MW")
        if hide_graphene:
            _selection = _selection.select_atoms("not resname GRA")
        view = nv.show_mdanalysis(_selection)
        view.add_unitcell()
        # solvent not loaded in later versions, so do this explicitly
        view.add_representation("ball+stick", "not resname pol")
        return view

    def compute_ion_counts(
        self, fname=None, restrained=True, correct_for_initial_salt_in_freshwater=False
    ):
        """
        Tracking ion numbers in each portion of the simulation box.
        The polymer is defined by the minimum and maximum z-coordinate of the polymer
        in each frame. Anything to the left is designated as saltwater and any space
        to the right is designated as freshwater. If restrained is True, the bounds of
        the polymer are computed once to save on computation time.

        Ion counts are computed per frame (given as time in ns) as a pandas dataframe, and
        (optionally) saved to a csv.

        If correct_for_initial_salt_in_freshwater is True, ions that are present in fresh water
        at t = 0 will be removed from the overall count. These ions are used to neutralise the
        membrane, and shouldn't contribute to overall ion migration.

        Returns:
            Pandas dataframe with columns:
            time, cl_salt, cl_mem, cl_fresh, na_salt, na_mem, na_fresh
        """
        # Polymer bounds constant with a restrained polymer
        min_z = self.polymer.positions[:, 2].min()
        max_z = self.polymer.positions[:, 2].max()
        counts = np.zeros((self.traj.n_frames, 7))
        for i, ts in enumerate(completion(self.traj)):
            if not restrained:  # need to recompute bounds
                min_z = self.polymer.positions[:, 2].min()
                max_z = self.polymer.positions[:, 2].max()
            cl_in_saltwater = self.universe.select_atoms(
                f"resname CL and prop z < {min_z}"
            )
            cl_in_membrane = self.universe.select_atoms(
                f"resname CL and prop z > {min_z} and prop z < {max_z}"
            )
            cl_in_freshwater = self.universe.select_atoms(
                f"resname CL and prop z > {max_z}"
            )
            na_in_saltwater = self.universe.select_atoms(
                f"resname NA and prop z < {min_z}"
            )
            na_in_membrane = self.universe.select_atoms(
                f"resname NA and prop z > {min_z} and prop z < {max_z}"
            )
            na_in_freshwater = self.universe.select_atoms(
                f"resname NA and prop z > {max_z}"
            )
            counts[i] = [
                ts.time,
                cl_in_saltwater.n_atoms,
                cl_in_membrane.n_atoms,
                cl_in_freshwater.n_atoms,
                na_in_saltwater.n_atoms,
                na_in_membrane.n_atoms,
                na_in_freshwater.n_atoms,
            ]
        counts[:, 0] *= PS_TO_NS
        if correct_for_initial_salt_in_freshwater:
            # Remove ions that are present in fresh water at t = 0
            # chloride
            counts[:, 3] -= counts[0, 3]
            counts[:, 3][counts[:, 3] < 0] = 0
            # sodium
            counts[:, 6] -= counts[0, 6]
            counts[:, 6][counts[:, 6] < 0] = 0
        df = pd.DataFrame(
            counts,
            columns=[
                "time_ns",
                "cl_salt",
                "cl_mem",
                "cl_fresh",
                "na_salt",
                "na_mem",
                "na_fresh",
            ],
        )
        if fname is not None:
            df.to_csv(fname, index=False)
        return df

    def compute_partial_densities(
        self,
        molecules=["pol", "SOL", "NA", "CL"],
        binwidth=0.25,
        dim="z",
        fname=None,
    ):
        """
        Compute the partial densities of molecules in the simulation box.
        If you wish to change which molecules are included, you can pass a list of
        resnames as the molecules argument.

        For a much quicker computation, do not unwrap the trajectory before calling this method.

        binwidth is the width of the bins in angstroms - density is computed on slices of the
        simulation box.

        Data is saved to a csv if fname is not None

        Returns:
            Pandas dataframe with columns:
            z-coord, density of each molecule passed in
        By default: coord, density of polymer, density of water, density of sodium, density of chloride

        """
        data = {}
        for res in molecules:
            _sel = self.universe.select_atoms(f"resname {res}")
            _lin = LinearDensity(_sel, binwidth=binwidth)
            _lin.run()
            data[res] = _lin.results[dim].pos

        data = pd.DataFrame(data)
        data.insert(0, f"{dim}_coord", np.arange(data.shape[0]) * binwidth)
        if fname is not None:
            data.to_csv(fname, index=False)
        return data
