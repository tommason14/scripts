import subprocess as sp
import sys
import pandas as pd
import MDAnalysis as mda
import MDAnalysis.transformations as trans
import nglview as nv
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

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
        rc={"mathtext.default": "regular"},
    )
    p = sns.relplot(
        data=tidy,
        x="time_ns",
        y="count",
        hue="region",
        col="ion",
        kind="line",
        height=3,
        aspect=1,
        palette="Set2",
        facet_kws={"sharey": False},
    )
    p.set(xlabel="Time (ns)", ylabel="Ion Count")
    sns.move_legend(
        p, loc="lower center", ncol=3, title=None
    )  # no direction param, so control orientation with number of columns
    p.set_titles(col_template="{col_name}")
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
        Unwrap the trajectory, and wrap molecules so that the centre of mass
        lies in the periodic box.
        For large trajectories, frame-by-frame unwrapping and analysis takes a long time,
        so gmx trjconv can be used, and the universe is re-created using the unwrapped trajectory.
        """
        if gmx:
            # Assumes simulation was created with the polymer first, then solvated afterwards - 
            # this means that the polymer is the group labelled 2
            _unwrapped = self.trajname.replace(".xtc", "_unwrapped.xtc")
            sp.call(f'printf "2\n0" | gmx trjconv -f {self.trajname} -s {self.coordname} -o {_unwrapped} -center -pbc mol', shell=True)
            self.trajname = _unwrapped
            self.load()
        else:
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

    def compute_ion_counts(self, fname=None):
        """
        Tracking ion numbers in each portion of the simulation box.
        The polymer is defined by the minimum and maximum z-coordinate of the polymer
        in each frame. Anything to the left is designated as saltwater and any space
        to the right is designated as freshwater.
        Pass in the residue name of the ion of interest.

        Ion counts are computed per frame (given as time in ns) as a pandas dataframe, and
        (optionally) saved to a csv.

        Returns:
            Pandas dataframe with columns:
            time, cl_salt, cl_mem, cl_fresh, na_salt, na_mem, na_fresh
        """
        # Polymer bounds constant with a restrained polymer
        min_z = self.polymer.positions[:, 2].min()
        max_z = self.polymer.positions[:, 2].max()
        counts = np.zeros((self.traj.n_frames, 7))
        for i, ts in enumerate(completion(self.traj)):
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
        # convert time from ps to ns
        counts[:, 0] /= 1000
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

    ### when using unwrapped traj, MDAnalysis and vmd are similar in time so might as well just use MDAnalysis
    # def compute_ion_counts(self):
    #     """
    #     Instead of numpy, use vmd to compute ion counts, then read in the csv into a pandas
    #     DataFrame.
    #     """
    #     # first unwrap traj with gmx trjconv - much faster than in vmd
    #     # assumes sim was created by placing polymer in centre of box first, so that the polymer
    #     # is the second group in the trjconv output
    #     self.unwrapped_trajname = "traj_unwrapped.xtc"
    #     os.system(
    #         f'echo -e "2\n0" | gmx trjconv -f {self.trajname} -s {self.coordname} -o {self.unwrapped_trajname} -center -pbc mol'
    #     )
    #     # vmd can't read tpr files, so convert to gro
    #     if self.coordname.endswith(".tpr"):
    #         self.grofile = self.coordname[:-4] + ".gro"
    #         self.universe.atoms.write(self.grofile)
    #     else:
    #         self.grofile = self.coordname
    #     os.system(
    #         f"vmd {self.grofile} {self.unwrapped_trajname} -e ~/.local/scripts/chem/vmd/ion_counts_unwrapped.tcl -dispdev text"
    #     )
    #     df = pd.read_csv("ion_counts.csv")
    #     times = np.arange(0, self.traj.n_frames + 1) * self.dt * PS_TO_NS
    #     df = df.drop(columns=["frame"])
    #     df.insert(0, "time_ns", times)
    #     return df
