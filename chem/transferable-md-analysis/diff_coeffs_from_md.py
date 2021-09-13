#!/usr/bin/env python3

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.core import groups
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import linregress


class EinsteinMSD(AnalysisBase):
    r"""Class to calculate Mean Squared Displacement by the Einstein relation.
    Taken from https://github.com/MDAnalysis/mdanalysis/blob/release-2.0.0-beta/package/MDAnalysis/analysis/msd.py
    Parameters
    ----------
    u : Universe or AtomGroup
        An MDAnalysis :class:`Universe` or :class:`AtomGroup`.
        Note that :class:`UpdatingAtomGroup` instances are not accepted.
    select : str
        A selection string. Defaults to "all" in which case
        all atoms are selected.
    msd_type : {'xyz', 'xy', 'yz', 'xz', 'x', 'y', 'z'}
        Desired dimensions to be included in the MSD. Defaults to 'xyz'.
    fft : bool
        If ``True``, uses a fast FFT based algorithm for computation of
        the MSD. Otherwise, use the simple "windowed" algorithm.
        The tidynamics package is required for `fft=True`.
        Defaults to ``True``.
    Attributes
    ----------
    dim_fac : int
        Dimensionality :math:`d` of the MSD.
    timeseries : :class:`numpy.ndarray`
        The averaged MSD over all the particles with respect to lag-time.
    msds_by_particle : :class:`numpy.ndarray`
        The MSD of each individual particle with respect to lag-time.
    ag : :class:`AtomGroup`
        The :class:`AtomGroup` resulting from your selection
    n_frames : int
        Number of frames included in the analysis.
    n_particles : int
        Number of particles MSD was calculated over.
    """

    def __init__(self, u, select="all", msd_type="xyz", fft=True, com=False, **kwargs):
        r"""
        Parameters
        ----------
        u : Universe or AtomGroup
            An MDAnalysis :class:`Universe` or :class:`AtomGroup`.
        select : str
            A selection string. Defaults to "all" in which case
            all atoms are selected.
        msd_type : {'xyz', 'xy', 'yz', 'xz', 'x', 'y', 'z'}
            Desired dimensions to be included in the MSD.
        fft : bool
            If ``True``, uses a fast FFT based algorithm for computation of
            the MSD. Otherwise, use the simple "windowed" algorithm.
            The tidynamics package is required for `fft=True`.
        com : bool
            If ``True``, compute MSDs for the centre of mass of each residue, 
            instead of all atoms
        """
        if isinstance(u, groups.UpdatingAtomGroup):
            raise TypeError("UpdatingAtomGroups are not valid for MSD " "computation")

        super(EinsteinMSD, self).__init__(u.universe.trajectory, **kwargs)

        # args
        self.select = select
        self.msd_type = msd_type
        self._parse_msd_type()
        self.fft = fft
        self.com = com

        # local
        self.ag = u.select_atoms(self.select)
        if self.com:
            self.n_particles = len(self.ag.residues.resids)  # 1 per molecule
        else:
            self.n_particles = len(self.ag)
        self._position_array = None

        # result
        self.msds_by_particle = None
        self.timeseries = None

    def _prepare(self):
        # self.n_frames only available here
        # these need to be zeroed prior to each run() call
        self.msds_by_particle = np.zeros((self.n_frames, self.n_particles))
        self._position_array = np.zeros((self.n_frames, self.n_particles, self.dim_fac))
        # self.timeseries not set here

    def _parse_msd_type(self):
        r""" Sets up the desired dimensionality of the MSD.
        """
        keys = {
            "x": [0],
            "y": [1],
            "z": [2],
            "xy": [0, 1],
            "xz": [0, 2],
            "yz": [1, 2],
            "xyz": [0, 1, 2],
        }

        self.msd_type = self.msd_type.lower()

        try:
            self._dim = keys[self.msd_type]
        except KeyError:
            raise ValueError(
                "invalid msd_type: {} specified, please specify one of xyz, "
                "xy, xz, yz, x, y, z".format(self.msd_type)
            )

        self.dim_fac = len(self._dim)

    def _single_frame(self):
        r""" Constructs array of positions for MSD calculation.
        """
        # shape of position array set here, use span in last dimension
        # from this point on
        if self.com:
            self._position_array[self._frame_index] = self.ag.center_of_mass(
                compound="residues"
            )[:, self._dim]
        else:
            self._position_array[self._frame_index] = self.ag.positions[:, self._dim]

    def _conclude(self):
        print(f"Computing {self.ag.atoms[0].resname} MSD")
        if self.fft:
            self._conclude_fft()
        else:
            self._conclude_simple()

    def _conclude_simple(self):
        r""" Calculates the MSD via the simple "windowed" algorithm.
        """
        lagtimes = np.arange(1, self.n_frames)
        positions = self._position_array.astype(np.float64)
        for lag in lagtimes:
            disp = positions[:-lag, :, :] - positions[lag:, :, :]
            sqdist = np.square(disp).sum(axis=-1)
            self.msds_by_particle[lag, :] = np.mean(sqdist, axis=0)
        self.timeseries = self.msds_by_particle.mean(axis=1)

    def _conclude_fft(self):  # with FFT, np.float64 bit prescision required.
        r""" Calculates the MSD via the FCA fast correlation algorithm.
        """
        try:
            import tidynamics
        except ImportError:
            raise ImportError(
                """ERROR --- tidynamics was not found!
                tidynamics is required to compute an FFT based MSD (default)
                try installing it using pip eg:
                    pip install tidynamics
                or set fft=False"""
            )

        positions = self._position_array.astype(np.float64)
        for n in range(self.n_particles):
            self.msds_by_particle[:, n] = tidynamics.msd(positions[:, n, :])
        self.timeseries = self.msds_by_particle.mean(axis=1)


def read_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute mean squared displacements for each residue, saving results to msd.csv. "
            "Note that PBC corrections to the diffusion coefficients printed might be required for "
            "agreement with experimental data."
        )
    )
    parser.add_argument(
        "-c",
        "--coords",
        help="Pass in a coordinate file (.gro/.tpr/.pdb)",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--trajectory",
        help="Pass in a trajectory file (.xtc/.trr)",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--start",
        help=(
            "Start time to compute diffusion coefficients. Default = 0.1, "
            "meaning that MSD is considered after 10 %% of the run"
        ),
        default=0.1,
        type=float,
    )
    parser.add_argument(
        "-e",
        "--end",
        help=(
            "Final time to compute diffusion coefficients. Default = 0.9, "
            "meaning that MSD values after 90 %% of the run are not considered"
        ),
        default=0.9,
        type=float,
    )
    parser.add_argument(
        "-dt",
        "--timestep",
        help="Interval between frames of the trajectory file, in fs. Default = 10000",
        default=10000,
        type=int,
    )
    parser.add_argument(
        "-p",
        "--plot",
        help=(
            "Plot MSD against time, saved to msd.png. "
            "The shaded area represents the region used to compute diffusion coefficients"
        ),
        action="store_true",
    )
    parser.add_argument(
        "-dim",
        "--dimensionality",
        help='Dimensions to consider for MSD, default="xyz"',
        default="xyz",
    )
    parser.add_argument(
        "-com",
        "--centre-of-mass",
        help=(
            "Compute MSDs and subsequent diffusion coefficients using the centre of mass of "
            "each residue instead of averaging over every atom."
        ),
        action="store_true",
    )
    parser.add_argument(
        "-r",
        "--replace",
        help=(
            "Pass in residue names you would like to change along with their replacement strings. "
            "For example, changing ch+ to ch and WAT to H2O by including: "
            "-r ch+ ch WAT H2O. "
            "LaTeX expressions are also carried through to the plot, if requested. This means you "
            "can write -r dhp 'DHP$^-$' and have a superscript minus in the plot legend. Avoid double "
            "quotes here otherwise bash/zsh variable expansion may occur"
        ),
        nargs="+",
    )
    return parser.parse_args()


def compute_msd(
    resname, universe, timestep=10000, dimensionality="xyz", fft=True, com=False
):
    """
    Compute mean squared displacements for a particular residue, according to the Einstein
    relation.
    The timestep is the interval between trajectory frames, in fs.
    The simulation is assumed to be 3D, otherwise change the dimensionality parameter.
    A Pandas DataFrame is returned.
    """
    FS_TO_NS = 1e-6
    MSD = EinsteinMSD(
        universe, select=f"resname {resname}", msd_type=dimensionality, fft=fft, com=com
    )
    MSD.run()

    msdisp = MSD.timeseries
    nframes = MSD.n_frames

    lagtimes = np.arange(nframes) * timestep * FS_TO_NS

    df = pd.DataFrame({"Time (ns)": lagtimes, "MSD": msdisp})
    df["resname"] = resname
    return df


def compute_diff_coeffs(df, start=0.1, end=0.9, dimensionality="xyz"):
    """
    Compute diffusion coefficients using MSD values returned by EinsteinMSD.
    Default is to use the middle 80 %, which is the default configuration of gmx msd
    and gives close agreement to TRAVIS in testing.
    """
    N = len(dimensionality)
    ANGSTROM2_PER_NS_TO_M2_PER_S = 1e-11

    maxtime = df["Time (ns)"].max()
    df = df[(start * maxtime <= df["Time (ns)"]) & (df["Time (ns)"] <= end * maxtime)]

    def diffcoeffs(group):
        regression = linregress(group["Time (ns)"], group["MSD"])
        return regression.slope / 2 * N * ANGSTROM2_PER_NS_TO_M2_PER_S

    return (
        df.groupby("resname")
        .apply(lambda g: diffcoeffs(g))
        .reset_index(name="D (m2/s)")
    )


def plot_msd(df, start=0.1, end=0.9):
    """
    Plot MSD against time, as computed by EinsteinMSD.
    Filled area represents the region used to compute 
    diffusion coefficients.
    """
    sns.set(style="white", palette="rainbow")
    plt.rcParams["mathtext.default"] = "regular"
    p = sns.lineplot(x="Time (ns)", y="MSD", hue="resname", data=df, ci=None)
    plt.axvspan(
        start * df["Time (ns)"].max(),
        end * df["Time (ns)"].max(),
        alpha=0.1,
        color="gray",
    )
    p.set_ylabel(r"MSD (${\AA}^2$)")
    p.legend().set_title(None)
    plt.tight_layout()  # y axis label sometimes cut off
    plt.savefig("msd.pdf", dpi=300)


def change_resnames(df, list_of_strings):
    """
    Pass in a list like ['ch+', 'Ch', 'WAT', 'SOL']
    to change residue names stored in df['resname'] 
    from ch+ to Ch and WAT to SOL.
    """
    original = list_of_strings[::2]
    replacements = list_of_strings[1::2]
    _repl = {i: j for i, j in zip(original, replacements)}
    df["resname"] = df["resname"].replace(_repl)
    return df


def check_args(args):
    _msg = (
        "Pass in a replacement for every residue name! "
        "The length of the list should be a multiple of 2"
    )
    if args.replace and len(args.replace) % 2 != 0:
        raise AttributeError(_msg)


def main():
    args = read_args()
    check_args(args)
    u = mda.Universe(args.coords, args.trajectory)
    msd = pd.concat(
        compute_msd(
            resname,
            u,
            timestep=args.timestep,
            dimensionality=args.dimensionality,
            fft=True,
            com=args.centre_of_mass,
        )
        for resname in np.unique(u.atoms.resnames)
    )
    if args.replace:
        msd = change_resnames(msd, args.replace)
    msd.to_csv("msd.csv", index=False)
    coeffs = compute_diff_coeffs(
        msd, start=args.start, end=args.end, dimensionality=args.dimensionality
    )
    coeffs.to_csv("diff_coeffs.csv", index=False)
    print(coeffs.to_string(index=False))
    if args.plot:
        plot_msd(msd, start=args.start, end=args.end)


if __name__ == "__main__":
    main()
