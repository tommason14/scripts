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
import sys


class EinsteinMSD(AnalysisBase):
    r"""Class to calculate Mean Squared Displacement by the Einstein relation.
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

    def __init__(self, u, select="all", msd_type="xyz", fft=True, **kwargs):
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
        """
        if isinstance(u, groups.UpdatingAtomGroup):
            raise TypeError("UpdatingAtomGroups are not valid for MSD " "computation")

        super(EinsteinMSD, self).__init__(u.universe.trajectory, **kwargs)

        # args
        self.select = select
        self.msd_type = msd_type
        self._parse_msd_type()
        self.fft = fft

        # local
        self.ag = u.select_atoms(self.select)
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
        self._position_array[self._frame_index] = self.ag.positions[:, self._dim]

    def _conclude(self):
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


################# END OF MSD #####################

parser = argparse.ArgumentParser()
parser.add_argument(
    "-c", "--coords", help="Pass in a coordinate file (.gro/.tpr/.pdb)", required=True
)
parser.add_argument(
    "-t", "--trajectory", help="Pass in a trajectory file (.xtc/.trr)", required=True
)
parser.add_argument(
    "-dt",
    "--timestep",
    help="Interval between frames of the trajectory file, in fs. Default = 10000",
    default=10000,
    type=int,
    required=True,
)
parser.add_argument(
    "-r",
    "--distance",
    help="Distance within which sodium ions are considered condensed to the polymer, in angstroms",
    type=float,
    required=True,
)

args = parser.parse_args()

u = mda.Universe(args.coords, args.trajectory)

# %%
######################
#  Diffusion coeffs  #
######################

# - msd with MDAnalysis
# - linear regression to find gradient

MSDc = EinsteinMSD(
    u, select=f"name Na and around {args.distance} name S*", msd_type="xyz", fft=True
)

MSDc.run()

msdisp_c = MSDc.timeseries
nframes = MSDc.n_frames

lagtimes = np.arange(nframes) * args.timestep

df = pd.DataFrame({"time": lagtimes, "MSD": msdisp_c,})
# time to ns
df["time"] /= 10 ** 6
df.to_csv("msd.csv", index=False)
sns.lineplot(x="time", y="MSD", data=df)
plt.xlabel("Time (ns)")
plt.ylabel(r"MSD (${\AA}^2$)")
plt.savefig("msd.png", dpi=300)

df = df[(df.time > 5) & (df.time < 9)]
df.to_csv("msd_fit.csv", index=False)
sns.lineplot(x="time", y="MSD", data=df)  # overlay onto original plot
plt.savefig("msd_fit.png", dpi=300)
lin = linregress(df["time"], df["MSD"])
grad = lin.slope
D = grad / 6  # 1/6 gradient = diff coeff
diff = f"Diffusion coeffs of sodiums within {args.distance} Å of sulfurs = {D /1e7 :<.4e} cm²/s"
print(diff)
with open("diff_coeff_5_to_9.txt", "w") as f:
    f.write(f'# {" ".join([sys.argv[0].split("/")[-1]] + sys.argv[1:])}\n')
    f.write(diff + "\n")
