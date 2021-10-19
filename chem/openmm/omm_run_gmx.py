#!/usr/bin/env python3
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from parmed import load_file
from mdtraj.reporters import XTCReporter
import ommhelper as oh
import sys
import argparse


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gro", help="Gromacs coordinate file, default = conf.gro", default="conf.gro"
    )
    parser.add_argument(
        "--top",
        help="Gromacs coordinate file, default = topol.top",
        default="topol.top",
    )
    parser.add_argument(
        "--thermostat",
        help="Choice of langevin or nose-hoover, default = langevin",
        default="langevin",
    )
    parser.add_argument(
        "-t",
        "--temp",
        help="Temperature in Kelvin, default = 300",
        default=300,
        type=float,
    )
    parser.add_argument(
        "-p",
        "--press",
        help="Pressure in bar. Default is to apply no pressure",
        type=float,
    )
    parser.add_argument(
        "-dt",
        "--timestep",
        help="Timestep in femtoseconds, default = 1",
        default=1,
        type=float,
    )
    parser.add_argument("--chk", help="Checkpoint file to restart simulation")
    parser.add_argument(
        "-n", "--steps", help="Number of steps to run simulation for", type=int
    )
    parser.add_argument(
        "-i",
        "--interval",
        help="Number of steps to save trajectory, default = 10000",
        default=10000,
        type=int,
    )
    parser.add_argument(
        "-s",
        "--scalecharge",
        help="Scale atomic charges by a multiplication factor",
        type=float,
    )
    parser.add_argument(
        "-r",
        "--restrain",
        help=(
            "Add position restraints to these residues. "
            "Pass in a list of residue names separated by spaces"
        ),
        nargs="+",
    )
    parser.add_argument(
        "-min",
        "--minimise",
        help="Minimise energy before simulation",
        action="store_true",
    )
    return parser.parse_args()


def gen_simulation(
    grofile,
    topfile,
    chk=None,
    thermostat="langevin",
    temp=300,
    press=None,
    timestep=1,
    interval=10000,
    charge_factor=None,
    restrain=None,
):
    # OpenMM parsers can't read virtual sites, so use parmed to be safe
    gro = load_file(grofile)
    top = load_file(topfile)
    top.box = gro.box[:]

    if charge_factor is not None:
        for atom in top.atoms:
            atom.charge *= charge_factor

    system = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=12 * angstroms,
        constraints=HBonds,
        rigidWater=True,
    )
    if press is not None:
        baro = MonteCarloBarostat(press * bar, temp * kelvin)
        system.addForce(baro)

    if restrain is not None:
        atoms_to_restrain = [
            a.idx for a in top.atoms for x in restrain if a.residue.name == x
        ]
        force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addGlobalParameter("k", 5.0 * kilocalories_per_mole / angstroms ** 2)
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        for i in atoms_to_restrain:
            force.addParticle(i, list(gro.positions[i]))
        system.addForce(force)

    if thermostat == "langevin":
        integrator = LangevinIntegrator(
            temp * kelvin, 1 / picosecond, timestep * femtoseconds
        )
    elif thermostat == "nose-hoover":
        integrator = LangevinIntegrator(
            temp * kelvin, 1 / picosecond, timestep * femtoseconds
        )
    else:
        raise AttributeError("Thermostat not supported - use langevin or nose-hoover")

    simulation = Simulation(top.topology, system, integrator)
    if chk is not None:
        simulation.loadCheckpoint(chk)
        simulation.currentStep = (
            round(
                simulation.context.getState().getTime().value_in_unit(picoseconds)
                / timestep
                / 10
            )
            * 10
        )
        simulation.context.setTime(simulation.currentStep * timestep)
    else:
        simulation.context.setPositions(gro.positions)
        simulation.context.setVelocitiesToTemperature(temp * kelvin)

    simulation.reporters.append(XTCReporter("traj.xtc", interval))
    simulation.reporters.append(oh.CheckpointReporter("cpt.cpt", interval))
    simulation.reporters.append(
        StateDataReporter(
            sys.stdout,
            1000,
            step=True,
            time=True,
            kineticEnergy=True,
            potentialEnergy=True,
            totalEnergy=True,
            temperature=True,
            density=True,
            volume=True,
            speed=True,
        )
    )
    return simulation


def main():
    args = arguments()
    sim = gen_simulation(
        grofile=args.gro,
        topfile=args.top,
        chk=args.chk,
        thermostat=args.thermostat,
        temp=args.temp,
        press=args.press,
        timestep=args.timestep,
        interval=args.interval,
        charge_factor=args.scalecharge,
        restrain=args.restrain,
    )

    if args.minimise:
        sim.minimizeEnergy()
        state = sim.context.getState(getPositions=True)
        with open("min.gro", "w") as f:
            oh.GroFile.writeFile(
                sim.topology,
                state.getPositions(asNumpy=True),
                state.getPeriodicBoxVectors(),
                f,
            )

    sim.step(args.steps)


if __name__ == "__main__":
    main()
