#!/usr/bin/env python3
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from parmed import load_file
from mdtraj.reporters import XTCReporter
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
    return parser.parse_args()


def gen_simulation(
    grofile, topfile, chk=None, temp=300, press=None, timestep=1, interval=10000
):
    # OpenMM parsers can't read virtual sites, so use parmed to be safe
    gro = load_file(grofile)
    top = load_file(topfile)
    top.box = gro.box[:]
    system = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=12 * angstroms,
        constraints=HBonds,
        rigidWater=True,
    )
    if press is not None:
        baro = MonteCarloBarostat(press * bar, temp * kelvin)
        system.addForce(baro)
    integrator = LangevinIntegrator(
        temp * kelvin, 1 / picosecond, timestep * femtoseconds
    )
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

    simulation.reporters.append(XTCReporter("output.xtc", interval))
    simulation.reporters.append(CheckpointReporter("output.chk", interval))
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
        temp=args.temp,
        press=args.press,
        timestep=args.timestep,
        interval=args.interval,
    )
    sim.step(args.steps)


if __name__ == "__main__":
    main()
