#!/usr/bin/env python3
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from parmed import load_file
from mdtraj.reporters import XTCReporter
import ommhelper as oh
import os
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
        "--cos",
        type=float,
        default=0,
        help="Cosine acceleration for periodic-perturbation viscosity calculations in units of nm/ps²",
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
    cos=0,
):
    # OpenMM parsers can't read virtual sites, so use parmed to be safe
    print("Creating simulation system...")
    gro = load_file(grofile)
    top = load_file(topfile)
    top.box = gro.box[:]

    if charge_factor is not None:
        print(f"Charge scaling by {charge_factor}...")
        for atom in top.atoms:
            atom.charge *= charge_factor

    system = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=12 * angstroms,
        constraints=HBonds,
        rigidWater=True,
    )
    if press is not None:
        print(f"Applying barostat at {press} bar...")
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

    print(f"Applying {thermostat} thermostat at {temp} K...")
    # use Zheng's integrator if cosine acceleration is desired
    if cos != 0:
        from velocityverletplugin import VVIntegrator

        integrator = VVIntegrator(
            temp * kelvin,
            10 / picoseconds,
            1 * kelvin,
            40 / picoseconds,
            timestep * femtoseconds,
        )
        integrator.setCosAcceleration(cos)
    if thermostat == "langevin":
        integrator = LangevinIntegrator(
            temp * kelvin, 1 / picosecond, timestep * femtoseconds
        )
    elif thermostat == "nose-hoover":
        integrator = NoseHooverIntegrator(
            temp * kelvin, 1 / picosecond, timestep * femtoseconds
        )
    else:
        raise AttributeError("Thermostat not supported - use langevin or nose-hoover")

    simulation = Simulation(top.topology, system, integrator)
    if chk is not None:
        print(f"Loading {chk}...")
        simulation.loadCheckpoint(chk)
        simulation.currentStep = (
            round(
                simulation.context.getState().getTime().value_in_unit(picoseconds)
                / (timestep / 1000)
                / 10
            )
            * 10
        )
        simulation.context.setTime(simulation.currentStep * (timestep / 1000))
    else:
        simulation.context.setPositions(gro.positions)
        simulation.context.setVelocitiesToTemperature(temp * kelvin)

    append_xtc = "traj.xtc" in os.listdir(".")
    simulation.reporters.append(XTCReporter("traj.xtc", interval, append=append_xtc))
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
        cos=args.cos,
    )

    if args.minimise:
        print("Minimising...")
        sim.minimizeEnergy()
        state = sim.context.getState(getPositions=True)
        with open("min.gro", "w") as f:
            oh.GroFile.writeFile(
                sim.topology,
                state.getPositions(asNumpy=True),
                state.getPeriodicBoxVectors(),
                f,
            )
    print("Running...")
    sim.step(args.steps)


if __name__ == "__main__":
    main()
