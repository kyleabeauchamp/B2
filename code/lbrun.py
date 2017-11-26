#!/usr/bin/env python
import numpy as np
import mdtraj as md

from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u

import openmoltools
import smirnoffmixture
from density_simulation_parameters import (EQUIL_FRICTION, EQUIL_TIMESTEP, PRESSURE, BAROSTAT_FREQUENCY,
                                           CUTOFF, STD_ERROR_TOLERANCE, N_STEPS, OUTPUT_DATA_FREQUENCY,
                                           OUTPUT_FREQUENCY, FRICTION, TIMESTEP, OUTPUT_FREQUENCY_EQUIL, N_EQUIL_STEPS)

from pymbar import timeseries as ts
import pandas as pd
import fire


def build_monomer(cas, mol2_filename, frcmod_filename):
    """Generate GAFF mol2 and frcmod files for each chemical."""
    print(cas)
    smiles_string = openmoltools.cirpy.resolve(cas, 'smiles')
    openmoltools.openeye.smiles_to_antechamber(smiles_string, mol2_filename, frcmod_filename)


def build_box(in_mol2, in_frcmod, out_pdb, n_monomers, out_inpcrd, out_prmtop, ffxml=None):
    """Build an initial box with packmol and use it to generate AMBER files."""
    # NB: I removed the mixure support here for simplicity.
    gaff_mol2_filenames = [in_mol2]
    frcmod_filenames = [in_frcmod]
    n_monomers = [n_monomers]
    print(gaff_mol2_filenames)
    print(n_monomers)
    packed_trj = openmoltools.packmol.pack_box([md.load(mol2) for mol2 in gaff_mol2_filenames], n_monomers)
    packed_trj.save(out_pdb)

    if ffxml is None:
        tleap_cmd = openmoltools.amber.build_mixture_prmtop(gaff_mol2_filenames, frcmod_filenames, out_pdb, out_prmtop, out_inpcrd)
    else:
        tleap_cmd = smirnoffmixture.build_mixture_prmtop(gaff_mol2_filenames, frcmod_filenames, out_pdb, out_prmtop, out_inpcrd, ffxml=ffxml)


def equilibrate(in_prmtop, in_inpcrd, out_dcd, out_pdb, temperature):
    temperature = temperature * u.kelvin  # TODO: recycle John's simtk.unit parser

    prmtop = app.AmberPrmtopFile(in_prmtop)
    inpcrd = app.AmberInpcrdFile(in_inpcrd)

    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=CUTOFF, constraints=app.HBonds)
    integrator = mm.LangevinIntegrator(temperature, EQUIL_FRICTION, EQUIL_TIMESTEP)
    system.addForce(mm.MonteCarloBarostat(PRESSURE, temperature, BAROSTAT_FREQUENCY))

    simulation = app.Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)

    state = simulation.context.getState(getEnergy=True)
    print(state.getPotentialEnergy())

    print('Minimizing.')
    simulation.minimizeEnergy()

    state = simulation.context.getState(getEnergy=True)
    print(state.getPotentialEnergy())

    simulation.context.setVelocitiesToTemperature(temperature)
    print('Equilibrating.')

    simulation.reporters.append(app.DCDReporter(out_dcd, OUTPUT_FREQUENCY_EQUIL))
    simulation.step(N_EQUIL_STEPS)

    # Re-write a better PDB with correct box sizes.
    traj = md.load(out_dcd, top=in_prmtop)[-1]
    traj.save(out_pdb)


def production(in_prmtop, in_pdb, out_dcd, out_csv, temperature):
    temperature = temperature * u.kelvin  # TODO: recycle John's simtk.unit parser

    prmtop = app.AmberPrmtopFile(in_prmtop)
    pdb = app.PDBFile(in_pdb)

    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=CUTOFF, constraints=app.HBonds)

    integrator = mm.LangevinIntegrator(temperature, FRICTION, TIMESTEP)
    system.addForce(mm.MonteCarloBarostat(PRESSURE, temperature, BAROSTAT_FREQUENCY))

    simulation = app.Simulation(prmtop.topology, system, integrator)

    simulation.context.setPositions(pdb.positions)
    simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())
    simulation.context.setVelocitiesToTemperature(temperature)

    print('Production.')

    simulation.reporters.append(app.DCDReporter(out_dcd, OUTPUT_FREQUENCY))
    simulation.reporters.append(app.StateDataReporter(out_csv, OUTPUT_DATA_FREQUENCY, step=True, potentialEnergy=True, temperature=True, density=True))

    converged = False
    while not converged:
        simulation.step(N_STEPS)
        d = pd.read_csv(out_csv, names=["step", "U", "Temperature", "Density"], skiprows=1)
        density_ts = np.array(d.Density)
        [t0, g, Neff] = ts.detectEquilibration(density_ts, nskip=1000)
        density_ts = density_ts[t0:]
        density_mean_stderr = density_ts.std() / np.sqrt(Neff)
        if density_mean_stderr < STD_ERROR_TOLERANCE:
            converged = True


if __name__ == "__main__":
    cmdline_dict = {
        "build_monomer": build_monomer,
        "build_box": build_box,
        "equilibrate": equilibrate,
        "production": production,
                    }
    fire.Fire(cmdline_dict)
