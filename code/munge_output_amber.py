#!/usr/bin/env python
import parmed
from openmoltools import cirpy
import mdtraj as md
import pymbar
import pandas as pd
import dipole_errorbars
import fire


def predict(in_prmtop, in_csv, in_dcd, out_csv, cas, temperature):
    num_bootstrap = 100
    fixed_block_length = 20  # 200 ps blocks for dielectric error bar block averaging.

    print(cas, temperature)

    traj = md.load(in_dcd, top=in_prmtop)

    if traj.unitcell_lengths is None:
        raise(ValueError("No unitcell lengths!"))

    rho = pd.read_csv(in_csv)["Density (g/mL)"].values * 1000.  # g / mL -> kg /m3
    initial_traj_length = len(traj)
    initial_density_length = len(rho)
    [t0, g, Neff] = pymbar.timeseries.detectEquilibration(rho)
    mu = rho[t0:].mean()
    sigma = rho[t0:].std() * Neff ** -0.5
    prmtop = parmed.load_file(in_prmtop)
    charges = prmtop.to_dataframe().charge.values
    temperature = float(temperature)
    traj = traj[t0 * len(traj) / len(rho):]
    dielectric = md.geometry.static_dielectric(traj, charges, temperature)
    dielectric_sigma_fixedblock = dipole_errorbars.bootstrap_old(traj, charges, temperature, fixed_block_length)[1]
    block_length = dipole_errorbars.find_block_size(traj, charges, temperature)
    dielectric_sigma = dipole_errorbars.bootstrap(traj, charges, temperature, block_length, num_bootstrap)
    formula = cirpy.resolve(cas, "formula")
    data = dict(cas=cas, temperature=temperature, n_trimmed=t0, inefficiency=g,
                initial_traj_length=initial_traj_length, initial_density_length=initial_density_length,
                density=mu, density_sigma=sigma, Neff=Neff, n_frames=traj.n_frames, dielectric=dielectric,
                dielectric_sigma=dielectric_sigma, dielectric_sigma_fixedblock=dielectric_sigma_fixedblock, block_length=block_length, formula=formula)

    data = pd.Series(data)

    data.to_csv(out_csv)


def merge(incsv, outcsv):
    print(incsv)
    csv_filenames = incsv.split(",")
    x = pd.DataFrame([pd.read_csv(filename) for filename in csv_filenames])
    x.to_csv(outcsv)


if __name__ == "__main__":
    commands = {"predict": predict, "merge": merge}
    fire.Fire(commands)
