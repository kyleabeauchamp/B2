#!/usr/bin/env python

import os
import subprocess
import openmoltools
import numpy as np
import pandas as pd
import fire


def runner(in_csv):
    X = pd.read_csv(in_csv)

    data = []
    for k, smiles in enumerate(X.smiles.unique()):
        filename = "./mol2/%d.mol2" % k
        if not os.path.exists(filename):
            oemol = openmoltools.openeye.smiles_to_oemol(smiles)
            oemol = openmoltools.openeye.generate_conformers(oemol, strictStereo=False)
            openmoltools.openeye.molecule_to_mol2(oemol, filename)
        checkmol = subprocess.check_output(["checkmol", filename])
        data.append(dict(smiles=smiles, checkmol=checkmol))

    data = pd.DataFrame(data)
    v = data.checkmol.map(lambda x: x.split("\n"))
    groups = np.unique(sum(v))
    groups = np.array([g for g in groups if g not in ["", "cation"]])

    for group in groups:
        data[group] = data.checkmol.map(lambda x: group in x)

    counts = data[groups].sum(0)
    counts = pd.DataFrame(counts)
    print(counts.to_latex())


if __name__ == "__main__":
    fire.Fire(runner)
