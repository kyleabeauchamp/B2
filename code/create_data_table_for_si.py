#!/usr/bin/env python

import numpy as np
import pandas as pd
import fire


def runner(in_filename, out_filename):
    """
    This table will be used for comparing the various error estimates in the SI.
    """
    experiments = ["Mass density, kg/m3", "Relative permittivity at zero frequency"]
    #experiments = ["Relative permittivity at zero frequency"]

    X = pd.read_csv(in_filename)

    # Precalculate sigma ** -2 for use in weighted variance calculation
    for e in experiments:
        key0 = e + "_std"
        key1 = e + "_std" + "_invsquare"
        X[key1] = X[key0] ** -2.

    data = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[experiments].mean().dropna()
    counts = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[experiments].count().ix[data.index]

    uncertainty_std = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[experiments].std().ix[data.index]

    # What is the groupwise average (and std) of author-reported uncertainties?
    # Currently, not particularly useful because the std is undefined for most cases due to having only a single author-reported uncertainty
    uncertainty_author_averaged = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[[e + "_std" for e in experiments]].mean().ix[data.index]
    uncertainty_author_median = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[[e + "_std" for e in experiments]].median().ix[data.index]
    uncertainty_author_std = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[[e + "_std" for e in experiments]].std().ix[data.index]

    # What is the weighted groupwise std of author-reported uncertainties?
    # V(y) = [\sum_k \sigma_k^{-2}]^{-1}
    # std(y) = [\sum_k \sigma_k^{-2}]^{-0.5}
    # Use precalculated inverse squared uncertainties
    uncertainty_author = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[[e + "_std" + "_invsquare" for e in experiments]].sum().ix[data.index] ** -0.5


    uncertainty_author_median.columns = uncertainty_std.columns  # Strip off the _std from column names for now.
    uncertainty_author.columns = uncertainty_std.columns  # Strip off the _std from column names for now.
    uncertainty_author_averaged.columns = uncertainty_std.columns  # Strip off the _std from column names for now.
    uncertainty_author_std.columns = uncertainty_std.columns  # Strip off the _std from column names for now.


    # Previously preferred the standard deviation, but now we pick the *larger* error estimate.
    #mask = (uncertainty_std.isnull() & (~uncertainty_author.isnull()))
    #uncertainty_bestguess = uncertainty_std.copy()
    #uncertainty_bestguess[mask] = uncertainty_author[mask]

    uncertainty_author_median.replace(np.nan, -np.inf, inplace=True)
    uncertainty_author_averaged.replace(np.nan, -np.inf, inplace=True)
    uncertainty_std.replace(np.nan, -np.inf, inplace=True)
    uncertainty_author.replace(np.nan, -np.inf, inplace=True)

    uncertainty_bestguess = uncertainty_std.copy()
    ind = uncertainty_author > uncertainty_std
    uncertainty_bestguess[ind] = uncertainty_author[ind]

    uncertainty_author_median.replace(-np.inf, np.nan, inplace=True)
    uncertainty_author_averaged.replace(-np.inf, np.nan, inplace=True)
    uncertainty_std.replace(-np.inf, np.nan, inplace=True)
    uncertainty_author.replace(-np.inf, np.nan, inplace=True)
    uncertainty_bestguess.replace(-np.inf, np.nan, inplace=True)

    for e in experiments:
        data[e + "_uncertainty_std"] = uncertainty_std[e]
        data[e + "_uncertainty_author_weighted"] = uncertainty_author[e]
        data[e + "_uncertainty_author_averaged"] = uncertainty_author_averaged[e]
        data[e + "_uncertainty_author_median"] = uncertainty_author_median[e]
        data[e + "_uncertainty_bestguess"] = uncertainty_bestguess[e]
        data[e + "_counts"] = counts[e]

    data.to_csv(out_filename)


if __name__ == "__main__":
    fire.Fire(runner)
