#!/usr/bin/env python
import numpy as np
import sklearn.metrics
import sklearn.cross_validation
import statsmodels.formula.api as sm
import simtk.unit as u
import polarizability
import matplotlib.pyplot as plt
import pandas as pd
import fire


def runner(expt_csv, pred_csv, dens_pdf, diff_pdf, diel_pdf, nocorr_pdf):
    FIGURE_SIZE = (6.5, 6.5)
    DPI = 1600

    expt = pd.read_csv(expt_csv)
    expt["temperature"] = expt["Temperature, K"]

    pred = pd.read_csv(pred_csv)
    pred["polcorr"] = pd.Series(dict((cas, polarizability.dielectric_correction_from_formula(formula, density * u.grams / u.milliliter)) for cas, (formula, density) in pred[["formula", "density"]].iterrows()))
    pred["corrected_dielectric"] = pred["polcorr"] + pred["dielectric"]

    expt = expt.set_index(["cas", "temperature"])  # Can't do this because of duplicates  # Should be fixed now, probably due to the CAS / name duplication issue found by Julie.
    #expt = expt.groupby(["cas", "temperature"]).mean()  # Fix a couple of duplicates, not sure how they got there.
    pred = pred.set_index(["cas", "temperature"])

    pred["expt_density"] = expt["Mass density, kg/m3"]
    pred["expt_dielectric"] = expt["Relative permittivity at zero frequency"]
    #pred["expt_density_std"] = expt["Mass density, kg/m3_std"]
    pred["expt_density_std"] = expt["Mass density, kg/m3_uncertainty_bestguess"]
    #pred["expt_dielectric_std"] = expt["Relative permittivity at zero frequency_std"]
    pred["expt_dielectric_std"] = expt["Relative permittivity at zero frequency_uncertainty_bestguess"]

    plt.figure(figsize=FIGURE_SIZE, dpi=DPI)

    for (formula, grp) in pred.groupby("formula"):
        x, y = grp["density"], grp["expt_density"]
        xerr = grp["density_sigma"]
        yerr = grp["expt_density_std"].replace(np.nan, 0.0)
        x = x / 1000.  # Convert kg / m3 to g / mL
        y = y / 1000.  # Convert kg / m3 to g / mL
        xerr = xerr / 1000.  # Convert kg / m3 to g / mL
        yerr = yerr / 1000.  # Convert kg / m3 to g / mL
        plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', label=formula)

    plt.plot([.600, 1.400], [.600, 1.400], 'k', linewidth=1)
    plt.xlim((.600, 1.400))
    plt.ylim((.600, 1.400))
    plt.xlabel("Predicted (GAFF)")
    plt.ylabel("Experiment (ThermoML)")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.draw()

    x, y = pred["density"], pred["expt_density"]
    relative_rms = (((x - y) / x)**2).mean() ** 0.5
    # cv = sklearn.cross_validation.Bootstrap(len(x), train_size=len(x) - 1, n_iter=100)
    # relative_rms_grid = np.array([(((x[ind] - y[ind]) / x[ind])**2).mean() ** 0.5 for ind, _ in cv])
    # relative_rms_err = relative_rms_grid.std()
    plt.title(r"Density [g cm$^{-3}$]")
    plt.savefig(dens_pdf, bbox_inches="tight")

    plt.figure(figsize=FIGURE_SIZE, dpi=DPI)
    for (formula, grp) in pred.groupby("formula"):
        x, y = grp["density"], grp["expt_density"]
        xerr = grp["density_sigma"]
        yerr = grp["expt_density_std"].replace(np.nan, 0.0)
        x = x / 1000.  # Convert kg / m3 to g / mL
        y = y / 1000.  # Convert kg / m3 to g / mL
        xerr = xerr / 1000.  # Convert kg / m3 to g / mL
        yerr = yerr / 1000.  # Convert kg / m3 to g / mL
        plt.errorbar(x - y, y, xerr=xerr, yerr=yerr, fmt='.', label=formula)

    plt.xlim((-0.1, 0.1))
    plt.ylim((.600, 1.400))
    plt.xlabel("Predicted - Experiment")
    plt.ylabel("Experiment (ThermoML)")
    plt.gca().set_aspect('auto', adjustable='box')
    plt.draw()

    x, y = pred["density"], pred["expt_density"]
    relative_rms = (((x - y) / x)**2).mean() ** 0.5
    # cv = sklearn.cross_validation.Bootstrap(len(x), train_size=len(x) - 1, n_iter=100)
    # relative_rms_grid = np.array([(((x[ind] - y[ind]) / x[ind]) ** 2).mean() ** 0.5 for ind, _ in cv])
    # relative_rms_err = relative_rms_grid.std()
    plt.title(r"Density [g cm$^{-3}$]")

    plt.savefig(diff_pdf, bbox_inches="tight")


    yerr = pred["expt_dielectric_std"].replace(np.nan, 0.0)
    xerr = pred["dielectric_sigma"].replace(np.nan, 0.0)

    plt.figure(figsize=FIGURE_SIZE, dpi=DPI)

    plt.xlabel("Predicted (GAFF)")
    plt.ylabel("Experiment (ThermoML)")
    plt.title("Inverse Static Dielectric Constant")

    #ticks = np.concatenate([np.arange(1, 10), 10 * np.arange(1, 10)])

    #xticks(ticks)
    #yticks(ticks)

    plt.plot([0.0, 1], [0.0, 1], 'k')  # Guide
    #xscale('log')
    #yscale('log')


    x, y = pred["dielectric"], pred["expt_dielectric"]
    ols_model = sm.OLS(y, x)
    ols_results = ols_model.fit()
    r2 = ols_results.rsquared
    #plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', label="GAFF (R^2 = %.3f)" % r2)
    #plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', label="GAFF")
    plt.errorbar(x ** -1, y ** -1, xerr=xerr * x ** -2, yerr=yerr * y ** -2, fmt='.', label="GAFF")  # Transform xerr and yerr for 1 / epsilon plot

    plt.xlim((0.0, 1))
    plt.ylim((0.0, 1))
    plt.legend(loc=0)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.draw()
    plt.savefig(nocorr_pdf, bbox_inches="tight")


    x, y = pred["corrected_dielectric"], pred["expt_dielectric"]
    ols_model = sm.OLS(y, x)
    ols_results = ols_model.fit()
    r2 = ols_results.rsquared
    #plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', label="Corrected (R^2 = %.3f)" % r2)
    #plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', label="Corrected")
    plt.errorbar(x ** -1, y ** -1, xerr=xerr * x ** -2, yerr=yerr * y ** -2, fmt='.', label="Corrected")  # Transform xerr and yerr for 1 / epsilon plot

    plt.xlim((0.0, 1.02))
    plt.ylim((0.0, 1.02))
    plt.legend(loc=0)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.draw()
    plt.savefig(diel_pdf, bbox_inches="tight")


if __name__ == "__main__":
    fire.Fire(runner)
