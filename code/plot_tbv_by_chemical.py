import statsmodels.formula.api as sm
import simtk.unit as u
import polarizability
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

sns.set_palette("bright")
sns.set_style("white")
sns.set(font_scale=1.2)
sns.set(style="ticks")

expt = pd.read_csv("./tables/data_with_metadata.csv")
expt["temperature"] = expt["Temperature, K"]


pred = pd.read_csv("./tables/predictions.csv")
pred["polcorr"] = pd.Series(dict((cas, polarizability.dielectric_correction_from_formula(formula, density * u.grams / u.milliliter)) for cas, (formula, density) in pred[["formula", "density"]].iterrows()))
pred["corrected_dielectric"] = pred["polcorr"] + pred["dielectric"]

expt = expt.set_index(["cas", "temperature"])  # Can't do this because of duplicates  # Should be fixed now, probably due to the CAS / name duplication issue found by Julie.
pred = pred.set_index(["cas", "temperature"])

pred["name"] = expt["components"]
pred["expt_density"] = expt["Mass density, kg/m3"]
pred["expt_dielectric"] = expt["Relative permittivity at zero frequency"]
pred["expt_density_std"] = expt["Mass density, kg/m3_uncertainty_bestguess"]
pred["expt_dielectric_std"] = expt["Relative permittivity at zero frequency_uncertainty_bestguess"]

pred = pred.reset_index()


# HACK TO PLOT DATA THAT LACKS ERRORBARS!
pred.expt_density_std = pred.expt_density_std.fillna(0.0)
names = pred.name.unique()
name_dict = {name:i for i, name in enumerate(names)}
pred["name_ind"] = pred.name.map(lambda x: name_dict[x])

num_groups = 5
num_cols = 3
DENSITY_BUFFER = 10. / 1000.  # Convert kg / m3 to g / mL
PANEL_SIZE = 4.0
pred["name_group"] = pred.name_ind % num_groups

for col in ["expt_density", "expt_density_std", "density", "density_sigma"]:
    pred[col] = pred[col] / 1000.  # Convert kg / m3 to g / mL

for name_group, pred_i in pred.groupby("name_group"):
    g = sns.FacetGrid(pred_i, col="name", col_wrap=num_cols, xlim=[270, 330], ylim=[pred.density.min() - DENSITY_BUFFER, pred.density.max() + DENSITY_BUFFER], size=PANEL_SIZE, sharex=False, sharey=False)
    g.map(plt.errorbar, "temperature", "expt_density", "expt_density_std", fmt='.', color='b', label="Expt")
    g.map(plt.errorbar, "temperature", "density", "density_sigma", fmt='.', color='g', label="MD")
    g.set_ylabels("Density [g / cm^3]")
    g.set_xlabels("Temperature [K]")
    legend(loc=4)
    g.set_titles(col_template="{col_name}")
    g.set_xticklabels(labels=np.arange(270, 340, 10))
    plt.draw()
    plt.savefig("./manuscript/figures/densities_versus_temperature_part%d.pdf" % name_group, bbox_inches="tight")


# HACK TO PLOT DATA THAT LACKS ERRORBARS!
pred.expt_dielectric_std = pred.expt_dielectric_std.fillna(0.0)
# NEED TO DO ERROR ESTIMATES ON EPSILON
pred["corrected_dielectric_sigma"] = pred.dielectric_sigma
#pred["corrected_dielectric_sigma"] = 0.0
#pred["dielectric_sigma"] = 0.0
# End HACKERY

# Add dataframe columns for inverse static dielectric and associated error bars
pred["inv_expt_dielectric"] = pred.expt_dielectric ** -1.
pred["inv_dielectric"] = pred.dielectric ** -1.
pred["inv_corrected_dielectric"] = pred.corrected_dielectric ** -1.
pred["inv_expt_dielectric_std"] = pred.expt_dielectric_std * pred.expt_dielectric ** -2.
pred["inv_dielectric_sigma"] = pred.dielectric_sigma * pred.dielectric ** -2.
pred["inv_corrected_dielectric_sigma"] = pred.corrected_dielectric_sigma * pred.corrected_dielectric ** -2.


for name_group, pred_i in pred.groupby("name_group"):
    g = sns.FacetGrid(pred_i, col="name", col_wrap=num_cols, xlim=[270, 330], ylim=[1.0, 150.], size=PANEL_SIZE, sharex=False, sharey=False)
    g.map(plt.errorbar, "temperature", "inv_expt_dielectric", "inv_expt_dielectric_std", fmt='.', color='b', label="Expt")
    g.map(plt.errorbar, "temperature", "inv_dielectric", "inv_dielectric_sigma", fmt='.', color='g', label="MD")
    g.map(plt.errorbar, "temperature", "inv_corrected_dielectric", "inv_corrected_dielectric_sigma", fmt='.', color='r', label="Corr.")
    g.set_ylabels("Inverse Dielectric Constant")
    g.set_xlabels("Temperature [K]")
    legend(loc=2)  # Upper Left
    g.set_titles(col_template="{col_name}")
    g.set(ylim=(0., 1.05))
    plt.draw()
    plt.savefig("./manuscript/figures/dielectric_versus_temperature_part%d.pdf" % name_group, bbox_inches="tight")
