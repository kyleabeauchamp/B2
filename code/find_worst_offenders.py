import sklearn.metrics, sklearn.cross_validation
import statsmodels.formula.api as sm
import simtk.unit as u
import polarizability
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

sns.set_palette("bright")
sns.set_style("whitegrid")
sns.set(font_scale=1.2)


expt0 = pd.read_csv("./tables/data_dielectric.csv")

expt = pd.read_csv("./tables/data_with_metadata.csv")
expt["temperature"] = expt["Temperature, K"]


pred = pd.read_csv("./tables/predictions.csv")
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


q = abs(pred.density - pred.expt_density)
q.sort()
cas = q.reset_index()[-20:].cas.unique()

expt0[expt0.cas.isin(cas)].components.unique()


pred.density.mean()
pred.density.std()
pred.density.std() / sqrt(len(pred))

pred.expt_density.mean()
pred.expt_density.std() / sqrt(len(pred))
