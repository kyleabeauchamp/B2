"""Extract druglike liquid density and dielectric data from ThermoML.
Also counts the number of measurements remaining 
"""
import pandas as pd
from thermopyl import thermoml_lib, cirpy

data = pd.read_hdf("./data.h5", 'data')

bad_filenames = ["./10.1016/j.fluid.2013.12.014.xml"]  # This file confirmed to have possible data entry errors.
data = data[~data.filename.isin(bad_filenames)]

experiments = ["Mass density, kg/m3", "Relative permittivity at zero frequency"]

ind_list = [data[exp].dropna().index for exp in experiments]
ind = reduce(lambda x,y: x.union(y), ind_list)
X = data.ix[ind]

name_to_formula = pd.read_hdf("./compound_name_to_formula.h5", 'data')
name_to_formula = name_to_formula.dropna()

X["n_components"] = X.components.apply(lambda x: len(x.split("__")))
X = X[X.n_components == 1]
X.dropna(axis=1, how='all', inplace=True)

counts_data = {}
counts_data["0.  Single Component"] = X.count()[experiments]

X["formula"] = X.components.apply(lambda chemical: name_to_formula[chemical])

heavy_atoms = ["N", "C", "O", "S", "Cl", "Br", "F"]
desired_atoms = ["H"] + heavy_atoms

X["n_atoms"] = X.formula.apply(lambda formula_string : thermoml_lib.count_atoms(formula_string))
X["n_heavy_atoms"] = X.formula.apply(lambda formula_string : thermoml_lib.count_atoms_in_set(formula_string, heavy_atoms))
X["n_desired_atoms"] = X.formula.apply(lambda formula_string : thermoml_lib.count_atoms_in_set(formula_string, desired_atoms))
X["n_other_atoms"] = X.n_atoms - X.n_desired_atoms

X = X[X.n_other_atoms == 0]

counts_data["1.  Druglike Elements"] = X.count()[experiments]

X = X[X.n_heavy_atoms > 0]
X = X[X.n_heavy_atoms <= 10]
X.dropna(axis=1, how='all', inplace=True)

counts_data["2.  Heavy Atoms"] = X.count()[experiments]

X["smiles"] = X.components.apply(lambda x: cirpy.resolve(x, "smiles"))  # This should be cached via sklearn.
X = X[X.smiles != None]
X = X.ix[X.smiles.dropna().index]
    
X["cas"] = X.components.apply(lambda x: thermoml_lib.get_first_entry(cirpy.resolve(x, "cas")))  # This should be cached via sklearn.
X = X[X.cas != None]
X = X.ix[X.cas.dropna().index]

# Neither names (components) nor smiles are unique.  Use CAS to ensure consistency.
cannonical_smiles_lookup = X.groupby("cas").smiles.first()
cannonical_components_lookup = X.groupby("cas").components.first()

X["smiles"] = X.cas.apply(lambda x: cannonical_smiles_lookup[x])
X["components"] = X.cas.apply(lambda x: cannonical_components_lookup[x])


X = X[X["Temperature, K"] > 270]
X = X[X["Temperature, K"] < 330]

counts_data["3.  Temperature"] = X.count()[experiments]

X = X[X["Pressure, kPa"] > 100.]
X = X[X["Pressure, kPa"] < 102.]

counts_data["4.  Pressure"] = X.count()[experiments]

X = X[X.phase == "Liquid"]

counts_data["5.  Liquid state"] = X.count()[experiments]


X.dropna(axis=1, how='all', inplace=True)



X["Pressure, kPa"] = 101.325  # Assume everything within range is comparable.  
X["Temperature, K"] = X["Temperature, K"].apply(lambda x: np.round(x, 1))  # Round at the 0.1 digit.  

X.to_csv("./tables/full_filtered_data.csv")


mu = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[experiments].mean()

counts_data["6.  Aggregate T, P"] = mu.count()[experiments]

counts_data = pd.DataFrame(counts_data).T

q = mu.reset_index()
q = q.ix[q[experiments].dropna().index]
q.to_csv("./tables/data_dielectric.csv")

counts_data.ix["7.  Density+Dielectric"] = len(q)

print counts_data.to_latex()
