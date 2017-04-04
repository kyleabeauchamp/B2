import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

sns.set_palette("bright")
sns.set_style("white")
sns.set(font_scale=1.2)
sns.set(style="ticks")

experiments = ["Mass density, kg/m3", "Relative permittivity at zero frequency"]

data = pd.read_csv("./tables/data_with_metadata.csv")

x = data["Mass density, kg/m3_uncertainty_author_weighted"]
y = data["Mass density, kg/m3_uncertainty_std"]
z = data["Mass density, kg/m3_uncertainty_author_median"]

x = x / 1000.  # Convert kg / m3 to g / mL
y = y / 1000.  # Convert kg / m3 to g / mL
z = z / 1000.  # Convert kg / m3 to g / mL

ind = x.dropna().index.intersection(y.dropna().index)
x, y = x[ind], y[ind]
z = z[ind]

M = max(x.max(), y.max())
plt.figure()
plt.plot(z, y, 'o', label="Median")
plt.plot(x, y, 'o', label="Weighted")
plt.plot([0, M], [0, M], 'k')
plt.legend(loc=0)

plt.xlabel("Author Uncertainty")
plt.ylabel("Standard Deviation of Measurements")
plt.title("Error Estimates: Density [g / cm^3]")

plt.xscale('log')
plt.yscale('log')

plt.savefig("./manuscript/figures/error_analysis_density.pdf", bbox_inches="tight")


x = data["Relative permittivity at zero frequency_uncertainty_author_weighted"]
y = data["Relative permittivity at zero frequency_uncertainty_std"]
z = data["Relative permittivity at zero frequency_uncertainty_author_median"]

ind = x.dropna().index.intersection(y.dropna().index)
x, y = x[ind], y[ind]
z = z[ind]

M = max(x.max(), y.max())
plt.figure()
plt.plot(z, y, 'o', label="Median")
plt.plot(x, y, 'o', label="Weighted")
plt.plot([0, M], [0, M], 'k')
plt.legend(loc=0)

plt.xlabel("Author Uncertainty")
plt.ylabel("Standard Deviation of Measurements")
plt.title("Error Estimates: Relative permittivity at zero frequency")

plt.xscale('log')
plt.yscale('log')

plt.savefig("./manuscript/figures/error_analysis_dielectric.pdf", bbox_inches="tight")


x = data["Mass density, kg/m3_uncertainty_author_weighted"]
y = data["Mass density, kg/m3_uncertainty_std"]
z = data["Mass density, kg/m3_uncertainty_author_median"]

x = x / 1000.  # Convert kg / m3 to g / mL
y = y / 1000.  # Convert kg / m3 to g / mL
z = z / 1000.  # Convert kg / m3 to g / mL

plt.figure()
plt.plot(z, 'o', label="Median(author)")  # plot this first so that weighted overwrites
plt.plot(x, 'o', label="Weighted(author)")
plt.plot(y, 'o', label="std()")

plt.xlabel("Measurment Index")
plt.ylabel("Error Estimate")
plt.title("Error Estimates: Density [g / cm^3]")
plt.legend(loc=0)

plt.yscale('log')

plt.savefig("./manuscript/figures/error_analysis_density_index.pdf", bbox_inches="tight")




x = data["Relative permittivity at zero frequency_uncertainty_author_weighted"]
y = data["Relative permittivity at zero frequency_uncertainty_std"]
z = data["Relative permittivity at zero frequency_uncertainty_author_median"]

plt.figure()
plt.plot(z, 'o', label="Median(author)")  # plot this first so that weighted overwrites
plt.plot(x, 'o', label="Weighted(author)")
plt.plot(y, 'o', label="std()")


plt.xlabel("Measurement Index")
plt.ylabel("Error Estimate")
plt.title("Error Estimates: Relative permittivity at zero frequency")
plt.legend(loc=0)

plt.yscale('log')

plt.savefig("./manuscript/figures/error_analysis_dielectric_index.pdf", bbox_inches="tight")
