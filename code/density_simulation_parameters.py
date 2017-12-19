"""
This file contains all simulation parameters used in this project.
"""
from simtk import unit as u


MOLECULES_PER_BOX = 1000  # smirnoff eats up 16GB of RAM with 1,000 molecules

CUTOFF = 0.95 * u.nanometers

PRESSURE = 1.0 * u.atmospheres
BAROSTAT_FREQUENCY = 25

FRICTION = 1.0 / u.picoseconds
EQUIL_FRICTION = 5.0 / u.picoseconds


EQUIL_TIMESTEP = 0.4 * u.femtoseconds
TIMESTEP = 1.0 * u.femtoseconds

FF_NAME = "SMIRNOFF" #Used for plot labels
#FF_NAME = "GAFF" #Used for plot labels

#N_STEPS = 1000000 # 1.0 ns (at a time)
#N_EQUIL_STEPS = 10000000
N_STEPS = 500000 # 1.0 ns (at a time)
N_EQUIL_STEPS = 400000

OUTPUT_FREQUENCY_EQUIL = 100000
OUTPUT_FREQUENCY = 5000 # 5ps
OUTPUT_DATA_FREQUENCY = 250 # 0.25ps

#STD_ERROR_TOLERANCE = 0.0002 # g/mL
STD_ERROR_TOLERANCE = 0.008 # g/mL


# DEBUG PARAMETERS
N_STEPS = 100000
N_EQUIL_STEPS = 50000
OUTPUT_FREQUENCY_EQUIL = 10000
STD_ERROR_TOLERANCE = 0.004 # g/mL
