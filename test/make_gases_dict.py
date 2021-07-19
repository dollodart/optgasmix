"""

This is not a module to reduce dependency on other packages. One could easily
process the input data to construct the Gas objects without using pandas.

"""

import pandas as pd
from optgasmix import lennard_jones_heat_capacities_data_file
from optgasmix.props import kB, mp, NA, cp
from optgasmix.gases import Gas

df = pd.read_csv(lennard_jones_heat_capacities_data_file)
df['B'] = df['B'] / 1e3 # constant factors (also defined in ogm)
df['C'] = df['C'] / 1e6
df['D'] = df['D'] / 1e-5
df['E'] = df['E'] / 1e9
df['sigma'] = df['sigma'] * 1.e-10 # Ang. -> m
df['mass'] = df['mol weight'] / NA / 1000 # g/mol -> kg
df['epsilon'] = df['eps/k'] * kB

gases = {}
for index, row in df.iterrows():
    row = row.to_dict()
    row['heat_capacity_calculator'] = lambda T: cp(T, row['A'],row['B'],row['C'],row['D'],row['E'])
    g = Gas(**row)
    gases[row['formula']] = g
