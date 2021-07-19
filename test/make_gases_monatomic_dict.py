from optgasmix import lennard_jones_data_file
from optgasmix.props import kB, mp, NA, cp
from optgasmix.gases import Gas
import pandas as pd
tdf = pd.read_csv(lennard_jones_data_file)
tdf['sigma'] = tdf['sigma'] * 1.e-10 # Ang. -> m
tdf['mass'] = tdf['mol weight'] / NA / 1000 # g/mol -> kg
tdf['epsilon'] = tdf['eps/k'] * kB

gases_monatomic = dict()
for index, row in tdf.iterrows():
    row = row.to_dict()
    row['name'] = row['gas']
    row['heat_capacity_calculator'] = lambda T: 5/2*kB
    g = Gas(**row)
    g.Tmin = 0.
    g.Tmax = 1.e9
    gases_monatomic[row['gas']] = g
