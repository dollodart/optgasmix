"""
Tests

Compare output toStandard conditions hydrogen 
molar mass: ~2 g/mol ~ 2 proton masses for every molecule

density: 0.09693 - 0.08078 kg/m^3 from Incropera

thermal conductivity: 157 - 183 mW/(m*K) according to Incropera's Heat Transfer, and 177.9 mW/(m*K) in BSL at 300 K
this is for 1 cm length 1 m^2 cross section and 100 K differnce (boiling water and ice) a heat rate of 0.183 kW
this is physically reasonable, though small: a vacuum cleaner draws something 10 Amps * 120 V = 1.2 kW
it is around 10 times less than the thermal conductivity of ice at 1.88 W/(m*K) from Incropera
which is suprisingly high given it is around 1000 times less dense than ice

viscosity: .009 cP = 9.e-6 Pa*s from Eng. toolbox, 89.6e-7 Pa*s = 8.96e-6 Pa*s from Incropera
this is around 1 % the viscosity of water

heat capacity: 3.46*R from the CRC, 28.6 J/mol*K = 3.43*R from Incropera
"""

import pandas as pd
from props import kB, mp, NA, mixrule, cp, Gas, GasMixture
from time import time

df = pd.read_csv('data/merged.csv')
df['B'] = df['B'] / 1e3 # constant factors
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

t0 = time()
for key in gases:
    g = gases[key]

    T = 298.15 # K
    P = 101325. # Pa
    grho = g.density(T, P)
    gcp = g.heat_capacity(T)
    gmu = g.viscosity(T)
    gcv = gcp - kB
    gkapprox = 5/2 * gmu * (gcv / g.mass) 
    print(key, 'm/proton mass', g.mass/mp)
    print(key, 'm/kg', g.mass)
    print(key, 'rho/rho_water', grho/1000) # water is 1000 kg/m^3 by definition of units
    print(key, 'rho/(kg/m^3)', grho) 
    print(key, 'Cp/R', gcp/kB) # fractions of gas constant
    print(key, 'Cp/(J/mol*K)', gcp/kB*8.314) 
    print(key, 'mu/mu_water', gmu/1e-3)
    print(key, 'mu/(Pa*s)', gmu)
    print(key, 'k/k_ice', gkapprox/1.88) # from Incropera, ice has 1.88 W/(m*K)
    print(key, 'k/(W/(m*K))', gkapprox)
    print(key, 'k/(W/(m*K)) (poly correction)', g.thermal_conductivity(T))

print('\n\n')
# example 1.4-2 from BSL Ed 2
mumix = mixrule( (1462e-7, 2031e-7, 1754e-7), (44.01,32.00,28.02), (0.133, 0.039, 0.828))
print(mumix, 1714e-7, mumix / 1714e-7) # compare to answer
assert abs(mumix - 1714e-7)/1714e-7 < 0.01

print('\n\n')
# different mixtures of O2 and N2
gs = [gases['O2'], gases['N2']]
gm = GasMixture(gs, None)
print('xO2 xN2 k/(W/m*k)')
for x in range(0, 11):
    x /= 10
    gm.compositions = (x, 1 - x)
    print(x, round(1-x, 3), round(gm.thermal_conductivity(T), 5))

print('\n\n')
# ternary mixtures
gs = [gases['O2'], gases['Cl2'], gases['N2']]
gm = GasMixture(gs, None)
from rand import gen_random_3
mx = 0
print('xO2','xCl2','xN2','therm. cond.')
tdict = dict()
scale = 20
for compositions in gen_random_3(scale):
    gm.compositions = compositions
    tc = gm.thermal_conductivity(T) # in SI units
#    print(compositions[0],compositions[1],compositions[2],tc)
    tdict[(int(scale*compositions[0]), int(scale*compositions[1]), int(scale*compositions[2]))] = tc
    if tc > mx:
        mxc = compositions
        mx = tc
print('optimum', mxc, mx)
print(f'total time {time() - t0:.3f}s')
