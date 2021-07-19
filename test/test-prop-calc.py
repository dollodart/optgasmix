"""
Tests property calculation valdidity by comparing outputs to hydrogen at
standard conditions. STP hydrogen has the following properties:

molar mass: ~2 g/mol ~ 2 proton masses for every molecule
density: 0.09693 - 0.08078 kg/m^3 from Incropera
thermal conductivity: 157 - 183 mW/(m*K) according to Incropera's Heat Transfer, and 177.9 mW/(m*K) in BSL at 300 K
    This is for 1 cm length 1 m^2 cross section and 100 K differnce (boiling
    water and ice) a heat rate of 0.183 kW.  This is physically reasonable heat
    rate, though small: a vacuum cleaner draws something 10 Amps * 120 V = 1.2
    kW.  The thermal conductivity is around 10 times less than that of ice at
    1.88 W/(m*K) from Incropera.  This is suprisingly high given hydrogen is
    around 1000 times less dense than ice, but clearly lattice vibrations are
    not as effective ballistic collisions.
viscosity: .009 cP = 9.e-6 Pa*s from Eng. toolbox, 89.6e-7 Pa*s = 8.96e-6 Pa*s from Incropera
    This is around 1 % the viscosity of water.
heat capacity: 3.46*R from the CRC, 28.6 J/mol*K = 3.43*R from Incropera
"""

from make_gases_dict import gases, kB, mp
from time import time

t0 = time()

# STP (according to some conventions)
T = 298.15 # K
P = 101325. # Pa

for key in gases:
    g = gases[key]
    #
    grho = g.density(T, P)
    gcp = g.heat_capacity(T)
    gmu = g.viscosity(T)
    gcv = gcp - kB
    gkapprox = 5/2 * gmu * (gcv / g.mass) 
    #
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

print(f'total time {1000*(time() - t0):.3f}ms')
