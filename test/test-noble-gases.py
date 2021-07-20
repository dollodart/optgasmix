from optgasmix.gases import Gas
from optgasmix.constants import NA, kB, pi
from make_gases_monatomic_dict import gases_monatomic

cp = lambda T: 5/2*kB

T = 798.15 # K
P = 101325. # Pa
gs = [gases_monatomic['He'], gases_monatomic['Ne'], 
      gases_monatomic['Ar'], gases_monatomic['Kr']]

# assume flat plate (Blasius solution)
a = 1/2
b = 1/3

fom_max = 0
y = []
for g in gs:
    gk = g.thermal_conductivity(T)
    gmu = g.viscosity(T)
    gcp = g.heat_capacity(T)
    gmm = g.mass
    fom = gk**(1-b) * gmu**(b-a) * gcp**(b) * gmm**(a-b)
    if fom > fom_max:
        fom_max = fom
    y.append(fom)

for i in range(len(gs)):
    print(gs[i].name, f'{round(y[i] / fom_max * 100):n}%')
