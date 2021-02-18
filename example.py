from gases import Gas, GasMixture
from constants import NA, kB

def shomates_cp(x, A, B, C, D, E):
    # x = T / 1000 in Kelvin
    # returns J/(mol*K)
    return A  + B * x + C * x**2 + D * x**3 + E/x**2

# between 298 and 6000 K
cp = lambda T: shomates_cp(T/1000, 20.78603,4.850683e-10,-1.582916e-10,1.525102e-11,3.196347e-11)/NA

# heat capacity data from NIST thermochemical 
he = Gas(name='Helium',
    mass=4.002602/NA/1000, 
    sigma=2.576e-10,
    epsilon=10.2*kB,
    heat_capacity_calculator=cp,
    Tmin=298,
    Tmax=6000)

def cp(T):
    x = T/1000
    if 100 < T <= 500:
        return shomates_cp(x, 28.98641, 1.853978, -9.647459, 16.63537, 0.000117)/NA
    elif 500 < T <= 2000:
        return shomates_cp(x, 19.60683, 19.88705, -8.598535, 1.367974, 0.527601)/NA
    elif 2000 < T <= 6000:
        return shomates_cp(x, 35.51872, 1.128728, -0.196103, 0.014662, -4.553760)/NA

n2 = Gas(name='Nitrogen',
        mass=28.0134/NA/1000,
        sigma=3.667e-10,
        epsilon=99.8*kB,
        heat_capacity_calculator=cp,
        Tmin=100,
        Tmax=6000)

# different mixtures of O2 and N2
T = 298.15
P = 101325.
gs = [he, n2]
gm = GasMixture(gs, (0.5,0.5))

# from Incropera, air at standard conditions should be around 1 kg/m^3
print('density (kg/m^3)')
for g in gs:
    print(g.name, g.density(T,P))
print('mixture', gm.density(T,P))

print('Prandtl number')
for g in gs:
    print(g.name, g.viscosity(T) * g.heat_capacity(T)/(g.thermal_conductivity(T)*g.mass))

pr = gm.viscosity(T) * gm.heat_capacity(T) / (gm.thermal_conductivity(T)*gm.mass()) # mass is geometrically averaged property?
print('mixture', pr) 

# flat plate (Blasius solution)
a = 1/2
b = 1/5
print(f'x{gs[0].name}=1-x{gs[1].name} fom')
ksolid = 30 # high temperature conductivity of Silicon, W/(m*K)

opt = -1
for x in range(0, 101):
    x /= 100
    gm.compositions = (x, 1 - x)
    gmm = gm.mass()
    gcp = gm.heat_capacity(T)
    gmu = gm.viscosity(T)
    gk = gm.thermal_conductivity(T)
    fom = gk**(1-b) * gmu**(b-a) * gcp**(b) * gmm**(a-b)
    grho = gm.density(T, P)

    re = grho*0.3*1./gmu
    pr = gmu * gcp / (gk * gmm)
    nu = 2/9* re**a * pr**b
    bi = nu*gk/ksolid
    assert bi < 0.1
    print(x, fom)
    if fom > opt:
        opt = fom
        optmix = x
print('-'*72)
print('optimum')
print(f'x{gs[0].name}={optmix*100:.1f}%, x{gs[1].name}={(1-optmix)*100:.1f}%')
