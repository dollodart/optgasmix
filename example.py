from props import gases, kB, mp, GasMixture

# different mixtures of O2 and N2
T = 298.15
P = 101325.
P = 1 # molecular flow regime
gs = [gases['H2'], gases['N2']]
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

#if pr < 100 and 0.01 < pr:
#    print('outside correlation ranges of prandtl',
#          'this will only work for high reynold\'s number,'
#          'but still must not be turbulent (<2300 for tubes)', sep='\n')
#    diam = 0.3 # 300 mm wafer
#    velocity = gm.density(T,P)/(diam*gm.viscosity(T))
#    print(gm.density(T,P))
#    print(gm.viscosity(T)) # viscosity is ~ 1e-5 Pa*s
#    print(f'velocity between {velocity*100} and {velocity*2300} m/s')


a = 1/2
b = 1/5
print(f'x{gs[0].name}=1-x{gs[1].name} fom')
ksolid = 30 # high temperature conductivity of Silicon, W/(m*K)
for x in range(0, 11):
    x /= 10
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
