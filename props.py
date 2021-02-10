import pandas as pd
import numpy as np
kB = 1.38e-23 # J / K
# note: 1 kB * NA = 8.314 J/(mol*K)
Cvhat = 3. / 2.
mp = 1.672e-27 # kg
NA = 6.022e23
# note (1 g/mol) / Na / 1000 ~ proton mass

def omega(kTovereps):
    """BSL appendix E, applies to viscosity and thermal conductivity."""
    a = 1.16145
    b = 0.14874
    c = 0.52487
    d = 0.77320
    e = 2.16178
    f = 2.43787
    x = kTovereps
    return a / x**b + c / np.exp(d * x) + e / np.exp(f * x)


def omega_D(kTovereps):
    """BSL appendix E, applies to diffusivity."""
    a = 1.06036
    b = 0.15610
    c = 0.19300
    d = 0.47635
    e = 1.03587
    f = 1.52996
    g = 1.76474
    h = 3.89411
    x = kTovereps
    return a / x**b + c / np.exp(d * x) + e / np.exp(f * x) + g / np.exp(h * x)

#for ktovereps in range(1, 1000):
#    ktovereps /= 10
#    print(omega(ktovereps))
#    print(omega_D(ktovereps))
#import sys; sys.exit()

def k(T, m, sigma, omega, Cp):
    """

    eqn (9.3-13) BSL
    only valid for monatomic gases! (Eucken formula needed for polyatomic gases)
    note ideal gases all obey Cp = Cv + kB as a thermodynamic law
    here mass specific heat capacity at constant volume is found
    from number specific heat capacity at constant pressure

    """
    Cvhat = (Cp - kB) / m
    return 25. / 32. * np.sqrt(np.pi * kB * m * T) / (np.pi * sigma**2 * omega) * Cvhat

def mu(T, m, sigma, omega):
    """
    eqn (1.4-14) BSL
    derived for monatomic molecules but satisfactory for polyatomic ones due to center of mass coordinate 
    being more important for momentum transfer
    than internal coordinates (the same is not true for heat)

    """
    return 5. / 16. * np.sqrt(np.pi * kB * m * T) / (np.pi * sigma**2 * omega)

def kpoly(T, Cp, m, mu):
    """
    equation (9.3-15) BSL
    note the heat capacity is mass-specific in BSL equation, here equation has been changed for molar specific input
    """
    return (Cp + 5/4 * kB ) * mu / m

def mixrule(mu, mm, x):
    """
    equation (1.4-15) and (9.3-17) of BSL
    """
    mumix = 0
    n = len(x)
    for a in range(n):
        phib = []
        for b in range(n):
            mmr = mm[a] / mm[b]
            mur = mu[a] / mu[b]
            res = (1/8)**(1/2)*(1 + mmr)**(-1/2)*(1 + mur**(1/2)*mmr**(1/4))**2
            phib.append(res)
        mumix += x[a]*mu[a]/sum(phib[i] * x[i] for i in range(n))
    return mumix


def cp(T, A, B, C, D, E):
    """Gas heat capacity, empirical equation from Koretsky.
    Number specific per molecule, units J/K. """
    return kB * (A + B * T + C * T**2 + D * T**(-2) + E * T**3)

class Gas:
    def __init__(self, name, mass, sigma, epsilon, A, B, C, D, E, *args, **kwargs):
        self.name = name
        self.mass = mass
        self.sigma = sigma
        self.epsilon = epsilon
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E
        self.Tmin = kwargs['Tmin']
        self.Tmax = kwargs['Tmax']
        self.memos = {}

    def _eval_w(self,temperature):
        factor = kB*temperature / self.epsilon
        w = omega(factor)
        return w

    def check_memo(self,temperature, quantity_name):
        if temperature in self.memos:
            if quantity_name in self.memos[temperature]:
                return self.memos[temperature][quantity_name]
        self.memos[temperature] = dict()
        return None

    def density(self, temperature, pressure):
        """
        Default: mass specific density
        """
        memo = self.check_memo(temperature, 'density')
        if memo is None:
            density = pressure * self.mass / (kB*temperature) 
            self.memos[temperature]['density'] = density
            return density
        return memo 

    def specific_volume(self, temperature, pressure):
        """
        Default: mole specific volume.
        """
        memo = self.check_memo(temperature, 'specific volume')
        if memo is None:
            specific_volume = kB*temperature / pressure
            self.memos[temperature]['specific_volume'] = specific_volume
            return specific_volume
        return memo

    def viscosity(self, temperature): # monatomic valid for polyatomic
        memo = self.check_memo(temperature, 'viscosity')
        if memo is None:
            viscosity = mu(temperature, self.mass, self.sigma, self._eval_w(temperature))
            self.memos[temperature]['viscosity'] = viscosity
            return viscosity
        return memo

    def thermal_conductivity_monatomic(self, temperature):
        memo = self.check_memo(temperature, 'thermal conductivity monatomic')
        if memo is None:
            thermal_conductivity_monatomic = mu(temperature, self.mass, self.sigma, self._eval_w(temperature))
            self.memos[temperature]['thermal conductivity monatomic'] = thermal_conductivity_monatomic
            return thermal_conductivity_monatomic
        return memo

    def heat_capacity(self, temperature):
        memo = self.check_memo(temperature, 'heat capacity')
        if memo is None:
            if temperature > self.Tmin and temperature < self.Tmax:
                heat_capacity = cp(temperature, self.A, self.B, self.C, self.D, self.E)
                self.memos[temperature]['heat capacity'] = heat_capacity
                return heat_capacity
            else:
                message = f"temperature {temperature} is outside of range [{self.Tmin}, {self.Tmax}]"
                raise Exception(message)
                return 0
        return memo


    def thermal_conductivity(self, temperature):
        memo = self.check_memo(temperature, 'thermal conductivity')
        if memo is None:
            thermal_conductivity = kpoly(temperature, self.heat_capacity(temperature), self.mass, self.viscosity(temperature))
            self.memos[temperature]['thermal conductivity'] = thermal_conductivity
            return thermal_conductivity
        return memo

class GasMixture:
    """
    TODO: check_memo viscosity calculations in the case temperature is held constant and mole fraction varied
    """
    def __init__(self, gases, compositions):
        self.gases = gases
        self._compositions = compositions

    @property
    def compositions(self):
        """
        Molar composition.
        """
        return self._compositions

    @compositions.setter
    def compositions(self, value):
        if abs(sum(value) - 1.0) > 0.01:
            raise Exception('mole fractions do not sum to 1')
        self._compositions = value


    def viscosity(self, temperature):
        mus = []
        mm = []
        for i in range(len(self.gases)):
            mm.append(self.gases[i].mass)
            mus.append(self.gases[i].viscosity(temperature))

        return mixrule(mus, mm, self.compositions)

    def thermal_conductivity(self, temperature):
        ts = []
        mm = []
        for i in range(len(self.gases)):
            mm.append(self.gases[i].mass)
            ts.append(self.gases[i].thermal_conductivity(temperature))

        return mixrule(ts, mm, self.compositions)

    def density(self, temperature, pressure):
        """Each gas acts independently, for thermodynamic properties."""
        weighted = sum(self.gases[i].mass * self.compositions[i] for i in range(len(self.gases)))
        return weighted*pressure/(kB*temperature)

    def heat_capacity(self, temperature):
        return sum(self.gases[i].heat_capacity(temperature) * self.compositions[i] for i in range(len(self.gases)))

    def mass(self):
        """Note: gas mixture doesn't have atomic property of mass, returning arithmetic mean"""
        return sum(self.gases[i].mass * self.compositions[i] for i in range(len(self.gases)))

    def __str__(self):
        st = ''
        for i in range(len(self.gases)):
            st += self.gases[i].name + ' ' + str(compositions[i]) + '\n'
        return st.strip('\n')

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
    g = Gas(**row)
    gases[row['formula']] = g
