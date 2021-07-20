from .constants import *
from .props import omega, mu, k, kpoly
from .mixing import mixrule, weighted_geometric_mean, weighted_arithmetic_mean, weighted_harmonic_mean

def memo_checker(func):
    """
    Decorator for memoization.
    """
    def inner_func(self, temperature, pressure=None): # will contain pressure
        if temperature in self.memos:
            if func in self.memos[temperature]:
                return self.memos[temperature][func]
                # you can use first class objects like functions as dict keys
        else:
            self.memos[temperature] = dict()
        if pressure is None:
            val = func(self, temperature)
        else:
            val = func(self, temperature, pressure)
        self.memos[temperature][func] = val
        return val
    return inner_func

class Gas:
    """
    Gas class. While temperature and pressure are state variables, they are not
    chemical identities or material properties, and so are inputs to methods
    rather than attributes.
    """
    def __init__(self, name, mass, sigma, epsilon, heat_capacity_calculator, *args, **kwargs):
        self.name = name
        self.mass = mass
        self.sigma = sigma
        self.epsilon = epsilon
        self.heat_capacity_calculator = heat_capacity_calculator

        self.Tmin = kwargs.get('Tmin', None)
        self.Tmax = kwargs.get('Tmax', None)
        self.memos = {}

    def _eval_w(self,temperature):
        factor = kB*temperature / self.epsilon
        w = omega(factor)
        return w

    @memo_checker
    def density(self, temperature, pressure):
        """
        Default: mass specific density
        """
        return pressure * self.mass / (kB*temperature) 

    @memo_checker
    def specific_volume(self, temperature, pressure):
        """
        Default: mole specific volume.
        """
        return kB*temperature / pressure

    @memo_checker
    def viscosity(self, temperature): # monatomic valid for polyatomic
        return mu(temperature, self.mass, self.sigma, self._eval_w(temperature))

    @memo_checker
    def thermal_conductivity_monatomic(self, temperature):
        return k(temperature, self.mass, self.sigma, self._eval_w(temperature))

    @memo_checker
    def heat_capacity(self, temperature):
        if temperature > self.Tmin and temperature < self.Tmax:
            heat_capacity = self.heat_capacity_calculator(temperature)
            return heat_capacity
        else:
            message = f"temperature {temperature} is outside of range [{self.Tmin}, {self.Tmax}]"
            raise Exception(message)

    @memo_checker
    def thermal_conductivity(self, temperature):
        return kpoly(temperature, self.heat_capacity(temperature), self.mass, self.viscosity(temperature))

class GasMixture:
    """
    Gas mixture, defined in terms of gases and their compositions as state
    variables. Temperature and pressure are input to methods rather than a state because
    it is a condition, rather than chemical identity or material property of the gas mixture.
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
        """

        Each gas acts independently for density. It's as if it occupies
        a fraction of the volume in proportion to its composition.

        Total Mass/Total Volume = sum(density_i*volume_i,i)/sum(volume) = sum(density_i*composition_i,i)/1

        This is an arithmetic mean.
        """

        return weighted_arithmetic_mean(tuple(g.mass*pressure/(kB*temperature) for g in self.gases), self.compositions)

    def heat_capacity(self, temperature):
        """

        Each gas acts independently for heat capacity. It's as if
        it occupies a fraction of the volume in proportion to its
        composition.

        Total Heat Absorbed/Temperature Delta = (sum(cp_i*deltaT*number_i,i)/deltaT)
        Total Heat Capacity = (Total Heat Absorbed/Temperature Delta)/Total Number 
        = sum(cp_i*number_i,i)/sum(number_i,i) = sum(cp_i*composition_i,i)/1

        This is an arithmetic mean.
        """
        return weighted_arithmetic_mean(tuple(g.heat_capacity(temperature) for g in self.gases), self.compositions)

    def mass(self):
        """
        Gas mixture don't have atomic property of mass, returning weighted geometric mean.
        """
        return weighted_geometric_mean(tuple(g.mass for g in self.gases), self.compositions)

    def __str__(self):
        st = ''
        for i in range(len(self.gases)):
            st += self.gases[i].name + ' ' + str(compositions[i]) + '\n'
        return st.strip('\n')
