import numpy as np
from .constants import *

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

def cp(T, A, B, C, D, E):
    """Gas heat capacity, empirical equation from Koretsky.
    Number specific per molecule, units J/K. """
    return kB * (A + B * T + C * T**2 + D * T**(-2) + E * T**3)
