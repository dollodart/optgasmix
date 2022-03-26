# Physics of the Optimum Gas for Heat Transfer

Gases are rather conductors of heat relative to liquieds. However, they can be
used to directly cool a product in many cases where a liquid heat transfer
needs a shell. Some applications in semiconductor processing require direct
cooling of devices by gases, and though limited by reactivity considerations,
there are still some considerations for choosing among eligible inert gases.
These applications are at standard pressure or less in vacuum systems and near
standard temperature or greater, well within the ideal gas range. The advantage
of gases is their transport properties are relatively simple to evaluate at low
density where the kinetic theory of gases is valid. In this particular case,
the flow is usually forced-convection, unconfined, and laminar, which is
further ideal, since transport estimates for turbulent flow are much harder to
make.

As BSL notes in section 1.4 and asks the reader to work out in
problem 1A.2 (though the nearest problem is actually 1A.4), the
mixture properties for viscosity (and also thermal conductivity) can
be highly non-linear. Moreover, it is not just that the behavior
is non-linear interpolating, but that the mixture property can be
more extreme than the property of any of the single components. For
example, they give data showing the viscosity of 75% hydrogen and 25%
dichlorodifluoromethane as 135.1 micropoise, compared to 124 micropoise
for dichlorodifluoromethane and 88.4 micropoise for hydrogen.

Hence the problem of optimizing the thermal conductivity is non-trivial.
Moreover in applications with combined mass and heat transfer,
optimizing the gas composition so the heat transfer coefficient is
maximum is further non-trivial, since the molar mass in addition to the
thermal conductivity is highly relevant for the momentum transfer. 

Since the effective heat transfer coefficient (Nusselt number) depends on
geometries, the optimum heat transfer fluid cannot even be given just
as a function of the intensive thermodynamic properties of the system
(e.g., what temperatures and pressures), but must include in principle
the geometry and flow velocity. However, the sensitivity of the solution
to these parameters may not be high.

These calculations should enable the optimum to be found, and can be
combined with cost coefficients to determine the optimum gas mixture for
heat transfer with process economics.

Fundamentally, what should be decisive is the molecular parameters of
interaction, notably the energetic parameter of the Lennard-Jones, relative to
the molecular weight. The semi-empirical mixing rule doesn't require
estimating interaction parameters between unlike species, and it is remarkable
it works as well as it does--for gases with differing polarities, it isn't
expected that the energies of interactions would follow such a general rule.

# Uniformly Sampling Composition Space

It is simple to uniformly sample composition space in the case of two
variables, since there is by the condition of fractions summing to 1 only one
degree of freedom, which can take on a linear space. However, for higher cases
that is not the case. A distribution which has a support satisfying the
constraint of composition, namely that all mole fractions sum to 1, is sampled
to quickly obtain sample points. This is not uniform, but near uniform when
sufficiently densely sampled.

# Related Projects

More than one chemical engineering libraries allow you to compute the transport
properties of gases, though they may not implement multicomponent mixing rules
and instead use simple linear weighting. One library in wide use and well-funded is
Cantera.

# Data sources
- Collision cross-sections and lennard-jones parameters: BSL Transport
  Phenomena, 2nd Edition
- Heat capacities: Koretsky, Chemical and Engineering Thermodynamics

# Related Literature
- Optimum composition of gas mixture in a novel chimney-based LED bulb
