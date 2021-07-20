from optgasmix import collision_cross_section_data_file
from optgasmix.props import omega, omega_D
import matplotlib.pyplot as plt
import numpy as np

lines = open(collision_cross_section_data_file, 'r').readlines()
vals = []
for line in lines[1:]:
    line = line.rstrip('\n')
    line = line.split(',')
    line = [float(val) for val in line]
    vals.append(line)

x = 10.**np.linspace(np.log10(0.3), np.log10(100), 100)
plt.loglog(x, omega(x), 'r-', label=r'$\Omega_k=\Omega_\mu$')
plt.loglog(x, omega_D(x), 'b-', label=r'$\Omega_D$')
xd, omega_table, omega_D_table = zip(*vals)
plt.loglog(xd, omega_table, 'rx')
plt.loglog(xd, omega_D_table, 'bx')

plt.xlabel(r'$k_B T/\epsilon$')
plt.ylabel(r'$\Omega$')
plt.legend()
plt.show()

# residual difference
# is nearly constant at high kBT/epsilon
plt.figure()
plt.plot(x, 2 * (omega(x) - omega_D(x)) / (omega(x) + omega_D(x)))
plt.xlabel(r'$k_B T/\epsilon$')
plt.ylabel(r'$2\frac{\Omega - \Omega_D}{\Omega + \Omega_D}$')
plt.show()
