from make_gases_dict import gases
from optgasmix.mixing import mixrule
from optgasmix.gases import GasMixture
from optgasmix.rand import gen_random_3
import ternary

T = 298.15 # K

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
mx = 0
tdict = dict()
scale = 20

print('xO2','xCl2','xN2','therm. cond.')
for compositions in gen_random_3(scale):
    gm.compositions = compositions
    tc = gm.thermal_conductivity(T) # in SI units
#    print(compositions[0],compositions[1],compositions[2],tc)
    tdict[(int(scale*compositions[0]), int(scale*compositions[1]), int(scale*compositions[2]))] = tc
    if tc > mx:
        mxc = compositions
        mx = tc
print('optimum', mxc, mx)

figure, tax = ternary.figure(scale=scale)
tax.boundary(linewidth=2.0)
tax.gridlines(color="black",multiple=5)
tax.gridlines(color="blue",multiple=1,linewidth=0.5)

fontsize = 20
tax.set_title("Ternary Gas Mixture Thermal Conductivity", fontsize=fontsize)
tax.left_axis_label(f"Fraction {gs[0].name}/%", fontsize=fontsize)
tax.right_axis_label(f"Fraction {gs[1].name}/%", fontsize=fontsize)
tax.bottom_axis_label(f"Fraction {gs[2].name}/%", fontsize=fontsize)

# set ticks for composition space
ticks = tuple(str(x) for x in range(0, 101, 5))
tax.ticks(ticks=ticks, axis='lbr', linewidth=1)

# remove matplotlib default ticks and axes
tax.clear_matplotlib_ticks()
tax.get_axes().axis('off')

# make heatmap
tax.heatmap(tdict, scale)

ternary.plt.show()
