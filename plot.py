from tests import tdict, scale, gs
import ternary

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
