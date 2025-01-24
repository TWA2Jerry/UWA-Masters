import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

bro = mpatches.Wedge((0, 0), 1.0, 0, 180, ec="none", alpha = 0.4)
dude = mpatches.Arc((0.0, 0.0), 0.2, 0.2, angle= 0.0, theta1 = 0.0, theta2 = 180.0)
fig, ax = plt.subplots()
ax.scatter(0.0, 0.0, color = 'k')
ax.add_artist(bro)
ax.add_artist(dude)
ax.set(
	aspect= 1,
	xlim = (-1.1, 1.1),
	ylim = (-1.1, 1.1)
)

ax.arrow(0., 0., 0., .5, head_width = 30*0.001, color = 'k')
ax.annotate(r"$\vec{v}_i$", (0.0, 0.), xytext = (0.1,0.5), size = 20)
ax.annotate(r"$\Phi$", (0.0, 0.0), xytext = (0.1, 0.1), size = 20) 
ax.set_axis_off()
plt.tight_layout(pad = 0.01, w_pad = 0.1)
#plt.show()
fig.savefig("../Records/stdod56/fov.pdf")
