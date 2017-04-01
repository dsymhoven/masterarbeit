import matplotlib.pyplot as plt
import numpy as np
import math


fig = plt.figure()
ax = plt.subplot(111)

x = np.arange(0,5,0.01)

f1 = 2 * np.arcsin(0.5 * np.sin(x / 2))
f2 = 0.5 * x
plt.plot(x,f1, label = 'numeric dispersion relation')
plt.plot(x,f2, label = 'analytic dispersion relation')

plt.plot([math.pi,math.pi],[0,3], 'k-', alpha = 0.7, linewidth = 0.5)
plt.plot([math.pi,math.pi],[2 * np.arcsin(0.5), math.pi/2], 'k-', alpha = 1.0, linewidth = 1.5)

plt.plot([math.pi], [ 2 * np.arcsin(0.5)], marker='x', markersize=5, color="red")
ax.legend()
plt.ylabel('$\omega h$')
plt.xlabel('$k_x\Delta x$')

# ax.annotate('local max', xy=(2, 1), xytext=(3, 1.5),
#             arrowprops=dict(facecolor='black', shrink=0.05),
#             )

filename = 'dispersionRelation'
fig.savefig("{}.png".format(filename), bbox_inches='tight', dpi = 300)
# close fig
plt.close(fig)
