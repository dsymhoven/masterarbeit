import numpy as np
import glob
import matplotlib.pyplot as plt

data = np.genfromtxt("Forces/dampingTermVSLorentzForce.txt")
time = data[:,0]
dampingForce = data[:,1]
lorentzForce = data[:,2]
ratio = data[:,3]

fig = plt.figure()
# plt.plot(time, dampingForce, label = "dampingForce")
# plt.plot(time, lorentzForce, label = "lorentzForce")
plt.plot(time, ratio, label = r"$\frac{F_{damping}}{F_{lorentz}}$")
plt.xlabel("Simulation time [s]")
plt.ylabel("Force [factor * N]")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

fig.savefig("Forces/dampingTermVSLorentzForce.png", bbox_inches='tight', dpi=300)

plt.close(fig)
