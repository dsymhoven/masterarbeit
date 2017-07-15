import numpy as np
import glob
import matplotlib.pyplot as plt

data = np.genfromtxt("Forces/dampingTermVSLorentzForce.txt")
time = data[:,0]
dampingForce = data[:,1]
lorentzForce = data[:,2]
ratio = data[:,3]
fontsize = 18

fig = plt.figure()
# plt.plot(time, dampingForce, label = "dampingForce")
# plt.plot(time, lorentzForce, label = "lorentzForce")
plt.plot(time, ratio)
plt.xlabel("Simulation time", fontsize = fontsize)
plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
plt.ylabel(r"$g^{\mu}$ / $F_{L}$", fontsize = fontsize)
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = fontsize)

fig.savefig("Forces/dampingTermVSLorentzForce.png", bbox_inches='tight', dpi=300)

plt.close(fig)
