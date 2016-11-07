import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import cm
import scipy.interpolate

N = 100

x, y, z = np.genfromtxt(r'fields.txt', unpack=True)

xi = np.linspace(x.min(), x.max(), N)
yi = np.linspace(y.min(), y.max(), N)
[X,Y]=np.meshgrid(xi,yi);
zi = griddata(x,y,z,X,Y, interp='linear');
#zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

fig = plt.figure()
plt.contourf(xi, yi, zi,100)
plt.colorbar()
plt.xlabel("X")
plt.ylabel("Y")
plt.show()