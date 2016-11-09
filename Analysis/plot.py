import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import cm
import scipy.interpolate

N = 100
files = glob.glob('fields/*.txt')

for file in files:
	x, y, z = np.genfromtxt(file, unpack=True)
	
	xi = np.linspace(x.min(), x.max(), N)
	yi = np.linspace(y.min(), y.max(), N)
	[X,Y]=np.meshgrid(xi,yi);
	#zi = griddata(x,y,z,X,Y, interp='linear');
	zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
	
	fig = plt.figure()
	plt.contourf(xi, yi, zi,100)
	plt.colorbar()
	plt.xlabel("X")
	plt.ylabel("Y")
	filename = file.replace(".txt", "")
	filename = file.replace("fields/", "")
	fig.savefig("png/" + "{}.png".format(filename), bbox_inches='tight')
	#plt.show()