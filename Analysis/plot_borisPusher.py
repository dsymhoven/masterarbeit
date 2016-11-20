import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import cm
import scipy.interpolate

N = 100
files = glob.glob('Particles/*.txt')

for file in files:
	# read data from text and save it into array data
	data = np.genfromtxt(file)
	# open figure
	fig = plt.figure()
	# creates contour plot of data array automatically. vmin and vamx sets values for colorbar
	plt.plot(data[0][1], data[0][2], marker='o', color = 'r')
	# set labels
	plt.xlabel("X")
	plt.ylabel("Y")
	# set axis 
	plt.xlim([0,32])
	plt.ylim([0,32])
	plt.xticks(np.arange(0,32 + 1,4))
	plt.yticks(np.arange(0,32 + 1,4))
	
	plt.grid(linestyle='-', color='red')
	# define filename for saving
	filename = file.replace(".txt", "")
	filename = filename.replace("Particles/", "")
	fig.savefig("png/" + "{}.png".format(filename), bbox_inches='tight')
	# close fig
	plt.close(fig)
	

