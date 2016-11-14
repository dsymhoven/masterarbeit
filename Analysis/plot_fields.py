import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import cm
import scipy.interpolate

N = 100
files = glob.glob('E_fields/*.txt')

for file in files:
	# read data from text and save it into array data
	data = np.genfromtxt(file)
	# open figure
	fig = plt.figure()
	# creates contour plot of data array automatically. vmin and vamx sets values for colorbar
	plt.imshow(data, vmin=0, vmax=1.0)
	plt.colorbar()
	# set labels
	plt.xlabel("X")
	plt.ylabel("Y")
	# define filename for saving
	filename = file.replace(".txt", "")
	filename = filename.replace("fields/", "")
	fig.savefig("png/" + "{}.png".format(filename), bbox_inches='tight')
	# close fig
	plt.close(fig)
