import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import fnmatch
import os

numberOfFieldFiles = len(fnmatch.filter(os.listdir('E_fields/'), '*.txt'))
gridParameters = np.genfromtxt('gridParameters.txt')

lengthOfSimulationBoxInX = gridParameters[6]
lengthOfSimulationBoxInY = gridParameters[7]

for i in range(0,numberOfFieldFiles):
	# read data from text and save it into array data
    field = np.genfromtxt('E_fields/E_field' + str(i) + '.txt')
	#data = np.genfromtxt(file)
    # open figure
    fig = plt.figure()
    # creates contour plot of data array automatically. vmin and vamx sets values for colorbar
    plt.imshow(field, aspect='auto', origin='lower', cmap = 'jet', extent=(0,lengthOfSimulationBoxInX,0,lengthOfSimulationBoxInY), vmin=0, vmax=1.0)
    plt.colorbar()
    # set labels
    plt.xlabel("X")
    plt.ylabel("Y")
    # define filename for saving
    filename = 'img' + str(i)
    fig.savefig("png/" + "{}.png".format(filename), bbox_inches='tight')
    # close fig
    plt.close(fig)
