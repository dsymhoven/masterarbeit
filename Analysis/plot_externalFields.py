import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import fnmatch
import os

numberOfFieldFiles = len(fnmatch.filter(os.listdir('E_extern/'), '*.txt'))
gridParameters = np.genfromtxt('gridParameters.txt')

lengthOfSimulationBoxInX = gridParameters[6]
lengthOfSimulationBoxInY = gridParameters[7]
lengthOfOneBoxInX = gridParameters[0] * gridParameters[3]
lengthOfOneBoxInY = gridParameters[1] * gridParameters[4]
numberOfBoxesInX =  lengthOfSimulationBoxInX / lengthOfOneBoxInX
numberOfBoxesInY =  lengthOfSimulationBoxInY / lengthOfOneBoxInY


for i in range(numberOfFieldFiles):
	fig = plt.figure()
	field = np.genfromtxt('E_extern/E_extern'+ str(i) +'.txt')
	# plot fields
	plt.imshow(field, aspect='auto', origin='lower', cmap = 'jet', extent=(0,lengthOfSimulationBoxInX,0,lengthOfSimulationBoxInY), vmin=0, vmax=1)
	plt.colorbar()
	# set labels
	plt.xlabel("X")
	plt.ylabel("Y")
	# set axis
	plt.xlim([0, lengthOfSimulationBoxInX])
	plt.ylim([0, lengthOfSimulationBoxInY])
	plt.xticks(np.arange(0, lengthOfSimulationBoxInX + 1, lengthOfOneBoxInX))
	plt.yticks(np.arange(0, lengthOfSimulationBoxInY + 1, lengthOfOneBoxInY))
	plt.grid(linestyle = "-", color='red')
	# define filename for saving
	filename = 'img' + str(i)
	fig.savefig("png/" + "{}.png".format(filename), bbox_inches='tight')
	# close fig
	plt.close(fig)
