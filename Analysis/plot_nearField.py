import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import fnmatch
import os

numberOfParticleFiles = len(fnmatch.filter(os.listdir('Particles/'), '*.txt'))

#files = glob.glob('Particles/*.txt')
gridParameters = np.genfromtxt('gridParameters.txt')
field = np.genfromtxt('E_fields/E_field0.txt')

lengthOfSimulationBoxInX = gridParameters[6]
lengthOfSimulationBoxInY = gridParameters[7]
lengthOfOneBoxInX = gridParameters[0] * gridParameters[3]
lengthOfOneBoxInY = gridParameters[1] * gridParameters[4] 
x=[]
y=[]

for i in range(numberOfParticleFiles):
	# read data from text and save it into array data
	data = np.genfromtxt('Particles/Particle'+ str(i) +'.txt')
	# open figure
	fig = plt.figure()
	# define variables
	x.append(data[0][1])
	y.append(data[0][2])
	# plot x and y value of particle as red dot
	plt.plot(x, y, color = 'r')
	plt.imshow(field, vmin=0, vmax=1.0)
	plt.colorbar()
	
	currentAxis = plt.gca()
	# define variables
	xMin = data[2][0]
	xMax = data[2][1]
	yMin = data[2][2]
	yMax = data[2][3]
	
	width = xMax - xMin
	heigth = yMax - yMin
	
	rect = patches.Rectangle(
        (xMin, yMin),	
        width,			# width
        heigth,			# heigth
        linewidth = 1,
        edgecolor = 'r',
        fill=False  # remove background
    )
	currentAxis.add_patch(rect)
    
	# set labels
	plt.xlabel("X")
	plt.ylabel("Y")
	# set axis 
	plt.xlim([0, lengthOfSimulationBoxInX])
	plt.ylim([0, lengthOfSimulationBoxInY])
	plt.xticks(np.arange(0, lengthOfSimulationBoxInX + 1, lengthOfOneBoxInX))
	plt.yticks(np.arange(0, lengthOfSimulationBoxInY + 1, lengthOfOneBoxInY))
	
	plt.grid(color='red')
	# define filename for saving
	filename = 'Particle' + str(i)
	fig.savefig("png/" + "{}.png".format(filename), bbox_inches='tight')
	# close fig
	plt.close(fig)
	

