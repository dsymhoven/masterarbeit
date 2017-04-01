import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import fnmatch
import os

simulationInfo = np.genfromtxt('simulationInfo.txt')
numberOfParticleFiles = len(fnmatch.filter(os.listdir('Particles/Particle0'), '*.txt'))
numberOfParticles = int(simulationInfo[0])
startTime = int(simulationInfo[1])

#files = glob.glob('Particles/*.txt')
gridParameters = np.genfromtxt('gridParameters.txt')
field = np.genfromtxt('E_fields/E_field0.txt')

lengthOfSimulationBoxInX = gridParameters[6]
lengthOfSimulationBoxInY = gridParameters[7]
lengthOfOneBoxInX = gridParameters[0] * gridParameters[3]
lengthOfOneBoxInY = gridParameters[1] * gridParameters[4]
X = np.zeros((numberOfParticles,1))
Y = np.zeros((numberOfParticles,1))
x=[]
y=[]

for i in range(startTime, startTime + numberOfParticleFiles):
	# open figure
	fig = plt.figure()
	for p in range(numberOfParticles):
		# read data from text and save it into array data
		data = np.genfromtxt('Particles/Particle'+ str(p) +'/Particle' + str(p) + '_' + str(i) + '.txt')
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
		# define variables
		x.append(data[0][1])
		y.append(data[0][2])
	X = np.c_[X,x]
	Y = np.c_[Y,y]
	x=[]
	y=[]
	if i == startTime or len(X[0]) > 40:
		X = np.delete(X,0,1)
		Y = np.delete(Y,0,1)
	for p in range(numberOfParticles):
		# plot x and y value of particle as red dot
		plt.plot(X[p], Y[p], color = 'r')

	# plot fields
	plt.imshow(field, aspect='auto', origin='lower', cmap = 'jet', extent=(0,lengthOfSimulationBoxInX,0,lengthOfSimulationBoxInY), vmin=0, vmax=0.1)
	plt.colorbar()
	# set labels
	plt.xlabel("X")
	plt.ylabel("Y")
	plt.text(14,15,'$e^{-}$',color='red')
	plt.text(22,23,'$e^{-}$',color='red')
	# set axis
	plt.xlim([0, lengthOfSimulationBoxInX])
	plt.ylim([0, lengthOfSimulationBoxInY])
	plt.xticks(np.arange(0, lengthOfSimulationBoxInX + 1, lengthOfOneBoxInX))
	plt.yticks(np.arange(0, lengthOfSimulationBoxInY + 1, lengthOfOneBoxInY))
	plt.grid(color='red')
	# define filename for saving
	filename = 'img' + str(i - startTime)
	fig.savefig("png/" + "{}.png".format(filename), bbox_inches='tight', dpi=300)
	# close fig
	plt.close(fig)
