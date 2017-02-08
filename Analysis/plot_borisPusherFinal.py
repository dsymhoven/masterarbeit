import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import fnmatch
import os

simulationInfo = np.genfromtxt('simulationInfo.txt')
numberOfParticles = int(simulationInfo[0])
startTime = int(simulationInfo[1])
numberOfParticleFiles = len(fnmatch.filter(os.listdir('Particles/Particle0/'), '*.txt'))
gridParameters = np.genfromtxt('gridParameters.txt')

lengthOfSimulationBoxInX = gridParameters[6]
lengthOfSimulationBoxInY = gridParameters[7]
lengthOfOneBoxInX = gridParameters[0] * gridParameters[3]
lengthOfOneBoxInY = gridParameters[1] * gridParameters[4]
numberOfBoxesInX =  lengthOfSimulationBoxInX / lengthOfOneBoxInX
numberOfBoxesInY =  lengthOfSimulationBoxInY / lengthOfOneBoxInY

X = np.zeros((numberOfParticles,1))
Y = np.zeros((numberOfParticles,1))
x = []
y = []
xnew = []
ynew = []

for i in range(startTime, startTime + numberOfParticleFiles):
	for p in range(numberOfParticles):
		# read data from text and save it into array data
		data = np.genfromtxt('Particles/Particle'+ str(p) +'/Particle' + str(p) + '_' + str(i) + '.txt')
		# define variables
		x.append(data[0][1])
		y.append(data[0][2])
	X = np.c_[X,x]
	Y = np.c_[Y,y]
	x=[]
	y=[]
	if i == startTime: #or len(X[0]) > 40:
		X = np.delete(X,0,1)
		Y = np.delete(Y,0,1)

# open figure
fig = plt.figure()
for p in range(numberOfParticles):
	# plot x and y value of particle as red dot
	plt.plot(X[p], Y[p], color = 'r')
# set labels
plt.xlabel("X")
plt.ylabel("Y")
# set axis
plt.axis('equal')
plt.xlim([0, lengthOfSimulationBoxInX])
plt.ylim([0, lengthOfSimulationBoxInY])
plt.xticks(np.arange(0, lengthOfSimulationBoxInX + 1, lengthOfOneBoxInX))
plt.yticks(np.arange(0, lengthOfSimulationBoxInY + 1, lengthOfOneBoxInY))
plt.grid(linestyle = "-", color='red')

# define filename for saving
filename = 'img'
fig.savefig("png/" + "{}.png".format(filename), bbox_inches='tight', dpi = 300)
# close fig
plt.close(fig)