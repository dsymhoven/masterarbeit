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

files = glob.glob('Particles/*.txt')
gridParameters = np.genfromtxt('gridParameters.txt')
field = np.genfromtxt('E_fields/E_field0.txt')

lengthOfSimulationBoxInX = gridParameters[6]
lengthOfSimulationBoxInY = gridParameters[7]
lengthOfOneBoxInX = gridParameters[0] * gridParameters[3]
lengthOfOneBoxInY = gridParameters[1] * gridParameters[4]
numberOfBoxesInX =  lengthOfSimulationBoxInX / lengthOfOneBoxInX
numberOfBoxesInY =  lengthOfSimulationBoxInY / lengthOfOneBoxInY
if gridParameters[9] >= 0.1:
	EMax = round(gridParameters[9],1)
if gridParameters[9] < 0.1:
	EMax = 0.05
if gridParameters[9] > 10:
	EMax = 0.5

x=[]
y=[]

# open figure only once because we want to plot several particles in one figure
fig = plt.figure()

for p in range(numberOfParticles):
	for i in range(startTime, startTime + numberOfParticleFiles):
		# read data from text and save it into array data
		data = np.genfromtxt('Particles/Particle'+ str(p) +'/Particle' + str(p) +'_' + str(i) + '.txt')
		# define variables
		x.append(data[0][1])
		y.append(data[0][2])
		# if len(x) > 40:
		# 	x.pop(0)
		# 	y.pop(0)
	# plot x and y value of particle as red dot
	plt.plot(x, y, color = 'r')
	# delete x,y array for next particle
	x=[]
	y=[]

# plot fields
plt.imshow(field, aspect='auto', origin='lower', cmap = 'jet', extent=(0,lengthOfSimulationBoxInX,0,lengthOfSimulationBoxInY), vmin=0, vmax=0.04)
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
filename = 'img'
fig.savefig("png/" + "{}.png".format(filename),bbox_inches='tight',dpi=300)
# close fig
plt.close(fig)
