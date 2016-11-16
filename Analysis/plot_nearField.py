import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as patches


files = glob.glob('Particles/*.txt')
gridParameters = np.genfromtxt('gridParameters.txt')
field = np.genfromtxt('E_fields/E_field0.txt')

lengthOfSimulationBoxInX = gridParameters[0]
lengthOfSimulationBoxInY = gridParameters[1]
lengthOfOneBoxInX = gridParameters[3]
lengthOfOneBoxInY = gridParameters[4] 


for file in files:
	# read data from text and save it into array data
	data = np.genfromtxt(file)
	# open figure
	fig = plt.figure()
	# define variables
	vx = data[0][1]
	vy = data[0][2]
	# plot x and y value of particle as red dot
	plt.plot(vx, vy, marker='o', color = 'r')
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
	
	#plt.grid(linestyle='-', color='red')
	# define filename for saving
	filename = file.replace(".txt", "")
	filename = filename.replace("Particles/", "")
	fig.savefig("png/" + "{}.png".format(filename), bbox_inches='tight')
	# close fig
	plt.close(fig)
	

