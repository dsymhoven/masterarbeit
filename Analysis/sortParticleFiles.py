import numpy as np
import os
import shutil
import glob

simulationInfo = np.genfromtxt('simulationInfo.txt')
numberOfParticles = int(simulationInfo[0])
for p in range(numberOfParticles):
	path = 'Particles/Particle' + str(p)
	try: 
	    os.makedirs(path)
	except OSError:
	    if not os.path.isdir(path):
	        raise
	particles = glob.glob('Particles/Particle'+str(p)+'_*.txt')
	for file in particles:
		shutil.copy(file,path)
		os.remove(file)
		