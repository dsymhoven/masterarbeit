import numpy as np
import os
import shutil
import glob

numberOfParticles = np.genfromtxt('numberOfParticles.txt')
for i in range(numberOfParticles):
	path = 'Particles/Particle' + str(i)
	try: 
	    os.makedirs(path)
	except OSError:
	    if not os.path.isdir(path):
	        raise
	particles = glob.glob('Particles/Particle'+str(i)+'_*.txt')
	for file in particles:
		shutil.copy(file,path)
		os.remove(file)
		
os.remove('numberOfParticles.txt')