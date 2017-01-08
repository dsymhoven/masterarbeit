import numpy as np
import os
import glob  
import shutil


pathToAnalysis = os.path.expanduser('~/Desktop/Projects/masterarbeit/Analysis')
pathToExecutable = os.path.expanduser('~/Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug/')
numberOfDirectories = 0

os.chdir(pathToAnalysis)
if not os.path.exists('initialFields/'):
    os.makedirs('initialFields/')
    
f = open("initialFields/numberOfDirectories.txt","w")
    
folders = glob.glob("initialFields/*")
numberOfDirectories = len(folders) - 1
f.write(str(numberOfDirectories))
f.close()

if not os.path.exists('initialFields/' + str(numberOfDirectories)):
    os.makedirs('initialFields/' + str(numberOfDirectories))
    
if os.path.exists(pathToExecutable + 'E_initialField.txt'):
	shutil.copy(pathToExecutable + 'E_initialField.txt','initialFields/' + str(numberOfDirectories))
	os.remove(pathToExecutable + 'E_initialField.txt')
if os.path.exists(pathToExecutable + 'H_initialField.txt'):
	shutil.copy(pathToExecutable + 'H_initialField.txt','initialFields/' + str(numberOfDirectories))
	os.remove(pathToExecutable + 'H_initialField.txt')
if os.path.exists(pathToExecutable + 'initialConditions.txt'):
	shutil.copy(pathToExecutable + 'initialConditions.txt','initialFields/' + str(numberOfDirectories))
	os.remove(pathToExecutable + 'initialConditions.txt')