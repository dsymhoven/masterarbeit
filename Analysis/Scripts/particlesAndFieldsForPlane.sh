#!/bin/bash
# define variables
PATHTOANALYSIS=~/Desktop/Projects/masterarbeit/Analysis
PATHTOEXECUTABLE=~/Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug

# check if file with name "B_field0.txt" exists in /Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug/
# if so, create new folder "B_Fields" in Analysis folder of masterarbeit directory and copy all B_field files to that folder
# analogue with E field and Particle files
# if gridParameters were written out, copy them into Analysis folder. folder
if [ -e B_field0.txt ]
then
	rm -rf $PATHTOANALYSIS/B_fields
	mkdir $PATHTOANALYSIS/B_fields
	mv $PATHTOEXECUTABLE/B_field*.txt $PATHTOANALYSIS/B_Fields
fi
if [ -e E_field0.txt ]
then
	rm -rf $PATHTOANALYSIS/E_fields
	mkdir -p $PATHTOANALYSIS/E_fields
	mv $PATHTOEXECUTABLE/E_field*.txt $PATHTOANALYSIS/E_Fields
fi
for file in Particle*; do
	if [[ -f $file ]]; then
		rm -rf $PATHTOANALYSIS/Particles
		mkdir -p $PATHTOANALYSIS/Particles
		mv $PATHTOEXECUTABLE/Particle*.txt $PATHTOANALYSIS/Particles
	fi
done
if [ -e gridParameters.txt ]
then
	mv $PATHTOEXECUTABLE/gridParameters.txt $PATHTOANALYSIS/
fi
if [ -e simulationInfo.txt ]
then
	mv $PATHTOEXECUTABLE/simulationInfo.txt $PATHTOANALYSIS/
fi
#  create png folder in Analysis folder. Neccessary for python script
mkdir -p $PATHTOANALYSIS/png
# change directory to Analysis folder in order for python script to work properly
cd $PATHTOANALYSIS
# execute python script
python2.7 sortParticleFiles.py
python2.7 plot_particlesAndFieldsForPlane.py

