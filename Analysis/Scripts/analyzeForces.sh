#!/bin/bash
# define variables
PATHTOANALYSIS=~/Desktop/Projects/masterarbeit/Analysis
PATHTOEXECUTABLE=~/Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug

if [ -e dampingTermVSLorentzForce.txt ]
then
	mkdir -p $PATHTOANALYSIS/Forces
	mv $PATHTOEXECUTABLE/dampingTermVSLorentzForce.txt $PATHTOANALYSIS/Forces
fi

# change directory to Analysis folder in order for python script to work properly
cd $PATHTOANALYSIS
# execute python script
python2.7 plot_analyzeForces.py
