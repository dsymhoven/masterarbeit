#!/bin/bash
# define variables
PATHTOANALYSIS=~/Desktop/Projects/masterarbeit/Analysis
PATHTOEXECUTABLE=~/Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug
# 1. create folder for all fields.txt files in Analysis folder of masterarbeit directory
mkdir -p $PATHTOANALYSIS/fields
# 2. copy created files from /Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug/ to masterarbeit folder
mv $PATHTOEXECUTABLE/*.txt $PATHTOANALYSIS/fields
# 3. create png folder in Analysis folder. Needed for python script
mkdir -p $PATHTOANALYSIS/png
#3. change directory to Analysis folder in order for python script to work properly
cd $PATHTOANALYSIS
#4. execute python script
python2.7 plot.py

