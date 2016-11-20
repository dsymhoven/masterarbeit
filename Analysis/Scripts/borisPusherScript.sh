#!/bin/bash
# define variables
PATHTOANALYSIS=~/Desktop/Projects/masterarbeit/Analysis
PATHTOEXECUTABLE=~/Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug
 
# check if file with name "B_field0.txt" exists in /Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug/
# if so, create new folder "B_Fields" in Analysis folder of masterarbeit directory and copy all B_field files to that folder
# analogue with E field and Particle files
if [ -e B_field0.txt ]
then
	mkdir -p $PATHTOANALYSIS/B_fields
	mv $PATHTOEXECUTABLE/B_field*.txt $PATHTOANALYSIS/B_Fields
fi
if [ -e E_field0.txt ]
then
	mkdir -p $PATHTOANALYSIS/E_fields
	mv $PATHTOEXECUTABLE/E_field*.txt $PATHTOANALYSIS/E_Fields
fi
if [ -e Particle0.txt ]
then
	mkdir -p $PATHTOANALYSIS/Particles
	mv $PATHTOEXECUTABLE/Particle*.txt $PATHTOANALYSIS/Particles
fi

# create png folder in Analysis folder. Neccessary for python script
mkdir -p $PATHTOANALYSIS/png

# change directory to Analysis folder in order for python script to work properly
cd $PATHTOANALYSIS
# execute python script
python2.7 plot_borisPusher.py

# chane to png folder
cd png/
# execute ffmpeg to create video file
~/../../opt/local/bin/ffmpeg -framerate 4 -start_number 0 -i Particle%d.png -vf "scale=720:trunc(ow/a/2)*2" -c:v libx264 -r 30 -pix_fmt yuv420p movie.mp4
# create new folder for Video file
cd ../
mkdir -p Movies/
# copy video into new folder
mv png/movie.mp4 Movies
# clean up
rm -rf png/