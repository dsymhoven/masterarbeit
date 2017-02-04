#!/bin/bash
# define variables
PATHTOANALYSIS=~/Desktop/Projects/masterarbeit/Analysis
PATHTOEXECUTABLE=~/Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug

# check if file with name "B_field0.txt" exists in /Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug/
# if so, create new folder "B_Fields" in Analysis folder of masterarbeit directory and copy all B_field files to that folder
# analogue with E field and Particle files
# if gridParameters were written out, copy them into Analysis folder. folder
for file in B_field*; do
	if [[ -f $file ]]; then
		rm -rf $PATHTOANALYSIS/B_fields
		mkdir -p $PATHTOANALYSIS/B_fields
		mv $PATHTOEXECUTABLE/B_field*.txt $PATHTOANALYSIS/B_Fields
	fi
done
for file in E_field*; do
	if [[ -f $file ]]; then
		rm -rf $PATHTOANALYSIS/E_fields
		mkdir -p $PATHTOANALYSIS/E_fields
		mv $PATHTOEXECUTABLE/E_field*.txt $PATHTOANALYSIS/E_Fields
	fi
done
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

# create png folder in Analysis folder. Neccessary for python script
mkdir -p $PATHTOANALYSIS/png

# change directory to Analysis folder in order for python script to work properly
cd $PATHTOANALYSIS
# execute python script
python2.7 sortParticleFiles.py
python2.7 plot_borisPusherFinal.py

# chane to png folder
#cd png/
# execute ffmpeg to create video file
#~/../../opt/local/bin/ffmpeg -framerate 4 -start_number 0 -i img%d.png -vf "scale=720:trunc(ow/a/2)*2" -c:v libx264 -r 30 -pix_fmt yuv420p borisPusher.mp4
# create new folder for Video file
#cd ../
#mkdir -p Movies/
# copy video into new folder
#mv png/borisPusher.mp4 Movies
