#!/bin/bash
# define variables
PATHTOANALYSIS=~/Desktop/Projects/masterarbeit/Analysis
PATHTOEXECUTABLE=~/Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug
# 1. create folders for all .txt files in Analysis folder of masterarbeit directory
  
# 2. copy created files from /Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug/ to masterarbeit folder
if [ -e B_field0.txt ]
then
	mkdir -p $PATHTOANALYSIS/E_fields
	mv $PATHTOEXECUTABLE/B_field*.txt $PATHTOANALYSIS/B_Fields
fi
if [ -e E_field0.txt ]
then
	mkdir -p $PATHTOANALYSIS/B_fields
	mv $PATHTOEXECUTABLE/E_field*.txt $PATHTOANALYSIS/E_Fields
fi
if [ -e Particle0.txt ]
then
	mkdir -p $PATHTOANALYSIS/Particles
	mv $PATHTOEXECUTABLE/Particle*.txt $PATHTOANALYSIS/Particles
fi
# 3. create png folder in Analysis folder. Needed for python script
mkdir -p $PATHTOANALYSIS/png
# 4. change directory to Analysis folder in order for python script to work properly
cd $PATHTOANALYSIS
# 5. execute python script
python2.7 plot_particles.py
# 6. chane to png folder
cd png/
# 7. execute ffmpeg to create video file
~/../../opt/local/bin/ffmpeg -framerate 4 -start_number 0 -i Particle%d.png -vf "scale=720:trunc(ow/a/2)*2" -c:v libx264 -r 30 -pix_fmt yuv420p movie.mp4
# 8. create new folder for Video file
cd ../
mkdir -p Movies/
# 10. copy video into new folder
mv png/movie.mp4 Movies
# clean up
rm -rf png/ #fields/