#!/bin/bash
# define variables
PATHTOANALYSIS=~/Desktop/Projects/masterarbeit/Analysis
PATHTOEXECUTABLE=~/Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug
# 1. create folder for all E_fields.txt files in Analysis folder of masterarbeit directory
for file in E_field*; do
	if [[ -f $file ]]; then
		rm -rf $PATHTOANALYSIS/E_fields
		mkdir -p $PATHTOANALYSIS/E_fields
		mv $PATHTOEXECUTABLE/E_field*.txt $PATHTOANALYSIS/E_Fields
	fi
done
if [ -e gridParameters.txt ]
then
	mv $PATHTOEXECUTABLE/gridParameters.txt $PATHTOANALYSIS/
fi
# 2. create png folder to store all image files
mkdir -p $PATHTOANALYSIS/png
# 4. change directory to Analysis folder in order for python script to work properly
cd $PATHTOANALYSIS
# 5. execute python script
python2.7 plot_maxwellPusher.py
# 6. chane to png folder
cd png/
# 7. execute ffmpeg to create video file
~/../../opt/local/bin/ffmpeg -framerate 4 -start_number 0 -i img%d.png -vf "scale=720:trunc(ow/a/2)*2" -c:v libx264 -r 30 -pix_fmt yuv420p maxwellPusher.mp4
# 8. create new folder for Video file
cd ../
mkdir -p Movies/
# 10. copy video into new folder
mv png/maxwellPusher.mp4 Movies
