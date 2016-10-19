//
//  grid.c
//  masterarbeit
//
//  Created by David Symhoven on 13.10.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#include "grid.h"

void calcualteNearFieldBoxes(double *x, int lengthOfSimulationBox, int numberOfGridPoints, double *edgeOfNearFieldBox){
    int i, j, k;
    double numberOfBoxesPerDimension = numberOfGridPoints / lengthOfSimulationBox;
    double lengthOfOneBox = lengthOfSimulationBox / numberOfBoxesPerDimension;
    
    i = x[1] / lengthOfOneBox;
    j = x[2] / lengthOfOneBox;
    k = x[3] / lengthOfOneBox;
    
    double xMin = (i - 1) * lengthOfOneBox;
    double xMax = (i + 2) * lengthOfOneBox;
    double yMin = (j - 1) * lengthOfOneBox;
    double yMax = (j + 2) * lengthOfOneBox;
    double zMin = (k - 1) * lengthOfOneBox;
    double zMax = (k + 2) * lengthOfOneBox;
    
    edgeOfNearFieldBox[0] = xMin;
    edgeOfNearFieldBox[1] = yMin;
    edgeOfNearFieldBox[2] = zMin;
    
    edgeOfNearFieldBox[3] = xMax;
    edgeOfNearFieldBox[4] = yMin;
    edgeOfNearFieldBox[5] = zMin;
    
    edgeOfNearFieldBox[6] = xMax;
    edgeOfNearFieldBox[7] = yMax;
    edgeOfNearFieldBox[8] = zMin;
    
    edgeOfNearFieldBox[9] = xMin;
    edgeOfNearFieldBox[10] = yMax;
    edgeOfNearFieldBox[11] = zMin;
    
    edgeOfNearFieldBox[12] = xMin;
    edgeOfNearFieldBox[13] = yMin;
    edgeOfNearFieldBox[14] = zMax;
    
    edgeOfNearFieldBox[15] = xMax;
    edgeOfNearFieldBox[16] = yMin;
    edgeOfNearFieldBox[17] = zMax;
    
    edgeOfNearFieldBox[18] = xMax;
    edgeOfNearFieldBox[19] = yMax;
    edgeOfNearFieldBox[20] = zMax;
    
    edgeOfNearFieldBox[21] = xMin;
    edgeOfNearFieldBox[22] = yMax;
    edgeOfNearFieldBox[23] = zMax;
    
    
}
