//
//  experiment.c
//  masterarbeit
//
//  Created by David Symhoven on 02.11.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#include "experiment.h"
#include "stdlib.h"
#include "math.h"
#include "stdio.h"
#include "grid.h"


void testMaxwellPusher(){
    
    
// ======================================================
#pragma mark: Initializations
// ======================================================
    Grid Grid;
    int numberOfGridPointsInX = 256;
    int numberOfGridPointsInY = 256;
    int numberOfGridPointsInZ = 256;
    int lengthOfSimulationBoxInX = 32;
    int lengthOfSimulationBoxInY = 32;
    int lengthOfSimulationBoxInZ = 32;
    
    initGrid(&Grid, numberOfGridPointsInX, numberOfGridPointsInY, numberOfGridPointsInZ, lengthOfSimulationBoxInX, lengthOfSimulationBoxInY, lengthOfSimulationBoxInZ);
    initFieldsOnGrid(&Grid);
    double dt = 0.5 * Grid.dx;
    double t = 0;
    double tEnd = 20;
    char filename[32] = "some";
    
// ======================================================
#pragma mark: Main Routine
// ======================================================
    for(int p = 0; p < tEnd / dt; p++){

        t += dt;
    }
    
    
    freeMemoryOnGrid(&Grid);
}
