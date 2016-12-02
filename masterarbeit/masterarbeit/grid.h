//
//  grid.h
//  masterarbeit
//
//  Created by David Symhoven on 13.10.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#ifndef grid_h
#define grid_h

#include <stdio.h>
#include "stdbool.h"

struct Grid {
    int numberOfBoxesInX;
    int numberOfBoxesInY;
    int numberOfBoxesInZ;
    
    int numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ;
    
    int numberOfGridPointsInX;
    int numberOfGridPointsInY;
    int numberOfGridPointsInZ;
    
    double lengthOfSimulationBoxInX;
    double lengthOfSimulationBoxInY;
    double lengthOfSimulationBoxInZ;
    
    double dx;
    double dy;
    double dz;
    
    double *E;
    double *B;
    
    double **Bz_im1;
    double **By_im1;
    double **Bx_jm1;
    double **Bz_jm1;
    double **Bx_km1;
    double **By_km1;
    double **Ey_ip1;
    double **Ez_ip1;
    double **Ez_jp1;
    double **Ex_jp1;
    double **Ex_kp1;
    double **Ey_kp1;
    
};

typedef struct Grid Grid;

void initGrid(Grid *Grid, double dx, double dy, double dz, int numberOfGridPointsForBoxInX, int numberOfGridPointsForBoxInY, int numberOfGridPointsForBoxInZ, int numberOfBoxesInX, int numberOfBoxesInY, int numberOfBoxesInZ);
void allocateMemoryOnGrid(Grid *Grid);
void freeMemoryOnGrid(Grid *Grid);
void initSamplePulseOnGrid(Grid *Grid);
void pushEFieldOnGrid(Grid *Grid, double dt);
void pushBFieldOnGrid(Grid *Grid, double dt);
void writeFieldsToFile(Grid *Grid, char *filename, int index, int planeForPlotting, bool plotE, bool plotB);
void writeGridParametersToFile(Grid *Grid);
void pushEFieldInsideBoxes(Grid *Grid, double dt);
void pushBFieldInsideBoxes(Grid *Grid, double dt);
void allocateFieldsOnBoxBorders(Grid *Grid);
#endif /* grid_h */
