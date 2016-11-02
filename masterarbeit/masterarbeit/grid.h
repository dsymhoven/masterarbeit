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

struct Grid {
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
    
};

typedef struct Grid Grid;

void initGrid(Grid *Grid, int numberOfGridPointsInX, int numberOfGridPointsInY, int numberOfGridPointsInZ, double lengthOfSimulationBoxInX, double lengthOfSimulationBoxInY, double lengthOfSimulationBoxInZ);
void allocateFieldsOnGrid(Grid *Grid);
void freeMemoryOnGrid(Grid *Grid);
void initFieldsOnGrid(Grid *Grid);
void pushEFieldOnGrid(Grid *Grid, double dt);
void PushBFieldOnGrid(Grid *Grid, double dt);
//void calcualteNearFieldBoxes(double x[4], int lengthOfSimulationBox, int numberOfGridPoints, double *edgeOfNearFieldBox);

#endif /* grid_h */
