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
    int numberOfGridPoints;
    double lengthOfSimulationBox;
    double *E;
    double *B;
    
};

typedef struct Grid Grid;

void initGrid(Grid Grid, int numberOfGridPoints, double lengthOfSimulationBox);
void allocateFieldsOn(Grid Grid);
//void calcualteNearFieldBoxes(double x[4], int lengthOfSimulationBox, int numberOfGridPoints, double *edgeOfNearFieldBox);

#endif /* grid_h */
