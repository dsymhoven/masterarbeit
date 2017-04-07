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
#include "particle.h"
#include "resolution.h"
#include "box.h"

struct Grid {
    int numberOfBoxesInX;
    int numberOfBoxesInY;
    int numberOfBoxesInZ;

    Box Box;
    Resolution Resolution;
    
    int numberOfGridPointsInX;
    int numberOfGridPointsInY;
    int numberOfGridPointsInZ;
    
    double lengthOfSimulationBoxInX;
    double lengthOfSimulationBoxInY;
    double lengthOfSimulationBoxInZ;
    
    double upmlLayerWidth;
    
    double EMax;
    double HMax;
    
    double *E;
    double *B;
    
    double *D;
    double *H;
    
    double **Hz_im1;
    double **Hy_im1;
    double **Hx_jm1;
    double **Hz_jm1;
    double **Hx_km1;
    double **Hy_km1;
    double **Ey_ip1;
    double **Ez_ip1;
    double **Ez_jp1;
    double **Ex_jp1;
    double **Ex_kp1;
    double **Ey_kp1;
    
    double *upml1E;
    double *upml2E;
    double *upml3E;
    double *upml4E;
    double *upml5E;
    double *upml6E;
    
    double *upml1H;
    double *upml2H;
    double *upml3H;
    double *upml4H;
    double *upml5H;
    double *upml6H;
    
};

typedef struct Grid Grid;
typedef struct Particle Particle;

void initGrid(Grid *Grid, Resolution *Resolution, Box *Box, const int numberOfBoxesInX, const int numberOfBoxesInY, const int numberOfBoxesInZ);
void allocateMemoryOnGrid(Grid *Grid);
void freeMemoryOnGrid(Grid *Grid);
void initSamplePulseOnGrid(Grid *Grid);
void pushEFieldOnGrid(Grid *Grid, const double dt);
void pushHFieldOnGrid(Grid *Grid, const double dt);
void externalPlaneWave(const double xParticle[4], const double tStart, double Eextern[3], double Bextern[3]);
void clearFieldsFromGrid(Grid *Grid);
#endif /* grid_h */
