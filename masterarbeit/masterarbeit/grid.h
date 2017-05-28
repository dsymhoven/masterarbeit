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

    Box Box;
    Resolution Resolution;
    
    int numberOfBoxesInX;
    int numberOfBoxesInY;
    int numberOfBoxesInZ;
    
    int numberOfGridPointsInX;
    int numberOfGridPointsInY;
    int numberOfGridPointsInZ;
    
    double lengthOfSimulationBoxInX;
    double lengthOfSimulationBoxInY;
    double lengthOfSimulationBoxInZ;
    
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
    
    double *upml1Ex;
    double *upml2Ex;
    double *upml3Ex;
    double *upml4Ex;
    double *upml5Ex;
    double *upml6Ex;
    
    double *upml1Ey;
    double *upml2Ey;
    double *upml3Ey;
    double *upml4Ey;
    double *upml5Ey;
    double *upml6Ey;
    
    double *upml1Ez;
    double *upml2Ez;
    double *upml3Ez;
    double *upml4Ez;
    double *upml5Ez;
    double *upml6Ez;
    
    double *upml1Hx;
    double *upml2Hx;
    double *upml3Hx;
    double *upml4Hx;
    double *upml5Hx;
    double *upml6Hx;
    
    double *upml1Hy;
    double *upml2Hy;
    double *upml3Hy;
    double *upml4Hy;
    double *upml5Hy;
    double *upml6Hy;
    
    double *upml1Hz;
    double *upml2Hz;
    double *upml3Hz;
    double *upml4Hz;
    double *upml5Hz;
    double *upml6Hz;
    
    double upmlLayerWidth;
    
    bool useUPML;
    
    
};

typedef struct Grid Grid;
typedef struct Particle Particle;

void initGrid(Grid *Grid, Resolution *Resolution, Box *Box, const int numberOfBoxesInX, const int numberOfBoxesInY, const int numberOfBoxesInZ, bool useUPML);
void freeMemoryOnGrid(Grid *Grid);
void initSamplePulseOnGrid(Grid *Grid);
void pushEFieldOnGrid(Grid *Grid, const double dt);
void pushHFieldOnGrid(Grid *Grid, const double dt);
void externalPlaneWave(const double xParticle[4], const double tStart, double Eextern[3], double Bextern[3]);
void externalPulse(const double x[4], const double tStart, double Eextern[3], double Hextern[3], Grid *Grid);
void clearFieldsFromGrid(Grid *Grid);
#endif /* grid_h */
