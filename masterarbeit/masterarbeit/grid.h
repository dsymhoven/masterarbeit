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

void initGrid(Grid *Grid, double dx, double dy, double dz, int numberOfGridPointsForBoxInX, int numberOfGridPointsForBoxInY, int numberOfGridPointsForBoxInZ, int numberOfBoxesInX, int numberOfBoxesInY, int numberOfBoxesInZ);
void allocateMemoryOnGrid(Grid *Grid);
void freeMemoryOnGrid(Grid *Grid);
void initSamplePulseOnGrid(Grid *Grid);
void pushEFieldOnGrid(Grid *Grid, double dt);
void pushBFieldOnGrid(Grid *Grid, double dt);
void writeFieldsToFile(Grid *Grid, char *filename, int index, int planeForPlotting, bool plotE, bool plotB);
void writeGridParametersToFile(Grid *Grid);
void pushEFieldInsideBoxes(Grid *Grid, double dt);
void pushHFieldInsideBoxes(Grid *Grid, double dt);
void allocateFieldsOnBoxBorders(Grid *Grid);
void adjustHFields(Grid *Grid, Particle *Particle, const double t);
void adjustHyz_im1(Grid *Grid, Particle *Particle, const int boxIndex, const int ib, const int jb, const int kb, const double t);
void adjustHxz_jm1(Grid *Grid, Particle *Particle, const int boxIndex, const int ib, const int jb, const int kb, const double t);
void adjustHxy_km1(Grid *Grid, Particle *Particle, const int boxIndex, const int ib, const int jb, const int kb, const double t);
void adjustEFields(Grid *Grid, Particle *Particle, const double t);
void adjustEyz_ip1(Grid *Grid, Particle *Particle, const int boxIndex, const int ib, const int jb, const int kb, const double t);
void adjustExz_jp1(Grid *Grid, Particle *Particle, const int boxIndex, const int ib, const int jb, const int kb, const double t);
void adjustExy_kp1(Grid *Grid, Particle *Particle, const int boxIndex, const int ib, const int jb, const int kb, const double t);
void pushHFieldAtBorders(Grid *Grid, double dt);
void pushEFieldAtBorders(Grid *Grid, double dt);
void setHFieldOnBorders(Grid *Grid);
void setEFieldOnBorders(Grid *Grid);
void pushHField(Grid *Grid, Particle *Particle, double t, double dt);
void pushEField(Grid *Grid, Particle *Particle, double t, double dt);
void allocateUPMLCoefficients(Grid *Grid);
#endif /* grid_h */
