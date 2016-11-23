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
#include "stdbool.h"
#include "particle.h"
#include "calculations.h"


void testMaxwellPusher(){
    
    
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    Grid Grid;
    double dx = 0.125;
    double dy = 0.125;
    double dz = 0.125;
    int numberOfGridPointsForBoxInX = 32;
    int numberOfGridPointsForBoxInY = 32;
    int numberOfGridPointsForBoxInZ = 32;
    int numberOfBoxesInX = 8;
    int numberOfBoxesInY = 8;
    int numberOfBoxesInZ = 8;
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateFieldsOnGrid(&Grid);
    initSamplePulseOnGrid(&Grid);
    
    double dt = 0.5 * Grid.dx;
    double t = 0;
    double tEnd = 3;
    char filename[32] = "some";
    int arrayLength = tEnd/dt;

    int planeForPlotting = 124;
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    
    for(int p = 0; p < arrayLength; p++){
        writeFieldsToFile(&Grid, filename, p, planeForPlotting, true, false);
        pushEFieldOnGrid(&Grid, dt);
        PushBFieldOnGrid(&Grid, dt);
        PushBFieldOnGrid(&Grid, dt);
        pushEFieldOnGrid(&Grid, dt);
        t += dt;
    }
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/maxwellPusherScript.sh");
    freeMemoryOnGrid(&Grid);
}


void testBorisPusher(){
    
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Particle Particle;
    double charge = 1.0;
    double mass = 1.0;
    
    double dt = 0.5 * 0.125;
    double t = 0;
    double tEnd = 3;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initParticle(&Particle, charge, mass, arrayLength);
    
    Particle.x[0] = 0;
    Particle.x[1] = 11.21;
    Particle.x[2] = 12.01;
    Particle.x[3] = 14.401;
    
    Particle.u[0] = 1.1;
    Particle.u[1] = 0.458;
    Particle.u[2] = 0;
    Particle.u[3] = 0;
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 1;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    for (int p = 0; p < tEnd / dt; p++){
        writeParticleToFile(&Particle, filename, p);
        updateVelocityWithBorisPusher(&Particle, Eextern, Bextern, dt);
        updateLocation(&Particle, dt);
        
        t += dt;
    }
    printf("executin bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/borisPusherScript.sh");
    freeMemoryOnParticle(&Particle, arrayLength);
}

void testNearFieldCalculation(){
    
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    double dx = 0.2;
    double dy = 0.2;
    double dz = 0.2;
    int numberOfGridPointsForBoxInX = 32;
    int numberOfGridPointsForBoxInY = 32;
    int numberOfGridPointsForBoxInZ = 32;
    int numberOfBoxesInX = 5;
    int numberOfBoxesInY = 5;
    int numberOfBoxesInZ = 5;
    
    Particle Particle;
    double charge = 1.0;
    double mass = 1.0;
    
    double dt = 0.2;
    double t = 0;
    double tEnd = 40;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateFieldsOnGrid(&Grid);
    writeGridParametersToFile(&Grid);
    initParticle(&Particle, charge, mass, arrayLength);
    
    Particle.x[0] = 0;
    Particle.x[1] = 11.21;
    Particle.x[2] = 12.01;
    Particle.x[3] = 14.401;
    
    Particle.u[1] = 0.458;
    Particle.u[2] = 0;
    Particle.u[3] = 0;
    Particle.u[0] = getGammaFromVelocityVector(&Particle);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0.2;
    
    int planeForPlotting = Particle.x[3] / dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    for (int p = 0; p < tEnd / dt; p++){
        writeFieldsToFile(&Grid, filename, p, planeForPlotting, true, false);
        getCurrentBoxIndexOfParticle(&Grid, &Particle);
        getEdgesOfNearFieldBox(&Grid, &Particle);
        writeParticleToFile(&Particle, filename, p);
        updateVelocityWithBorisPusher(&Particle, Eextern, Bextern, dt);
        updateLocation(&Particle, dt);
        t += dt;
    }
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/nearFieldScript.sh");
    freeMemoryOnParticle(&Particle, arrayLength);
    freeMemoryOnGrid(&Grid);
}

void testLWFieldCalculationForPlane(){
    
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    double dx = 0.125;
    double dy = 0.125;
    double dz = 0.125;
    int numberOfGridPointsForBoxInX = 32;
    int numberOfGridPointsForBoxInY = 32;
    int numberOfGridPointsForBoxInZ = 32;
    int numberOfBoxesInX = 5;
    int numberOfBoxesInY = 5;
    int numberOfBoxesInZ = 5;
    
    Particle Particle;
    double charge = 1.0;
    double mass = 1.0;
    
    double dt = 0.5 * dx;
    double t = 0;
    double tEnd = 6;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateFieldsOnGrid(&Grid);
    writeGridParametersToFile(&Grid);
    initParticle(&Particle, charge, mass, arrayLength);

    
    Particle.x[0] = 0;
    Particle.x[1] = 11.21;
    Particle.x[2] = 12.01;
    Particle.x[3] = 14.401;
    
    Particle.u[1] = 0.458;
    Particle.u[2] = 0;
    Particle.u[3] = 0;
    Particle.u[0] = getGammaFromVelocityVector(&Particle);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 1;
    
    int planeForPlotting = Particle.x[3] / dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    for (int p = 0; p < tEnd / dt; p++){
        writeParticleToFile(&Particle, filename, p);
        addCurrentStateToParticleHistory(&Particle, p);
        updateVelocityWithBorisPusher(&Particle, Eextern, Bextern, dt);
        updateLocation(&Particle, dt);
        t += dt;
    }

    calcLWFieldsForPlane(&Grid, &Particle, t, planeForPlotting);
    writeFieldsToFile(&Grid, filename, 0, planeForPlotting, true, false);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFieldsForPlane.sh");
    freeMemoryOnParticle(&Particle, arrayLength);
    freeMemoryOnGrid(&Grid);
}

