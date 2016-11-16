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
    int numberOfGridPointsInX = 256;
    int numberOfGridPointsInY = 256;
    int numberOfGridPointsInZ = 256;
    int lengthOfSimulationBoxInX = 32;
    int lengthOfSimulationBoxInY = 32;
    int lengthOfSimulationBoxInZ = 32;
    
    initGrid(&Grid, numberOfGridPointsInX, numberOfGridPointsInY, numberOfGridPointsInZ, lengthOfSimulationBoxInX, lengthOfSimulationBoxInY, lengthOfSimulationBoxInZ);
    allocateFieldsOnGrid(&Grid);
    initSamplePulseOnGrid(&Grid);
    
    double dt = 0.5 * Grid.dx;
    double t = 0;
    double tEnd = 3;
    char filename[32] = "some";
    int arrayLength = tEnd/dt;

    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    for(int p = 0; p < arrayLength; p++){
        writeFieldsToFile(&Grid, filename, p, true, false);
        
        pushEFieldOnGrid(&Grid, dt);
        PushBFieldOnGrid(&Grid, dt);
        PushBFieldOnGrid(&Grid, dt);
        pushEFieldOnGrid(&Grid, dt);
        t += dt;
    }
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/script.sh");
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
    double gamma = 2.1;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initParticle(&Particle, charge, mass, arrayLength);
    
    Particle.x[0] = 0;
    Particle.x[1] = 10;
    Particle.x[2] = 14;
    Particle.x[3] = 10;
    
    Particle.u[0] = gamma;
    Particle.u[1] = 1.5;
    Particle.u[2] = 1.5;
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
    system("~/Desktop/Projects/masterarbeit/Analysis/particleScript.sh");
    freeMemoryOnParticle(&Particle, arrayLength);
}

void testNearFieldCalculation(){
    
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    int numberOfGridPointsInX = 256;
    int numberOfGridPointsInY = 256;
    int numberOfGridPointsInZ = 256;
    int lengthOfSimulationBoxInX = 16;
    int lengthOfSimulationBoxInY = 16;
    int lengthOfSimulationBoxInZ = 16;
    
    Particle Particle;
    double charge = 1.0;
    double mass = 1.0;
    
    
    double dt = 0.5 * 0.125;
    double t = 0;
    double tEnd = 3;
    double gamma = 2.1;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initGrid(&Grid, numberOfGridPointsInX, numberOfGridPointsInY, numberOfGridPointsInZ, lengthOfSimulationBoxInX, lengthOfSimulationBoxInY, lengthOfSimulationBoxInZ);
    allocateFieldsOnGrid(&Grid);
    writeGridParametersToFile(&Grid);
    initParticle(&Particle, charge, mass, arrayLength);
    
    Particle.x[0] = 0;
    Particle.x[1] = 10;
    Particle.x[2] = 14;
    Particle.x[3] = 10;
    
    Particle.u[0] = gamma;
    Particle.u[1] = 1.5;
    Particle.u[2] = 1.5;
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
        writeFieldsToFile(&Grid, filename, p, true, false);
        writeParticleToFile(&Particle, filename, p);
        updateVelocityWithBorisPusher(&Particle, Eextern, Bextern, dt);
        updateLocation(&Particle, dt);
        
        t += dt;
    }
    printf("executin bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/postProcessing.sh");
    freeMemoryOnParticle(&Particle, arrayLength);
}
