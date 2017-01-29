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
#include "lwFields.h"


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
    allocateMemoryOnGrid(&Grid);
    initSamplePulseOnGrid(&Grid);
    writeGridParametersToFile(&Grid);
    
    double dt = 0.5 * Grid.dx;
    double t = 0;
    double tEnd = 3;
    char filename[32] = "some";
    int arrayLength = tEnd/dt;

    int planeForPlotting = 124;
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    
    for(int step = 0; step < arrayLength; step++){
        printf("step %d of %d\n", step, arrayLength);
        writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
        pushEFieldOnGrid(&Grid, dt);
        pushHFieldOnGrid(&Grid, dt);
        pushHFieldOnGrid(&Grid, dt);
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
    
    double dt = 0.5 * dx;
    double t = 0;
    double tEnd = 3;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle = &Particles[0];
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateMemoryOnGrid(&Grid);
    writeGridParametersToFile(&Grid);
    initParticles(Particles, numberOfParticles, arrayLength);
    
    Particle->mass = 1;
    Particle->charge = 1;
    Particle->x[0] = 0;
    Particle->x[1] = 11.21;
    Particle->x[2] = 12.01;
    Particle->x[3] = 14.401;
    
    Particle->u[0] = 1.1;
    Particle->u[1] = 0.458;
    Particle->u[2] = 0;
    Particle->u[3] = 0;
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 1;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    writeSimulationInfoToFile(numberOfParticles, t);
    for (int step = 0; step < tEnd / dt; step++){
        printf("step %d of %d\n", step, arrayLength);
        writeParticlesToFile(Particles, numberOfParticles, filename, step);
        for(int p = 0; p < numberOfParticles; p++){
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(Particle, &Grid, dt);
        }
        t += dt;
    }
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/borisPusherScript.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
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
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle = &Particles[0];;
    
    double dt = 0.5 * dx;
    double t = 0;
    double tEnd = 6;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateMemoryOnGrid(&Grid);
    writeGridParametersToFile(&Grid);
    initParticles(Particles, numberOfParticles, arrayLength);

    Particle->mass = 1;
    Particle->charge = 1;
    Particle->x[0] = 0;
    Particle->x[1] = 11.21;
    Particle->x[2] = 12.01;
    Particle->x[3] = 14.401;
    
    Particle->u[1] = 0.458;
    Particle->u[2] = 0;
    Particle->u[3] = 0;
    Particle->u[0] = getGammaFromVelocityVector(Particle->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 1;
    
    int planeForPlotting = Particle->x[3] / dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    writeSimulationInfoToFile(numberOfParticles, t);
    for (int step = 0; step < tEnd / dt; step++){
        printf("step %d of %d\n", step, (int)(tEnd / dt));
        writeParticlesToFile(Particles, numberOfParticles, filename, step);
        for(int p = 0; p < numberOfParticles; p++){
            addCurrentStateToParticleHistory(&Particles[p], step);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
        }
        
        t += dt;
    }
    

    calcLWFieldsForPlane(&Grid, Particles, numberOfParticles, t, planeForPlotting);
    writeFieldsToFile(&Grid, filename, 0, planeForPlotting, true, false);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFieldsForPlane.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
}


void testNearAndFarFields(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    double dx = 0.2;
    double dy = 0.2;
    double dz = 0.2;
    int numberOfGridPointsForBoxInX = 20;
    int numberOfGridPointsForBoxInY = 20;
    int numberOfGridPointsForBoxInZ = 20;
    int numberOfBoxesInX = 8;
    int numberOfBoxesInY = 8;
    int numberOfBoxesInZ = 8;
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle = &Particles[0];
    
    double dt = 0.5 * dx;
    double t = 0;
    double tEnd = 10;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateMemoryOnGrid(&Grid);
    initParticles(Particles, numberOfParticles, arrayLength);
    calcUPMLCoefficients(&Grid);
    
    Particle->mass = 1;
    Particle->charge = 1;
    Particle->x[0] = 0;
    Particle->x[1] = 14.21;
    Particle->x[2] = 14.01;
    Particle->x[3] = 14.401;
    
    Particle->u[1] = 0.458;
    Particle->u[2] = 0;
    Particle->u[3] = 0;
    Particle->u[0] = getGammaFromVelocityVector(Particle->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 1;
    
    int planeForPlotting = Particle->x[3] / dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    writeSimulationInfoToFile(numberOfParticles, t);
    
    for (int step = 0; step < tEnd / dt; step++){
        printf("step %d of %d\n", step, arrayLength);
        writeParticlesToFile(Particles, numberOfParticles, filename, step);
        writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);

        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        
        for(int p = 0; p < numberOfParticles; p++){
            addCurrentStateToParticleHistory(&Particles[p], step);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
        }
        
        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        
        t += dt;
    }
    clearFieldsFromGrid(&Grid);
    calcLWFieldsForPlaneWithNearField(&Grid, Particles, numberOfParticles, t, planeForPlotting);
    writeFieldsToFile(&Grid, filename, (int)(tEnd/dt), planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFields.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
 
}


void testUPML(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    double dx = 0.2;
    double dy = 0.2;
    double dz = 0.2;
    int numberOfGridPointsForBoxInX = 16;
    int numberOfGridPointsForBoxInY = 16;
    int numberOfGridPointsForBoxInZ = 16;
    int numberOfBoxesInX = 8;
    int numberOfBoxesInY = 8;
    int numberOfBoxesInZ = 8;
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle = &Particles[0];
    
    double dt = 0.5 * dx;
    double t = 0;
    double tEnd = 20;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateMemoryOnGrid(&Grid);
    calcUPMLCoefficients(&Grid);
    initParticles(Particles, numberOfParticles, arrayLength);
    
    Particle->mass = 1;
    Particle->charge = 1;
    Particle->x[0] = 0;
    Particle->x[1] = 12.21;
    Particle->x[2] = 12.01;
    Particle->x[3] = 12.401;
    
    Particle->u[1] = 0.458;
    Particle->u[2] = 0;
    Particle->u[3] = 0;
    Particle->u[0] = getGammaFromVelocityVector(Particle->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 1;
    
    int planeForPlotting = Particle->x[3] / dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    writeSimulationInfoToFile(numberOfParticles, t);
    for (int step = 0; step < tEnd / dt; step++){
        printf("step %d of %d\n", step, arrayLength);
        writeParticlesToFile(Particles, numberOfParticles, filename, step);
        writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
        
        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        
        for(int p = 0; p < numberOfParticles; p++){
            addCurrentStateToParticleHistory(&Particles[p], step);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
        }
        
        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        t += dt;
    }
    clearFieldsFromGrid(&Grid);
    calcLWFieldsForPlaneWithNearField(&Grid, Particles, numberOfParticles, t, planeForPlotting);
    writeFieldsToFile(&Grid, filename, (int)(tEnd/dt), planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFields.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
    
}

void testNearFieldUpdate(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    double dx = 0.2;
    double dy = 0.2;
    double dz = 0.2;
    int numberOfGridPointsForBoxInX = 20;
    int numberOfGridPointsForBoxInY = 20;
    int numberOfGridPointsForBoxInZ = 20;
    int numberOfBoxesInX = 8;
    int numberOfBoxesInY = 8;
    int numberOfBoxesInZ = 8;
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle = &Particles[0];
    
    double dt = 0.5 * dx;
    double t = 0;
    double tEnd = 20;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateMemoryOnGrid(&Grid);
    calcUPMLCoefficients(&Grid);
    initParticles(Particles, numberOfParticles, arrayLength);
    
    Particle->mass = 1;
    Particle->charge = 1;
    Particle->x[0] = 0;
    Particle->x[1] = 11.41;
    Particle->x[2] = 12.21;
    Particle->x[3] = 11.401;
    
    Particle->u[1] = 0.458;
    Particle->u[2] = 0;
    Particle->u[3] = 0;
    Particle->u[0] = getGammaFromVelocityVector(Particle->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0.2;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 1;
    
    int planeForPlotting = Particle->x[3] / dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    writeSimulationInfoToFile(numberOfParticles, t);
    for (int step = 0; step < tEnd / dt; step++){
        printf("step %d of %d\n", step, arrayLength);
        writeParticlesToFile(Particles, numberOfParticles, filename, step);
        writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
        
        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        
        for(int p = 0; p < numberOfParticles; p++){
            addCurrentStateToParticleHistory(&Particles[p], step);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
            updateNearField(&Grid, &Particles[p], t);
        }
        
        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        
        t += dt;
    }
    clearFieldsFromGrid(&Grid);
    calcLWFieldsForPlaneWithNearField(&Grid, Particles, numberOfParticles, t, planeForPlotting);
    writeFieldsToFile(&Grid, filename, (int)(tEnd/dt), planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFields.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
    
}


void testMultipleParticles(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    double dx = 0.2;
    double dy = 0.2;
    double dz = 0.2;
    int numberOfGridPointsForBoxInX = 20;
    int numberOfGridPointsForBoxInY = 20;
    int numberOfGridPointsForBoxInZ = 20;
    int numberOfBoxesInX = 5;
    int numberOfBoxesInY = 5;
    int numberOfBoxesInZ = 5;
    
    int numberOfParticles = 2;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    Particle *Particle2 = &Particles[1];
    
    double dt = 0.5 * dx;
    double t = 0;
    double tEnd = 15;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateMemoryOnGrid(&Grid);
    calcUPMLCoefficients(&Grid);
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle1->mass = 1;
    Particle1->charge = 1;
    Particle1->x[0] = 0;
    Particle1->x[1] = 6.41;
    Particle1->x[2] = 11.60;
    Particle1->x[3] = 11.401;
    
    Particle1->u[1] = 1.0;
    Particle1->u[2] = 0;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Particle2->mass = 1;
    Particle2->charge = 1;
    Particle2->x[0] = 0;
    Particle2->x[1] = 14.41;
    Particle2->x[2] = 10.00;
    Particle2->x[3] = 11.401;
    
    Particle2->u[1] = -1.0;
    Particle2->u[2] = 0;
    Particle2->u[3] = 0;
    Particle2->u[0] = getGammaFromVelocityVector(Particle2->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0;
    
    int planeForPlotting = Particle1->x[3] / dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    
    
    for (int step = 0; step < tEnd / dt; step++){
        printf("step %d of %d\n", step, arrayLength);
        writeParticlesToFile(Particles, numberOfParticles, filename, step);
        writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
        
        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        
        for(int p = 0; p < numberOfParticles; p++){
            addCurrentStateToParticleHistory(&Particles[p], step);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
            updateNearField(&Grid, &Particles[p], t);
        }
        
        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        t += dt;
    }
    //calcLWFieldsForPlane(&Grid, Particles, numberOfParticles, t, planeForPlotting);
    //writeFieldsToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFields.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
    
}

void testHistoryBeforeSimulation(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    double dx = 0.125;
    double dy = 0.125;
    double dz = 0.125;
    int numberOfGridPointsForBoxInX = 24;
    int numberOfGridPointsForBoxInY = 24;
    int numberOfGridPointsForBoxInZ = 24;
    int numberOfBoxesInX = 8;
    int numberOfBoxesInY = 8;
    int numberOfBoxesInZ = 8;
    
    int numberOfParticles = 2;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    Particle *Particle2 = &Particles[1];
    
    double dt = 0.5 * dx;
    double t = 8;
    double tEnd = 12;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateMemoryOnGrid(&Grid);
    calcUPMLCoefficients(&Grid);
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle1->mass = 1;
    Particle1->charge = 1;
    
    Particle1->x[0] = 0;
    Particle1->x[1] = 16.81;
    Particle1->x[2] = 16.20;
    Particle1->x[3] = 11.401;
    
    Particle1->u[1] = 0.458;
    Particle1->u[2] = 0;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Particle2->mass = 1;
    Particle2->charge = 1;
    Particle2->x[0] = 0;
    Particle2->x[1] = 8.81;
    Particle2->x[2] = 7.20;
    Particle2->x[3] = 11.401;
    
    Particle2->u[1] = 0.458;
    Particle2->u[2] = 0;
    Particle2->u[3] = 0;
    Particle2->u[0] = getGammaFromVelocityVector(Particle2->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 1;
    
    int planeForPlotting = Particle1->x[3] / dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    extendParticleHistory(Particles, &Grid, numberOfParticles, Eextern, Bextern, dt, t);
    calcFieldsOnGridWithoutNearField(Particles, &Grid, numberOfParticles, t);
    writeSimulationInfoToFile(numberOfParticles, t / dt);
    for (int step = t / dt; step < tEnd / dt; step++){
        printf("step %d of %f\n", step, tEnd / dt);
        writeParticlesToFile(Particles, numberOfParticles, filename, step);
        writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
        
        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        
        for(int p = 0; p < numberOfParticles; p++){
            addCurrentStateToParticleHistory(&Particles[p], step);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
            updateNearField(&Grid, &Particles[p], t);
        }

        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        t += dt;
   }
//    calcLWFieldsForPlaneWithNearField(&Grid, Particles, numberOfParticles, t, planeForPlotting);
//    writeFieldsToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFields.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
    
}

void testScattering(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    double dx = 0.125;
    double dy = 0.125;
    double dz = 0.125;
    int numberOfGridPointsForBoxInX = 20;
    int numberOfGridPointsForBoxInY = 20;
    int numberOfGridPointsForBoxInZ = 20;
    int numberOfBoxesInX = 8;
    int numberOfBoxesInY = 8;
    int numberOfBoxesInZ = 8;
    
    int numberOfParticles = 2;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    Particle *Particle2 = &Particles[1];
    
    double dt = 0.5 * dx;
    double t = 10;
    double tEnd = 24;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateMemoryOnGrid(&Grid);
    calcUPMLCoefficients(&Grid);
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle1->mass = 1;
    Particle1->charge = 1;
    
    Particle1->x[0] = 0;
    Particle1->x[1] = 14.1;
    Particle1->x[2] = 10.10;
    Particle1->x[3] = 11.401;
    
    Particle1->u[1] = -0.458;
    Particle1->u[2] = 0;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Particle2->mass = 1;
    Particle2->charge = 1;
    Particle2->x[0] = 0;
    Particle2->x[1] = 7.81;
    Particle2->x[2] = 9.90;
    Particle2->x[3] = 11.401;
    
    Particle2->u[1] = 0.458;
    Particle2->u[2] = 0;
    Particle2->u[3] = 0;
    Particle2->u[0] = getGammaFromVelocityVector(Particle2->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0;
    
    int planeForPlotting = Particle1->x[3] / dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    
    extendParticleHistory(Particles, &Grid, numberOfParticles, Eextern, Bextern, dt, t);
    writeSimulationInfoToFile(numberOfParticles, t / dt);
//    if(!readInitialFieldFromFileIfExists(&Grid, Particles, numberOfParticles, t, Eextern, Bextern)){
//        calcFieldsOnGridWithoutNearField(Particles, &Grid, numberOfParticles, t);
//        writeFieldsFromCompleteGridToFile(&Grid);
//        writeInitialConditionsToFile(&Grid, Particles, numberOfParticles, t, tEnd, Eextern, Bextern);
//        system("python2.7 ~/Desktop/Projects/masterarbeit/Analysis/initialFields.py");
//    }
    
    for (int step = t / dt; step < tEnd / dt; step++){
        printf("step %d of %d\n", step, (int)(tEnd / dt));
        writeParticlesToFile(Particles, numberOfParticles, filename, step);
//        writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
//        
//        pushEField(&Grid, Particles, numberOfParticles, t, dt);
//        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        
        for(int p = 0; p < numberOfParticles; p++){
            addCurrentStateToParticleHistory(&Particles[p], step);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
//            updateNearField(&Grid, &Particles[p], t);
        }
//        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
//        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        t += dt;
    }
//    clearFieldsFromGrid(&Grid);
    calcLWFieldsForPlaneWithNearField(&Grid, Particles, numberOfParticles, t, planeForPlotting);
    writeFieldsToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFieldsForPlane.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);

    
}

void extendHistoryAndCalcLWFieldsIndependantly(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    double dx = 0.125;
    double dy = 0.125;
    double dz = 0.125;
    int numberOfGridPointsForBoxInX = 20;
    int numberOfGridPointsForBoxInY = 20;
    int numberOfGridPointsForBoxInZ = 20;
    int numberOfBoxesInX = 8;
    int numberOfBoxesInY = 8;
    int numberOfBoxesInZ = 8;
    
    int numberOfParticles = 2;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    Particle *Particle2 = &Particles[1];
    
    double dt = 0.5 * dx;
    double t = 30;
    double tEnd = 44;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateMemoryOnGrid(&Grid);
    calcUPMLCoefficients(&Grid);
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle1->mass = 1;
    Particle1->charge = 1;
    
    Particle1->x[0] = 0;
    Particle1->x[1] = 14.1;
    Particle1->x[2] = 10.10;
    Particle1->x[3] = 11.401;
    
    Particle1->u[1] = -0.458;
    Particle1->u[2] = 0;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Particle2->mass = 1;
    Particle2->charge = 1;
    Particle2->x[0] = 0;
    Particle2->x[1] = 7.81;
    Particle2->x[2] = 9.90;
    Particle2->x[3] = 11.401;
    
    Particle2->u[1] = 0.458;
    Particle2->u[2] = 0;
    Particle2->u[3] = 0;
    Particle2->u[0] = getGammaFromVelocityVector(Particle2->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0;
    
    int planeForPlotting = Particle1->x[3] / dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    
    extendParticleHistory(Particles, &Grid, numberOfParticles, Eextern, Bextern, dt, t);
    writeSimulationInfoToFile(numberOfParticles, t / dt);
    if(!readInitialFieldFromFileIfExists(&Grid, Particles, numberOfParticles, t, Eextern, Bextern)){
        calcFieldsOnGridWithoutNearField(Particles, &Grid, numberOfParticles, t);
        writeFieldsFromCompleteGridToFile(&Grid);
        writeInitialConditionsToFile(&Grid, Particles, numberOfParticles, t, tEnd, Eextern, Bextern);
        system("python2.7 ~/Desktop/Projects/masterarbeit/Analysis/initialFields.py");
    }
    
    for (int step = t / dt; step < tEnd / dt; step++){
        printf("step %d of %d\n", step, (int)(tEnd / dt));
        writeParticlesToFile(Particles, numberOfParticles, filename, step);
        writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
        
        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        
        for(int p = 0; p < numberOfParticles; p++){
            addCurrentStateToParticleHistory(&Particles[p], step);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
            updateNearField(&Grid, &Particles[p], t);
        }
        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        t += dt;
    }
    clearFieldsFromGrid(&Grid);
    calcLWFieldsForPlaneWithNearField(&Grid, Particles, numberOfParticles, t, planeForPlotting);
    writeFieldsToFile(&Grid, filename, (int)(tEnd / dt), planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFields.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
    
}


void electronScatteringSmallGrid_init9(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    double dx = 0.125;
    double dy = 0.125;
    double dz = 0.125;
    int numberOfGridPointsForBoxInX = 20;
    int numberOfGridPointsForBoxInY = 20;
    int numberOfGridPointsForBoxInZ = 20;
    int numberOfBoxesInX = 8;
    int numberOfBoxesInY = 8;
    int numberOfBoxesInZ = 8;
    
    int numberOfParticles = 2;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    Particle *Particle2 = &Particles[1];
    
    double dt = 0.5 * dx;
    double t = 20;
    double tEnd = 32;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateMemoryOnGrid(&Grid);
    calcUPMLCoefficients(&Grid);
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle1->mass = 1;
    Particle1->charge = 1;
    
    Particle1->x[0] = 0;
    Particle1->x[1] = 14.1;
    Particle1->x[2] = 10.10;
    Particle1->x[3] = 11.401;
    
    Particle1->u[1] = -0.458;
    Particle1->u[2] = 0;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Particle2->mass = 1;
    Particle2->charge = 1;
    Particle2->x[0] = 0;
    Particle2->x[1] = 7.81;
    Particle2->x[2] = 9.90;
    Particle2->x[3] = 11.401;
    
    Particle2->u[1] = 0.458;
    Particle2->u[2] = 0;
    Particle2->u[3] = 0;
    Particle2->u[0] = getGammaFromVelocityVector(Particle2->u);
    
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0;
    
    int planeForPlotting = Particle1->x[3] / dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    
    extendParticleHistory(Particles, &Grid, numberOfParticles, Eextern, Bextern, dt, t);
    writeSimulationInfoToFile(numberOfParticles, t / dt);
//    if(!readInitialFieldFromFileIfExists(&Grid, Particles, numberOfParticles, t, Eextern, Bextern)){
//        calcFieldsOnGridWithoutNearField(Particles, &Grid, numberOfParticles, t);
//        writeFieldsFromCompleteGridToFile(&Grid);
//        writeInitialConditionsToFile(&Grid, Particles, numberOfParticles, t, tEnd, Eextern, Bextern);
//        system("python2.7 ~/Desktop/Projects/masterarbeit/Analysis/initialFields.py");
//    }
    
    for (int step = t / dt; step < tEnd / dt; step++){
        printf("step %d of %d\n", step, (int)(tEnd / dt));
        writeParticlesToFile(Particles, numberOfParticles, filename, step);
//        writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
//        
//        pushEField(&Grid, Particles, numberOfParticles, t, dt);
//        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        
        for(int p = 0; p < numberOfParticles; p++){
            addCurrentStateToParticleHistory(&Particles[p], step);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
//            updateNearField(&Grid, &Particles[p], t);
       }
//        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
//        pushEField(&Grid, Particles, numberOfParticles, t, dt);


        t += dt;
    }
//    clearFieldsFromGrid(&Grid);
    calcLWFieldsForPlane(&Grid, Particles, numberOfParticles, t, planeForPlotting);
    writeFieldsToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFieldsForPlane.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
    
    
}

void electronScatteringLargeGrid_init8(){
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
    int numberOfBoxesInX = 9;
    int numberOfBoxesInY = 9;
    int numberOfBoxesInZ = 9;
    
    int numberOfParticles = 2;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    Particle *Particle2 = &Particles[1];
    
    double dt = 0.5 * dx;
    double t = 40;
    double tEnd = 40;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
    initGrid(&Grid, dx, dy, dz, numberOfGridPointsForBoxInX, numberOfGridPointsForBoxInY, numberOfGridPointsForBoxInZ, numberOfBoxesInX, numberOfBoxesInY, numberOfBoxesInZ);
    allocateMemoryOnGrid(&Grid);
    calcUPMLCoefficients(&Grid);
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle1->mass = 1;
    Particle1->charge = 1;
    
    Particle1->x[0] = 0;
    Particle1->x[1] = 26.1;
    Particle1->x[2] = 18.1;
    Particle1->x[3] = 16.401;
    
    Particle1->u[1] = -0.458;
    Particle1->u[2] = 0;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Particle2->mass = 1;
    Particle2->charge = 1;
    Particle2->x[0] = 0;
    Particle2->x[1] = 10.81;
    Particle2->x[2] = 16.1;
    Particle2->x[3] = 16.401;
    
    Particle2->u[1] = 0.458;
    Particle2->u[2] = 0;
    Particle2->u[3] = 0;
    Particle2->u[0] = getGammaFromVelocityVector(Particle2->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0;
    
    int planeForPlotting = Particle1->x[3] / dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    
    extendParticleHistory(Particles, &Grid, numberOfParticles, Eextern, Bextern, dt, t);
    writeSimulationInfoToFile(numberOfParticles, t / dt);
    if(!readInitialFieldFromFileIfExists(&Grid, Particles, numberOfParticles, t, Eextern, Bextern)){
        calcFieldsOnGridWithoutNearField(Particles, &Grid, numberOfParticles, t);
        writeFieldsFromCompleteGridToFile(&Grid);
        writeInitialConditionsToFile(&Grid, Particles, numberOfParticles, t, tEnd, Eextern, Bextern);
        system("python2.7 ~/Desktop/Projects/masterarbeit/Analysis/initialFields.py");
    }
    
    for (int step = t / dt; step < tEnd / dt; step++){
        printf("step %d of %d\n", step, (int)(tEnd / dt));
        //        writeParticlesToFile(Particles, numberOfParticles, filename, step);
        //        writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
        //
        //        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        //        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        
        //        for(int p = 0; p < numberOfParticles; p++){
        //            addCurrentStateToParticleHistory(&Particles[p], step);
        //            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
        //            updateLocation(&Particles[p], &Grid, dt);
        //            updateNearField(&Grid, &Particles[p], t);
        //       }
        //        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        //        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        
        
        t += dt;
    }
    //clearFieldsFromGrid(&Grid);
    //calcLWFieldsForPlane(&Grid, Particles, numberOfParticles, t, planeForPlotting);
    writeFieldsToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFieldsForPlane.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
    
}




