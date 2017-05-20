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
#include "serializers.h"


void testMaxwellPusher(){
    
    
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    Grid Grid;
    Resolution Resolution;
    Box Box;
 
    
    initResolution(&Resolution, 0.2, 0.2, 0.2);
    initBox(&Box, 20, 20, 20);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, false);
    
    initSamplePulseOnGrid(&Grid);
    writeGridParametersToFile(&Grid);
    
    double dt = 0.5 * Resolution.dx;
    double t = 0;
    double tEnd = 10;
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
    Resolution Resolution;
    Box Box;
 
    
    initResolution(&Resolution, 0.2, 0.2, 0.2);
    initBox(&Box, 20, 20, 20);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, false);
    
    
    double dt = 0.5 * Resolution.dx;
    double t = 0;
    double tEnd = 40;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle = &Particles[0];
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    Particle->mass = 1;
    Particle->charge = 1;
    
    Particle->x[0] = 0;
    Particle->x[1] = 14.1;
    Particle->x[2] = 18.1;
    Particle->x[3] = 14.1;
    
    Particle->u[1] = 0.958;
    Particle->u[2] = 0;
    Particle->u[3] = 0;
    Particle->u[0] = getGammaFromVelocityVector(Particle->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0.2;
    
    int planeForPlotting = Particle->x[3] / Resolution.dz;
    
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
    writeFieldsToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFieldsForPlane.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
     
}

void testNearFieldCalculation(){
    
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    Resolution Resolution;
    Box Box;
 
    
    initResolution(&Resolution, 0.2, 0.2, 0.2);
    initBox(&Box, 20, 20, 20);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, false);
    
    int numberOfParticles = 2;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    Particle *Particle2 = &Particles[1];
    
    double dt = 0.2;
    double t = 0;
    double tEnd = 1;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle1->mass = 1;
    Particle1->charge = 1;
    Particle1->x[0] = 0;
    Particle1->x[1] = 14.21;
    Particle1->x[2] = 14.01;
    Particle1->x[3] = 14.401;
    
    Particle1->u[1] = 0.458;
    Particle1->u[2] = 0;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Particle2->mass = 1;
    Particle2->charge = 1;
    Particle2->x[0] = 0;
    Particle2->x[1] = 22.21;
    Particle2->x[2] = 22.01;
    Particle2->x[3] = 14.401;
    
    Particle2->u[1] = 0.458;
    Particle2->u[2] = 0;
    Particle2->u[3] = 0;
    Particle2->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0.2;
    
    int planeForPlotting = Particle1->x[3] / Resolution.dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    writeSimulationInfoToFile(numberOfParticles, t);
    for (int step = 0; step < tEnd / dt; step++){
        printf("step %d of %d\n", step, (int)(tEnd / dt));
        writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
        writeParticlesToFile(Particles, numberOfParticles, filename, step);
        for (int p = 0; p < numberOfParticles; p++){
            getCurrentBoxIndexArrayOfParticle(&Grid, &Particles[p]);
            getEdgesOfNearFieldBox(&Grid, &Particles[p]);
            addCurrentStateToParticleHistory(&Particles[p], step);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
        }
        t += dt;
    }
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/nearFieldScript.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
     
}


void testLWFieldCalculationForPlane(){
    
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    Resolution Resolution;
    Box Box;
 
    
    initResolution(&Resolution, 0.1, 0.1, 0.1);
    initBox(&Box, 30, 30, 30);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, false);
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle = &Particles[0];
    
    double dt = 0.5 * Resolution.dx;
    double t = 0;
    double tEnd = 7;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    Particle->mass = 1;
    Particle->charge = 1;
    Particle->x[0] = 0;
    Particle->x[1] = 11.1;
    Particle->x[2] = 15.4;
    Particle->x[3] = 10.1;
    
    Particle->u[1] = 0.658;
    Particle->u[2] = 0;
    Particle->u[3] = 0;
    Particle->u[0] = getGammaFromVelocityVector(Particle->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0.2;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 1;
    
    int planeForPlotting = Particle->x[3] / Resolution.dz;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    writeSimulationInfoToFile(numberOfParticles, t);
    writeGridParametersToFile(&Grid);
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
    
    
    calcLWFieldsForPlaneWithNearField(&Grid, Particles, numberOfParticles, t, planeForPlotting);
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
    Resolution Resolution;
    Box Box;
 

    
    initResolution(&Resolution, 0.2, 0.2, 0.2);
    initBox(&Box, 22, 22, 22);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, false);
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle = &Particles[0];
    
    double dt = 0.5 * Resolution.dx;
    double t = 0;
    double tEnd = 12;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    Particle->mass = 1;
    Particle->charge = 1;
    Particle->x[0] = 0;
    Particle->x[1] = 16.21;
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
    
    int planeForPlotting = Particle->x[3] / Resolution.dz;
    
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
    Resolution Resolution;
    Box Box;
 
    
    initResolution(&Resolution, 0.2, 0.2, 0.2);
    initBox(&Box, 16, 16, 16);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, true);
    
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle = &Particles[0];
    
    double dt = 0.5 * Resolution.dx;
    double t = 0;
    double tEnd = 20;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    Particle->mass = 1;
    Particle->charge = 1;
    Particle->x[0] = 0;
    Particle->x[1] = 11.0;
    Particle->x[2] = 11.8;
    Particle->x[3] = 14.8;
    
    Particle->u[1] = -0.2;
    Particle->u[2] = 0;
    Particle->u[3] = 0;
    Particle->u[0] = getGammaFromVelocityVector(Particle->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0.2;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 1;
    
    int planeForPlotting = Particle->x[3] / Resolution.dz;
    
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
    Resolution Resolution;
    Box Box;
 
    
    initResolution(&Resolution, 0.1, 0.1, 0.1);
    initBox(&Box, 30, 30, 30);
      
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, true);
     
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle = &Particles[0];
    
    double dt = 0.5 * Resolution.dx;
    double t = 0;
    double tEnd = 12;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    Particle->mass = 1;
    Particle->charge = 1;
    Particle->x[0] = 0;
    Particle->x[1] = 11.1;
    Particle->x[2] = 15.5;
    Particle->x[3] = 10.1;
    
    Particle->u[1] = 0.658;
    Particle->u[2] = 0;
    Particle->u[3] = 0;
    Particle->u[0] = getGammaFromVelocityVector(Particle->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0.2;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 1;
    
    int planeForPlotting = Particle->x[3] / Resolution.dz;
    
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
    Resolution Resolution;
    Box Box;
 
    
    initResolution(&Resolution, 0.1, 0.1, 0.1);
    initBox(&Box, 32, 32, 32);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, true);
     
    
    int numberOfParticles = 2;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    Particle *Particle2 = &Particles[1];
    
    double dt = 0.5 * Resolution.dx;
    double t = 0;
    double tEnd = 12;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = tEnd / dt;
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle1->mass = 1;
    Particle1->charge = 1;
    Particle1->x[0] = 0;
    Particle1->x[1] = 4.0;
    Particle1->x[2] = 4.5;
    Particle1->x[3] = 6.8;
    
    Particle1->u[1] = 0.458;
    Particle1->u[2] = 0;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Particle2->mass = 1;
    Particle2->charge = 1;
    Particle2->x[0] = 0;
    Particle2->x[1] = 8.5;
    Particle2->x[2] = 8.0;
    Particle2->x[3] = 6.8;
    
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
    
    int planeForPlotting = Particle1->x[3] / Resolution.dz;
    
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

void testHistoryBeforeSimulation(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    Resolution Resolution;
    Box Box;
 
    

    initResolution(&Resolution, 0.1, 0.1, 0.1);
    initBox(&Box, 30, 30, 30);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, true);
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle = &Particles[0];
    
    double dt = 0.5 * Resolution.dx;
    double t = 12;
    double tEnd = 13;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle->mass = 1;
    Particle->charge = 1;
    Particle->x[0] = 0;
    Particle->x[1] = 11.1;
    Particle->x[2] = 15.5;
    Particle->x[3] = 10.1;
    
    Particle->u[1] = 0.658;
    Particle->u[2] = 0;
    Particle->u[3] = 0;
    Particle->u[0] = getGammaFromVelocityVector(Particle->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 1;
    
    int planeForPlotting = Particle->x[3] / Resolution.dz;
    
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
        printf("step %d of %f\n", step, tEnd / dt);
        writeParticlesToFile(Particles, numberOfParticles, filename, step);
//        writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
//        
//        pushEField(&Grid, Particles, numberOfParticles, t, dt);
//        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
//        
        for(int p = 0; p < numberOfParticles; p++){
            addCurrentStateToParticleHistory(&Particles[p], step);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
//            updateNearField(&Grid, &Particles[p], t);
        }
//
//        pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
//        pushEField(&Grid, Particles, numberOfParticles, t, dt);
        t += dt;
    }
//    clearFieldsFromGrid(&Grid);
//    calcLWFieldsForPlaneWithNearField(&Grid, Particles, numberOfParticles, t, planeForPlotting);
//    writeFieldsToFile(&Grid, filename, (int)(tEnd/dt), planeForPlotting, true, false);
    writeFieldsToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFieldsForPlane.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
     
    
}

void testScattering(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    Resolution Resolution;
    Box Box;
 
    
    initResolution(&Resolution, 0.2, 0.2, 0.2);
    initBox(&Box, 16, 16, 16);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, true);
     
    
    int numberOfParticles = 2;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    Particle *Particle2 = &Particles[1];
    
    double dt = 0.5 * Resolution.dx;
    double t = 0;
    double tEnd = 20;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle1->mass = 1;
    Particle1->charge = 1;
    
    Particle1->x[0] = 0;
    Particle1->x[1] = 3.0;
    Particle1->x[2] = 11.0;
    Particle1->x[3] = 14.8;
    
    Particle1->u[1] = 0.958;
    Particle1->u[2] = 0;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Particle2->mass = 1;
    Particle2->charge = 1;
    Particle2->x[0] = 0;
    Particle2->x[1] = 25.0;
    Particle2->x[2] = 12.0;
    Particle2->x[3] = 14.8;
    
    Particle2->u[1] = -0.958;
    Particle2->u[2] = 0;
    Particle2->u[3] = 0;
    Particle2->u[0] = getGammaFromVelocityVector(Particle2->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0;
    
    int planeForPlotting = Particle1->x[3] / Resolution.dz;
    
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
    writeFieldsToFile(&Grid, filename, (int)(tEnd/dt), planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFields.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
     
    
    
}

void extendHistoryAndCalcLWFieldsIndependantly(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    Resolution Resolution;
    Box Box;
 
    
    initResolution(&Resolution, 0.125, 0.125, 0.125);
    initBox(&Box, 20, 20, 20);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, true);
     
    
    int numberOfParticles = 2;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    Particle *Particle2 = &Particles[1];
    
    double dt = 0.5 * Resolution.dx;
    double t = 30;
    double tEnd = 44;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
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
    
    int planeForPlotting = Particle1->x[3] / Resolution.dz;
    
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
    Resolution Resolution;
    Box Box;
 
    
    initResolution(&Resolution, 0.125, 0.125, 0.125);
    initBox(&Box, 20, 20, 20);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, true);
     
    
    int numberOfParticles = 2;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    Particle *Particle2 = &Particles[1];
    
    double dt = 0.5 * Resolution.dx;
    double t = 20;
    double tEnd = 32;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
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
    
    int planeForPlotting = Particle1->x[3] / Resolution.dz;
    
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
    Resolution Resolution;
    Box Box;
 
    
    initResolution(&Resolution, 0.1, 0.1, 0.1);
    initBox(&Box, 32, 32, 32);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, true);
     
    
    int numberOfParticles = 2;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    Particle *Particle2 = &Particles[1];
    
    double dt = 0.5 * Resolution.dx;
    double t = 30;
    double tEnd = 50;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle1->mass = 1;
    Particle1->charge = 1;
    
    Particle1->x[0] = 0;
    Particle1->x[1] = 7;
    Particle1->x[2] = 11;
    Particle1->x[3] = 12.8;
    
    Particle1->u[1] = 0.979;
    Particle1->u[2] = 0;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Particle2->mass = 1;
    Particle2->charge = 1;
    
    Particle2->x[0] = 0;
    Particle2->x[1] = 20;
    Particle2->x[2] = 12;
    Particle2->x[3] = 12.8;
    
    Particle2->u[1] = -0.979;
    Particle2->u[2] = 0;
    Particle2->u[3] = 0;
    Particle2->u[0] = getGammaFromVelocityVector(Particle2->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0;
    
    int planeForPlotting = Particle1->x[3] / Resolution.dz;
    
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

void testTimeDependentExternalFields(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    Resolution Resolution;
    Box Box;
 
    
    initResolution(&Resolution, pow(10, 7), pow(10, 7), pow(10, 7));
    initBox(&Box, 32, 32, 32);
    initGrid(&Grid, &Resolution, &Box, 22, 10, 10, true);
     
    
    
    double dt = 0.5 * Resolution.dx;
    double t = 0;
    double tStart = t;
    double tEnd = 600 * pow(10, 7);
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    
    int planeForPlotting = Grid.numberOfGridPointsInZ * Resolution.dx / 2;
    
    // ======================================================
#pragma mark: Main Routine
    // ======================================================
    
    for (int step = t / dt; step < tEnd / dt; step++){
        printf("step %d of %d\n", step, (int)(tEnd / dt));
        
        t += dt;
    }
    writeExternalFieldsToFile(&Grid, Eextern, Bextern, t, tStart, filename, 0, planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/externalFields.sh");
    freeMemoryOnGrid(&Grid);
     
    
}

void scatteringInEMWave_analytic(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    Resolution Resolution;
    Box Box;
    Forces Forces;
 
    
    initResolution(&Resolution, 10, 10, 10);
    initBox(&Box, 20, 20, 20);
    initGrid(&Grid, &Resolution, &Box, 25, 25, 25, true);
    initForces(&Forces);
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    
    double dt = 0.5 * Resolution.dx;
    double t = 4000;
    double tStart = t;
    double tEnd = 7500;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle1->mass = 1;
    Particle1->charge = 1;
    
    Particle1->x[0] = 0;
    Particle1->x[1] = 3600;
    Particle1->x[2] = 100;
    Particle1->x[3] = 1600;
    
    Particle1->u[1] = -1.6;
    Particle1->u[2] = 1.6;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0;
    
    int planeForPlotting = Particle1->x[3] / Resolution.dz;
    
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
        //                writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
        //
        //                pushEField(&Grid, Particles, numberOfParticles, t, dt);
        //                pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        
        
        for(int p = 0; p < numberOfParticles; p++){
            analyzeForces(&Particles[p], &Forces, Eextern, Bextern, t);
            addCurrentStateToParticleHistory(&Particles[p], step);
            externalPlaneWave(Particles[p].x, tStart, Eextern, Bextern);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
            //                    updateNearField(&Grid, &Particles[p], t);
        }
        //                pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        //                pushEField(&Grid, Particles, numberOfParticles, t, dt);
        
        
        t += dt;
    }
    clearFieldsFromGrid(&Grid);
    calcLWFieldsForPlaneWithNearField(&Grid, Particles, numberOfParticles, t, planeForPlotting);
    writeFieldComponentsForFourierAnalysisToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeFieldsToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeExternalFieldsToFile(&Grid, Eextern, Bextern, t, tStart, filename, 0, planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/fourierAnalysis.sh");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/externalFields.sh");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/analyzeForces.sh");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFieldsForPlane.sh");
    
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
     
    
    
}

void scatteringInEMWave_simulation(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    Resolution Resolution;
    Box Box;
    
    
    initResolution(&Resolution, 0.1, 0.1, 0.1);
    initBox(&Box, 32, 32, 32);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, true);
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    
    double dt = 0.5 * Resolution.dx;
    double t = 40;
    double tStart = t;
    double tEnd = 75;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle1->mass = 1;
    Particle1->charge = 1;
    
    Particle1->x[0] = 0;
    Particle1->x[1] = 19;
    Particle1->x[2] = 1;
    Particle1->x[3] = 16.401;
    
    Particle1->u[1] = -0.4;
    Particle1->u[2] = 0.4;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0;
    
    int planeForPlotting = Particle1->x[3] / Resolution.dz;
    
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
        //                writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
        //
        //                pushEField(&Grid, Particles, numberOfParticles, t, dt);
        //                pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        
        
        for(int p = 0; p < numberOfParticles; p++){
            addCurrentStateToParticleHistory(&Particles[p], step);
            externalPlaneWave(Particles[p].x, tStart, Eextern, Bextern);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
            //                    updateNearField(&Grid, &Particles[p], t);
        }
        //                pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        //                pushEField(&Grid, Particles, numberOfParticles, t, dt);
        
        
        t += dt;
    }
    clearFieldsFromGrid(&Grid);
    calcLWFieldsForPlaneWithNearField(&Grid, Particles, numberOfParticles, t, planeForPlotting);
    //    writeFieldComponentsForFourierAnalysisToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeFieldsToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeExternalFieldsToFile(&Grid, Eextern, Bextern, t, tStart, filename, 0, planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    //    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/fourierAnalysis.sh");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/externalFields.sh");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFieldsForPlane.sh");
    
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
    
    
    
}

void testRadiationDampingVSLorentzForce(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    Resolution Resolution;
    Box Box;
    Forces Forces;
 
    
    initResolution(&Resolution, 0.1, 0.1, 0.1);
    initBox(&Box, 32, 32, 32);
    initGrid(&Grid, &Resolution, &Box, 8, 8, 8, true);
    initForces(&Forces);
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    
    double dt = 0.5 * Resolution.dx;
    double t = 0;
    double tStart = t;
    double tEnd = 50;
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    Particle1->mass = 1;
    Particle1->charge = 1;
    
    Particle1->x[0] = 0;
    Particle1->x[1] = 13;
    Particle1->x[2] = 17;
    Particle1->x[3] = 16;
    
    Particle1->u[1] = 0.458;
    Particle1->u[2] = 0;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0.1;
    
    int planeForPlotting = Particle1->x[3] / Resolution.dz;
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
            analyzeForces(&Particles[p], &Forces, Eextern, Bextern, t);
            addCurrentStateToParticleHistory(&Particles[p], step);
//            externalPlaneWave(Particles[p].x, tStart, Eextern, Bextern);
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
//    writeFieldComponentsForFourierAnalysisToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeFieldsToFile(&Grid, filename, 0, planeForPlotting, true, false);
//    writeExternalFieldsToFile(&Grid, Eextern, Bextern, t, tStart, filename, 0, planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
//    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/fourierAnalysis.sh");
//    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/externalFields.sh");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFieldsForPlane.sh");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/analyzeForces.sh");
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
     

}

void scatteringInEMWave_analytic_largeScale(){
    // ======================================================
#pragma mark: Initializations
    // ======================================================
    
    Grid Grid;
    Resolution Resolution;
    Box Box;
    Forces Forces;
    
    
    initResolution(&Resolution, pow(10, 7), pow(10, 7), pow(10, 7));
    initBox(&Box, 32, 32, 32);
    initGrid(&Grid, &Resolution, &Box, 22, 10, 10, true);
    initForces(&Forces);
    
    int numberOfParticles = 1;
    
    Particle Particles[numberOfParticles];
    Particle *Particle1 = &Particles[0];
    
    double dt = 0.5 * Resolution.dx;
    double t = 0;
    double tStart = t;
    double tEnd = 700 * pow(10, 7);
    
    char filename[32] = "some";
    double Eextern[3];
    double Bextern[3];
    int arrayLength = (tEnd - t) / dt;
    
    initParticles(Particles, numberOfParticles, arrayLength);
    
    
    Particle1->mass = pow(10, 10);
    Particle1->charge = pow(10, 10);
    
    Particle1->x[0] = 0;
    Particle1->x[1] = 6.0 * pow(10, 9);
    Particle1->x[2] = 1.3 * pow(10, 9);;
    Particle1->x[3] = 1.6 * pow(10, 9);
    
    Particle1->u[1] = -0.2;
    Particle1->u[2] = 0.1;
    Particle1->u[3] = 0;
    Particle1->u[0] = getGammaFromVelocityVector(Particle1->u);
    
    Eextern[0] = 0;
    Eextern[1] = 0;
    Eextern[2] = 0;
    
    Bextern[0] = 0;
    Bextern[1] = 0;
    Bextern[2] = 0;
    
    int planeForPlotting = Particle1->x[3] / Resolution.dz;
    
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
        //                writeFieldsToFile(&Grid, filename, step, planeForPlotting, true, false);
        //
        //                pushEField(&Grid, Particles, numberOfParticles, t, dt);
        //                pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        
        
        for(int p = 0; p < numberOfParticles; p++){
            analyzeForces(&Particles[p], &Forces, Eextern, Bextern, t);
            addCurrentStateToParticleHistory(&Particles[p], step);
            externalPulse(Particles[p].x, tStart, Eextern, Bextern, &Grid);
            updateVelocityWithBorisPusher(Particles, &Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], &Grid, dt);
            //                    updateNearField(&Grid, &Particles[p], t);
        }
        //                pushHField(&Grid, Particles, numberOfParticles, t + dt / 2., dt);
        //                pushEField(&Grid, Particles, numberOfParticles, t, dt);
        
        
        t += dt;
    }
    clearFieldsFromGrid(&Grid);
    calcLWFieldsForPlaneWithNearField(&Grid, Particles, numberOfParticles, t, planeForPlotting);
    writeFieldComponentsForFourierAnalysisToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeFieldsToFile(&Grid, filename, 0, planeForPlotting, true, false);
    writeExternalFieldsToFile(&Grid, Eextern, Bextern, t, tStart, filename, 0, planeForPlotting, true, false);
    writeGridParametersToFile(&Grid);
    printf("executing bash-script ...\n");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/fourierAnalysis.sh");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/externalFields.sh");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/analyzeForces.sh");
    system("~/Desktop/Projects/masterarbeit/Analysis/Scripts/particlesAndFieldsForPlane.sh");
    
    freeMemoryOnParticles(Particles, numberOfParticles);
    freeMemoryOnGrid(&Grid);
    
    
    
}






