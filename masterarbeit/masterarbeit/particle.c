//
//  particle.c
//  masterarbeit
//
//  Created by David Symhoven on 10.11.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#include "particle.h"
#include "stdlib.h"

void initParticle(Particle *Particle, double const charge, double const mass, int const arrayLength){
    for (int i = 0; i < 4; i++){
        Particle->x[i] = 0;
        Particle->u[i] = 0;
    }
    Particle->charge = charge;
    Particle->mass = mass;
    double **trajectoryHistory = Particle->trajectoryHistory;
    double **velocityHistory = Particle->velocityHistory;
    
    // initialize trajectory[tEnd/dt][4]. First allocate array of pointers and then allocate for each pointer another pointer of length 4
    trajectoryHistory = (double **) malloc(arrayLength * sizeof(double *));
    if (trajectoryHistory == NULL){
        printf("ERROR: Could not allocate memory for trajectoryHistory[]");
    }
    else{
        for (int i = 0; i < arrayLength; i++){
            trajectoryHistory[i] = (double *) malloc(4 * sizeof(double));
            if (trajectoryHistory[i] == NULL){
                printf("ERROR: Could not allocate memory for trajectoryHistory[][]");
            }
        }
    }
    
    // same for velocity
    velocityHistory = (double **) malloc(arrayLength * sizeof(double *));
    if (velocityHistory == NULL){
        printf("ERROR: Could not allocate memory for velocityHistory[]");
    }
    else{
        for (int i = 0; i < arrayLength; i++){
            velocityHistory[i] = (double *) malloc(4 * sizeof(double));
            if (velocityHistory[i] == NULL){
                printf("ERROR: Could not allocate memory for velocityHistory[][]");
            }
        }
    }
}
