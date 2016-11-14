//
//  particle.c
//  masterarbeit
//
//  Created by David Symhoven on 10.11.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#include "particle.h"
#include "stdlib.h"
#include "string.h"

void initParticle(Particle *Particle, double const charge, double const mass, int const arrayLength){
    printf("initializing Particle ...\n");
    for (int i = 0; i < 4; i++){
        Particle->x[i] = 0;
        Particle->u[i] = 0;
    }
    Particle->charge = charge;
    Particle->mass = mass;
    allocateParticleHistories(Particle, arrayLength);
}

void allocateParticleHistories(Particle *Particle, int const arrayLength){
    printf("allocating memory for Particle histories\n");
    
    // initialize xHistory[tEnd/dt][4]. First allocate array of pointers and then allocate for each pointer another pointer of length 4
    Particle->xHistory = (double **) malloc(arrayLength * sizeof(double *));
    if (Particle->xHistory == NULL){
        printf("ERROR: Could not allocate memory for xHistory[]\n");
    }
    else{
        for (int i = 0; i < arrayLength; i++){
            Particle->xHistory[i] = (double *) malloc(4 * sizeof(double));
            if (Particle->xHistory[i] == NULL){
                printf("ERROR: Could not allocate memory for xHistory[][]\n");
            }
        }
    }
    
    // same for velocity
    Particle->uHistory = (double **) malloc(arrayLength * sizeof(double *));
    if (Particle->uHistory == NULL){
        printf("ERROR: Could not allocate memory for uHistory[]\n");
    }
    else{
        for (int i = 0; i < arrayLength; i++){
            Particle->uHistory[i] = (double *) malloc(4 * sizeof(double));
            if (Particle->uHistory[i] == NULL){
                printf("ERROR: Could not allocate memory for uHistory[][]\n");
            }
        }
    }

    
}

///@brief writes current position and velocity vectors to file. First line is four vector x. Second line is four vector u.
void writeParticleToFile(Particle *Particle, char *filename, int index){
    printf("Writing Particle to file ...\n");
    FILE *fid = NULL;
    sprintf(filename, "Particle%d", index);
    strcat(filename, ".txt");
    fid = fopen(filename,"w");
    if (fid == NULL){
        printf("ERROR: Could not open file Particle!\n");
    }
    
    for (int i = 0; i < 4; i++){
        fprintf(fid, "%f\t", Particle->x[i]);
    }
    fprintf(fid, "\n");
    
    for (int i = 0; i < 4; i++){
        fprintf(fid, "%f\t", Particle->u[i]);
    }
    fclose(fid);
}

///@brief  method for releasing all previously allocated memory in struct Particle. Put all free() invokations in here
void freeMemoryOnParticle(Particle *Particle, int const arrayLength){
    printf("releasing allocated memory in Particle...\n");
    for (int i = 0; i < arrayLength; i++){
        free(Particle->xHistory[i]);
        free(Particle->uHistory[i]);
    }
    free(Particle->xHistory);
    free(Particle->uHistory);
    
    
}
