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

/// @brief initializes all properties of struct Particle. x and u are initilized with 0
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

///@brief allocates memory for trajectroy and velocity history.
///@remark Since x and u are both four vectors we need to allocate an array of arrays. We need "arrayLength = tEnd / dt" many arrays with length 4.
///@param arrayLength number of simulation steps. Meaning tEnd / dt
///@throws ERROR: Could not allocate memory for one of the arrays
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
///@param filename pointer to a char. Gets modified inside the method
///@param index outer loop index. Is used to name the output file
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

void getCurrentBoxIndexOfParticle(Grid *Grid, Particle *Particle, int currentBoxIndexArray[3]){
    currentBoxIndexArray[0] = Particle->x[1] / Grid->boxLengthInX;
    currentBoxIndexArray[1] = Particle->x[2] / Grid->boxLengthInY;
    currentBoxIndexArray[2] = Particle->x[3] / Grid->boxLengthInZ;
    
}

void getEdgesOfNearFieldBox(Grid *Grid, int currentBoxIndexArray[3], int sizeOfNearFieldBox, int edgesOfNearFieldBox[6]){
    int xMin = (currentBoxIndexArray[0] - sizeOfNearFieldBox) * Grid->boxLengthInX;
    int xMax = (currentBoxIndexArray[0] + sizeOfNearFieldBox + 1) * Grid->boxLengthInX;
    int yMin = (currentBoxIndexArray[1] - sizeOfNearFieldBox) * Grid->boxLengthInY;
    int yMax = (currentBoxIndexArray[1] - sizeOfNearFieldBox + 1) * Grid->boxLengthInY;
    int zMin = (currentBoxIndexArray[2] - sizeOfNearFieldBox) * Grid->boxLengthInZ;
    int zMax = (currentBoxIndexArray[2] - sizeOfNearFieldBox + 1) * Grid->boxLengthInZ;
    
    edgesOfNearFieldBox[0] = xMin;
    edgesOfNearFieldBox[1] = xMax;
    edgesOfNearFieldBox[2] = yMin;
    edgesOfNearFieldBox[3] = yMax;
    edgesOfNearFieldBox[4] = zMin;
    edgesOfNearFieldBox[5] = zMax;
    
}
