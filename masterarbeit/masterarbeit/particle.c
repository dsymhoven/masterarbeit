//
//  particle.c
//  masterarbeit
//
//  Created by David Symhoven on 10.11.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#include "particle.h"
#include "stdlib.h"
#include "math.h"
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
    Particle->lengthOfHistoryArray = arrayLength;
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

///@brief writes current position, velocity vectors and near field info to file. First line is four vector x. Second line is four vector u. Third line is xMin up to yMax of near field box
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
    fprintf(fid, "\n");
    for (int i = 0; i < 4; i++){
        fprintf(fid, "%f\t", Particle->edgesOfNearFieldBox[i]);
    }
    fprintf(fid, "\n");
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

///@brief saves the current box index for x,y and z in currentBoxIndexArray of Particle struct
void getCurrentBoxIndexOfParticle(Grid *Grid, Particle *Particle){
    Particle->currentBoxIndexArray[0] = Particle->x[1] / (Grid->dx * Grid->numberOfGridPointsForBoxInX);
    Particle->currentBoxIndexArray[1] = Particle->x[2] / (Grid->dy * Grid->numberOfGridPointsForBoxInY);
    Particle->currentBoxIndexArray[2] = Particle->x[3] / (Grid->dz * Grid->numberOfGridPointsForBoxInZ);
    
}

///@brief saves min and max values of x,y and z of the near field box in edgesOfNearFieldArray of Particle struct
void getEdgesOfNearFieldBox(Grid *Grid, Particle *Particle){
    double xMin = (Particle->currentBoxIndexArray[0] - 1) * (Grid->dx * Grid->numberOfGridPointsForBoxInX);
    double xMax = (Particle->currentBoxIndexArray[0] + 1 + 1) * (Grid->dx * Grid->numberOfGridPointsForBoxInX);
    double yMin = (Particle->currentBoxIndexArray[1] - 1) * (Grid->dy * Grid->numberOfGridPointsForBoxInY);
    double yMax = (Particle->currentBoxIndexArray[1] + 1 + 1) * (Grid->dy * Grid->numberOfGridPointsForBoxInY);
    double zMin = (Particle->currentBoxIndexArray[2] - 1) * (Grid->dz * Grid->numberOfGridPointsForBoxInZ);
    double zMax = (Particle->currentBoxIndexArray[2] + 1 + 1) * (Grid->dz * Grid->numberOfGridPointsForBoxInZ);
    
    Particle->edgesOfNearFieldBox[0] = xMin;
    Particle->edgesOfNearFieldBox[1] = xMax;
    Particle->edgesOfNearFieldBox[2] = yMin;
    Particle->edgesOfNearFieldBox[3] = yMax;
    Particle->edgesOfNearFieldBox[4] = zMin;
    Particle->edgesOfNearFieldBox[5] = zMax;
    
}
///@brief adds current position and velocity information (the entire four vectors) to respective history array of Particle struct
///@param index outer loop index to indicate the current time step.
void addCurrentStateToParticleHistory(Particle *Particle, int index){
    for (int i = 0; i < 4; i++){
        Particle->xHistory[index][i] = Particle->x[i];
        Particle->uHistory[index][i] = Particle->u[i];
    }
}

///@brief calculates gamma from spatial components of given velocity vector via gamma = sqrt(1 + u^2)
///@remark Don't use this method before spatial components of u are initialized
double getGammaFromVelocityVector(Particle *Particle){
    double result = 0.0;
    for (int i = 1; i < 4; i++){
        result += Particle->u[i] * Particle->u[i];
    }
    return sqrt(1 + result);

}
