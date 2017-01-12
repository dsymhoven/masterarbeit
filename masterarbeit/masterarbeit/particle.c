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

///@brief initializes all particles in Particles array with initParticle() method
///@param Particles pointer to Particle struct containing all particles
///@param numberOfParticles number of particles
///@param arrayLength length of particle history array
void initParticles(Particle *Particles, int const numberOfParticles, int const arrayLength){
    for(int i = 0; i < numberOfParticles; i++){
        initParticle(&Particles[i], arrayLength);
    }
}


/// @brief initializes all properties of struct Particle. x and u are initilized with 0
///@param Particle pointer to Particl struct
///@param arrayLength length of particle history array
void initParticle(Particle *Particle, int const arrayLength){
    printf("initializing Particle ...\n");
    for (int i = 0; i < 4; i++){
        Particle->x[i] = 0;
        Particle->u[i] = 0;
    }
    Particle->lengthOfHistoryArray = arrayLength;
    Particle->currentHistoryLength = 0;
    allocateParticleHistories(Particle, arrayLength);
}

///@brief allocates memory for trajectroy and velocity history.
///@remark Since x and u are both four vectors we need to allocate an array of arrays. We need "arrayLength = tEnd / dt" many arrays with length 4.
///@param Particle pointer to Particl struct
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
///@param Particles pointer to Particle struct containing all particles
///@param numberOfParticles number of particles
void writeParticlesToFile(Particle *Particles, int numberOfParticles, char *filename, int index){
    for(int p = 0; p < numberOfParticles; p++){
        
        printf("Writing Particle%d to file ...\n", p);
        FILE *fid = NULL;
        sprintf(filename, "Particle%d_%d", p, index);
        strcat(filename, ".txt");
        fid = fopen(filename,"w");
        if (fid == NULL){
            printf("ERROR: Could not open file Particle!\n");
        }
        
        for (int i = 0; i < 4; i++){
            fprintf(fid, "%f\t", Particles[p].x[i]);
        }
        fprintf(fid, "\n");
        
        for (int i = 0; i < 4; i++){
            fprintf(fid, "%f\t", Particles[p].u[i]);
        }
        fprintf(fid, "\n");
        for (int i = 0; i < 4; i++){
            fprintf(fid, "%f\t", Particles[p].edgesOfNearFieldBox[i]);
        }
        fprintf(fid, "\n");
        fclose(fid);
    }
    
}


///@brief method for releasing all previously allocated memory in struct Particle. Put all free() invokations in here
///@param Particles pointer to Particle struct containing all particles
///@param numberOfParticles number of particles
void freeMemoryOnParticles(Particle *Particles, int numberOfParticles){
    for (int p = 0; p < numberOfParticles; p++){
        printf("releasing allocated memory in Particle%d...\n", p);
        int arrayLength = Particles[p].lengthOfHistoryArray;
        for (int i = 0; i < arrayLength; i++){
            free(Particles[p].xHistory[i]);
            Particles[p].xHistory[i] = NULL;
            free(Particles[p].uHistory[i]);
            Particles[p].uHistory[i] = NULL;
        }
        free(Particles[p].xHistory);
        Particles[p].xHistory = NULL;
        
        free(Particles[p].uHistory);
        Particles[p].uHistory = NULL;
    }
    
}

///@brief saves the current box index for x,y and z in currentBoxIndexArray of Particle struct
///@param Grid pointer to Grid struct
///@param Particle pointer to Particle struct
void getCurrentBoxIndexArrayOfParticle(Grid *Grid, Particle *Particle){
    Particle->currentBoxIndexArray[0] = Particle->x[1] / (Grid->dx * Grid->numberOfGridPointsForBoxInX);
    Particle->currentBoxIndexArray[1] = Particle->x[2] / (Grid->dy * Grid->numberOfGridPointsForBoxInY);
    Particle->currentBoxIndexArray[2] = Particle->x[3] / (Grid->dz * Grid->numberOfGridPointsForBoxInZ);
    
}

/** @brief saves min and max values of x,y and z of the near field box in edgesOfNearFieldArray of Particle struct
 @code edgesOfNearFieldBox[0] = xMin;
 edgesOfNearFieldBox[1] = xMax;
 edgesOfNearFieldBox[2] = yMin;
 edgesOfNearFieldBox[3] = yMax;
 edgesOfNearFieldBox[4] = zMin;
 edgesOfNearFieldBox[5] = zMax;
 @endcode
 @param Grid pointer to Grid struct
 @param Particle pointer to Particle struct
 */
void getEdgesOfNearFieldBox(Grid *Grid, Particle *Particle){
    
    getCurrentBoxIndexArrayOfParticle(Grid, Particle);
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
///@param Particle instance of Particle struct
///@param index outer loop index to indicate the current time step.
void addCurrentStateToParticleHistory(Particle *Particle,  int index){
    
    for (int i = 0; i < 4; i++){
        Particle->xHistory[index][i] = Particle->x[i];
        Particle->uHistory[index][i] = Particle->u[i];
    }
    Particle->currentHistoryLength = index;

}

///@brief calculates gamma from spatial components of given velocity vector via gamma = sqrt(1 + u^2)
///@remark Don't use this method before spatial components of u are initialized
double getGammaFromVelocityVector(double u[4]){
    double result = 0.0;
    for (int i = 1; i < 4; i++){
        result += u[i] * u[i];
    }
    return sqrt(1 + result);
    
}

///@brief writes numberOfParticles and startTime of simulation to file. This is needed for python scripts, if simulation does not start at t = 0.
///@param numberOfParticles number of particles
///@param startTime start time of simulation
void writeSimulationInfoToFile(int numberOfParticles, int startTime){
    printf("Writing simulationInfo to file ...\n");
    
    FILE *fid = fopen("simulationInfo.txt","w");
    fprintf(fid, "%d %d", numberOfParticles, startTime);
    fclose(fid);
}

