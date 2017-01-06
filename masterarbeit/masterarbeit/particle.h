//
//  particle.h
//  masterarbeit
//
//  Created by David Symhoven on 10.11.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#ifndef particle_h
#define particle_h

#include <stdio.h>
#include "grid.h"

struct Particle {
    double charge;
    double mass;
    double x[4];
    double u[4];
    double **xHistory;
    double **uHistory;
    double edgesOfNearFieldBox[6];
    int currentBoxIndexArray[3];
    int lengthOfHistoryArray;
    int currentHistoryLength;
    bool didChangeBox;
    int boxIndicesOfNearFieldBoxesBeforePush[27];
    int boxIndicesOfNearFieldBoxesAfterPush[27];
    
};

typedef struct Particle Particle;
typedef struct Grid Grid;

void initParticle(Particle *Particle, int const arrayLength);
void initParticles(Particle *Particles, int const numberOfParticles, int const arrayLength);
void writeParticlesToFile(Particle *Particles, int numberOfParticles, char *filename, int index);
void freeMemoryOnParticles(Particle *Particles, int const numberOfParticles);
void allocateParticleHistories(Particle *Particle, int const arrayLength);
void getCurrentBoxIndexArrayOfParticle(Grid *Grid, Particle *Particle);
void getEdgesOfNearFieldBox(Grid *Grid, Particle *Particle);
void addCurrentStateToParticleHistory(Particle *Particle, int index);
double getGammaFromVelocityVector(double u[4]);
void writeSimulationInfoToFile(int numberOfParticles, int startTime);
#endif /* particle_h */
