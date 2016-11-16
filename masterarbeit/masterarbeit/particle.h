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
    
};

typedef struct Particle Particle;

void initParticle(Particle *Particle, double const charge, double const mass, int const arrayLength);
void writeParticleToFile(Particle *Particle, char *filename, int index);
void freeMemoryOnParticle(Particle *Particle, int const arrayLength);
void allocateParticleHistories(Particle *Particle, int const arrayLength);
void getCurrentBoxIndexOfParticle(Grid *Grid, Particle *Particle, int currentBoxIndexArray[3]);
void getEdgesOfNearFieldBox(Grid *Grid, int currentBoxIndexArray[3], int sizeOfNearFieldBox, int edgesOfNearFieldBox[6]);
#endif /* particle_h */
