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

struct Particle {
    double charge;
    double mass;
    double x[4];
    double u[4];
    double **trajectoryHistory;
    double **velocityHistory;
    
};

typedef struct Particle Particle;

void initParticle(Particle *Particle, double const charge, double const mass, int const arrayLength);
#endif /* particle_h */
