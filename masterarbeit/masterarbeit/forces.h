//
//  forces.h
//  masterarbeit
//
//  Created by David Symhoven on 07.04.17.
//  Copyright © 2017 David Symhoven. All rights reserved.
//

#ifndef forces_h
#define forces_h

#include <stdio.h>
#include "particle.h"



struct Forces{
    double lorentz[3];
    double damping[3];
};

typedef struct Forces Forces;

void initForces(Forces *Forces);
void analyzeForces(Particle *Particle, Forces *Forces, double Eextern[3], double Bextern[3], double t);

#endif /* forces_h */
