//
//  forces.c
//  masterarbeit
//
//  Created by David Symhoven on 07.04.17.
//  Copyright © 2017 David Symhoven. All rights reserved.
//

#include "forces.h"
#include "particle.h"
#include "calculations.h"
#include "serializers.h"

void initForces(Forces *Forces){
    for(int i = 0; i < 3; i++){
        Forces->damping[i] = 0;
        Forces->lorentz[i] = 0;
    }
    
}

void analyzeForces(Particle *Particle, Forces *Forces, double Eextern[3], double Bextern[3], double t){
    
    calcRadiationDamping(Eextern, Bextern, Particle->u, Forces->damping);
    calculateLorentzForce(Particle, Bextern, Forces->lorentz);
    writeForcesToFile(Forces, t);
    
}
