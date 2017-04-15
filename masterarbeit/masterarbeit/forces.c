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

///@brief initializes dampingForce array and lorentzFroce array with 0.
///@param Forces pointer to Forces struct
void initForces(Forces *Forces){
    for(int i = 0; i < 3; i++){
        Forces->damping[i] = 0;
        Forces->lorentz[i] = 0;
    }
    
}

///@brief analyzes lorentz and damping force by comparing them at every time step.
///@param Particle pointer to Particle struct
///@param Forces pointer to Forces struct
///@param Eextern vector containing external E-field components
///@param Bextern vector containing external B-field components
///@param t simulation time
void analyzeForces(Particle *Particle, Forces *Forces, double Eextern[3], double Bextern[3], double t){
    
    calcRadiationDamping(Eextern, Bextern, Particle->u, Forces->damping);
    calculateLorentzForce(Particle, Bextern, Forces->lorentz);
    writeForcesToFile(Forces, t);
    
}
