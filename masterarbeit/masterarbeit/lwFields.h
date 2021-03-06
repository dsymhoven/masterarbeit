//
//  lwFields.h
//  masterarbeit
//
//  Created by David Symhoven on 16.12.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#ifndef lwFields_h
#define lwFields_h

#include <stdio.h>
#include "grid.h"
#include "particle.h"


void adjustHFields(Grid *Grid, Particle *Particles, int numberOfParticles, const double t);
void adjustHyz_im1(Grid *Grid, Particle *Particles, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t);
void adjustHxz_jm1(Grid *Grid, Particle *Particles, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t);
void adjustHxy_km1(Grid *Grid, Particle *Particles, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t);
void adjustEFields(Grid *Grid, Particle *Particles, int numberOfParticles, const double t);
void adjustEyz_ip1(Grid *Grid, Particle *Particles, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t);
void adjustExz_jp1(Grid *Grid, Particle *Particles, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t);
void adjustExy_kp1(Grid *Grid, Particle *Particles, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t);
void pushHFieldAtBorders(Grid *Grid, double dt);
void pushEFieldAtBorders(Grid *Grid, double dt);
void setHFieldOnBorders(Grid *Grid);
void setEFieldOnBorders(Grid *Grid);
void pushHField(Grid *Grid, Particle *Particles, int numberOfParticles, double t, double dt);
void pushEField(Grid *Grid, Particle *Particles, int numberOfParticles, double t, double dt);
void pushEFieldInsideBoxes(Grid *Grid, double dt);
void pushHFieldInsideBoxes(Grid *Grid, double dt);
void addLWField(Grid *Grid, Particle *Particle, double *destination, double xObserver[4], int component);
void subLWField(Grid *Grid, Particle *Particle, double *destination, double xObserver[4], int component);
void addLWFieldsInBox(Grid *Grid, Particle *Particle, int boxIndex, double t);
void subLWFieldsInBox(Grid *Grid, Particle *Particle, int boxIndex, double t);
void updateFieldsForParticlePush(Particle *Particle, Grid *Grid, double Eextern[3], double Bextern[3], double E[3], double B[3]);
#endif /* lwFields_h */
