//
//  serializers.h
//  masterarbeit
//
//  Created by David Symhoven on 01.04.17.
//  Copyright © 2017 David Symhoven. All rights reserved.
//

#ifndef serializers_h
#define serializers_h

#include <stdio.h>
#include "stdbool.h"
#include "grid.h"
#include "particle.h"
#include "forces.h"

#endif /* serializers_h */

void writeSimulationInfoToFile(int numberOfParticles, int startTime);
void writeFieldsToFile(Grid *Grid, char *filename, int index, int planeForPlotting, bool plotE, bool plotB);
void writeFieldComponentsForFourierAnalysisToFile(Grid *Grid, char *filename, int index, int planeForPlotting, bool plotE, bool plotB);
void writeFieldsFromCompleteGridToFile(Grid *Grid);
void writeParticlesToFile(Particle *Particles, int numberOfParticles, char *filename, int index);
void writeForcesToFile(Forces *Forces, double t);
void writeGridParametersToFile(Grid *Grid);
void writeExternalFieldsToFile(Grid *Grid, double Eextern[3], double Bextern[3], const double t, char *filename, const int index, const int planeForPlotting, bool plotE, bool plotB);
void writeSimulationInfoToFile(int numberOfParticles, int startTime);
void writeInitialConditionsToFile(Grid *Grid, Particle *Particles, int numberOfParticles, double t, double tEnd, double Eextern[3], double Bextern[3]);
