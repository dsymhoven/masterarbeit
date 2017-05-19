//
//  experiment.h
//  masterarbeit
//
//  Created by David Symhoven on 02.11.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#ifndef experiment_h
#define experiment_h

#include <stdio.h>

void testMaxwellPusher();
void testBorisPusher();
void testNearFieldCalculation();
void testLWFieldCalculationForPlane();
void testNearAndFarFields();
void testUPML();
void testNearFieldUpdate();
void testMultipleParticles();
void testHistoryBeforeSimulation();
void testScattering();
void extendHistoryAndCalcLWFieldsIndependantly();
void electronScatteringSmallGrid_init9();
void electronScatteringLargeGrid_init8();
void testTimeDependentExternalFields();
void scatteringInEMWave_analytic();
void scatteringInEMWave_simulation();
void testRadiationDampingVSLorentzForce();
void scatteringInEMWave_analytic_largeScale();
#endif /* experiment_h */
