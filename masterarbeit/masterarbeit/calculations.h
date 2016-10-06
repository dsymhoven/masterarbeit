/* =============================================
HEADER FÜR EIGENE FUNKTIONEN
================================================ */
#import <stdbool.h>

#ifndef CALCULATIONS
#define CALCULATIONS

void borisPusher(double *u, double *E, double *B, double dt, double chargeOverMass);
void crossProduct(double *a, double *b, double *result);
bool isInsideBackwardLightcone(double *xParticle, double *xObserver);
bool isInsideForwardLightcone(double *xParticle, double *xObserver);
void minkowskiDifference(double *x, double *y, double *result);
double minkowskiProduct(double *x, double *y);
double calculateLambdaForLinearInterpolation(double *xInside, double *xOutside, double *xObserver);
#endif
