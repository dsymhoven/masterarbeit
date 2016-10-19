/* =============================================
HEADER FÜR EIGENE FUNKTIONEN
================================================ */
#import <stdbool.h>

#ifndef CALCULATIONS
#define CALCULATIONS

void borisPusher(double *u, double *E, double *B, double dt, double chargeOverMass);
void updateLocation(double *u, double *x, double dt);
void crossProduct(double *a, double *b, double *result);
bool isInsideBackwardLightcone(double *xParticle, double *xObserver);
bool isInsideForwardLightcone(double *xParticle, double *xObserver);
void vectorDifference(double *x, double *y, double *result);
double vectorProduct(double *x, double *y);
double calculateLambdaForLinearInterpolation(double *xInside, double *xOutside, double *xObserver);
double calculateDistance(double *x, double *y);
void calculateLienardWiechertParameters(double gamma_sq, double R_sq, double R, double *n, double *beta, double *xParticle, double *u, double *xObserver);
void calculateBetaDot(double *betaDot, double *u1, double *u2, double dt);
void calcuateLienardWiechertFields(double gamma_sq, double R_sq, double R, double *n, double *beta, double *beta_dot, double charge, double *E, double *B);

#endif
