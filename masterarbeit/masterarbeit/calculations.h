/* =============================================
HEADER FÜR EIGENE FUNKTIONEN
================================================ */
#import <stdbool.h>
#import "particle.h"

#ifndef CALCULATIONS
#define CALCULATIONS

void updateVelocityWithBorisPusher(Particle *Particle, double *Eextern, double *Bextern, double dt);
void updateLocation(Particle *Particle, double dt);
void crossProduct(double a[3], double b[3], double result[3]);
bool isInsideBackwardLightcone(double xParticle[4], double xObserver[4]);
bool isInsideForwardLightcone(double xParticle[4], double xObserver[4]);
void vectorDifference(double x[4], double y[4], double result[4]);
double vectorProduct(double x[4], double y[4]);
double calculateLambdaForLinearInterpolation(double xInside[4], double xOutside[4], double xObserver[4]);
double calculateDistance(double *x, double *y);
void calculateLienardWiechertParameters(double xParticle[4], double xObserver[4], double u[4], double *gamma_sq, double *R_sq, double *R, double n[3], double beta[3]);
void calculateBetaDot(double *uOld, double *uNew, double dt, double betaDot[3]);
void calcuateLienardWiechertFields(double gamma_sq, double R_sq, double R, double *n, double *beta, double *beta_dot, double charge, double *E, double *B);
void calculateIntersectionPoint(double xInside[4], double xOutside[4], double xObserver[4], double intersectionPoint[4]);
void calculateBeta(double xOld[4], double xNew[4], double beta[3]);
#endif
