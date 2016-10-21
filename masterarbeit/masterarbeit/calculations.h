/* =============================================
HEADER FÜR EIGENE FUNKTIONEN
================================================ */
#import <stdbool.h>

#ifndef CALCULATIONS
#define CALCULATIONS

void borisPusher(double *u, double *E, double *B, double dt, double chargeOverMass);
void updateLocation(double *u, double *x, double dt);
void crossProduct(double a[3], double b[3], double result[3]);
bool isInsideBackwardLightcone(double xParticle[4], double xObserver[4]);
bool isInsideForwardLightcone(double xParticle[4], double xObserver[4]);
void vectorDifference(double x[4], double y[4], double result[4]);
double vectorProduct(double x[4], double y[4]);
double calculateLambdaForLinearInterpolation(double xInside[4], double xOutside[4], double xObserver[4]);
double calculateDistance(double *x, double *y);
void calculateLienardWiechertParameters(double gamma_sq, double R_sq, double R, double *n, double *beta, double *xParticle, double *u, double *xObserver);
void calculateBetaDot(double *betaDot, double *u1, double *u2, double dt);
void calcuateLienardWiechertFields(double gamma_sq, double R_sq, double R, double *n, double *beta, double *beta_dot, double charge, double *E, double *B);
void calculateIntersectionPoint(double xInside[4], double xOutside[4], double xObserver[4], double intersectionPoint[4]);
#endif
