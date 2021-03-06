/* =============================================
HEADER FÜR EIGENE FUNKTIONEN
================================================ */
#import <stdbool.h>
#import "particle.h"

#ifndef CALCULATIONS
#define CALCULATIONS

void updateVelocityWithBorisPusher(Particle *Particles, Grid *Grid, int numberOfParticles, int particleIndex, double Eextern[3], double Bextern[3], double dt);
void updateLocation(Particle *Particle, Grid *Grid, double dt);
void crossProduct(double a[3], double b[3], double result[3]);
void scaleVector(double x[3], double factor);
bool isInsideBackwardLightcone(double xParticle[4], double xObserver[4]);
bool isInsideForwardLightcone(double xParticle[4], double xObserver[4]);
void vectorDifference(double x[4], double y[4], double result[4]);
double vectorProduct(double x[4], double y[4]);
double calculateLambdaForLinearInterpolation(double xInside[4], double xOutside[4], double xObserver[4]);
double calculateDistance(double *x, double *y);
void calculateLienardWiechertParameters(double xParticle[4], double xObserver[4], double u[4], double *gamma_sq, double *R_sq, double *R, double n[3], double beta[3]);
void calculateBetaDot(double *uOld, double *uNew, double dt, double betaDot[3]);
void calcuateLienardWiechertFields(double gamma_sq, double R_sq, double R, double *n, double *beta, double *beta_dot, double charge, double *E, double *B);
void calculateIntersectionPoint(double xInside[4], double xOutside[4], double uInside[4], double uOutside[4], double xObserver[4], double intersectionPoint[4], double velocityAtIntersectionPoint[4]);
void calcLWFieldsEverywhereOnGrid(Grid *Grid, Particle *Particle, int timeStep);
void calculateBeta(double xOld[4], double xNew[4], double beta[3]);
void calcLWFieldsOnGrid(Grid *Grid, Particle *Particle, double t);
void calcLWFieldsForPlane(Grid *Grid, Particle *Particles, int numberOfParticles, double t, int planeForPlotting);
void calcLWFieldsForPlaneWithNearField(Grid *Grid, Particle *Particles, int numberOfParticles, double t, int planeForPlotting);
int calcCurrentBoxIndexOfParticle(Particle *Particle, Grid *Grid);
void calcBoxIndizesOfNextNeighbourBoxes(Grid *Grid, Particle *Particle, int boxIndizesOfNextNeighbourBoxes[27]);
bool boxIsInNearFieldOfParticle(Grid *Grid, Particle *Particle, int boxIndex);
int calcBoxIndexIm1(Grid *Grid, const int boxIndex);
int calcBoxIndexJm1(Grid *Grid, const int boxIndex);
int calcBoxIndexKm1(Grid *Grid, const int boxIndex);
int calcBoxIndexIp1(Grid *Grid, const int boxIndex);
int calcBoxIndexJp1(Grid *Grid, const int boxIndex);
int calcBoxIndexKp1(Grid *Grid, const int boxIndex);
void updateNearField(Grid *Grid, Particle *Particle, double t);
void calcInteractionWithOtherParticles(Particle *Particles, Grid *Grid, int numberOfParticles, int particleIndex, double E[3], double B[3]);
void calcFieldsOnGridWithoutNearField(Particle *Particles, Grid *Grid, int numberOfParticles, double t);
void extendParticleHistory(Particle *Particles, Grid *Grid, int numberOfParticles, double Eextern[3], double Bextern[3], double dt, double t);
bool readInitialFieldFromFileIfExists(Grid *Grid, Particle *Particles, int numberOfParticles, double t, double Eextern[3], double Bextern[3]);
void interpolateFields(Grid *Grid, Particle *Particle, double E[3], double B[3]);
void calcGridIndizesNextNeighboursEx(Grid *Grid, Particle *Particle, int gridIndizesNextNeighbours[8]);
void calculateLorentzForce(Particle *Particle, double Bextern[3], double lorentzForce[3]);
void calcRadiationDamping(double Eextern[3], double Bextern[3], double u[4], double dampingTerm[4]);
void calcUPMLCoefficients(Grid *Grid);
void calcLWFieldsForPlaneInNearFieldRegion(Grid *Grid, Particle *Particles, int numberOfParticles, double t, int planeForPlotting);
#endif
