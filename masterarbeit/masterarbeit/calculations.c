/* =============================================
IMPLEMENTATION OF CALCULATIONS HEADER
================================================ */

#import <math.h>
#import <stdbool.h>
#import <stdio.h>
#include "calculations.h"

/**
calculates the cross product of two vectors "a" and "b" and saves the result in array "result"
 
 - remark:
 array "result" is automatically available in outer scope
 */
void crossProduct(double *a, double *b, double *result)
{
	result[0] = a[1]*b[2]-a[2]*b[1];
	result[1] = a[2]*b[0]-a[0]*b[2];
	result[2] = a[0]*b[1]-a[1]*b[0];
}

/**
returns true if position of particle is inside the backward lightcone of the observer, otherwise false.

 - remark: two conditions have to be fulfilled. Since the information needs time to travel from particle position to observer position, time difference has to be positive: xObserver[0]-xParticle[0] > 0.
 And Since c > v <=> (ct)^2 > x^2
 */
bool isInsideBackwardLightcone(double xParticle[4], double xObserver[4]) {
    double dt;
    double dxsq;
    double xObserverMinusxParticle[4];
    
    vectorDifference(xObserver, xParticle, xObserverMinusxParticle);
    dt = xObserverMinusxParticle[0];
    
    dxsq = vectorProduct(xObserverMinusxParticle, xObserverMinusxParticle);
    
    return (dt > 0 && dt*dt > dxsq);
}

/**
 returns true if position of particle is inside the forward lightcone of the observer, otherwise false.
 
 - remark: two conditions have to be fulfilled. Since the information needs time to travel from particle position to observer position, time difference has to be negative: xObserver[0]-xParticle[0] < 0.
 And Since c > v <=> (ct)^2 > x^2
 */
bool isInsideForwardLightcone(double xParticle[4], double xObserver[4]) {
    
    double dt;
    double dxsq;
    double xObserverMinusxParticle[4];
    
    vectorDifference(xObserver, xParticle, xObserverMinusxParticle);
    dt = xObserverMinusxParticle[0];
    
    dxsq = vectorProduct(xObserverMinusxParticle, xObserverMinusxParticle);
    
    return (dt < 0 && dt*dt > dxsq);
}


/**
 calcualtes the vector difference of spatial dimensions of the two four vectors and saves the result in the array passed as argument.
 */
void vectorDifference(double x[4], double y[4], double result[4]){
    for (int i = 0; i < 4; i++){
        result[i] = x[i] - y[i];
    }
}


/**
 returns the vector product in three dimensions.
 */
double vectorProduct(double x[4], double y[4]){
    
    return x[1]*y[1] + x[2]*y[2] + x[3]*y[3];
}

/**
 returns the spatial distance between two four vectors.
 */
double calculateDistance(double *x, double *y){
    double difference[4];
    vectorDifference(x, y, difference);
    return sqrt(difference[1]*difference[1] + difference[2]*difference[2] + difference[3]*difference[3]);
}

/**
 returns parameter lambda for the linear interpolation between particle trajectory and backwards lightcone at observation point.
 
 In order to calculate the Liénard-Wiechart fields at the observation point the intersection of the trajectory and the backward lightcone at the observation point is needed. Since the calculated trajectory is discrete we need to interpolate between the first point outside the lightcone and the last point inside the lightcone.
 */
double calculateLambdaForLinearInterpolation(double xInside[4], double xOutside[4], double xObserver[4]) {
    
    double a, b, c;
    double lambda;
    double xInsideMinusxOutside[4];
    double xObserverMinusxOutside[4];
    
    vectorDifference(xInside, xOutside, xInsideMinusxOutside);
    vectorDifference(xObserver, xOutside, xObserverMinusxOutside);
    
    a = vectorProduct(xInsideMinusxOutside, xInsideMinusxOutside);
    b = 2.0 * vectorProduct(xInsideMinusxOutside, xObserverMinusxOutside);
    c = vectorProduct(xOutside, xOutside) + vectorProduct(xObserver, xObserver) - 2.0 * vectorProduct(xOutside, xObserver);
    
    lambda = (-b+sqrt(b*b-4.0*a*c))/(2.0*a);
    return lambda;
}
/**
 updates velocity with boris method. 
 
 First obtain “uMinus” by adding half acceleration to the initial velocity.
 Then obtain "uPlus" by performing rotation with "uPrime" and internally computed assisting values.
 Finally, add another half acceleration.
 */
void borisPusher(double *u, double *E, double *B, double dt, double chargeOverMass){

    int dimension = 3;
    double uPrime[3];
    double t[3];
    double s[3];
    double uMinusCrossT[3];
    double uPrimeCrossS[3];
    double absoluteSquareValueOfT;
    double uModifiedForCrossProduct[3];
    
    
    // assisting values
    for (int i = 0; i < dimension; i++){
        t[i] = chargeOverMass * B[i] * dt * 0.5;
    }
    
    absoluteSquareValueOfT = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
    
    for(int i = 0; i < dimension; i++){
        s[i] = 2 * t[i] / (1 + absoluteSquareValueOfT);
    }
    
    // implementation of boris method
    // first obtain “uMinus” by adding half acceleration to the initial velocity
    for (int i = 0; i < dimension; i++){
        u[i+1] +=  chargeOverMass * E[i] * dt * 0.5;
    }
    
    // Since u is a four vector rewrite it to a three dimensional vector to make it usabale for crossProduct
    for(int i = 0; i < dimension; i++){
        uModifiedForCrossProduct[i] = u[i+1];
    }
    
    crossProduct(uModifiedForCrossProduct, t, uMinusCrossT);
    
    for (int i = 0; i < dimension; i++){
        uPrime[i] = u[i+1] + uMinusCrossT[i];
    }
    
    crossProduct(uPrime, s, uPrimeCrossS);
    
    // then obtain "uPlus" by performing rotation with "uPrime" and "s"
    for (int i = 0; i < dimension; i++){
        u[i+1] +=  uPrimeCrossS[i];
    }
    
    // finally, add another half acceleration,
    for (int i = 0; i < dimension; i++){
        u[i+1] += chargeOverMass * E[i] * dt * 0.5;
    }
    u[0] += dt;
}

/**
 updates location x for particle with velocity u using x = v * dt
 time component x[0] gets updated as well.
 */
void updateLocation(double *u, double *x, double dt){
    int dimension = 3;
    
    for(int i = 0; i < dimension; i++){
        x[i+1] += u[i+1] * dt;
    }
    x[0] += dt;
}

void calculateLienardWiechertParameters(double gamma_sq, double R_sq, double R, double *n, double *beta, double *xParticle, double *u, double *xObserver) {
    n[0] = xObserver[1]-xParticle[1];
    n[1] = xObserver[2]-xParticle[2];
    n[2] = xObserver[3]-xParticle[3];
    
    R_sq = n[0]*n[0]+n[1]*n[1]+n[2]*n[2];
    R = sqrt(R_sq);
    
    n[0] /= R;
    n[1] /= R;
    n[2] /= R;
    
    gamma_sq = u[0]*u[0];
    
    beta[0] = u[1]/u[0];
    beta[1] = u[2]/u[0];
    beta[2] = u[3]/u[0];
}


void calculateBetaDot(double *betaDot, double *u1, double *u2,
                  double dt) {
    betaDot[0] = (u2[1]/u2[0]-u1[1]/u1[0])/dt;
    betaDot[1] = (u2[2]/u2[0]-u1[2]/u1[0])/dt;
    betaDot[2] = (u2[3]/u2[0]-u1[3]/u1[0])/dt;
}

void calculateBeta(double *beta, double *x1, double *x2) {
    double dt = x2[0]-x1[0];
    beta[0] = (x2[1]-x1[1])/dt;
    beta[1] = (x2[2]-x1[2])/dt;
    beta[2] = (x2[3]-x1[3])/dt;
}


void calcuateLienardWiechertFields(double gamma_sq, double R_sq, double R, double *n,
                    double *beta, double *beta_dot, double charge, double *E,
                    double *B) {
    double denominator1;
    double denominator2;
    double oneMinusBetaN;
    double oneMinuBetaNCubed;
    double nMinusBeta[3], nCrossnMinusBetaCrossBetaDot[3];
    double nMinusBetaCrossBetaDot[3];
    
    if (R_sq==0||R==0) {
        E[0] = 0;
        E[1] = 0;
        E[2] = 0;
        B[0] = 0;
        B[1] = 0;
        B[2] = 0;
        return;
    }
    
    oneMinusBetaN = 1.0-(beta[0]*n[0]+beta[1]*n[1]+beta[2]*n[2]);
    oneMinuBetaNCubed = oneMinusBetaN*oneMinusBetaN*oneMinusBetaN;
    denominator1 = 1.0/(gamma_sq*oneMinuBetaNCubed*R_sq);
    denominator2 = 1.0/(oneMinuBetaNCubed*R);
    nMinusBeta[0] = n[0]-beta[0];
    nMinusBeta[1] = n[1]-beta[1];
    nMinusBeta[2] = n[2]-beta[2];
    
    crossProduct(nMinusBeta, beta_dot, nMinusBetaCrossBetaDot);
    crossProduct(n, nMinusBetaCrossBetaDot, nCrossnMinusBetaCrossBetaDot);
    
    E[0] = charge*(denominator1*nMinusBeta[0]+denominator2*nCrossnMinusBetaCrossBetaDot[0]);
    E[1] = charge*(denominator1*nMinusBeta[1]+denominator2*nCrossnMinusBetaCrossBetaDot[1]);
    E[2] = charge*(denominator1*nMinusBeta[2]+denominator2*nCrossnMinusBetaCrossBetaDot[2]);
    
    crossProduct(n, E, B);
}


