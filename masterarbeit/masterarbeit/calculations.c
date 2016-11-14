/* =============================================
IMPLEMENTATION OF CALCULATIONS HEADER
================================================ */

#import <math.h>
#import <stdbool.h>
#import <stdio.h>
#include "calculations.h"
#include "particle.h"
/**
@brief calculates the cross product of two vectors "a" and "b" and saves the result in array "result"
@param a first vector of cross product
@param b second vector of cross product
@param result the result vector which  is automatically available in outer scope

 */
void crossProduct(double a[3], double b[3], double result[3]){
	result[0] = a[1]*b[2]-a[2]*b[1];
	result[1] = a[2]*b[0]-a[0]*b[2];
	result[2] = a[0]*b[1]-a[1]*b[0];
}

/**
@brief true if position of particle is inside the backward lightcone of the observer, otherwise false.
@remark two conditions have to be fulfilled. Since the information needs time to travel from particle position to observer position, time difference has to be positive: 
@code xObserver[0]-xParticle[0] > 0
@endcode
And Since
@code c > v <=> (ct)^2 > x^2
@endcode
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
@brief returns true if position of particle is inside the forward lightcone of the observer, otherwise false.
@remark two conditions have to be fulfilled. Since the information needs time to travel from particle position to observer position, time difference has to be negative:
@code xObserver[0]-xParticle[0] < 0
@endcode
And Since
@code c > v <=> (ct)^2 > x^2
@endcode
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
@brief calcualtes the vector difference of spatial dimensions of the two four vectors and saves the result in the array passed as argument.
@param x vector you want to substract from
@param y vector to substract from x
@param result the resulting vector of substracting y from x
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
 returns the minkowski product with metric (+,-,-,-)
 */
double minkowskiProduct(double x[4], double y[4]){
    return x[0] * y[0] - x[1] * y[1] - x[2] * y[2] - x[3] * y[3];
}

/**
 returns parameter lambda for the linear interpolation between particle trajectory and backwards lightcone at observation point.
 
 In order to calculate the Liénard-Wiechart fields at the observation point the intersection of the trajectory and the backward lightcone at the observation point is needed. Since the calculated trajectory is discrete we need to interpolate between the first point outside the lightcone and the last point inside the lightcone.
 */
double calculateLambdaForLinearInterpolation(double xInside[4], double xOutside[4], double xObserver[4]) {
    
    double a, b, c;
    double lambda;
    double xInsideMinusxOutside[4];
    double xOutsideMinusxObserver[4];
    
    vectorDifference(xInside, xOutside, xInsideMinusxOutside);
    vectorDifference(xOutside, xObserver, xOutsideMinusxObserver);
    
    a = minkowskiProduct(xInsideMinusxOutside, xInsideMinusxOutside);
    b = 2.0 * minkowskiProduct(xInsideMinusxOutside, xOutsideMinusxObserver);
    c = minkowskiProduct(xOutside, xOutside) + minkowskiProduct(xObserver, xObserver) - 2.0 * minkowskiProduct(xOutside, xObserver);

    lambda = (-b+sqrt(b*b-4.0*a*c))/(2.0*a);
    return lambda;
}

/**
 calculates the intersection point of the particle trajectory and the backward lightcone of the observation point with linear interpolation.
 */
void calculateIntersectionPoint(double xInside[4], double xOutside[4], double xObserver[4], double intersectionPoint[4]){
    double lambda = calculateLambdaForLinearInterpolation(xInside, xOutside, xObserver);
    for(int i = 0; i < 4; i++){
        intersectionPoint[i] = xOutside[i] + lambda * (xInside[i] - xOutside[i]);
    }
    
}

/**
 updates velocity with boris method. 
 
 First obtain “uMinus” by adding half acceleration to the initial velocity.
 Then obtain "uPlus" by performing rotation with "uPrime" and internally computed assisting values.
 Finally, add another half acceleration.
 */
void updateVelocityWithBorisPusher(Particle *Particle, double *Eextern, double *Bextern, double dt){

    int dimension = 3;
    double uPrime[3];
    double t[3];
    double s[3];
    double uMinusCrossT[3];
    double uPrimeCrossS[3];
    double absoluteSquareValueOfT;
    double uModifiedForCrossProduct[3];
    double chargeOverMass = Particle->charge / Particle->mass;
    
    // assisting values
    for (int i = 0; i < dimension; i++){
        t[i] = chargeOverMass * Bextern[i] * dt * 0.5;
    }
    
    absoluteSquareValueOfT = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
    
    for(int i = 0; i < dimension; i++){
        s[i] = 2 * t[i] / (1 + absoluteSquareValueOfT);
    }
    
    // implementation of boris method
    // first obtain “uMinus” by adding half acceleration to the initial velocity
    for (int i = 0; i < dimension; i++){
        Particle->u[i+1] +=  chargeOverMass * Eextern[i] * dt * 0.5;
    }
    
    // Since u is a four vector rewrite it to a three dimensional vector to make it usabale for crossProduct
    for(int i = 0; i < dimension; i++){
        uModifiedForCrossProduct[i] = Particle->u[i+1];
    }
    
    crossProduct(uModifiedForCrossProduct, t, uMinusCrossT);
    
    for (int i = 0; i < dimension; i++){
        uPrime[i] = Particle->u[i+1] + uMinusCrossT[i];
    }
    
    crossProduct(uPrime, s, uPrimeCrossS);
    
    // then obtain "uPlus" by performing rotation with "uPrime" and "s"
    for (int i = 0; i < dimension; i++){
        Particle->u[i+1] +=  uPrimeCrossS[i];
    }
    
    // finally, add another half acceleration,
    for (int i = 0; i < dimension; i++){
        Particle->u[i+1] += chargeOverMass * Eextern[i] * dt * 0.5;
    }
    Particle->u[0] += dt;
}

/**
 updates location x for particle with velocity u using x = v * dt
 time component x[0] gets updated as well.
 */
void updateLocation(Particle *Particle, double dt){
    int dimension = 3;
    
    for(int i = 0; i < dimension; i++){
        Particle->x[i+1] += Particle->u[i+1] * dt;
    }
    Particle->x[0] += dt;
}

/**
@brief calcualtes Liénard Wiechert Parameters gamma_sq, R_sq, R, n, beta from xParticle, u and xObserver
@remark in order to have the values of R_sq, R and gamma_sq available in the outer scope, we need to work with pointers like we do when we passing arrays as function arguments
@param xParticle position vector of particle
@param xObserver observation point where you want to calcualte the lw fields
@param u velocity vector the your particle
@param gamma_sq gets calculated
@param R_sq gets calculated
@param R gets calculated
@param n gets calculated
@param beta gets calculated
 */
void calculateLienardWiechertParameters(double xParticle[4], double xObserver[4], double u[4], double *gamma_sq, double *R_sq, double *R, double n[3], double beta[3]) {
    n[0] = xObserver[1]-xParticle[1];
    n[1] = xObserver[2]-xParticle[2];
    n[2] = xObserver[3]-xParticle[3];
    
    *R_sq = n[0]*n[0]+n[1]*n[1]+n[2]*n[2];
    *R = sqrt(*R_sq);
    
    n[0] /= *R;
    n[1] /= *R;
    n[2] /= *R;
    
    *gamma_sq = u[0]*u[0];
    
    beta[0] = u[1]/u[0];
    beta[1] = u[2]/u[0];
    beta[2] = u[3]/u[0];
}


/**
 @brief calcultes time derivative in each component of beta. 
 @param uOld velocity vector of particle before it gets pushed
 @param uNew actual velocity vector of particle
 @param dt step size for calculating the derivative
 @param betaDot the result vector containing the derivatives
 */
void calculateBetaDot(double uOld[4], double uNew[4], double dt, double betaDot[3]) {
    betaDot[0] = (uNew[1]/uNew[0]-uOld[1]/uOld[0])/dt;
    betaDot[1] = (uNew[2]/uNew[0]-uOld[2]/uOld[0])/dt;
    betaDot[2] = (uNew[3]/uNew[0]-uOld[3]/uOld[0])/dt;
}

/**
@brief calcuates acceleration of particle
@param xOld position of particle before it was pushed via borisPusher
@param xNew actual position of particle
@param beta acceleration vector of spatial components
*/
void calculateBeta(double xOld[4], double xNew[4], double beta[3]) {
    double dt = xNew[0]-xOld[0];
    beta[0] = (xNew[1]-xOld[1])/dt;
    beta[1] = (xNew[2]-xOld[2])/dt;
    beta[2] = (xNew[3]-xOld[3])/dt;
}


void calcuateLienardWiechertFields(double gamma_sq, double R_sq, double R, double *n,
                    double *beta, double *beta_dot, double charge, double *E,
                    double *B) {
    double denominator1;
    double denominator2;
    double oneMinusBetaN;
    double oneMinuBetaNCubed;
    double nMinusBeta[3];
    double nCrossnMinusBetaCrossBetaDot[3];
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


