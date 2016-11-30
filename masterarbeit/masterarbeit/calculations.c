/* =============================================
 IMPLEMENTATION OF CALCULATIONS HEADER
 ================================================ */

#import <math.h>
#import <stdbool.h>
#import <stdio.h>
#include "calculations.h"
#include "particle.h"
#include "string.h"
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
void calculateIntersectionPoint(double xInside[4], double xOutside[4], double uInside[4], double uOutside[4], double xObserver[4], double intersectionPoint[4], double velocityAtIntersectionPoint[4]){
    double lambda = calculateLambdaForLinearInterpolation(xInside, xOutside, xObserver);
    for(int i = 0; i < 4; i++){
        intersectionPoint[i] = xOutside[i] + lambda * (xInside[i] - xOutside[i]);
        velocityAtIntersectionPoint[i] = uOutside[i] + lambda * (uInside[i] - uOutside[i]);
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
    
    Particle->u[0] = getGammaFromVelocityVector(Particle->u);
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
///@brief calcualtes Liénard-Wiechert Fields everywhere on the grid. Each point on the grid is considered to be the observation point, where the zeroth component is the current simulation time and also the time at which we want to calculate the fields. Therefore loop through the particle xhistory vector up to the current simulation time and search for the pair of positions which fulfil the condition that the old position is inside and the new position is outside the backward lightcone of the observation point. Once that pair of positions is found the rest of the xHistory vector can be skipped because all following points will be outside as well.
///@param timeStep outer loop index indicating the current simulation time
void calcLWFieldsEverywhereOnGrid(Grid *Grid, Particle *Particle, int timeStep){
    printf("Calculating LW Fields on grid\n");
    double xObserver[4] = {0};
    double beta[3] = {0};
    double intersectionPoint[4] = {0};
    double velocityAtIntersectionPoint[4] = {0};
    double gamma_sq;
    double R_sq;
    double R;
    double n[3] = {0};
    double betaDot[3] = {0};
    double dt = 0.5 * Grid->dx;
    double E[3] = {0};
    double B[3] = {0};
    
    int nx = Grid->numberOfGridPointsInX;
    int ny = Grid->numberOfGridPointsInY;
    int nz = Grid->numberOfGridPointsInZ;
    
    for(int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            for(int k = 115; k < nz; k++){
                
                xObserver[0] = timeStep * dt;
                xObserver[1] = (i + 0.5)* Grid->dx;
                xObserver[2] = (j + 0.5) * Grid->dy;
                xObserver[3] = (k + 0.5) * Grid->dz;
                
                for (int index = 0; index < timeStep; index ++){
                    if(isInsideBackwardLightcone(Particle->xHistory[index], xObserver) && !isInsideBackwardLightcone(Particle->xHistory[index+1], xObserver)){
                        calculateIntersectionPoint(Particle->xHistory[index], Particle->xHistory[index+1], Particle->uHistory[index], Particle->uHistory[index+1], xObserver, intersectionPoint, velocityAtIntersectionPoint);
                        calculateBeta(Particle->xHistory[index], Particle->xHistory[index+1], beta);
                        calculateLienardWiechertParameters(intersectionPoint, xObserver, velocityAtIntersectionPoint, &gamma_sq, &R_sq, &R, n, beta);
                        calculateBetaDot(Particle->uHistory[index], Particle->uHistory[index+1], dt, betaDot);
                        calcuateLienardWiechertFields(gamma_sq, R_sq, R, n, beta, betaDot, Particle->charge, E, B);
                        Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 0] = E[0];
                        Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1] = E[1];
                        Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2] = E[2];
                        Grid->B[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 0] = B[0];
                        Grid->B[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1] = B[1];
                        Grid->B[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2] = B[2];
                        //printf("%d %d %d\n", i,j,k);
                        break;
                    }
                }
            }
        }
    }
    
    
}

///@brief calculates LW fields
void calcLWFieldsOnGrid(Grid *Grid, Particle *Particle, double t){
    int totalNumberOfBoxes = Grid->numberOfBoxesInX * Grid->numberOfBoxesInY * Grid->numberOfBoxesInZ;
    for (int boxIndex = 0; boxIndex < totalNumberOfBoxes; boxIndex++){
        addLWFieldsInBox(Grid, Particle, boxIndex, t);
    }
}

void addLWFieldsInBox(Grid *Grid, Particle *Particle, int boxIndex, double t){
    printf("calculating LW fields in box %d\n", boxIndex);
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    double dx = Grid->dx;
    double dy = Grid->dy;
    double dz = Grid->dz;
    
    int ib = boxIndex / (numberOfBoxesInY * numberOfBoxesInZ);
    int jb = (boxIndex - ( ib * numberOfBoxesInY * numberOfBoxesInZ)) / numberOfBoxesInY;
    int kb = (boxIndex - ( ib * numberOfBoxesInY * numberOfBoxesInZ)) - jb * numberOfBoxesInY;
    
    double xObserver[4] = {0};
    xObserver[0] = t;
    
    int lowerLeftGridIndexInBox = ib * numberOfGridPointsForBoxInX * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + jb * numberOfGridPointsForBoxInY * numberOfGridPointsInZ * 3 + kb * numberOfGridPointsForBoxInZ * 3;
    
    for (int id = 0; id < numberOfGridPointsForBoxInX; id++)
        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++)
            for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++)
            {
                int gridIndexInBox = lowerLeftGridIndexInBox + 3 * kd + 3 * jd * numberOfGridPointsInZ + 3 * id * numberOfGridPointsInZ * numberOfGridPointsInY;
                
                xObserver[1] = (ib * numberOfGridPointsForBoxInX + id) * dx;
                xObserver[2] = (jb * numberOfGridPointsForBoxInY + jd) * dy;
                xObserver[3] = (kb * numberOfGridPointsForBoxInZ + kd) * dz;
                
                AddLWField(Grid, xObserver, 3, gridIndexInBox, Particle);
                AddLWField(Grid, xObserver, 4, gridIndexInBox, Particle);
                AddLWField(Grid, xObserver, 5, gridIndexInBox, Particle);
                
                AddLWField(Grid, xObserver, 0, gridIndexInBox, Particle);
                AddLWField(Grid, xObserver, 1, gridIndexInBox, Particle);
                AddLWField(Grid, xObserver, 2, gridIndexInBox, Particle);
                
            }
}

void AddLWField(Grid *Grid, double xObserver[4], int component, int gridIndexInBox, Particle *Particle){
    double xObserverCopy[4];
    memcpy(xObserverCopy, xObserver, 4 * sizeof(double));
    switch (component)
    {
        case 0:
            xObserverCopy[1] += 0.5 * (Grid->dx);
            break;
        case 1:
            xObserverCopy[2] += 0.5 * (Grid->dy);
            break;
        case 2:
            xObserverCopy[3] += 0.5 * (Grid->dz);
            break;
        case 3:
            xObserverCopy[2] += 0.5 * (Grid->dy);
            xObserverCopy[3] += 0.5 * (Grid->dz);
            break;
        case 4:
            xObserverCopy[3] += 0.5 * (Grid->dz);
            xObserverCopy[1] += 0.5 * (Grid->dx);
            break;
        case 5:
            xObserverCopy[1] += 0.5 * (Grid->dx);
            xObserverCopy[2] += 0.5 * (Grid->dy);
            break;
    }
    
    
    double beta[3] = {0};
    double intersectionPoint[4] = {0};
    double velocityAtIntersectionPoint[4] = {0};
    double gamma_sq;
    double R_sq;
    double R;
    double n[3] = {0};
    double betaDot[3] = {0};
    double dt = 0.5 * Grid->dx;
    double E[3] = {0};
    double B[3] = {0};
    int currentHistoryLength = Particle->currentHistoryLength;
    
    
    for (int index = 0; index < currentHistoryLength - 1; index ++){
        if(isInsideBackwardLightcone(Particle->xHistory[index], xObserverCopy) && !isInsideBackwardLightcone(Particle->xHistory[index+1], xObserverCopy)){
            calculateIntersectionPoint(Particle->xHistory[index], Particle->xHistory[index+1], Particle->uHistory[index], Particle->uHistory[index+1], xObserverCopy, intersectionPoint, velocityAtIntersectionPoint);
            calculateBeta(Particle->xHistory[index], Particle->xHistory[index+1], beta);
            calculateLienardWiechertParameters(intersectionPoint, xObserverCopy, velocityAtIntersectionPoint, &gamma_sq, &R_sq, &R, n, beta);
            calculateBetaDot(Particle->uHistory[index], Particle->uHistory[index+1], dt, betaDot);
            calcuateLienardWiechertFields(gamma_sq, R_sq, R, n, beta, betaDot, Particle->charge, E, B);
            break;
        }
    }
    
    switch (component)
    {
        case 0:
            Grid->E[gridIndexInBox] += E[0];
            break;
        case 1:
            Grid->E[gridIndexInBox + 1] += E[1];
            break;
        case 2:
            Grid->E[gridIndexInBox + 2] += E[2];
            break;
        case 3:
            Grid->B[gridIndexInBox] += B[0];
            break;
        case 4:
            Grid->B[gridIndexInBox + 1] += B[1];
            break;
        case 5:
            Grid->B[gridIndexInBox + 2] += B[2];
            break;
    }
    
}

void calcLWFieldsForPlane(Grid *Grid, Particle *Particle, double t, int planeForPlotting){
    printf("Calculating LW Fields on plane %d, heigth:%f\n", planeForPlotting, planeForPlotting * Grid->dz);

    double xObserver[4];
    int k = planeForPlotting;
    int nx = Grid->numberOfGridPointsInX;
    int ny = Grid->numberOfGridPointsInY;
    int nz = Grid->numberOfGridPointsInZ;
    
    for(int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            
            
            xObserver[0] = t;
            xObserver[1] = (i)* Grid->dx;
            xObserver[2] = (j) * Grid->dy;
            xObserver[3] = (k) * Grid->dz;
            
            double gridIndexInBox = 3 * ny * nz * i + 3 * nz * j + 3 * k;
            
            AddLWField(Grid, xObserver, 3, gridIndexInBox, Particle);
            AddLWField(Grid, xObserver, 4, gridIndexInBox, Particle);
            AddLWField(Grid, xObserver, 5, gridIndexInBox, Particle);
            
            AddLWField(Grid, xObserver, 0, gridIndexInBox, Particle);
            AddLWField(Grid, xObserver, 1, gridIndexInBox, Particle);
            AddLWField(Grid, xObserver, 2, gridIndexInBox, Particle);
            
            
        }
    }
    
    
}



