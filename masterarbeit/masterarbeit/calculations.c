/* =============================================
 IMPLEMENTATION OF CALCULATIONS HEADER
 ================================================ */

#import <math.h>
#import <stdbool.h>
#import <stdio.h>
#include "calculations.h"
#include "string.h"
#include "lwFields.h"

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


void scaleVector(double x[3], double factor){
    for (int i = 0; i < 3; i++){
        x[i] *= factor;
    }
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
void updateVelocityWithBorisPusher(Particle *Particles, Grid *Grid, int numberOfParticles, int particleIndex, double Eextern[3], double Bextern[3], double dt){
    
    Particle *Particle = &Particles[particleIndex];
    
    int dimension = 3;
    double uPrime[3];
    double t[3];
    double s[3];
    double uMinusCrossT[3];
    double uPrimeCrossS[3];
    double E[3];
    double B[3];
    double absoluteSquareValueOfT;
    double uModifiedForCrossProduct[3];
    double chargeOverMass = Particle->charge / Particle->mass;
    
    updateFieldsForParticlePush(Particle, Grid, Eextern, Bextern, E, B);
    calcInteractionWithOtherParticles(Particles, Grid, numberOfParticles, particleIndex, E, B);
    
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
        Particle->u[i+1] +=  chargeOverMass * E[i] * dt * 0.5;
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
        Particle->u[i+1] += chargeOverMass * E[i] * dt * 0.5;
    }
    
    Particle->u[0] = getGammaFromVelocityVector(Particle->u);
    
}



///@brief updates location x for particle with velocity u using x = v * dt time component x[0] gets updated as well.
///@remark when particle is pushed the old and new box indeces are calculted. When this index changes the Paticle property "didChangeBox" is set accordingly.

void updateLocation(Particle *Particle, Grid *Grid, double dt){
    int dimension = 3;
    
    int oldBoxIndexOfParticle = calcCurrentBoxIndexOfParticle(Particle, Grid);
    calcBoxIndizesOfNextNeighbourBoxes(Grid, Particle, Particle->boxIndicesOfNearFieldBoxesBeforePush);
    for(int i = 0; i < dimension; i++){
        Particle->x[i+1] += Particle->u[i+1] * dt;
    }
    int newBoxIndexOfParticle = calcCurrentBoxIndexOfParticle(Particle, Grid);
    calcBoxIndizesOfNextNeighbourBoxes(Grid, Particle, Particle->boxIndicesOfNearFieldBoxesAfterPush);
    
    if(oldBoxIndexOfParticle != newBoxIndexOfParticle){
        Particle->didChangeBox = true;
    }
    else{
        Particle->didChangeBox = false;
    }
    Particle->x[0] += dt;
    
}

///@brief updates location x for particles with velocity u using x = v * dt time component x[0] gets updated as well.
///@remark when particle is pushed the old and new box indeces are calculted. When this index changes the Paticle property "didChangeBox" is set accordingly.
void updateLocationForParticles(Particle *Particles, int numberOfParticles, Grid *Grid, double dt){
    for(int p = 0; p < numberOfParticles; p++){
        updateLocation(&Particles[p], Grid, dt);
    }
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
                        Grid->H[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 0] = B[0];
                        Grid->H[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1] = B[1];
                        Grid->H[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2] = B[2];
                        //printf("%d %d %d\n", i,j,k);
                        break;
                    }
                }
            }
        }
    }
    
    
}

///@brief calculates LW fields on complete grid by looping through every box and calculating LW fields in there.
///@param Grid pointer to an instance of a Grid struct.
///@param Particle pointer to an instance of a Particle struct
///@param t time at which the fields shall be calcualted. This is also the zeroth component of xObserver
void calcLWFieldsOnGrid(Grid *Grid, Particle *Particle, double t){
    int totalNumberOfBoxes = Grid->numberOfBoxesInX * Grid->numberOfBoxesInY * Grid->numberOfBoxesInZ;
    for (int boxIndex = 0; boxIndex < totalNumberOfBoxes; boxIndex++){
        addLWFieldsInBox(Grid, Particle, boxIndex, t);
    }
}


///@remark calcuates LW fields at simulation time t only on specyfied plane with index "planeForPlotting".
///@param Grid pointer to an instance of a Grid struct.
///@param Particles pointer to  Particle struct array
///@param t time at which the fields shall be calcualted. This is also the zeroth component of xObserver
///@param planeForPlotting number of the plane in which fields shall be calcualted. planeForPlotting = x[3] / dz
void calcLWFieldsForPlane(Grid *Grid, Particle *Particles, int numberOfParticles, double t, int planeForPlotting){
    for (int p = 0; p < numberOfParticles; p++){
        printf("Calculating LW Fields on plane %d, heigth:%f for Particle%d\n", planeForPlotting, planeForPlotting * Grid->dz, p);
        
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
                
                int gridIndexInBox = 3 * ny * nz * i + 3 * nz * j + 3 * k;
                
                addLWField(Grid, &Particles[p], &Grid->H[gridIndexInBox], xObserver, 3);
                addLWField(Grid, &Particles[p], &Grid->H[gridIndexInBox + 1], xObserver, 4);
                addLWField(Grid, &Particles[p], &Grid->H[gridIndexInBox + 2], xObserver, 5);
                
                addLWField(Grid, &Particles[p], &Grid->E[gridIndexInBox], xObserver, 0);
                addLWField(Grid, &Particles[p], &Grid->E[gridIndexInBox + 1], xObserver, 1);
                addLWField(Grid, &Particles[p], &Grid->E[gridIndexInBox + 2], xObserver, 2);
                
                
            }
        }
    }
    
}


void calcLWFieldsForPlaneWithNearField(Grid *Grid, Particle *Particle, double t, int planeForPlotting){
    printf("Calculating LW Fields on plane %d, heigth:%f\n", planeForPlotting, planeForPlotting * Grid->dz);
    
    double xObserver[4];
    int k = planeForPlotting;
    
    int nx = Grid->numberOfGridPointsInX;
    int ny = Grid->numberOfGridPointsInY;
    int nz = Grid->numberOfGridPointsInZ;
    
    getEdgesOfNearFieldBox(Grid, Particle);
    
    for(int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            
            
            xObserver[0] = t;
            xObserver[1] = (i)* Grid->dx;
            xObserver[2] = (j) * Grid->dy;
            xObserver[3] = (k) * Grid->dz;
            
            if(xObserver[1] >= Particle->edgesOfNearFieldBox[0] && xObserver[1] <= Particle->edgesOfNearFieldBox[1] && xObserver[2] >= Particle->edgesOfNearFieldBox[2] && xObserver[2] <= Particle->edgesOfNearFieldBox[3] && xObserver[3] >= Particle->edgesOfNearFieldBox[4] && xObserver[3] <= Particle->edgesOfNearFieldBox[5]){
                continue;
            }
            
            int gridIndexInBox = 3 * ny * nz * i + 3 * nz * j + 3 * k;
            
            addLWField(Grid, Particle, &Grid->H[gridIndexInBox], xObserver, 3);
            addLWField(Grid, Particle, &Grid->H[gridIndexInBox + 1], xObserver, 4);
            addLWField(Grid, Particle, &Grid->H[gridIndexInBox + 2], xObserver, 5);
            
            addLWField(Grid, Particle, &Grid->E[gridIndexInBox], xObserver, 0);
            addLWField(Grid, Particle, &Grid->E[gridIndexInBox + 1], xObserver, 1);
            addLWField(Grid, Particle, &Grid->E[gridIndexInBox + 2], xObserver, 2);
            
            
        }
    }
    
    
}

///@brief calcualtes the box index in which the particle is currently located.
///@param Grid pointer to an instance of a Grid struct.
///@param Particle pointer to an instance of a Particle struct
///@return index of box where particle is currenty located.
///@remark index for box is counted the same way as for E and B field array. First z component, then y and then x.
int calcCurrentBoxIndexOfParticle(Particle *Particle, Grid *Grid){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    double dx = Grid->dx;
    double dy = Grid->dy;
    double dz = Grid->dz;
    
    int ib = Particle->x[1] / dx / numberOfGridPointsForBoxInX;
    int jb = Particle->x[2] / dy / numberOfGridPointsForBoxInY;
    int kb = Particle->x[3] / dz / numberOfGridPointsForBoxInZ;
    
    return ib * numberOfBoxesInY * numberOfBoxesInZ + jb * numberOfBoxesInZ + kb;
}

///@brief calcualtes the indices of next neighbour boxes from current particle position and saves them in boxIndizesOfNextNeighbourBoxes[27] array. If index is smaller than zero or larger than the maximal number of boxes then index is set to -1.
///@param Grid pointer to an instance of a Grid struct.
///@param Particle pointer to an instance of a Particle struct
///@param boxIndizesOfNextNeighbourBoxes array in which the indicies will be stored.
void calcBoxIndizesOfNextNeighbourBoxes(Grid *Grid, Particle *Particle, int boxIndizesOfNextNeighbourBoxes[27]){
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    
    int currentBoxIndexOfParticle = calcCurrentBoxIndexOfParticle(Particle, Grid);
    int indexOfNextNeighbourBox = 0;
    int maxBoxIndex = numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ;
    int index = 0;
    
    for (int ii = -1; ii <= 1; ii++){
        for(int ij = -1; ij <= 1; ij++){
            for(int ik = -1; ik <= 1; ik++){
                indexOfNextNeighbourBox = currentBoxIndexOfParticle + (ik + numberOfBoxesInZ * ij + numberOfBoxesInZ * numberOfBoxesInY * ii);
                if (indexOfNextNeighbourBox < 0 || indexOfNextNeighbourBox > maxBoxIndex){
                    boxIndizesOfNextNeighbourBoxes[index] = -1;
                }
                else{
                    boxIndizesOfNextNeighbourBoxes[index] = indexOfNextNeighbourBox;
                }
                index++;
            }
        }
        
    }
}

///@param Grid pointer to an instance of a Grid struct.
///@param Particle pointer to an instance of a Particle struct
///@param boxIndex current box index from outter loop
///@return true if box with Index "boxIndex" is in near field of particle or false if not.
bool boxIsInNearFieldOfParticle(Grid *Grid, Particle *Particle, int boxIndex){
    int boxIndizesOfNextNeighbourBoxes[27] = {0};
    calcBoxIndizesOfNextNeighbourBoxes(Grid, Particle, boxIndizesOfNextNeighbourBoxes);
    
    bool isInNF = false;
    
    for (int i = 0; i < 27; i++){
        if(boxIndizesOfNextNeighbourBoxes[i] == boxIndex){
            isInNF = true;
            break;
        }
    }
    return isInNF;
}

int calcLowerLeftGridIndexInBox(Grid *Grid, Particle *Particle){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    getCurrentBoxIndexArrayOfParticle(Grid, Particle);
    int ib = Particle->currentBoxIndexArray[0];
    int jb = Particle->currentBoxIndexArray[1];
    int kb = Particle->currentBoxIndexArray[2];
    
    return ib * numberOfGridPointsForBoxInX * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + jb * numberOfGridPointsForBoxInY * numberOfGridPointsInZ * 3 + kb * numberOfGridPointsForBoxInZ * 3;
}

///@param Grid pointer to Grid struct.
///@param boxIndex current boxIndex of particle
///@returns the boxIndex to the left of the current box where particle is located.
int calcBoxIndexIm1(Grid *Grid, const int boxIndex){
    return boxIndex - (Grid->numberOfBoxesInZ) * (Grid->numberOfBoxesInY);
}

///@param Grid pointer to Grid struct.
///@param boxIndex current boxIndex of particle
///@returns the boxIndex in front of the current box where particle is located.
int calcBoxIndexJm1(Grid *Grid, const int boxIndex){
    return boxIndex - (Grid->numberOfBoxesInZ);
}

///@param Grid pointer to Grid struct.
///@param boxIndex current boxIndex of particle
///@returns the boxIndex below of the current box where particle is located.
int calcBoxIndexKm1(Grid *Grid, const int boxIndex){
    return boxIndex - 1;
}

///@param Grid pointer to Grid struct.
///@param boxIndex current boxIndex of particle
///@returns the boxIndex to the right of the current box where particle is located.
int calcBoxIndexIp1(Grid *Grid, const int boxIndex){
    return boxIndex + (Grid->numberOfBoxesInZ) * (Grid->numberOfBoxesInY);
}

///@param Grid pointer to Grid struct.
///@param boxIndex current boxIndex of particle
///@returns the boxIndex behind the current box where particle is located.
int calcBoxIndexJp1(Grid *Grid, const int boxIndex){
    return boxIndex + (Grid->numberOfBoxesInZ);
}

///@param Grid pointer to Grid struct.
///@param boxIndex current boxIndex of particle
///@returns the boxIndex above the current box where particle is located.
int calcBoxIndexKp1(Grid *Grid, const int boxIndex){
    return boxIndex + 1;
}

///@brief calculates the UPML Coefficients according to Taflove page 301.
///@param Grid pointer to Grid struct
///@remark The idea is to introduce a medium at the edge of the grid which shall absrob all incident waves, regardless of angle, polarisation and impedance, so that there are no reflections at the border of the grid. A quite longish calculation yields the coefficients implmeneted below. Those coefficients need to be implemented into the Maxwell Equations so that in the inner part of the grid everything is calculated via Yee algorithm (sigma = 0, kappa = 1). At the edge however we want to absorb the incident waves quite strongly, so conductiity sigma and loss kappa need to change. One questions still remains unanswered: How can we calculate optimal values for sigmaMax and kappaMax ?
void calcUPMLCoefficients(Grid *Grid){
    printf("calculating UPML coefficients ...\n");
    
    double sigmaMax = 18.0;
    double kappaMax = 1.0;
    double m = 3.5;
    int upmlLayerWidth = Grid->upmlLayerWidth;
    
    double sigma = 0.0;
    double kappa = 0.0;
    double depth = 0.0;
    
    double dt = Grid->dx * 0.25;
    int numberOfGridPointsInX = Grid->numberOfGridPointsInX;
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    for (int j = 0; j < numberOfGridPointsInY; j++){
        if(j < upmlLayerWidth){
            depth = upmlLayerWidth - j;
            sigma = pow(depth / upmlLayerWidth, m) * sigmaMax;
            kappa = 1 + (kappaMax - 1) * pow(depth / upmlLayerWidth, m);
        }
        
        else if(j >= numberOfGridPointsInY - upmlLayerWidth - 1){
            depth = j - numberOfGridPointsInY + upmlLayerWidth + 1;
            sigma = pow(depth / upmlLayerWidth, m) * sigmaMax;
            kappa = 1 + (kappaMax - 1) * pow(depth / upmlLayerWidth, m);
        }
        else{
            sigma = 0.0;
            kappa = 1.0;
        }
        Grid->upml1E[j] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml2E[j] = (2 * dt) / (2 * kappa + sigma * dt);
        
    }
    
    for (int k = 0; k < numberOfGridPointsInZ; k++){
        if(k < upmlLayerWidth){
            depth = upmlLayerWidth - k;
            sigma = pow(depth / upmlLayerWidth, m) * sigmaMax;
            kappa = 1 + (kappaMax - 1) * pow(depth / upmlLayerWidth, m);
        }
        
        else if(k >= numberOfGridPointsInZ - upmlLayerWidth - 1){
            depth = k - numberOfGridPointsInZ + upmlLayerWidth + 1;
            sigma = pow(depth / upmlLayerWidth, m) * sigmaMax;
            kappa = 1 + (kappaMax - 1) * pow(depth / upmlLayerWidth, m);
        }
        else{
            sigma = 0.0;
            kappa = 1.0;
        }
        Grid->upml3E[k] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml4E[k] = 1.0 / (2 * kappa + sigma * dt);
        
    }
    
    for (int i = 0; i < numberOfGridPointsInX; i++){
        if(i < upmlLayerWidth){
            depth = upmlLayerWidth - i - 0.5;
            sigma = pow(depth / upmlLayerWidth, m) * sigmaMax;
            kappa = 1 + (kappaMax - 1) * pow(depth / upmlLayerWidth, m);
        }
        
        else if(i >= numberOfGridPointsInX - upmlLayerWidth){
            depth = i - numberOfGridPointsInX + upmlLayerWidth + 0.5;
            sigma = pow(depth / upmlLayerWidth, m) * sigmaMax;
            kappa = 1 + (kappaMax - 1) * pow(depth / upmlLayerWidth, m);
        }
        else{
            sigma = 0.0;
            kappa = 1.0;
        }
        Grid->upml5E[i] = (2 * kappa + sigma * dt);
        Grid->upml6E[i] = (2 * kappa - sigma * dt);
    }
    
    for (int j = 0; j < numberOfGridPointsInY; j++){
        if(j < upmlLayerWidth){
            depth = upmlLayerWidth - j - 0.5;
            sigma = pow(depth / upmlLayerWidth, m) * sigmaMax;
            kappa = 1 + (kappaMax - 1) * pow(depth / upmlLayerWidth, m);
        }
        
        else if(j >= numberOfGridPointsInY - upmlLayerWidth){
            depth = j - numberOfGridPointsInY + upmlLayerWidth + 0.5;
            sigma = pow(depth / upmlLayerWidth, m) * sigmaMax;
            kappa = 1 + (kappaMax - 1) * pow(depth / upmlLayerWidth, m);
        }
        else{
            sigma = 0.0;
            kappa = 1.0;
        }
        Grid->upml1H[j] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml2H[j] = (2 * dt) / (2 * kappa + sigma * dt);
        
    }
    
    for (int k = 0; k < numberOfGridPointsInZ; k++){
        if(k < upmlLayerWidth){
            depth = upmlLayerWidth - k - 0.5;
            sigma = pow(depth / upmlLayerWidth, m) * sigmaMax;
            kappa = 1 + (kappaMax - 1) * pow(depth / upmlLayerWidth, m);
        }
        
        else if(k >= numberOfGridPointsInZ - upmlLayerWidth){
            depth = k - numberOfGridPointsInZ + upmlLayerWidth + 0.5;
            sigma = pow(depth / upmlLayerWidth, m) * sigmaMax;
            kappa = 1 + (kappaMax - 1) * pow(depth / upmlLayerWidth, m);
        }
        else{
            sigma = 0.0;
            kappa = 1.0;
        }
        Grid->upml3H[k] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml4H[k] = 1.0 / (2 * kappa + sigma * dt);
        
    }
    
    
    for (int i = 0; i < numberOfGridPointsInX; i++){
        if(i < upmlLayerWidth){
            depth = upmlLayerWidth - i;
            sigma = pow(depth / upmlLayerWidth, m) * sigmaMax;
            kappa = 1 + (kappaMax - 1) * pow(depth / upmlLayerWidth, m);
        }
        
        else if(i >= numberOfGridPointsInX - upmlLayerWidth - 1){
            depth = i - numberOfGridPointsInX + upmlLayerWidth + 1;
            sigma = pow(depth / upmlLayerWidth, m) * sigmaMax;
            kappa = 1 + (kappaMax - 1) * pow(depth / upmlLayerWidth, m);
        }
        else{
            sigma = 0.0;
            kappa = 1.0;
        }
        Grid->upml5H[i] = (2 * kappa + sigma * dt);
        Grid->upml6H[i] = (2 * kappa - sigma * dt);
        
    }
}

void updateNearField(Grid *Grid, Particle *Particle, double t){
    
    if (Particle->didChangeBox == true){
        printf("updating NearField ...\n");
        
        for(int i = 0; i < 27; i++){
            if(Particle->boxIndicesOfNearFieldBoxesBeforePush[i] == -1 || Particle->boxIndicesOfNearFieldBoxesAfterPush[i] == -1){
                continue;
            }
            if(!boxIsInNearFieldOfParticle(Grid, Particle, Particle->boxIndicesOfNearFieldBoxesBeforePush[i])){
                addLWFieldsInBox(Grid, Particle, Particle->boxIndicesOfNearFieldBoxesBeforePush[i], t);
            }
            
            bool wasInNFBefore = false;
            for (int j = 0; j < 27; j++){
                if ( Particle->boxIndicesOfNearFieldBoxesAfterPush[i] == Particle->boxIndicesOfNearFieldBoxesBeforePush[j] )
                    wasInNFBefore = true;
            }
            if (wasInNFBefore == false){
                subLWFieldsInBox(Grid, Particle, Particle->boxIndicesOfNearFieldBoxesAfterPush[i], t);
            }
        }
    }
    
}

void calcInteractionWithOtherParticles(Particle *Particles, Grid *Grid, int numberOfParticles, int particleIndex, double E[3], double B[3]){
    int currentBoxIndexOfOtherParticle = -1;
    
    for (int p = 0; p < numberOfParticles; p++){
        if(p != particleIndex){
            currentBoxIndexOfOtherParticle = calcCurrentBoxIndexOfParticle(&Particles[p], Grid);
            if(boxIsInNearFieldOfParticle(Grid, &Particles[particleIndex], currentBoxIndexOfOtherParticle)){
                
                
                double beta[3] = {0};
                double intersectionPoint[4] = {0};
                double velocityAtIntersectionPoint[4] = {0};
                double gamma_sq;
                double R_sq;
                double R;
                double n[3] = {0};
                double betaDot[3] = {0};
                double dt = 0.5 * Grid->dx;
                double E_LW[3] = {0};
                double B_LW[3] = {0};
                
                int currentHistoryLength = Particles[p].currentHistoryLength;
                
                for (int index = 0; index < currentHistoryLength - 1; index ++){
                    if(isInsideBackwardLightcone(Particles[p].xHistory[index], Particles[particleIndex].x) && !isInsideBackwardLightcone(Particles[p].xHistory[index+1], Particles[particleIndex].x)){
                        printf("calcuating interaction with Particle%d ...\n", p);
                        calculateIntersectionPoint(Particles[p].xHistory[index], Particles[p].xHistory[index+1], Particles[p].uHistory[index], Particles[p].uHistory[index+1], Particles[particleIndex].x, intersectionPoint, velocityAtIntersectionPoint);
                        calculateBeta(Particles[p].xHistory[index], Particles[p].xHistory[index+1], beta);
                        calculateLienardWiechertParameters(intersectionPoint, Particles[particleIndex].x, velocityAtIntersectionPoint, &gamma_sq, &R_sq, &R, n, beta);
                        calculateBetaDot(Particles[p].uHistory[index], Particles[p].uHistory[index+1], dt, betaDot);
                        calcuateLienardWiechertFields(gamma_sq, R_sq, R, n, beta, betaDot, Particles[p].charge, E_LW, B_LW);
                        for (int i = 0; i < 3; i++){
                            E[i] += E_LW[i];
                            B[i] += B_LW[i];
                        }
                        break;
                    }
                }
            }
        }
    }
}


