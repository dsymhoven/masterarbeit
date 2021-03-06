/* =============================================
 IMPLEMENTATION OF CALCULATIONS HEADER
 ================================================ */

#import <math.h>
#import <stdbool.h>
#import <stdio.h>
#import <stdlib.h>
#include "calculations.h"
#include "string.h"
#include "lwFields.h"
#include "unistd.h"

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
 @param xParticle four  position vector of particle
 @param xObserver four  position vector of observation point
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
 @param xParticle four  position vector of particle
 @param xObserver four  position vector of observation point
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
 @brief returns the spatial vector product in three dimensions of two four vectors
 @param x first four vector
 @param y second four vector
 @returns x[1]*y[1] + x[2]*y[2] + x[3]*y[3]
 
 */
double vectorProduct(double x[4], double y[4]){
    
    return x[1]*y[1] + x[2]*y[2] + x[3]*y[3];
}


/**
 @brief scales a given vector x with factor
 @param x three vector to be scaled
 @param factor the factor to scale x with
 */
void scaleVector(double x[3], double factor){
    for (int i = 0; i < 3; i++){
        x[i] *= factor;
    }
}

/**
 @brief calculates the spatial distance between two four vectors.
 @param x first four vector
 @param y second four vector
 @returns the spatial distance between two four vectors.
 */
double calculateDistance(double x[4], double y[4]){
    double difference[4];
    vectorDifference(x, y, difference);
    return sqrt(difference[1]*difference[1] + difference[2]*difference[2] + difference[3]*difference[3]);
}

/**
 @brief calculates the minkowski product with metric (+,-,-,-)
 @param x first four vector
 @param y second four vector
 @returns minkowski product with metric (+,-,-,-)
 */
double minkowskiProduct(double x[4], double y[4]){
    return x[0] * y[0] - x[1] * y[1] - x[2] * y[2] - x[3] * y[3];
}

/**
@brief returns parameter lambda for the linear interpolation between particle trajectory and backwards lightcone at observation point.
 
 In order to calculate the Liénard-Wiechart fields at the observation point the intersection of the trajectory and the backward lightcone at the observation point is needed. Since the calculated trajectory is discrete we need to interpolate between the first point outside the lightcone and the last point inside the lightcone.
 @param xInside last four position vector of particle inside the lightcone
 @param xOutside first four position vector of particle outside the lightcone
 @param xObserver four  position vector of observation point
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
 @param xInside last four position vector of particle inside the lightcone
 @param xOutside first four position vector of particle outside the lightcone
 @param uInside last four velocity vector of particle inside the lightcone
 @param uOutside first four velocity vector of particle outside the lightcone
 @param xObserver four  position vector of observation point
 @param intersectionPoint four vector of intersection point
 @param velocityAtIntersectionPoint four velocity vector at intersection point
 */
void calculateIntersectionPoint(double xInside[4], double xOutside[4], double uInside[4], double uOutside[4], double xObserver[4], double intersectionPoint[4], double velocityAtIntersectionPoint[4]){
    double lambda = calculateLambdaForLinearInterpolation(xInside, xOutside, xObserver);
    for(int i = 0; i < 4; i++){
        intersectionPoint[i] = xOutside[i] + lambda * (xInside[i] - xOutside[i]);
        velocityAtIntersectionPoint[i] = uOutside[i] + lambda * (uInside[i] - uOutside[i]);
    }
    
}

/**
@brief updates velocity with boris method.
 
 First obtain “uMinus” by adding half acceleration to the initial velocity.
 Then obtain "uPlus" by performing rotation with "uPrime" and internally computed assisting values.
 Finally, add another half acceleration.
 
 @param Particles struct containing all particles
 @param Grid instance of Grid struct
 @param numberOfParticles number of Particles
 @param Eextern vector containing external E-Field components
 @param Bextern vector containing external B-Field components
 @param dt time increment
 */
void updateVelocityWithBorisPusher(Particle *Particles, Grid *Grid, int numberOfParticles, int particleIndex, double Eextern[3], double Bextern[3], double dt){
    
    Particle *Particle = &Particles[particleIndex];
    
    int dimension = 3;
    double uPrime[3];
    double t[3];
    double s[3];
    double uMinusCrossT[3];
    double uPrimeCrossS[3];
    double E[3]={0};
    double B[3]={0};
    double absoluteSquareValueOfT;
    double uModifiedForCrossProduct[3];
    double chargeOverMass = Particle->charge / Particle->mass;
        
    if (particleIsInsideSimulationArea(Grid, Particle)){
        updateFieldsForParticlePush(Particle, Grid, Eextern, Bextern, E, B);
        calcInteractionWithOtherParticles(Particles, Grid, numberOfParticles, particleIndex, E, B);
    }
    
    
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
///@param Particle instance of particle struct
///@param Grid instance of Grid struct
///@param dt time increment
void updateLocation(Particle *Particle, Grid *Grid, double dt){
    int dimension = 3;
    
    int oldBoxIndexOfParticle = calcCurrentBoxIndexOfParticle(Particle, Grid);
    calcBoxIndizesOfNextNeighbourBoxes(Grid, Particle, Particle->boxIndicesOfNearFieldBoxesBeforePush);
    
    for(int i = 0; i < dimension; i++){
        Particle->x[i+1] += (Particle->u[i+1] * dt) / Particle->u[0];
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
///@param Particles struct containing all particles
///@param Grid instance of Grid struct
///@param numberOfParticles number of Particles
//void updateLocationForParticles(Grid *Grid, Particle *Particles, int numberOfParticles,  double dt){
//    for(int p = 0; p < numberOfParticles; p++){
//        updateLocation(&Particles[p], Grid, dt);
//    }
//}


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
 @param uNew current velocity vector of particle
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
 @param xNew current position of particle
 @param beta acceleration vector of spatial components
 */
void calculateBeta(double xOld[4], double xNew[4], double beta[3]) {
    double dt = xNew[0]-xOld[0];
    beta[0] = (xNew[1]-xOld[1])/dt;
    beta[1] = (xNew[2]-xOld[2])/dt;
    beta[2] = (xNew[3]-xOld[3])/dt;
}

///@brief calculates LW fields from formula
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
///@param numberOfParticles number of particles
///@param t time at which the fields shall be calcualted. This is also the zeroth component of xObserver
///@param planeForPlotting number of the plane in which fields shall be calcualted. planeForPlotting = x[3] / dz
void calcLWFieldsForPlane(Grid *Grid, Particle *Particles, int numberOfParticles, double t, int planeForPlotting){
    for (int p = 0; p < numberOfParticles; p++){
        printf("Calculating LW Fields on plane %d, heigth:%f for Particle%d\n", planeForPlotting, planeForPlotting * Grid->Resolution.dz, p);
        
        double xObserver[4];
        int k = planeForPlotting;
        int nx = Grid->numberOfGridPointsInX;
        int ny = Grid->numberOfGridPointsInY;
        int nz = Grid->numberOfGridPointsInZ;
        
        for(int i = 0; i < nx; i++){
            for (int j = 0; j < ny; j++){
                
                
                xObserver[0] = t;
                xObserver[1] = (i)* Grid->Resolution.dx;
                xObserver[2] = (j) * Grid->Resolution.dy;
                xObserver[3] = (k) * Grid->Resolution.dz;
                
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

///@remark calcuates LW fields at simulation time t only on specyfied plane with index "planeForPlotting" under consideration of particles nearField.
///@param Grid pointer to an instance of a Grid struct.
///@param Particles pointer to  Particle struct array
///@param numberOfParticles number of particles
///@param t time at which the fields shall be calcualted. This is also the zeroth component of xObserver
///@param planeForPlotting number of the plane in which fields shall be calcualted. planeForPlotting = x[3] / dz
void calcLWFieldsForPlaneWithNearField(Grid *Grid, Particle *Particles, int numberOfParticles, double t, int planeForPlotting){
    
    for (int p = 0; p < numberOfParticles; p++){
        printf("Calculating LW Fields on plane %d, heigth:%f\n", planeForPlotting, planeForPlotting * Grid->Resolution.dz);
        Particle *Particle = &Particles[p];
        
        double xObserver[4];
        int k = planeForPlotting;
        
        int nx = Grid->numberOfGridPointsInX;
        int ny = Grid->numberOfGridPointsInY;
        int nz = Grid->numberOfGridPointsInZ;
        
        getEdgesOfNearFieldBox(Grid, Particle);
        
        for(int i = 0; i < nx; i++){
            for (int j = 0; j < ny; j++){
                
                
                xObserver[0] = t;
                xObserver[1] = (i) * Grid->Resolution.dx;
                xObserver[2] = (j) * Grid->Resolution.dy;
                xObserver[3] = (k) * Grid->Resolution.dz;
                
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
    
    
}

///@remark calcuates LW fields at simulation time t only on specyfied plane with index "planeForPlotting" and only within the near field region
///@param Grid pointer to an instance of a Grid struct.
///@param Particles pointer to  Particle struct array
///@param numberOfParticles number of particles
///@param t time at which the fields shall be calcualted. This is also the zeroth component of xObserver
///@param planeForPlotting number of the plane in which fields shall be calcualted. planeForPlotting = x[3] / dz
void calcLWFieldsForPlaneInNearFieldRegion(Grid *Grid, Particle *Particles, int numberOfParticles, double t, int planeForPlotting){
    
    for (int p = 0; p < numberOfParticles; p++){
        printf("Calculating LW Fields on plane %d, heigth:%f\n", planeForPlotting, planeForPlotting * Grid->Resolution.dz);
        Particle *Particle = &Particles[p];
        
        double xObserver[4];
        int k = planeForPlotting;
        
        int nx = Grid->numberOfGridPointsInX;
        int ny = Grid->numberOfGridPointsInY;
        int nz = Grid->numberOfGridPointsInZ;
        
        getEdgesOfNearFieldBox(Grid, Particle);
        
        for(int i = 0; i < nx; i++){
            for (int j = 0; j < ny; j++){
                
                
                xObserver[0] = t;
                xObserver[1] = (i) * Grid->Resolution.dx;
                xObserver[2] = (j) * Grid->Resolution.dy;
                xObserver[3] = (k) * Grid->Resolution.dz;
                
                if(xObserver[1] >= Particle->edgesOfNearFieldBox[0] && xObserver[1] < Particle->edgesOfNearFieldBox[1] && xObserver[2] >= Particle->edgesOfNearFieldBox[2] && xObserver[2] < Particle->edgesOfNearFieldBox[3] && xObserver[3] >= Particle->edgesOfNearFieldBox[4] && xObserver[3] < Particle->edgesOfNearFieldBox[5]){
                    int gridIndexInBox = 3 * ny * nz * i + 3 * nz * j + 3 * k;
                    
                    addLWField(Grid, Particle, &Grid->H[gridIndexInBox], xObserver, 3);
                    addLWField(Grid, Particle, &Grid->H[gridIndexInBox + 1], xObserver, 4);
                    addLWField(Grid, Particle, &Grid->H[gridIndexInBox + 2], xObserver, 5);
                    
                    addLWField(Grid, Particle, &Grid->E[gridIndexInBox], xObserver, 0);
                    addLWField(Grid, Particle, &Grid->E[gridIndexInBox + 1], xObserver, 1);
                    addLWField(Grid, Particle, &Grid->E[gridIndexInBox + 2], xObserver, 2);
                }
                else{
                    continue;
                }
            }
        }
    }
}


///@brief calcualtes the box index in which the particle is currently located.
///@param Grid pointer to an instance of a Grid struct.
///@param Particle pointer to an instance of a Particle struct
///@return index of box where particle is currenty located.
///@remark index for box is counted the same way as for E and B field array. First z component, then y and then x.
int calcCurrentBoxIndexOfParticle(Particle *Particle, Grid *Grid){
    int numberOfGridPointsForBoxInX = Grid->Box.numberOfGridPointsInX;
    int numberOfGridPointsForBoxInY = Grid->Box.numberOfGridPointsInY;
    int numberOfGridPointsForBoxInZ = Grid->Box.numberOfGridPointsInZ;
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    double dx = Grid->Resolution.dx;
    double dy = Grid->Resolution.dy;
    double dz = Grid->Resolution.dz;
    
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

///@brief calcualtes the lower left grid index of the box the particle is currently located in.
///@Grid instance of Grid struct
///@Particle instance of Particle struct
int calcLowerLeftGridIndexInBox(Grid *Grid, Particle *Particle){
    int numberOfGridPointsForBoxInX = Grid->Box.numberOfGridPointsInX;
    int numberOfGridPointsForBoxInY = Grid->Box.numberOfGridPointsInY;
    int numberOfGridPointsForBoxInZ = Grid->Box.numberOfGridPointsInZ;
    
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
    
    double m = 3.5;
    double sigmaMax = 0.8 * (m + 1) / Grid->Resolution.dx;
    double kappaMax = 1.0;
    
    int upmlLayerWidth = Grid->upmlLayerWidth;
    
    double sigma = 0.0;
    double kappa = 0.0;
    double depth = 0.0;
    
    double dt = Grid->Resolution.dx * 0.25;
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
        Grid->upml1Ex[j] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml2Ex[j] = (2 * dt) / (2 * kappa + sigma * dt);
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
        Grid->upml1Ey[k] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml2Ey[k] = (2 * dt) / (2 * kappa + sigma * dt);
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
        Grid->upml1Ez[i] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml2Ez[i] = (2 * dt) / (2 * kappa + sigma * dt);
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
        Grid->upml3Ex[k] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml4Ex[k] = 1.0 / (2 * kappa + sigma * dt);
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
        Grid->upml3Ey[i] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml4Ey[i] = 1.0 / (2 * kappa + sigma * dt);
    }
    
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
        Grid->upml3Ez[j] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml4Ez[j] = 1.0 / (2 * kappa + sigma * dt);
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
        Grid->upml5Ex[i] = (2 * kappa + sigma * dt);
        Grid->upml6Ex[i] = (2 * kappa - sigma * dt);
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
        Grid->upml5Ey[j] = (2 * kappa + sigma * dt);
        Grid->upml6Ey[j] = (2 * kappa - sigma * dt);
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
        Grid->upml5Ez[k] = (2 * kappa + sigma * dt);
        Grid->upml6Ez[k] = (2 * kappa - sigma * dt);
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
        Grid->upml1Hx[j] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml2Hx[j] = (2 * dt) / (2 * kappa + sigma * dt);
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
        Grid->upml1Hy[k] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml2Hy[k] = (2 * dt) / (2 * kappa + sigma * dt);
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
        Grid->upml1Hz[i] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml2Hz[i] = (2 * dt) / (2 * kappa + sigma * dt);
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
        Grid->upml3Hx[k] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml4Hx[k] = 1.0 / (2 * kappa + sigma * dt);
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
        Grid->upml3Hy[i] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml4Hy[i] = 1.0 / (2 * kappa + sigma * dt);
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
        Grid->upml3Hz[j] = (2 * kappa - sigma * dt) / (2 * kappa + sigma * dt);
        Grid->upml4Hz[j] = 1.0 / (2 * kappa + sigma * dt);
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
        Grid->upml5Hx[i] = (2 * kappa + sigma * dt);
        Grid->upml6Hx[i] = (2 * kappa - sigma * dt);
    }
    
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
        Grid->upml5Hy[j] = (2 * kappa + sigma * dt);
        Grid->upml6Hy[j] = (2 * kappa - sigma * dt);
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
        Grid->upml5Hz[k] = (2 * kappa + sigma * dt);
        Grid->upml6Hz[k] = (2 * kappa - sigma * dt);
    }
}

///@brief updates the near field region of a particle. If the particle moves into another box the near field region changes and therefore LW fields need to be added in those boxes which are not in the near field any longer and subtracted in those which were just added to the near field region.
///@param Grid instance of a Grid struct
///@param Particle instance of a Particle struct
///@param t current simulation time
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

///@brief calcualtes the field interactions in the near field region of a particle with other particles. If a particle enters the near field region of another particle the interaction due to LW fields need to be taken into account. Therfore calculte the LW fields analytically each time step for as long as the particle is inside the near field region of the other particle.
///@param Particles struct containing all particles
///@param Grid instance of a Grid struct
///@param numberOfParticles number of particles
///@param particleIndex outer loop index for current particle for which the interaction shall be calcualted
///@param E vector containing all E field components. At this point E contains not just the external field components but also the field contributions from other particles propagated on the grid inside the near field region of Particle with particleIndex
///@param B same as E
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
                double dt = 0.5 * Grid->Resolution.dx;
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

///@brief performs a reverse Boris Push. Uses same methods like standard BorisPusher but uses -E und -B, so that particle moves in the opposite direction. The particles velocites are also reversed.
///@param Particles struct containing all particles
///@param Grid instance of a Grid struct
///@param numberOfParticles number of particles
///@param Eextern vector containing external E-Field components
///@param Bextern vector containing external B-Field components
///@param dt time increment
///@param t current simulation time
void reverseBorisPush(Particle *Particles, Grid *Grid, int numberOfParticles, double Eextern[3], double Bextern[3], double dt, double t){
//    scaleVector(Eextern, -1.0);
    scaleVector(Bextern, -1.0);
    
    for (int p = 0; p < numberOfParticles; p++){
        for (int i = 0; i < 3; i++){
            Particles[p].u[i+1] *= -1;
        }
    }

    for (int step = 0; step < t / dt; step++){
        for(int p = 0; p < numberOfParticles; p++){
            addCurrentStateToParticleHistory(&Particles[p], step);
            updateVelocityWithBorisPusher(Particles, Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], Grid, dt);
        }
    }
    
}

///@brief After a reverse borisBush was performed the initial conditions have to be reset in order to start the actual simulation at the initial conditions the user has set in the beginning.
///@param Particles struct containing all particles
///@param numberOfParticles number of particles
///@param Eextern vector containing external E-Field components
///@param Bextern vector containing external B-Field components
void resetInitialConditions(Particle *Particles, int numberOfParticles, double Eextern[3], double Bextern[3]){
//    scaleVector(Eextern, -1.0);
    scaleVector(Bextern, -1.0);
    
    for (int p = 0; p < numberOfParticles; p++){
        for (int i = 0; i < 3; i++){
            Particles[p].u[i+1] *= -1;
        }
        Particles[p].x[0] = 0;
        Particles[p].currentHistoryLength = 0;
    }
}

///@brief reallocates the size of particle xHistory and uHistory array by t / dt. If simulation starts at t > 0 then the history array needs to be extended by t / dt entries.
///@param Particles struct containing all particles
///@param numberOfParticles number of particles
///@param dt time increment
///@param t current simulation time
void reallocateParticleHistory(Particle *Particles, int numberOfParticles, double dt, double t){
    for(int p = 0; p < numberOfParticles; p++){
        Particles[p].xHistory = (double**) realloc(Particles[p].xHistory, (Particles[p].lengthOfHistoryArray + t / dt) * sizeof(double*));
        if(Particles[p].xHistory == NULL){
            printf("ERROR: Could not reallocate memory for particle xHistory!\n");
            return;
        }
        for(int i = Particles[p].lengthOfHistoryArray; i < Particles[p].lengthOfHistoryArray + t / dt; i++){
            Particles[p].xHistory[i] = (double*) malloc(4 * sizeof(double));
            if(Particles[p].xHistory[i] == NULL){
                printf("ERROR: Could not allocate extra memory in particle xHistory!\n");
                return;
            }
        }
        
        Particles[p].uHistory = (double**) realloc(Particles[p].uHistory, (Particles[p].lengthOfHistoryArray + t / dt) * sizeof(double*));
        if(Particles[p].uHistory == NULL){
            printf("ERROR: Could not reallocate memory for particle uHistory!\n");
            return;
        }
        for(int i = Particles[p].lengthOfHistoryArray; i < Particles[p].lengthOfHistoryArray + t / dt; i++){
            Particles[p].uHistory[i] = (double*) malloc(4 * sizeof(double));
            if(Particles[p].uHistory[i] == NULL){
                printf("ERROR: Could not allocate extra memory in particle uHistory!\n");
                return;
            }
        }
        Particles[p].lengthOfHistoryArray += t / dt;

    }
}

///@brief extends and calculates particle history up to time t. Example: Let's say t = 8 and tEnd = 12 then the actual simulation starts at t = 8 but all LW fields have been calculated up to this time. Therefore particle histories got expanded in such a way that particle positions are know from t = 0 up to t = 8, where the position at t = 8 equals the inital position specified by the user in the beginning. This is realized by a reverse borisPush with inverted velocity and inverted external fields. Now the particle position is known at t = -8. At this point the current velocity is reverted again and zeroth component is set to 0. The particle can now be pushed to t = 8 again, but now all positions are known and LW fields can be calculated.
///@param Particles struct containing all particles
///@param Grid instance of a Grid struct
///@param numberOfParticles number of particles
///@param Eextern vector containing external E-Field components
///@param Bextern vector containing external B-Field components
///@param dt time increment
///@param t current simulation time
void extendParticleHistory(Particle *Particles, Grid *Grid, int numberOfParticles, double Eextern[3], double Bextern[3], double dt, double t){
    reallocateParticleHistory(Particles, numberOfParticles, dt, t);
    reverseBorisPush(Particles, Grid, numberOfParticles, Eextern, Bextern, dt, t);
    resetInitialConditions(Particles, numberOfParticles, Eextern, Bextern);
    
    for (int step = 0; step < t / dt; step++){
        for(int p = 0; p < numberOfParticles; p++){
            addCurrentStateToParticleHistory(&Particles[p], step);
            updateVelocityWithBorisPusher(Particles, Grid, numberOfParticles, p, Eextern, Bextern, dt);
            updateLocation(&Particles[p], Grid, dt);
        }
    }
}

///@brief calcualtes LW fields on complete grid by looping through all boxes. In the near field region of each particle no fields are calculated.
///@param Particles struct containing all particles
///@param Grid instance of a Grid struct
///@param numberOfParticles number of particles
///@param t current simulation time
void calcFieldsOnGridWithoutNearField(Particle *Particles, Grid *Grid, int numberOfParticles, double t){
  
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    for (int p = 0; p < numberOfParticles; p++){
        for (int boxIndex = 0; boxIndex < numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ; boxIndex++){
            if(!boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndex)){
                addLWFieldsInBox(Grid, &Particles[p], boxIndex, t);
            }
            
        }
    }
}

///@brief checks if two toubles are equal. Due to floating point arithmetic errors "==" is a very error prone way to go. If difference between two doubles is smaller then 10^-12 we consider them equal and true will be returned.
///@param a first double
///@param b second double
///@returns bool true if |a - b| < 10^-12, false otherwise
bool double_equals(double a, double b){
    return fabs(a - b) < pow(10,-12);
}

///@brief checks if initalFields have already been calcualted and reads them into Grid struct if so. Changes working directory to "~/Desktop/Projects/masterarbeit/Analysis/initialFields/" and reads in "numberOfDirectories.txt". Loops through all those directories with name "0", "1", "2", ... reads in "initialConditions.txt" into coresponding variables and checks if they are equal. Check "writeInitialConditionsToFile()" for parameter ordering. If at any point one parameter doesn't match then break and immediately return false.
bool readInitialFieldFromFileIfExists(Grid *Grid, Particle *Particles, int numberOfParticles, double t, double Eextern[3], double Bextern[3]){
    
    double dx;
    double dy;
    double dz;
    int numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ;
    int numberOfBoxesInX;
    int numberOfBoxesInY;
    int numberOfBoxesInZ;
    int numberOfParticlesTest;
    double tTest;
    double x0, u0;
    double x1, u1;
    double x2, u2;
    double x3, u3;
    double E0, E1, E2;
    double B0, B1, B2;
    bool doesExist = false;
    
    int numberOfDirectories;
    char file[32] = "some";
    char command[128] = "../../../../../../../../Desktop/Projects/masterarbeit/Analysis/initialFields/";
    FILE *fid;
    FILE *fid2;
    // here were are in /initialFields
    chdir(command);
    fid = fopen("numberOfDirectories.txt", "r");
    if(fid == NULL){
        return doesExist;
    }
    fscanf(fid, "%d\n", &numberOfDirectories);

    
    for (int i = 0; i <= numberOfDirectories; i++){
        sprintf(file, "%d", i);
        // here were are in /'file'
        chdir(file);
        fid = fopen("initialConditions.txt", "r");
        // here were are in /initialFields
        chdir("../");
        // first time without any initalFields there is no "initalConditions.txt" file yet.
        if (fid == NULL){
            continue;
        }
        fscanf(fid, "%lf %lf %lf %d %d %d %d %d %d %d %lf", &dx, &dy, &dz, &numberOfGridPointsForBoxInX, &numberOfGridPointsForBoxInY, &numberOfGridPointsForBoxInZ, &numberOfBoxesInX, &numberOfBoxesInY, &numberOfBoxesInZ, &numberOfParticlesTest, &tTest);
        
        if (double_equals(Grid->Resolution.dx, dx) && double_equals(Grid->Resolution.dy, dy) && double_equals(Grid->Resolution.dz, dz) && Grid->Box.numberOfGridPointsInX == numberOfGridPointsForBoxInX && Grid->Box.numberOfGridPointsInY == numberOfGridPointsForBoxInY && Grid->Box.numberOfGridPointsInZ == numberOfGridPointsForBoxInZ && Grid->numberOfBoxesInX == numberOfBoxesInX && Grid->numberOfBoxesInY == numberOfBoxesInY && Grid->numberOfBoxesInZ == numberOfBoxesInZ && numberOfParticles == numberOfParticlesTest && double_equals(t, tTest)){
            
            for(int p = 0; p < numberOfParticles; p++){
                fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf", &x0, &x1, &x2, &x3, &u0, &u1, &u2, &u3);
                if(double_equals(Particles[p].x[0], x0) && double_equals(Particles[p].x[1], x1) && double_equals(Particles[p].x[2], x2) && double_equals(Particles[p].x[3], x3) && double_equals(Particles[p].u[0], u0) && double_equals(Particles[p].u[1], u1) && double_equals(Particles[p].u[2], u2) && double_equals(Particles[p].u[3], u3)){
                    doesExist = true;
                }
                else{
                    doesExist = false;
                    continue;
                }
            }
            // only continue when all parameters for all particles match
            if(doesExist){
                fscanf(fid, "%lf %lf %lf %lf %lf %lf", &E0, &E1, &E2, &B0, &B1, &B2);
                if (double_equals(Eextern[0], E0) && double_equals(Eextern[1], E1) && double_equals(Eextern[2], E2) && double_equals(Bextern[0], B0) && double_equals(Bextern[1], B1) && double_equals(Bextern[2], B2)){
                    doesExist = true;
                    // here were are in /'file'
                    chdir(file);
                    fid = fopen("E_initialField.txt", "r");
                    if (fid == NULL){
                        printf("ERROR: Could not read from E_initialField.txt!\n");
                        break;
                    }
                    fid2 = fopen("H_initialField.txt", "r");
                    if (fid2 == NULL){
                        printf("ERROR: Could not read from H_initialField.txt!\n");
                        break;
                    }
                    // here were are in /initialFields
                    chdir("../");
                    printf("initial Field does already exist! Reading in ...\n");
                    for (int i = 0; i < Grid->numberOfGridPointsInX * Grid->numberOfGridPointsInY * Grid->numberOfGridPointsInZ * 3; i++){
                        fscanf(fid, "%lf", &Grid->E[i]);
                        fscanf(fid2, "%lf", &Grid->H[i]);
                    }
                    fclose(fid2);
                    // if an initial field was found don't continue in other directories.
                    break;
                }
                else{
                    doesExist = false;
                    break;
                }
            }
        }
        else{
            doesExist = false;
            continue;
        }

    }
    fclose(fid);
    // here were are back to executable directory
    chdir("../../../../../Library/Developer/Xcode/DerivedData/masterarbeit-gqaflsvzotvzahddhxaxmsryswvq/Build/Products/Debug/");
    return doesExist;
}


///@brief calcualtes the surrounding eigth grid points respective to the current particle location and saves them into gridIndizesNextNeighbours array.
///@param Grid Pointer to Grid Struct
///@param particleGridIndexInX number of grid point in x
///@param particleGridIndexInY number of grid point in y
///@param particleGridIndexInZ number of grid point in z
///@param gridIndizesNextNeighbours array containing the next neighbour indizes
void calcGridIndizesNextNeighbours(Grid *Grid, int const particleGridIndexInX, int const particleGridIndexInY, int const particleGridIndexInZ, int gridIndizesNextNeighbours[8]){
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    gridIndizesNextNeighbours[0] = particleGridIndexInX * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + particleGridIndexInY * numberOfGridPointsInZ * 3 + particleGridIndexInZ * 3;
    gridIndizesNextNeighbours[1] = particleGridIndexInX * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + (particleGridIndexInY + 1) * numberOfGridPointsInZ * 3 + particleGridIndexInZ * 3;
    gridIndizesNextNeighbours[2] = particleGridIndexInX * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + particleGridIndexInY * numberOfGridPointsInZ * 3 + (particleGridIndexInZ + 1) * 3;
    gridIndizesNextNeighbours[3] = particleGridIndexInX * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + (particleGridIndexInY + 1) * numberOfGridPointsInZ * 3 + (particleGridIndexInZ + 1) * 3;
    gridIndizesNextNeighbours[4] = (particleGridIndexInX + 1) * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + particleGridIndexInY * numberOfGridPointsInZ * 3 + particleGridIndexInZ * 3;
    gridIndizesNextNeighbours[5] = (particleGridIndexInX + 1) * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + (particleGridIndexInY + 1) * numberOfGridPointsInZ * 3 + particleGridIndexInZ * 3;
    gridIndizesNextNeighbours[6] = (particleGridIndexInX + 1) * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + particleGridIndexInY * numberOfGridPointsInZ * 3 + (particleGridIndexInZ + 1) * 3;
    gridIndizesNextNeighbours[7] = (particleGridIndexInX + 1) * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + (particleGridIndexInY + 1) * numberOfGridPointsInZ * 3 + (particleGridIndexInZ + 1) * 3;
}


///@brief performs trilinear interpolation of values A,B,C,D,E,F,G,H at interpolationPoint interpolation parameters u,v,w for each dimension respectively.
///@param interpolationPoint vector containing the position you wants to interpolate at
///@param A value at the lower left grid point in the front
///@param B value at the upper left grid point in the front
///@param C value at the lower left grid point in the back
///@param D value at the upper left grid point in the back
///@param E value at the lower right grid point in the front
///@param F value at the upper right grid point in the front
///@param G value at the lower right grid point in the back
///@param H value at the upper right grid point in the back
///@param u interpolation paramter for x in [0,1].
///@param v interpolation paramter for y in [0,1].
///@param w interpolation paramter for z in [0,1].
///@code u = (x(Particle) - x(A))/(x(A) - x(E));
/// v = (y(Particle) - y(A))/(y(A) - y(C));
/// w = (z(Particle) - z(A))/(z(A) - z(B));
double trilinearInterpolation(double interpolationPoint[3], double A, double B, double C, double D, double E, double F, double G, double H, double u, double v, double w){
    double P1, P2, P3, P4, Pu, Po;

    P1 = A + u * ( B - A );
    P2 = C + u * ( D - C );
    Pu = P1 + v * ( P2 - P1 );
    
    P3 = E + u * ( F - E );
    P4 = G + u * ( H - G );
    Po = P3 + v * ( P4 - P3 );
    
    return Pu + w * ( Po - Pu);
    
}

///@brief the particle is usually not located exactly at a grid point, but inbetween eight sourrounding grid points. Those eight sourrounding grid points are defined as box here. This box is divided in eight partial boxes where the index is built similiar to the big grid boxes. 0 is the lower left corner in the front and 7 is the upper right corner in the back where we first increase z, then y and finally x. The seperation is necessary for the staggered grid to apply properly, because the location of the respective field components which we want to interpolate from differ for all components.
///@param Particle Pointer to Particle Struct
///@param Grid Pointer to Grid Struct
///returns partialBoxIndex index from 0 to 7 where 0 denotes the front lower left part of the box and 7 the back upper right part.
int getPartialBoxIndex(Particle *Particle, Grid *Grid){
    int ip, jp, kp;
    int partialBoxIndex = -1;
    
    ip = Particle->x[1] / Grid->Resolution.dx;
    jp = Particle->x[2] / Grid->Resolution.dy;
    kp = Particle->x[3] / Grid->Resolution.dz;
    
    if(fmod(Particle->x[1], Grid->Resolution.dx) <= 0.5 * Grid->Resolution.dx && fmod(Particle->x[2], Grid->Resolution.dy) <= 0.5 * Grid->Resolution.dy && fmod(Particle->x[3], Grid->Resolution.dz) <= 0.5 * Grid->Resolution.dz){
        partialBoxIndex = 0;
    }
    else if(fmod(Particle->x[1], Grid->Resolution.dx) <= 0.5 * Grid->Resolution.dx && fmod(Particle->x[2], Grid->Resolution.dy) <= 0.5 * Grid->Resolution.dy && fmod(Particle->x[3], Grid->Resolution.dz) > 0.5 * Grid->Resolution.dz){
        partialBoxIndex = 1;
    }
    else if(fmod(Particle->x[1], Grid->Resolution.dx) <= 0.5 * Grid->Resolution.dx && fmod(Particle->x[2], Grid->Resolution.dy) > 0.5 * Grid->Resolution.dy && fmod(Particle->x[3], Grid->Resolution.dz) <= 0.5 * Grid->Resolution.dz){
        partialBoxIndex = 2;
    }
    else if(fmod(Particle->x[1], Grid->Resolution.dx) <= 0.5 * Grid->Resolution.dx && fmod(Particle->x[2], Grid->Resolution.dy) > 0.5 * Grid->Resolution.dy && fmod(Particle->x[3], Grid->Resolution.dz) > 0.5 * Grid->Resolution.dz){
        partialBoxIndex = 3;
    }

    else if(fmod(Particle->x[1], Grid->Resolution.dx) > 0.5 * Grid->Resolution.dx && fmod(Particle->x[2], Grid->Resolution.dy) <= 0.5 * Grid->Resolution.dy && fmod(Particle->x[3], Grid->Resolution.dz) <= 0.5 * Grid->Resolution.dz){
        partialBoxIndex = 4;
    }
    else if(fmod(Particle->x[1], Grid->Resolution.dx) > 0.5 * Grid->Resolution.dx && fmod(Particle->x[2], Grid->Resolution.dy) <= 0.5 * Grid->Resolution.dy && fmod(Particle->x[3], Grid->Resolution.dz) > 0.5 * Grid->Resolution.dz){
        partialBoxIndex = 5;
    }
    else if(fmod(Particle->x[1], Grid->Resolution.dx) > 0.5 * Grid->Resolution.dx && fmod(Particle->x[2], Grid->Resolution.dy) > 0.5 * Grid->Resolution.dy && fmod(Particle->x[3], Grid->Resolution.dz) <= 0.5 * Grid->Resolution.dz){
        partialBoxIndex = 6;
    }
    else if(fmod(Particle->x[1], Grid->Resolution.dx) > 0.5 * Grid->Resolution.dx && fmod(Particle->x[2], Grid->Resolution.dy) > 0.5 * Grid->Resolution.dy && fmod(Particle->x[3], Grid->Resolution.dz) > 0.5 * Grid->Resolution.dz){
        partialBoxIndex = 7;
    }

    return partialBoxIndex;
    
}

///@brief trilinear interpolation of LW fields at particle position. In order to calculate the interaction of the particle with the far fields of other particles properly, we need to interpolate the field values. Field values only exist at the grid points, whereas the particle can be at any point on the grid. This method uses trilinear interpolation with the field values on the sorrounding eight grid points. We also consider the staggered grid. See Taflove page 59 Fig 3.1 for details and consider that we use another orientation of the coordinate system. Due to the staggered grid, the eight grid points, forming the box to interpolate in, differ for all six components Ex, Ey, Ez, Hx, Hy, Hz.
///@param Grid Pointer to Grid Struct
///@param Particle Pointer to Particle Struct
///@param E vector where interpolated E fields at particle position are stored
///@param B vector where interpolated H fields at particle position are stored
void interpolateFields(Grid *Grid, Particle *Particle, double E[3], double B[3]){
    int ip, jp, kp;
    double u,v,w;
    double x0=0.0,x1=0.0,y0=0.0,y1=0.0,z0=0.0,z1=0.0;
    int gridIndizesNextNeighbours[8];
    double interpolationPoint[3];
    int partialBoxIndex;
    
    ip = Particle->x[1] / Grid->Resolution.dx;
    jp = Particle->x[2] / Grid->Resolution.dy;
    kp = Particle->x[3] / Grid->Resolution.dz;
    
    
    if (ip < Grid->numberOfGridPointsInX - 1 && jp < Grid->numberOfGridPointsInY - 1 && kp < Grid->numberOfGridPointsInZ - 1 && ip > 0 && jp > 0 && kp > 0){
        
        interpolationPoint[0] = Particle->x[1];
        interpolationPoint[1] = Particle->x[2];
        interpolationPoint[2] = Particle->x[3];

        partialBoxIndex = getPartialBoxIndex(Particle, Grid);

        // MARK: Hx
        switch (partialBoxIndex) {
            case 0:
                calcGridIndizesNextNeighbours(Grid, ip, jp - 1, kp - 1, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 0.5) * Grid->Resolution.dy;
                
                z0 = (kp - 0.5) * Grid->Resolution.dz;
                z1 = (kp + 0.5) * Grid->Resolution.dz;
                break;
            case 1:
                calcGridIndizesNextNeighbours(Grid, ip, jp - 1, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 0.5) * Grid->Resolution.dy;
                
                z0 = (kp + 0.5) * Grid->Resolution.dz;
                z1 = (kp + 1.5) * Grid->Resolution.dz;
                break;
            case 2:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp - 1, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp + 0.5) * Grid->Resolution.dy;
                y1 = (jp + 1.5) * Grid->Resolution.dy;
                
                z0 = (kp - 0.5) * Grid->Resolution.dz;
                z1 = (kp + 0.5) * Grid->Resolution.dz;
                break;
            case 3:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp + 0.5) * Grid->Resolution.dy;
                y1 = (jp + 1.5) * Grid->Resolution.dy;
                
                z0 = (kp + 0.5) * Grid->Resolution.dz;
                z1 = (kp + 1.5) * Grid->Resolution.dz;
                break;
            case 4:
                calcGridIndizesNextNeighbours(Grid, ip, jp - 1, kp - 1, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 0.5) * Grid->Resolution.dy;
                
                z0 = (kp - 0.5) * Grid->Resolution.dz;
                z1 = (kp + 0.5) * Grid->Resolution.dz;
                break;
            case 5:
                calcGridIndizesNextNeighbours(Grid, ip, jp - 1, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 0.5) * Grid->Resolution.dy;
                
                z0 = (kp + 0.5) * Grid->Resolution.dz;
                z1 = (kp + 1.5) * Grid->Resolution.dz;
                break;
            case 6:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp - 1, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp + 0.5) * Grid->Resolution.dy;
                y1 = (jp + 1.5) * Grid->Resolution.dy;
                
                z0 = (kp - 0.5) * Grid->Resolution.dz;
                z1 = (kp + 0.5) * Grid->Resolution.dz;
                break;
            case 7:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp + 0.5) * Grid->Resolution.dy;
                y1 = (jp + 1.5) * Grid->Resolution.dy;
                
                z0 = (kp + 0.5) * Grid->Resolution.dz;
                z1 = (kp + 1.5) * Grid->Resolution.dz;
                break;
            default:
                break;
        }

        u = (Particle->x[1] - x0)/(x1 - x0);
        v = (Particle->x[2] - y0)/(y1 - y0);
        w = (Particle->x[3] - z0)/(z1 - z0);
        
        
        B[0] = trilinearInterpolation(interpolationPoint, Grid->H[gridIndizesNextNeighbours[0]], Grid->H[gridIndizesNextNeighbours[1]], Grid->H[gridIndizesNextNeighbours[2]], Grid->H[gridIndizesNextNeighbours[3]], Grid->H[gridIndizesNextNeighbours[4]], Grid->H[gridIndizesNextNeighbours[5]], Grid->H[gridIndizesNextNeighbours[6]], Grid->H[gridIndizesNextNeighbours[7]], u, v, w);
    
        // MARK: Hy
        switch (partialBoxIndex) {
            case 0:
                calcGridIndizesNextNeighbours(Grid, ip - 1, jp, kp - 1, gridIndizesNextNeighbours);
                x0 = (ip - 0.5) * Grid->Resolution.dx;
                x1 = (ip + 0.5) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp - 0.5) * Grid->Resolution.dz;
                z1 = (kp + 0.5) * Grid->Resolution.dz;
                break;
            case 1:
                calcGridIndizesNextNeighbours(Grid, ip - 1, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip - 0.5) * Grid->Resolution.dx;
                x1 = (ip + 0.5) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp + 0.5) * Grid->Resolution.dz;
                z1 = (kp + 1.5) * Grid->Resolution.dz;
                break;
            case 2:
                calcGridIndizesNextNeighbours(Grid, ip - 1, jp, kp - 1, gridIndizesNextNeighbours);
                x0 = (ip - 0.5) * Grid->Resolution.dx;
                x1 = (ip + 0.5) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp - 0.5) * Grid->Resolution.dz;
                z1 = (kp + 0.5) * Grid->Resolution.dz;
                break;
            case 3:
                calcGridIndizesNextNeighbours(Grid, ip - 1, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip - 0.5) * Grid->Resolution.dx;
                x1 = (ip + 0.5) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp + 0.5) * Grid->Resolution.dz;
                z1 = (kp + 1.5) * Grid->Resolution.dz;
                break;
            case 4:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp - 1, gridIndizesNextNeighbours);
                x0 = (ip + 0.5) * Grid->Resolution.dx;
                x1 = (ip + 1.5) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp - 0.5) * Grid->Resolution.dz;
                z1 = (kp + 0.5) * Grid->Resolution.dz;
                break;
            case 5:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip + 0.5) * Grid->Resolution.dx;
                x1 = (ip + 1.5) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp + 0.5) * Grid->Resolution.dz;
                z1 = (kp + 1.5) * Grid->Resolution.dz;
                break;
            case 6:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp - 1, gridIndizesNextNeighbours);
                x0 = (ip + 0.5) * Grid->Resolution.dx;
                x1 = (ip + 1.5) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp - 0.5) * Grid->Resolution.dz;
                z1 = (kp + 0.5) * Grid->Resolution.dz;
                break;
            case 7:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip + 0.5) * Grid->Resolution.dx;
                x1 = (ip + 1.5) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp + 0.5) * Grid->Resolution.dz;
                z1 = (kp + 1.5) * Grid->Resolution.dz;
                break;
            default:
                break;
        }
        

        u = (Particle->x[1] - x0)/(x1 - x0);
        v = (Particle->x[2] - y0)/(y1 - y0);
        w = (Particle->x[3] - z0)/(z1 - z0);
        
        
        B[1] = trilinearInterpolation(interpolationPoint, Grid->H[gridIndizesNextNeighbours[0] + 1], Grid->H[gridIndizesNextNeighbours[1] + 1], Grid->H[gridIndizesNextNeighbours[2] + 1], Grid->H[gridIndizesNextNeighbours[3] + 1], Grid->H[gridIndizesNextNeighbours[4] + 1], Grid->H[gridIndizesNextNeighbours[5] + 1], Grid->H[gridIndizesNextNeighbours[6] + 1], Grid->H[gridIndizesNextNeighbours[7] + 1], u, v, w);
        
        // MARK: Hz
        switch (partialBoxIndex) {
            case 0:
                calcGridIndizesNextNeighbours(Grid, ip - 1, jp - 1, kp, gridIndizesNextNeighbours);
                x0 = (ip - 0.5) * Grid->Resolution.dx;
                x1 = (ip + 0.5) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 0.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 1:
                calcGridIndizesNextNeighbours(Grid, ip - 1, jp - 1, kp, gridIndizesNextNeighbours);
                x0 = (ip - 0.5) * Grid->Resolution.dx;
                x1 = (ip + 0.5) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 0.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 2:
                calcGridIndizesNextNeighbours(Grid, ip - 1, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip - 0.5) * Grid->Resolution.dx;
                x1 = (ip + 0.5) * Grid->Resolution.dx;
                
                y0 = (jp + 0.5) * Grid->Resolution.dy;
                y1 = (jp + 1.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 3:
                calcGridIndizesNextNeighbours(Grid, ip - 1, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip - 0.5) * Grid->Resolution.dx;
                x1 = (ip + 0.5) * Grid->Resolution.dx;
                
                y0 = (jp + 0.5) * Grid->Resolution.dy;
                y1 = (jp + 1.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 4:
                calcGridIndizesNextNeighbours(Grid, ip, jp - 1, kp, gridIndizesNextNeighbours);
                x0 = (ip + 0.5) * Grid->Resolution.dx;
                x1 = (ip + 1.5) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 0.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 5:
                calcGridIndizesNextNeighbours(Grid, ip, jp - 1, kp, gridIndizesNextNeighbours);
                x0 = (ip + 0.5) * Grid->Resolution.dx;
                x1 = (ip + 1.5) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 0.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 6:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip + 0.5) * Grid->Resolution.dx;
                x1 = (ip + 1.5) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 1.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 7:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip + 0.5) * Grid->Resolution.dx;
                x1 = (ip + 1.5) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 1.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            default:
                break;
        }
        
        
        u = (Particle->x[1] - x0)/(x1 - x0);
        v = (Particle->x[2] - y0)/(y1 - y0);
        w = (Particle->x[3] - z0)/(z1 - z0);
        
        
        B[2] = trilinearInterpolation(interpolationPoint, Grid->H[gridIndizesNextNeighbours[0] + 2], Grid->H[gridIndizesNextNeighbours[1] + 2], Grid->H[gridIndizesNextNeighbours[2] + 2], Grid->H[gridIndizesNextNeighbours[3] + 2], Grid->H[gridIndizesNextNeighbours[4] + 2], Grid->H[gridIndizesNextNeighbours[5] + 2], Grid->H[gridIndizesNextNeighbours[6] + 2], Grid->H[gridIndizesNextNeighbours[7] + 2], u, v, w);

    
        // MARK: Ex
        switch (partialBoxIndex) {
            case 0:
                calcGridIndizesNextNeighbours(Grid, ip - 1, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip - 0.5) * Grid->Resolution.dx;
                x1 = (ip + 0.5) * Grid->Resolution.dx;
                
                y0 = jp  * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 1:
                calcGridIndizesNextNeighbours(Grid, ip - 1, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip - 0.5) * Grid->Resolution.dx;
                x1 = (ip + 0.5) * Grid->Resolution.dx;
                
                y0 = jp  * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 2:
                calcGridIndizesNextNeighbours(Grid, ip - 1, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip - 0.5) * Grid->Resolution.dx;
                x1 = (ip + 0.5) * Grid->Resolution.dx;
                
                y0 = jp  * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 3:
                calcGridIndizesNextNeighbours(Grid, ip - 1, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip - 0.5) * Grid->Resolution.dx;
                x1 = (ip + 0.5) * Grid->Resolution.dx;
                
                y0 = jp  * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 4:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip + 0.5) * Grid->Resolution.dx;
                x1 = (ip + 1.5) * Grid->Resolution.dx;
                
                y0 = jp  * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 5:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip + 0.5) * Grid->Resolution.dx;
                x1 = (ip + 1.5) * Grid->Resolution.dx;
                
                y0 = jp  * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 6:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip + 0.5) * Grid->Resolution.dx;
                x1 = (ip + 1.5) * Grid->Resolution.dx;
                
                y0 = jp  * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 7:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = (ip + 0.5) * Grid->Resolution.dx;
                x1 = (ip + 1.5) * Grid->Resolution.dx;
                
                y0 = jp  * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            default:
                break;
        }
        
        
        u = (Particle->x[1] - x0)/(x1 - x0);
        v = (Particle->x[2] - y0)/(y1 - y0);
        w = (Particle->x[3] - z0)/(z1 - z0);
        
        
        E[0] = trilinearInterpolation(interpolationPoint, Grid->E[gridIndizesNextNeighbours[0]], Grid->E[gridIndizesNextNeighbours[1]], Grid->E[gridIndizesNextNeighbours[2]], Grid->E[gridIndizesNextNeighbours[3]], Grid->E[gridIndizesNextNeighbours[4]], Grid->E[gridIndizesNextNeighbours[5]], Grid->E[gridIndizesNextNeighbours[6]], Grid->E[gridIndizesNextNeighbours[7]], u, v, w);
        
        // MARK: Ey
        switch (partialBoxIndex) {
            case 0:
                calcGridIndizesNextNeighbours(Grid, ip, jp - 1, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 0.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 1:
                calcGridIndizesNextNeighbours(Grid, ip, jp - 1, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 0.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 2:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp + 0.5) * Grid->Resolution.dy;
                y1 = (jp + 1.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 3:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp + 0.5) * Grid->Resolution.dy;
                y1 = (jp + 1.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 4:
                calcGridIndizesNextNeighbours(Grid, ip, jp - 1, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 0.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 5:
                calcGridIndizesNextNeighbours(Grid, ip, jp - 1, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp - 0.5) * Grid->Resolution.dy;
                y1 = (jp + 0.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 6:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp + 0.5) * Grid->Resolution.dy;
                y1 = (jp + 1.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            case 7:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = (jp + 0.5) * Grid->Resolution.dy;
                y1 = (jp + 1.5) * Grid->Resolution.dy;
                
                z0 = kp * Grid->Resolution.dz;
                z1 = (kp + 1) * Grid->Resolution.dz;
                break;
            default:
                break;
        }
        
        
        u = (Particle->x[1] - x0)/(x1 - x0);
        v = (Particle->x[2] - y0)/(y1 - y0);
        w = (Particle->x[3] - z0)/(z1 - z0);
        
        
        E[1] = trilinearInterpolation(interpolationPoint, Grid->E[gridIndizesNextNeighbours[0] + 1], Grid->E[gridIndizesNextNeighbours[1] + 1], Grid->E[gridIndizesNextNeighbours[2] + 1], Grid->E[gridIndizesNextNeighbours[3] + 1], Grid->E[gridIndizesNextNeighbours[4] + 1], Grid->E[gridIndizesNextNeighbours[5] + 1], Grid->E[gridIndizesNextNeighbours[6] + 1], Grid->E[gridIndizesNextNeighbours[7] + 1], u, v, w);
        
        // MARK: Ez
        switch (partialBoxIndex) {
            case 0:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp - 1, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp - 0.5) * Grid->Resolution.dz;
                z1 = (kp + 0.5) * Grid->Resolution.dz;
                break;
            case 1:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp + 0.5) * Grid->Resolution.dz;
                z1 = (kp + 1.5) * Grid->Resolution.dz;
                break;
            case 2:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp - 1, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp - 0.5) * Grid->Resolution.dz;
                z1 = (kp + 0.5) * Grid->Resolution.dz;
                break;
            case 3:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp + 0.5) * Grid->Resolution.dz;
                z1 = (kp + 1.5) * Grid->Resolution.dz;
                break;
            case 4:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp - 1, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp - 0.5) * Grid->Resolution.dz;
                z1 = (kp + 0.5) * Grid->Resolution.dz;
                break;
            case 5:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp + 0.5) * Grid->Resolution.dz;
                z1 = (kp + 1.5) * Grid->Resolution.dz;
                break;
            case 6:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp - 1, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp - 0.5) * Grid->Resolution.dz;
                z1 = (kp + 0.5) * Grid->Resolution.dz;
                break;
            case 7:
                calcGridIndizesNextNeighbours(Grid, ip, jp, kp, gridIndizesNextNeighbours);
                x0 = ip * Grid->Resolution.dx;
                x1 = (ip + 1) * Grid->Resolution.dx;
                
                y0 = jp * Grid->Resolution.dy;
                y1 = (jp + 1) * Grid->Resolution.dy;
                
                z0 = (kp + 0.5) * Grid->Resolution.dz;
                z1 = (kp + 1.5) * Grid->Resolution.dz;
                break;
            default:
                break;
        }
        
        
        u = (Particle->x[1] - x0)/(x1 - x0);
        v = (Particle->x[2] - y0)/(y1 - y0);
        w = (Particle->x[3] - z0)/(z1 - z0);
        
        
        E[2] = trilinearInterpolation(interpolationPoint, Grid->E[gridIndizesNextNeighbours[0] + 2], Grid->E[gridIndizesNextNeighbours[1] + 2], Grid->E[gridIndizesNextNeighbours[2] + 2], Grid->E[gridIndizesNextNeighbours[3] + 2], Grid->E[gridIndizesNextNeighbours[4] + 2], Grid->E[gridIndizesNextNeighbours[5] + 2], Grid->E[gridIndizesNextNeighbours[6] + 2], Grid->E[gridIndizesNextNeighbours[7] + 2], u, v, w);
    }
    
    


}
///@brief simple tensor - vector contraction. In this case it's similiar to a matrix - vector multiplication. Each tensor element of the first row is multiplied with the respective component of the vector you are contracting with.
///@param tensor 4x4 tensor you want to contract
///@param vector four vector you are contracing with
///@param result four vector containing the result of the contraction
void contractTensorWithFourVector(double tensor[4][4], double vector[4], double result[4]){
    
    for (int i = 0; i < 4; i++){
        result[i] = minkowskiProduct(tensor[i], vector);
    }
    
}

///@brief sets the field values of field tensor. For simplicity reasons we use the mixed notation (one index up, one down)
///@param tensor 4x4 tensor where the field values are getting stored in
///@param Eextern vector containing components of external electric field
///@param Bextern vector containing components of external magnetic field
void updateFieldTensor(double tensor[4][4], double Eextern[3], double Bextern[3]){
    tensor[0][0] = 0;
    tensor[0][1] = Eextern[0];
    tensor[0][2] = Eextern[1];
    tensor[0][3] = Eextern[2];
    tensor[1][0] = Eextern[1];
    tensor[1][1] = 0;
    tensor[1][2] = Bextern[2];
    tensor[1][3] = -Bextern[1];
    tensor[2][0] = Eextern[1];
    tensor[2][1] = -Bextern[2];
    tensor[2][2] = 0;
    tensor[2][3] = Bextern[0];
    tensor[3][0] = Eextern[2];
    tensor[3][1] = Bextern[1];
    tensor[3][2] = -Bextern[0];
    tensor[3][3] = 0;
}

///@brief calculates radiation damping with Landau-Lifschitz-Equation. See equation D.15 in dissertation from Christian Herzing. For simplicity we use the mixed notation for the field tensor (one Index up, one down). First term of equation D.15 is neglectible, as was proven in "Radiation Reaction Effects on Radiation Pressure Acceleration".
///@param Eextern vector containing components of external electric field
///@param Bextern vector containing components of external magnetic field
///@param u velocity vector
///@param dampingTerm four vector containing the result of the calculation
void calcRadiationDamping(double Eextern[3], double Bextern[3], double u[4], double dampingTerm[3]){
    double fieldTensor[4][4] = {0};
    double Fu[4] = {0};
    double FFu[4] = {0};
    
    updateFieldTensor(fieldTensor, Eextern, Bextern);
    contractTensorWithFourVector(fieldTensor, u, Fu);
    contractTensorWithFourVector(fieldTensor, Fu, FFu);
    
    for (int i = 0; i < 3; i++){
        dampingTerm[i] = 2. / 3. * (FFu[i+1] + u[i+1] * minkowskiProduct(Fu, Fu));
    }
    
}

void calculateLorentzForce(Particle *Particle, double Bextern[3], double lorentzForce[3]){
    double chargeOverMass = Particle->charge / Particle->mass;
    double u[3] = {0};
    double uCrossB[3] = {0};
    
    for(int i = 0; i < 3; i++){
        u[i] = Particle->u[i+1];
    }
    
    crossProduct(u, Bextern, uCrossB);
    
    for(int i = 0; i < 3; i++){
        lorentzForce[i] = uCrossB[i] * chargeOverMass;
    }
    
}





