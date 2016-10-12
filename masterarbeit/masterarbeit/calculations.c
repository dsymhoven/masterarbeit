/* =============================================
IMPLEMENTATION OF CALCULATIONS HEADER
================================================ */

#import <math.h>
#import <stdbool.h>

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
bool isInsideBackwardLightcone(double *xParticle, double *xObserver) {
    double dt;
    double dxsq;
    
    dt = xObserver[0]-xParticle[0];
    dxsq = (xObserver[1]-xParticle[1])*(xObserver[1]-xParticle[1])+(xObserver[2]-xParticle[2])*(xObserver[2]-xParticle[2])
    +(xObserver[3]-xParticle[3])*(xObserver[3]-xParticle[3]);
    
    return (dt>0&&dt*dt>dxsq);
}

/**
 returns true if position of particle is inside the forward lightcone of the observer, otherwise false.
 
 - remark: two conditions have to be fulfilled. Since the information needs time to travel from particle position to observer position, time difference has to be negative: xObserver[0]-xParticle[0] < 0.
 And Since c > v <=> (ct)^2 > x^2
 */
bool isInsideForwardLightcone(double *xParticle, double *xObserver) {
    
    double dt;
    double dxsq;
    
    dt = xObserver[0]-xParticle[0];
    dxsq = (xObserver[1]-xParticle[1])*(xObserver[1]-xParticle[1])+(xObserver[2]-xParticle[2])*(xObserver[2]-xParticle[2])
    +(xObserver[3]-xParticle[3])*(xObserver[3]-xParticle[3]);
    
    return (dt<0&&dt*dt>dxsq);
}


/**
 calcualtes the difference between two three dimensional vectors and saves the result in the array passed as argument.
 */
void vectorDifference(double *x, double *y, double *result){
    for (int i = 0; i < 3; i++){
        result[i] = x[i] - y[i];
    }
}


/**
 returns the vector product in three dimensions.
 */
double vectorProduct(double *x, double *y){
    
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

/**
 returns parameter lambda for the linear interpolation between particle trajectory and backwards lightcone at observation point.
 
 In order to calculate the Liénard-Wiechart fields at the observation point the intersection of the trajectory and the backward lightcone at the observation point is needed. Since the calculated trajectory is discrete we need to interpolate between the first point outside the lightcone and the last point inside the lightcone.
 */
double calculateLambdaForLinearInterpolation(double *xInside, double *xOutside, double *xObserver) {
    
    double a, b, c;
    double lambda;
    double xInsideMinusxOutside[3];
    double xObserverMinusxOutside[3];
    
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
    
}

void updateLocation(double *u, double *x, double dt){
    int dimension = 3;
    
    for(int i = 0; i < dimension; i++){
        x[i+1] += u[i+1] * dt;
    }
    
}







