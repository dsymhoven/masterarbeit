/* =============================================
IMPLEMENTATION OF CALCULATIONS HEADER
================================================ */

#import <math.h>

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

void borisPusher(double *u, double *E, double *B, double dt, double chargeOverMass){

    int dimension = 3;
    double uPrime[3];
    double t[3];
    double s[3];
    double uMinusCrossT[3];
    double uMinus[3];
    double uPlus[3];
    double uPrimeCrossS[3];
    double absoluteSquareValueOfT;
    
    for (int i = 0; i < dimension; i++){
        uMinus[i] = u[i] + chargeOverMass * E[i] * dt * 0.5;
    }
    
    for (int i = 0; i < dimension; i++){
        t[i] = chargeOverMass * B[i] * dt * 0.5;
    }
    
    crossProduct(uMinus, t, uMinusCrossT);
    
    for (int i = 0; i < dimension; i++){
        uPrime[i] = uMinus[i] + uMinusCrossT[i];
    }
    
    absoluteSquareValueOfT = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
    
    for(int i = 0; i < dimension; i++){
        s[i] = 2 * t[i] / (1 + absoluteSquareValueOfT);
    }
    
    
    crossProduct(uPrime, s, uPrimeCrossS);
    
    for (int i = 0; i < dimension; i++){
        uPlus[i] = uMinus[i] + uPrimeCrossS[i];
    }
    
    for (int i = 0; i < dimension; i++){
        u[i] = uPlus[i] + chargeOverMass * E[i] * dt * 0.5;
    }
    
}







