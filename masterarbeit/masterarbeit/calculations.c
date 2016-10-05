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

    double uPrime[3];
    double t[3];
    double s[3];
    double uCrossT[3];
    double uPrimeCrossS[3];
    double dtNew = dt / 2.0;
    double absoluteSquareValueOfT;
    
    u[1] += chargeOverMass * E[0] * dtNew;
    u[2] += chargeOverMass * E[1] * dtNew;
    u[3] += chargeOverMass * E[2] * dtNew;
    
    t[0] = chargeOverMass * B[0] * dtNew;
    t[1] = chargeOverMass * B[1] * dtNew;
    t[2] = chargeOverMass * B[2] * dtNew;
    
    crossProduct(u, t, uCrossT);
    
    uPrime[0] = u[1] + uCrossT[0];
    uPrime[1] = u[2] + uCrossT[1];
    uPrime[2] = u[3] + uCrossT[2];
    
    absoluteSquareValueOfT = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
    s[0] = 2 * t[0] / (1 + absoluteSquareValueOfT);
    s[1] = 2 * t[1] / (1 + absoluteSquareValueOfT);
    s[2] = 2 * t[2] / (1 + absoluteSquareValueOfT);
    
    crossProduct(uPrime, s, uPrimeCrossS);
    
    u[1] += uPrimeCrossS[0];
    u[2] += uPrimeCrossS[1];
    u[3] += uPrimeCrossS[2];
    
    u[1] += chargeOverMass * E[0] * dtNew;
    u[2] += chargeOverMass * E[1] * dtNew;
    u[3] += chargeOverMass * E[2] * dtNew;
    
}

void push_u_boris(double *u, double charge_over_mass, double dt,
                  double *E,double *B) {
    double us[3];
    double ucrossb[3];
    double f1, f2;
    
    // Definition p- in Arbeit.
    u[1] += (dt*charge_over_mass/2.0)*E[0];
    u[2] += (dt*charge_over_mass/2.0)*E[1];
    u[3] += (dt*charge_over_mass/2.0)*E[2];
    // entspricht gamma in der Arbeit. c wird wohl 1 sein. m kürzt sich mit dem m aus p
    u[0] = sqrt(1.0+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
    
    // u ist jetzt u-
    
    f1 = (dt*charge_over_mass)/(2.0*u[0]);
    f2 = (2.0*f1)/(1.0+f1*f1*(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]));
    
    crossProduct(u+1, B, ucrossb);
    us[0] = u[1]+f1*ucrossb[0];
    us[1] = u[2]+f1*ucrossb[1];
    us[2] = u[3]+f1*ucrossb[2];
    
    // us ist u_schlange und entspricht dem p Schlange in der Arbeit
    
    crossProduct(us, B, ucrossb);
    u[1] += f2*ucrossb[0];
    u[2] += f2*ucrossb[1];
    u[3] += f2*ucrossb[2];
    
    // u ist jetzt u+
    u[1] += (dt*charge_over_mass/2.0)*E[0];
    u[2] += (dt*charge_over_mass/2.0)*E[1];
    u[3] += (dt*charge_over_mass/2.0)*E[2];
    //	u[0] = sqrt(1.0+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
}








