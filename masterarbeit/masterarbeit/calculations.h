/* =============================================
HEADER FÜR EIGENE FUNKTIONEN
================================================ */

#ifndef CALCULATIONS
#define CALCULATIONS

void borisPusher(double *u, double *E, double *B, double dt, double chargeOverMass);
void crossProduct(double *a, double *b, double *result);
void push_u_boris(double *u, double charge_over_mass, double dt,
                  double *E,double *B);
#endif
