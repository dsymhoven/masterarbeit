
// Coulomb Interaction of electron with proton

// x'' = -1/4 M_pi e_0 x / (x^2 + y^2)^(3/2)
// y'' = -1/4 M_pi e_0 y / (x^2 + y^2)^(3/2)

// x(0) = 0.5, y(0) = 0, vx(0) = 0, vy(0) = 1.0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, const char * argv[]) {
    
    // Parameter definitions and allocations
    double h = 0.01;
    double kx1,kx2,kx3,kx4,ky1,ky2,ky3,ky4;
    double t = 0;
    double GM = 1.0;
    double e = 1;
    double tEnd = 10000;
    double x;
    double vx;
    double y;
    double vy;
    double Ekin,Epot;

    
    // open file identifier
    FILE *fid = fopen("solution.txt","w");
    
    // boundary condition
    x = 0.5;
    vx = 0;
    y = 0;
    vy = 1.63;
    
    // 4th order Runge Kutta
    while(t <= tEnd){
        x = x + vx * h;
        y = y + vy * h;
       	kx1 = h * -GM * e * e * x / pow(x*x + y*y, 1.5);
       	kx2 = h * -GM * e * e * (x + 0.5 * kx1) / pow((x + 0.5 * kx1)*(x + 0.5 * kx1) + y*y, 1.5);
       	kx3 = h * -GM * e * e * (x + 0.5 * kx2) / pow((x + 0.5 * kx2)*(x + 0.5 * kx2) + y*y, 1.5);
       	kx4 = h * -GM * e * e * (x + kx3) / pow((x + kx3)*(x + kx3) + y*y, 1.5);
       	ky1 = h * -GM * e * e * y / pow(x*x + y*y, 1.5);
       	ky2 = h * -GM * e * e * (y + 0.5 * ky1) / pow((y + 0.5 * ky1)*(y + 0.5 * ky1) + x*x, 1.5);
       	ky3 = h * -GM * e * e * (y + 0.5 * ky2) / pow((y + 0.5 * ky2)*(y + 0.5 * ky2) + x*x, 1.5);
       	ky4 = h * -GM * e * e * (y + ky3) / pow((y + ky3)*(y + ky3) + x*x, 1.5);  
       	vx = vx + 1.0/6.0 * (kx1 + 2.0*kx2 + 2.0*kx3 + kx4);
       	vy = vy + 1.0/6.0 * (ky1 + 2.0*ky2 + 2.0*ky3 + ky4);
       	
       	Ekin = 0.5 * (vx*vx + vy*vy);
 		Epot = -GM / sqrt(x*x + y*y);
       	
        fprintf(fid, "%f %f %f %f %f %f\n", t, x, vx, y, vy, Ekin+Epot);
      	t += 1;
        
    }

    fclose(fid);

    return 0;
}

