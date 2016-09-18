
// x'' = -1/4 M_pi e_0 x / (x^2 + y^2)^(3/2)
// y'' = -1/4 M_pi e_0 y / (x^2 + y^2)^(3/2)

// x(0) = 0.5, y(0) = 0, vx(0) = 0, vy(0) = 1.0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, const char * argv[]) {
    
    // Parameter definitions and allocations
    double h = 0.001;
    double e = 1;
    double GM = 1;
    double t = 0;
    double tEnd = 700;
    double x,xn;
    double vx;
    double y,yn;
    double vy;
    double Ekin,Epot;

    
    // open file identifier
    FILE *fid = fopen("solution_lf.txt","w");
    
    // boundary condition
    x = 0.5;
    vx = 0;
    y = 0;
    vy = 0.2;
    
    // leap frog
    while(t <= tEnd){
    
        xn = x + h * vx - h * h / 2.0 * GM * e * e * x / pow(x*x + y*y, 1.5);
        yn = y + h * vy - h * h / 2.0 * GM * e * e * y / pow(x*x + y*y, 1.5); 
        vx = vx + h / 2.0 * -GM * e * e * ( xn / pow(xn*xn + yn*yn, 1.5) + x / pow(x*x + y*y, 1.5));
        vy = vy + h / 2.0 * -GM * e * e * ( yn / pow(xn*xn + yn*yn, 1.5) + y / pow(x*x + y*y, 1.5)); 
 	
 		x = xn;
 		y = yn;	
       	
       	Ekin = 0.5 * (vx*vx + vy*vy);
 		Epot = -GM / sqrt(x*x + y*y);
 		
        fprintf(fid, "%f %f %f %f %f %f\n", t, x, vx, y, vy, Ekin+Epot);
      	t += 1;
        
    }

    fclose(fid);

    return 0;
}

