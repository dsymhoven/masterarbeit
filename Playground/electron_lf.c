
// x'' = -1/4 M_pi e_0 x1 / (x^2 + y^2)^(3/2)
// y'' = -1/4 M_pi e_0 y1 / (x^2 + y^2)^(3/2)

// x(0) = 0.5, y(0) = 0, vx(0) = 0, vy(0) = 1.0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, const char * argv[]) {
    
    // Parameter definitions and allocations
    double h = 0.01;
    double e1 = 1;
    double e2 = 1;
    double GM = 1;
    double t = 0;
    double tEnd = 10000;
    double x1,x1n,y1,y1n;
    double vx1,vy1;
    double x2,x2n,y2,y2n;
    double vx2, vy2;
    double Ekin,Epot;

    
    // open file identifier
    FILE *fid = fopen("solution_lf.txt","w");
    
    // boundary condition
    x1 = 0.5;
    y1 = 0.5;
    vx1 = -0.9;
    vy1 = 0.0;
    
    x2 = 0.0;
    y2 = 0.0;
    vx2 = 0.9;
    vy2 = 0.0;
    
    // leap frog
    while(t <= tEnd){
    
    	// EoM for first particle
        x1n = x1 + h * vx1 - h * h / 2.0 * GM * e1 * e2 * (x1 - x2) / pow((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2), 1.5);
        y1n = y1 + h * vy1 - h * h / 2.0 * GM * e1 * e2 * (y1 - y2) / pow((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2), 1.5);
        vx1 = vx1 + h / 2.0 * -GM * e1 * e2 * ( (x1n - x2n) / pow((x1n - x2n)*(x1n - x2n) + (y1n - y2n)*(y1n - y2n), 1.5) + (x1 - x2) / pow((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2), 1.5));
        vy1 = vy1 + h / 2.0 * -GM * e1 * e2 * ( (y1n - y2n)/ pow((x1n - x2n)*(x1n - x2n) + (y1n - y2n)*(y1n - y2n), 1.5) + (y1 - y2) / pow((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2), 1.5));
	
       	
       	// EoM for second particle
        x2n = x2 + h * vx2 - h * h / 2.0 * GM * e1 * e2 * (x2 - x1) / pow((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1), 1.5);
        y2n = y2 + h * vy2 - h * h / 2.0 * GM * e1 * e2 * (y2 - y1) / pow((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1), 1.5); 
        vx2 = vx2 + h / 2.0 * -GM * e1 * e2 * ( (x2n - x1n) / pow((x2n - x1n)*(x2n - x1n) + (y2n - y1n)*(y2n - y1n), 1.5) + (x2 - x1) / pow((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1), 1.5));
        vy2 = vy2 + h / 2.0 * -GM * e1 * e2 * ( (y2n - y1n) / pow((x2n - x1n)*(x2n - x1n) + (y2n - y1n)*(y2n - y1n), 1.5) + (y2 - y1) / pow((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1), 1.5));
 	
 	 	
 		x1 = x1n;
 		y1 = y1n;
 		
 		x2 = x2n;
 		y2 = y2n;
 		
       	Ekin = 0.5 * (vx2*vx2 + vy1*vy1);
 		Epot = -GM / sqrt(x2*x2 + y1*y1);
 		
        fprintf(fid, "%f %f %f %f %f %f %f %f %f %f\n", t, x1, y1, x2, y2, vx1, vy1, vx2, vy2, Ekin+Epot);
      	t += 1;
        
    }

    fclose(fid);

    return 0;
}

