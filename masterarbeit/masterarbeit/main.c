//
//  main.c
//  masterarbeit
//
//  Created by David Symhoven on 17.09.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include "calculations.h"

int main(int argc, const char * argv[]) {

    double B[3] = {0,0,1};
    double E[3] = {0,0,0};
    double u[4] = {0,0,1,0};
    double dt = 0.01;
    double t = 0;
    double tEnd = 15;
    double chargeOverMass = 1;
    
    FILE *fid = fopen("borisPusher.txt", "w");
    
    while(t < tEnd){
        fprintf(fid, "%f %f %f %f\n", t, u[1], u[2], u[3]);
        //push_u_boris(u, chargeOverMass, dt, E, B);
        borisPusher(u, E, B, -0.5*dt, chargeOverMass);
        t += dt;
        
    }
    
    
    
    fclose(fid);
    return 0;
}
