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
    double u[3] = {0,1,0};
    double dt = 0.01;
    double t = 0;
    double tEnd = 100;
    double chargeOverMass = 1;
    
    FILE *fid = fopen("borisPusher.txt", "w");
    
    while(t < tEnd){
        fprintf(fid, "%f %f %f %f\n", t, u[0], u[1], u[2]);
        borisPusher(u, E, B, -0.5*dt, chargeOverMass);
        t += dt;
    }

    
    
    fclose(fid);
    return 0;
}
