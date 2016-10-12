//
//  main.c
//  masterarbeit
//
//  Created by David Symhoven on 17.09.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calculations.h"

int main(int argc, const char * argv[]) {

    int nx = 256, ny = 256, nz = 256;
    double dx = 0.1, dy = 0.1, dz = 0.1;
    int lx = 32, ly = 32, lz = 32;
    double x[4], u[4];
    double gamma = 1.1;
    double dir = M_PI/4.0;
    double dt = 0.5 * dx;
    double charge = 1;
    double mass = 1;
    double t = 0;
    double tEnd = 10;
    
    double *E = (double *) malloc(3 * sizeof(double));
    double *B = (double *) malloc(3 * sizeof(double));
    
    if( B == NULL || E == NULL){
        printf("ERROR: Speicher für E und B konnte nicht reserviert werden\n");
        return -1;
    }
    else{
        printf("Speicher für E und B wurde reserviert.\n");
    }
    
    for (int i = 0; i < 3; i++){
        E[i] = 0;
        B[i] = 0;
    }
    
    double Eextern[3] = {0,0,0};
    double Bextern[3] = {0,0,1};
    
    x[0] = 0;
    x[1] = 10;
    x[2] = 10;
    x[3] = 10;
    
    u[0] = gamma;
    u[1] = cos(dir)*sqrt(gamma*gamma-1.0);
    u[2] = sin(dir)*sqrt(gamma*gamma-1.0);
    u[3] = 0;
    
    
    
    FILE *fid = fopen("test.txt","w");
    while(t < tEnd){
        borisPusher(u, Eextern, Bextern, dt, charge/mass);
        updateLocation(u, x, dt);
        t += dt;
        fprintf(fid, "%f %f %f %f %f\n", t, u[1], u[2], x[1], x[2]);
    }
    
    
    fclose(fid);
    free(E);
    free(B);
    

    return 0;
}
