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
#include <string.h>
#include "calculations.h"
#include "grid.h"


int main(int argc, const char * argv[]) {


#pragma mark: Initialisierungen
    int numberOfGridPoints = 256;
    int lengthOfSimulationBox = 32;
    double dx = (double)lengthOfSimulationBox / numberOfGridPoints;
    double dy = dx;
    double dz = dx;
    double gamma = 2.1;
    double dir = M_PI/4.0;
    double dt = 0.5 * dx;
    double charge = 1;
    double mass = 1;
    double t = 0;
    double tEnd = 20;
    
    double Eextern[3] = {0,0,0};
    double Bextern[3] = {0,0,0.5};

    char filename[16] = "some";
    

#pragma mark: Allocations
    double x[4], u[4];
    double xOld[4];
    double xObserver[4];
    double edgeOfNearFieldBox[24];
    
    double *E = (double *) malloc(3 * numberOfGridPoints * numberOfGridPoints * numberOfGridPoints * sizeof(double));
    double *B = (double *) malloc(3 * numberOfGridPoints * numberOfGridPoints * numberOfGridPoints * sizeof(double));
    
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

    
    x[0] = 0;
    x[1] = 10;
    x[2] = 14;
    x[3] = 10;
    
    u[0] = gamma;
    u[1] = cos(dir)*sqrt(gamma*gamma-1.0);
    u[2] = sin(dir)*sqrt(gamma*gamma-1.0);
    u[3] = 0;
    
    // Vorgehen: Startgeschwindigkeit und Position ist gegeben.
    // 1. Teilchen Push aufgrund äußeren Feldes mit borisPusher und updateLocation
    // 2. Zwischenspeichern der alten und der neuen Teilchenposition, um bestimmen zu können, ob Lichtkegel des Beobachtungspunktes durchquert wurde
    // 3. Loop über alle Gitterpunkte, um LW Fields zu berechnen.
    

    FILE *fid2 = fopen("nearFieldBox.txt","w");
    FILE *fid3 = fopen("completeTrajectory.txt","w");
    
    for(int j = 0; j < tEnd / dt; j++){
        // Write in File for Video
        if(j % 5 == 0){
            sprintf(filename, "%d", j);
            strcat(filename, ".txt");
            FILE *fid = fopen(filename,"w");
            fprintf(fid, "%f %f\n", x[1], x[2]);
            fclose(fid);
        }
        
        fprintf(fid3, "%f %f %f %f %f\n", t, u[1], u[2], x[1], x[2]);
        fprintf(fid2, "%f\n", edgeOfNearFieldBox[0]);

        
        borisPusher(u, Eextern, Bextern, dt, charge/mass);
        updateLocation(u, x, dt);
        calcualteNearFieldBoxes(x, lengthOfSimulationBox, numberOfGridPoints, edgeOfNearFieldBox);
        t += dt;

    }
    
    
    fclose(fid2);
    fclose(fid3);
    free(E);
    free(B);
    

    return 0;
}
