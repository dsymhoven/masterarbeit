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

    char filename[32] = "some";

#pragma mark: Allocations
    double x[4], u[4];
    double xOld[4];
    double uOld[4];
    double uNew[4];
    double xNew[4];
    double xObserver[4];
    double edgeOfNearFieldBox[24];
    double xTrajectory[320][4] = {0};
//    double *E = (double *) malloc(3 * numberOfGridPoints * numberOfGridPoints * numberOfGridPoints * sizeof(double));
//    double *B = (double *) malloc(3 * numberOfGridPoints * numberOfGridPoints * numberOfGridPoints * sizeof(double));
//    
//    if( B == NULL || E == NULL){
//        printf("ERROR: Speicher für E und B konnte nicht reserviert werden\n");
//        return -1;
//    }
//    else{
//        printf("Speicher für E und B wurde reserviert.\n");
//    }
//    
//    for (int i = 0; i < 3; i++){
//        E[i] = 0;
//        B[i] = 0;
//    }

    
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
    

    FILE *fid3 = fopen("completeTrajectory.txt","w");
    
    for(int p = 0; p < tEnd / dt; p++){
        calcualteNearFieldBoxes(x, lengthOfSimulationBox, numberOfGridPoints, edgeOfNearFieldBox);
        
        // Write in File for Video
        if(p % 5 == 0){
            sprintf(filename, "%d", p);
            strcat(filename, ".txt");
            FILE *fid = fopen(filename,"w");
            sprintf(filename, "nearFieldBox%d", p);
            strcat(filename, ".txt");
            FILE *fid2 = fopen(filename,"w");
            
            fprintf(fid, "%f %f\n", x[1], x[2]);
            fprintf(fid2, "%f %f %f %f %f %f %f %f\n", edgeOfNearFieldBox[0], edgeOfNearFieldBox[1], edgeOfNearFieldBox[3], edgeOfNearFieldBox[4], edgeOfNearFieldBox[6], edgeOfNearFieldBox[7], edgeOfNearFieldBox[9], edgeOfNearFieldBox[10]);
            
            fclose(fid);
            fclose(fid2);
        }
        
        fprintf(fid3, "%f %f %f %f %f\n", t, u[1], u[2], x[1], x[2]);

//        for (int i = 0; i < 4; i++){
//            xOld[i] = x[i];
//            uOld[i] = u[i];
//        }
        
        // save position of particle for lienard wiechert fields calculation
//        for(int i = 0; i < 4; i++){
//            xTrajectory[p][i] = x[i];
//        }
//        
        xOld[0] = 0;
        xOld[1] = 10;
        xOld[2] = 14;
        xOld[3] = 10;
        
        xNew[0] = 0.0625;
        xNew[1] = 10.084120;
        xNew[2] = 14.079020;
        xNew[3] = 10.0;
        
        uOld[0] = gamma;
        uOld[1] = cos(dir)*sqrt(gamma*gamma-1.0);
        uOld[2] = sin(dir)*sqrt(gamma*gamma-1.0);
        uOld[3] = 0;
        
        uNew[0] = gamma + dt;
        uNew[1] = 1.3459;
        uNew[2] = 1.2643;
        uNew[3] = 0.0;
        
        xObserver[0] = t;
        xObserver[1] = edgeOfNearFieldBox[21];
        xObserver[2] = edgeOfNearFieldBox[22];
        xObserver[3] = 11;//edgeOfNearFieldBox[23];
        
        double lengthOfNearFieldBox = edgeOfNearFieldBox[18] - edgeOfNearFieldBox[21];
        double numberOfGridPointsOfNearFieldBox = lengthOfNearFieldBox / dx;
        double intersectionPoint[4];
        double gamma_sq = 0;
        double R_sq;
        double R;
        double n[3];
        double beta[3];
        double betaDot[3];
        double E[3];
        double B[3];

        // vorgehen:
        // 1. teilchen push
        // 2. loop über nahfeldbox als beobachtungspunkte
        // 3. check ob alte Position des Teilchens außerhalb des Lichtkegels des Beobachtungspunktes liegt und neue Position innerhalb
        // 4. Schnittpunkt mit Lichtkegel berechnen --> intersection Point berücksichtigt bereits die endliche Ausbreitungsgschwindigkeit.
        // 5. intersection Point nutzen um LW Felder am Beobachtungsort zu berechnen
        
        xObserver[0] = t;
        for(int i = 0; i < numberOfGridPointsOfNearFieldBox; i++){
            if(!isInsideBackwardLightcone(xOld, xObserver) && isInsideBackwardLightcone(xNew, xObserver)){
                calculateIntersectionPoint(xOld, xNew, xObserver, intersectionPoint);
                calculateBeta(xOld, xNew, beta);
                calculateLienardWiechertParameters(intersectionPoint, xObserver, u, &gamma_sq, &R_sq, &R, n, beta);
                calculateBetaDot(uOld, uNew, dt, betaDot);
                calcuateLienardWiechertFields(gamma_sq, R_sq, R, n, beta, betaDot, charge, E, B);
                printf("%f E = %f in x = %f\n", t, E[0]*E[0] + E[1]*E[1] + E[2]*E[2], xObserver[1]);
            }
            xObserver[1] += dx;
        }
        
        borisPusher(u, Eextern, Bextern, dt, charge/mass);
        updateLocation(u, x, dt);
        t += dt;

    }
    
    fclose(fid3);
//    free(E);
//    free(B);
    

    return 0;
}
