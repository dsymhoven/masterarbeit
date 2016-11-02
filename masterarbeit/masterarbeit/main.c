//
//  main.c
//  masterarbeit
//
//  Created by David Symhoven on 17.09.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "experiment.h"


int main(int argc, const char * argv[]) {
    
    experiment1();
    
    return 0;
}

//// ======================================================
//    #pragma mark: Initializations
//// ======================================================
//    int numberOfGridPoints = 256;
//    int lengthOfSimulationBox = 32;
//    double dx = (double)lengthOfSimulationBox / numberOfGridPoints;
//    double dy = dx;
//    double dz = dx;
//    double gamma = 2.1;
//    double dir = M_PI/4.0;
//    double dt = 0.5 * dx;
//    double charge = 1;
//    double mass = 1;
//    double t = 0;
//    double tEnd = 20;
//    
//    double Eextern[3] = {0,0,0};
//    double Bextern[3] = {0,0,0.5};
//    
//    double x[4], u[4];
//    double xObserver[4];
//    double edgeOfNearFieldBox[24];
//    
//    char filename[32] = "some";
//    
//    
//// ======================================================
//    #pragma mark: Allocations
//// ======================================================
//
//    
//    // initialize trajectory[tEnd/dt][4]. First allocate array of pointers and then allocate for each pointer another pointer of length 4
//    double **trajectory = (double **) malloc(tEnd / dt * sizeof(double *));
//    if (trajectory == NULL){
//        printf("ERROR: Could not allocate memory for trajectory[]");
//    }
//    else{
//        for (int i = 0; i < tEnd / dt; i++){
//            trajectory[i] = (double *) malloc(4 * sizeof(double));
//            if (trajectory[i] == NULL){
//                printf("ERROR: Could not allocate memory for trajectory[][]");
//            }
//        }
//    }
//    
//    // same for velocity
//    double **velocity = (double **) malloc(tEnd / dt * sizeof(double *));
//    if (velocity == NULL){
//        printf("ERROR: Could not allocate memory for velocity[]");
//    }
//    else{
//        for (int i = 0; i < tEnd / dt; i++){
//            velocity[i] = (double *) malloc(4 * sizeof(double));
//            if (velocity[i] == NULL){
//                printf("ERROR: Could not allocate memory for velocity[][]");
//            }
//        }
//    }
//    
//    
////    double *E = (double *) malloc(3 * numberOfGridPoints * numberOfGridPoints * numberOfGridPoints * sizeof(double));
////    double *B = (double *) malloc(3 * numberOfGridPoints * numberOfGridPoints * numberOfGridPoints * sizeof(double));
////    
////    if( B == NULL || E == NULL){
////        printf("ERROR: Speicher für E und B konnte nicht reserviert werden\n");
////        return -1;
////    }
////    else{
////        printf("Speicher für E und B wurde reserviert.\n");
////    }
////    
////    for (int i = 0; i < 3; i++){
////        E[i] = 0;
////        B[i] = 0;
////    }
//
//    
//// ======================================================
//    # pragma mark: initial conditions
//// ======================================================
//    x[0] = 0;
//    x[1] = 10;
//    x[2] = 14;
//    x[3] = 10;
//    
//    u[0] = gamma;
//    u[1] = cos(dir)*sqrt(gamma*gamma-1.0);
//    u[2] = sin(dir)*sqrt(gamma*gamma-1.0);
//    u[3] = 0;
//    
//    
//    
//
//// ======================================================
//    #pragma mark: main routine
//// ======================================================
//    FILE *fid3 = fopen("completeTrajectory.txt","w");
//    
//    for(int p = 0; p < tEnd / dt; p++){
//        calcualteNearFieldBoxes(x, lengthOfSimulationBox, numberOfGridPoints, edgeOfNearFieldBox);
//        
//        // Write in File for Video
//        if(p % 5 == 0){
//            sprintf(filename, "%d", p);
//            strcat(filename, ".txt");
//            FILE *fid = fopen(filename,"w");
//            sprintf(filename, "nearFieldBox%d", p);
//            strcat(filename, ".txt");
//            FILE *fid2 = fopen(filename,"w");
//            
//            fprintf(fid, "%f %f\n", x[1], x[2]);
//            fprintf(fid2, "%f %f %f %f %f %f %f %f\n", edgeOfNearFieldBox[0], edgeOfNearFieldBox[1], edgeOfNearFieldBox[3], edgeOfNearFieldBox[4], edgeOfNearFieldBox[6], edgeOfNearFieldBox[7], edgeOfNearFieldBox[9], edgeOfNearFieldBox[10]);
//            
//            fclose(fid);
//            fclose(fid2);
//        }
//        
//        fprintf(fid3, "%f %f %f %f %f\n", t, u[1], u[2], x[1], x[2]);
//
//        
//        // save position and velocity of particle for lw field calculations
//        for(int i = 0; i < 4; i++){
//            trajectory[p][i] = x[i];
//            velocity[p][i] = u[i];
//        }
//        
//        
//        xObserver[0] = t;
//        xObserver[1] = edgeOfNearFieldBox[21];
//        xObserver[2] = edgeOfNearFieldBox[22];
//        xObserver[3] = 11; //same height as particle, so shortest distance to travel
//        
//        double lengthOfNearFieldBox = edgeOfNearFieldBox[18] - edgeOfNearFieldBox[21];
//        double numberOfGridPointsOfNearFieldBox = lengthOfNearFieldBox / dx;
//        double intersectionPoint[4];
//        double gamma_sq = 0;
//        double R_sq;
//        double R;
//        double n[3];
//        double beta[3];
//        double betaDot[3];
//        double E[3];
//        double B[3];
//
//        // vorgehen:
//        // 1. teilchen push
//        // 2. loop über nahfeldbox als beobachtungspunkte
//        // 3. check ob alte Position des Teilchens außerhalb des Lichtkegels des Beobachtungspunktes liegt und neue Position innerhalb
//        // 4. Schnittpunkt mit Lichtkegel berechnen --> intersection Point berücksichtigt bereits die endliche Ausbreitungsgschwindigkeit.
//        // 5. intersection Point nutzen um LW Felder am Beobachtungsort zu berechnen
//        
//        
//        // for each observation point loop through the entire trajectory and check which points are inside or outisde the lightcone respectively.
//        // but if first two points of trajectory do not cross the lightcone, then don't bother with the following trajectory points. They can't cross the lightcone either.
//        for(int i = 0; i < numberOfGridPointsOfNearFieldBox; i++){
//            for (int j = 0; j < p; j++){
//                // old position needs to be outside lightcone and new position inside lightcone.
//                // outisde the lightcone means that the information would be traveling faster than c, inside less than c
//                // on the lightcone information travels with c
//                if(!isInsideBackwardLightcone(trajectory[j], xObserver) && isInsideBackwardLightcone(trajectory[j+1], xObserver)){
//                    calculateIntersectionPoint(trajectory[j], trajectory[j+1], xObserver, intersectionPoint);
//                    calculateBeta(trajectory[j], trajectory[j+1], beta);
//                    calculateLienardWiechertParameters(intersectionPoint, xObserver, u, &gamma_sq, &R_sq, &R, n, beta);
//                    calculateBetaDot(velocity[j], velocity[j+1], dt, betaDot);
//                    calcuateLienardWiechertFields(gamma_sq, R_sq, R, n, beta, betaDot, charge, E, B);
//                    printf("%f E = %f in x = %f\n", t, E[0]*E[0] + E[1]*E[1] + E[2]*E[2], xObserver[1]);
//                }
//                else{
//                    break;
//                }
//            }
//            xObserver[1] += dx;
//        }
//        
//        borisPusher(u, Eextern, Bextern, dt, charge/mass);
//        updateLocation(u, x, dt);
//        t += dt;
//
//    }
//    
//    fclose(fid3);
////    free(E);
////    free(B);
//    
//

