//
//  grid.c
//  masterarbeit
//
//  Created by David Symhoven on 13.10.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#include "grid.h"
#include "stdlib.h"



void initGrid(Grid Grid, int numberOfGridPoints, double lengthOfSimulationBox){
    printf("initializing Grid ...\n");
    Grid.numberOfGridPoints = numberOfGridPoints;
    Grid.lengthOfSimulationBox = lengthOfSimulationBox;
    allocateFieldsOn(Grid);
    
}

void allocateFieldsOn(Grid Grid){
    printf("allocating Fields ... \n");
    int arrayLength = 3 * Grid.numberOfGridPoints * Grid.numberOfGridPoints * Grid.numberOfGridPoints;
    Grid.E = (double *) malloc(arrayLength * sizeof(double));
    Grid.B = (double *) malloc(arrayLength * sizeof(double));
    
    if( Grid.B == NULL || Grid.E == NULL){
        printf("ERROR: allocation for E and B failed!\n");
    }
    else{
        for (int i = 0; i < arrayLength; i++){
            Grid.E[i] = 0;
            Grid.B[i] = 0;
        }
    }
}
//void calcualteNearFieldBoxes(double x[4], int lengthOfSimulationBox, int numberOfGridPoints, double *edgeOfNearFieldBox){
////    int i, j, k;
////    int numberOfBoxesPerDimension = numberOfGridPoints / lengthOfSimulationBox;
////    double lengthOfOneBox = lengthOfSimulationBox / numberOfBoxesPerDimension;
////    int numberOfGridPointsInOneBox = numberOfGridPoints / numberOfBoxesPerDimension;
////    int numberOfGridPointsInNearFieldBox = numberOfGridPointsInOneBox * 3;
//
////    edgeOfNearFieldBox = (double ***) malloc(numberOfGridPointsInNearFieldBox * sizeof(double **));
////    if (edgeOfNearFieldBox == NULL){
////        printf("An error occured while allocation edgeOfNearFieldBox");
////    }
////    else{
////        for (int i = 0; i < numberOfGridPointsInNearFieldBox; i++){
////            edgeOfNearFieldBox[i] = (double **) malloc(numberOfGridPointsInNearFieldBox * sizeof(double *));
////            if (edgeOfNearFieldBox[i] == NULL){
////                printf("An error occured while allocation edgeOfNearFieldBox at position i");
////            }
////            else{
////                for (int j = 0; j < numberOfGridPointsInNearFieldBox; j++){
////                    edgeOfNearFieldBox[i][j] = (double *) malloc(numberOfGridPointsInNearFieldBox * sizeof(double));
////                    if (edgeOfNearFieldBox[i][j] == NULL){
////                        printf("An error occured while allocation edgeOfNearFieldBox at position i,j");
////                    }
////                    else{
////                        printf("edgeOfNearFieldBox allocated succesfully");
////                    }
////                }
////
////            }
////        }
////    }
//    
//    i = x[1] / lengthOfOneBox;
//    j = x[2] / lengthOfOneBox;
//    k = x[3] / lengthOfOneBox;
//    
//    double xMin = (i - 1) * lengthOfOneBox;
//    double xMax = (i + 2) * lengthOfOneBox;
//    double yMin = (j - 1) * lengthOfOneBox;
//    double yMax = (j + 2) * lengthOfOneBox;
//    double zMin = (k - 1) * lengthOfOneBox;
//    double zMax = (k + 2) * lengthOfOneBox;
//    
////    for (int i = 0; i < numberOfGridPointsInNearFieldBox; i++){
////        for (int j = 0; i < numberOfGridPointsInNearFieldBox; i++){
////            for (int k = 0; i < numberOfGridPointsInNearFieldBox; i++){
////                edgeOfNearFieldBox[i][j][k] = xMin;
////                printf("%f\n", edgeOfNearFieldBox[i][j][k]);
////            }
////        }
////    }
//    
//    edgeOfNearFieldBox[0] = xMin;
//    edgeOfNearFieldBox[1] = yMin;
//    edgeOfNearFieldBox[2] = zMin;
//    
//    edgeOfNearFieldBox[3] = xMax;
//    edgeOfNearFieldBox[4] = yMin;
//    edgeOfNearFieldBox[5] = zMin;
//    
//    edgeOfNearFieldBox[6] = xMax;
//    edgeOfNearFieldBox[7] = yMax;
//    edgeOfNearFieldBox[8] = zMin;
//    
//    edgeOfNearFieldBox[9] = xMin;
//    edgeOfNearFieldBox[10] = yMax;
//    edgeOfNearFieldBox[11] = zMin;
//    
//    edgeOfNearFieldBox[12] = xMin;
//    edgeOfNearFieldBox[13] = yMin;
//    edgeOfNearFieldBox[14] = zMax;
//    
//    edgeOfNearFieldBox[15] = xMax;
//    edgeOfNearFieldBox[16] = yMin;
//    edgeOfNearFieldBox[17] = zMax;
//    
//    edgeOfNearFieldBox[18] = xMax;
//    edgeOfNearFieldBox[19] = yMax;
//    edgeOfNearFieldBox[20] = zMax;
//    
//    edgeOfNearFieldBox[21] = xMin;
//    edgeOfNearFieldBox[22] = yMax;
//    edgeOfNearFieldBox[23] = zMax;
//    
//    
//}
