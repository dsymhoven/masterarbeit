//
//  grid.c
//  masterarbeit
//
//  Created by David Symhoven on 13.10.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#include "grid.h"
#include "stdlib.h"



void initGrid(Grid *Grid, int numberOfGridPointsInX, int numberOfGridPointsInY, int numberOfGridPointsInZ, double lengthOfSimulationBoxInX, double lengthOfSimulationBoxInY, double lengthOfSimulationBoxInZ){
    printf("initializing Grid ...\n");
    Grid->numberOfGridPointsInX = numberOfGridPointsInX;
    Grid->numberOfGridPointsInY = numberOfGridPointsInY;
    Grid->numberOfGridPointsInZ = numberOfGridPointsInZ;
    Grid->lengthOfSimulationBoxInX = lengthOfSimulationBoxInX;
    Grid->lengthOfSimulationBoxInY = lengthOfSimulationBoxInY;
    Grid->lengthOfSimulationBoxInZ = lengthOfSimulationBoxInZ;
    
    Grid->dx = lengthOfSimulationBoxInX / numberOfGridPointsInX;
    Grid->dy = lengthOfSimulationBoxInY / numberOfGridPointsInY;
    Grid->dz = lengthOfSimulationBoxInZ / numberOfGridPointsInZ;
    
    allocateFieldsOnGrid(Grid);
    
}

void allocateFieldsOnGrid(Grid *Grid){
    printf("allocating Fields ... \n");
    int numberOfGridPointsInX = Grid->numberOfGridPointsInX;
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    int arrayLength = 3 * numberOfGridPointsInX * numberOfGridPointsInY * numberOfGridPointsInZ;
    Grid->E = (double *) malloc(arrayLength * sizeof(double));
    Grid->B = (double *) malloc(arrayLength * sizeof(double));
    
    if( Grid->B == NULL || Grid->E == NULL){
        printf("ERROR: allocation for E and B failed!\n");
    }
    else{
        for (int i = 0; i < arrayLength; i++){
            Grid->E[i] = 0;
            Grid->B[i] = 0;
        }
    }
}

void freeMemoryOnGrid(Grid *Grid){
    printf("releasing allocated memory ...\n");
    free(Grid->B);
    free(Grid->E);
}

void pushEFieldOnGrid(Grid *Grid, double dt){
    double Bx_ijk;
    double By_ijk;
    double Bz_ijk;
    double Bz_ijm1k;
    double By_ijkm1;
    double Bx_ijkm1;
    double Bz_im1jk;
    double By_im1jk;
    double Bx_ijm1k;
    
    double cnx = .5 * dt / Grid->dx;
    double cny = .5 * dt / Grid->dy;
    double cnz = .5 * dt / Grid->dz;
    
    
    int nx = Grid->numberOfGridPointsInX;
    int ny = Grid->numberOfGridPointsInY;
    int nz = Grid->numberOfGridPointsInZ;
    
    int i, j, k;
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            for (k = 0; k < nz; k++)
            {

                Bx_ijk = Grid->B[3 * nz * ny * i + 3 * nz * j + 3 * k + 0];
                By_ijk = Grid->B[3 * nz * ny * i + 3 * nz * j + 3 * k + 1];
                Bz_ijk = Grid->B[3 * nz * ny * i + 3 * nz * j + 3 * k + 2];
                
                Bz_ijm1k = Grid->B[3 * nz * ny * i + 3 * nz * (j - 1) + 3 * k + 2];
                By_ijkm1 = Grid->B[3 * nz * ny * i + 3 * nz * j + 3 * (k - 1) + 1];
                Bx_ijkm1 = Grid->B[3 * nz * ny * i + 3 * nz * j + 3 * (k - 1) + 0];
                Bz_im1jk = Grid->B[3 * nz * ny * (i - 1) + 3 * nz * j + 3 * k + 2];
                By_im1jk = Grid->B[3 * nz * ny * (i - 1) + 3 * nz * (j) + 3 * (k) + 1];
                Bx_ijm1k = Grid->B[3 * nz * ny * (i) + 3 * nz * (j - 1) + 3 * (k) + 0];
                
                
                double val_x = cny * (Bz_ijk - Bz_ijm1k) - cnz * (By_ijk - By_ijkm1);
                double val_y = cnz * (Bx_ijk - Bx_ijkm1) - cnx * (Bz_ijk - Bz_im1jk);
                double val_z = cnx * (By_ijk - By_im1jk) - cny * (Bx_ijk - Bx_ijm1k);
                    
                Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 0] += val_x;
                Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1] += val_y;
                Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2] += val_z;
                
            }
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
