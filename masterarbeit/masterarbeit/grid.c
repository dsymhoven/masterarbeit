//
//  grid.c
//  masterarbeit
//
//  Created by David Symhoven on 13.10.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#include "grid.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "particle.h"


/// @brief initializes all properties of struct Grid.
void initGrid(Grid *Grid, double dx, double dy, double dz, int numberOfGridPointsForBoxInX, int numberOfGridPointsForBoxInY, int numberOfGridPointsForBoxInZ, int numberOfBoxesInX, int numberOfBoxesInY, int numberOfBoxesInZ){
    
    printf("initializing Grid ...\n");
    Grid->dx = dx;
    Grid->dy = dy;
    Grid->dz = dz;
    
    Grid->numberOfBoxesInX = numberOfBoxesInX;
    Grid->numberOfBoxesInY = numberOfBoxesInY;
    Grid->numberOfBoxesInZ = numberOfBoxesInZ;
    
    Grid->numberOfGridPointsForBoxInX = numberOfGridPointsForBoxInX;
    Grid->numberOfGridPointsForBoxInY = numberOfGridPointsForBoxInY;
    Grid->numberOfGridPointsForBoxInZ = numberOfGridPointsForBoxInZ;
    
    Grid->numberOfGridPointsInX = numberOfGridPointsForBoxInX * numberOfBoxesInX;
    Grid->numberOfGridPointsInY = numberOfGridPointsForBoxInY * numberOfBoxesInY;
    Grid->numberOfGridPointsInZ = numberOfGridPointsForBoxInZ * numberOfBoxesInZ;
    
    Grid->lengthOfSimulationBoxInX = Grid->numberOfGridPointsInX * dx;
    Grid->lengthOfSimulationBoxInY = Grid->numberOfGridPointsInY * dy;
    Grid->lengthOfSimulationBoxInZ = Grid->numberOfGridPointsInZ * dz;
}
/// @brief Allocation of E and B field array. Inits E and B with 0 by default
/// @remark arrayLength = 3 * numberOfGridPointsInX * numberOfGridPointsInY * numberOfGridPointsInZ. Factor 3 because we need x,y and z components on each grid point
/// @throws ERROR: allocation for E and B failed! if memory couldn't be allocated
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

///@brief method for releasing all previously allocated memory in struct Grid. Put all free() invokations in here
void freeMemoryOnGrid(Grid *Grid){
    printf("releasing allocated memory in Grid...\n");
    free(Grid->B);
    free(Grid->E);
}

///@brief maxwellPusher for E field.
///@remark E field is calculated via negative curl. Therefore value of B on the left side of E on the grid is required. Thus start i,j,k with 1
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
    
    double cnx = 0.5 * dt / Grid->dx;
    double cny = 0.5 * dt / Grid->dy;
    double cnz = 0.5 * dt / Grid->dz;
    
    
    int nx = Grid->numberOfGridPointsInX;
    int ny = Grid->numberOfGridPointsInY;
    int nz = Grid->numberOfGridPointsInZ;
    
    int i, j, k;
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            for (k = 1; k < nz - 1; k++)
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

///@brief inits a sample E and B field onto the grid for testing purposes.
void initSamplePulseOnGrid(Grid *Grid){
    
    int ny = Grid->numberOfGridPointsInY;
    int nz = Grid->numberOfGridPointsInZ;
    double mu = 128;
    
    for (int i = 32*2; i < 32*6; i++){
        for (int j = 32*2; j < 32*6; j++){
            for (int k = 32*2; k < 32*6; k++){
                Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2] =  exp(-(pow((i - mu) * Grid->dx, 2) + pow((j - mu) * Grid->dy, 2) + pow((k - mu) * Grid->dz, 2)));
                Grid->B[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1] =  exp(-(pow((i - mu) * Grid->dx, 2) + pow((j - mu) * Grid->dy, 2) + pow((k - mu) * Grid->dz, 2)));
            }
        }
    }
    //    for (int i = 16 * 4; i < 16 * 5; i++)
    //    {
    //        for (int j = 1; j < ny; j++)
    //        {
    //            for (int l = 1; l < nz; l++)
    //            {
    //                Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (l) + 2] = sin(i * 2 * M_PI / 32)*sin(M_PI*j/256);
    //
    //                Grid->B[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (l) + 1] = sin(i * 2 * M_PI / 32)*sin(M_PI*j/256);
    //            }
    //        }
    //
    //    }
}

///@brief maxwellPusher for B field.
///@remark B field is calculated via positive curl. Therefore value of E on the right side of B on the grid is required. Thus stop i,j,k with at n - 1 where n denotes the numberOfGridPoints
void PushBFieldOnGrid(Grid *Grid, double dt){
    double Ex_ijk;
    double Ey_ijk;
    double Ez_ijk;
    double Ey_ijkp1;
    double Ez_ijp1k;
    double Ez_ip1jk;
    double Ex_ijkp1;
    double Ex_ijp1k;
    double Ey_ip1jk;
    
    double cnx = 0.5 * dt / Grid->dx;
    double cny = 0.5 * dt / Grid->dy;
    double cnz = 0.5 * dt / Grid->dz;
    
    
    int nx = Grid->numberOfGridPointsInX;
    int ny = Grid->numberOfGridPointsInY;
    int nz = Grid->numberOfGridPointsInZ;
    
    int i, j, k;
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            for (k = 1; k < nz - 1; k++)
            {
                
                Ey_ijkp1 = Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k + 1) + 1];
                Ey_ijk = Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1];
                Ez_ijp1k = Grid->E[3 * nz * ny * (i) + 3 * nz * (j + 1) + 3 * (k) + 2];
                Ez_ijk = Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2];
                Ez_ip1jk = Grid->E[3 * nz * ny * (i + 1) + 3 * nz * (j) + 3 * (k) + 2];
                Ex_ijkp1 = Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k + 1) + 0];
                Ex_ijk = Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 0];
                Ex_ijp1k = Grid->E[3 * nz * ny * (i) + 3 * nz * (j + 1) + 3 * (k) + 0];
                Ey_ip1jk = Grid->E[3 * nz * ny * (i + 1) + 3 * nz * (j) + 3 * (k) + 1];
                
                
                double add_x = cnz * (Ey_ijkp1 - Ey_ijk) - cny * (Ez_ijp1k - Ez_ijk);
                double add_y = cnx * (Ez_ip1jk - Ez_ijk) - cnz * (Ex_ijkp1 - Ex_ijk);
                double add_z = cny * (Ex_ijp1k - Ex_ijk) - cnx * (Ey_ip1jk - Ey_ijk);
                
                Grid->B[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 0] += add_x;
                Grid->B[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1] += add_y;
                Grid->B[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2] += add_z;
                
            }
        }
    }
}

///@brief loops through the entire E and B array and writes |B|^2 and |E|^2 to seperate files. File is structured similiar to the grid.
///@param filename pointer to a char. Gets modified inside the method
///@param index outer loop index. Is used to name the output file
///@param plotE set to true, if you want to write E-field to file
///@param plotB set to true, if you want to write B-field to file
///@throws ERROR: Could not open file for E or B field
void writeFieldsToFile(Grid *Grid, char *filename, int index, bool plotE, bool plotB){
    printf("Writing fields to file ...\n");
    FILE *fid = NULL;
    FILE *fid2 = NULL;
    
    if (plotE){
        sprintf(filename, "E_field%d", index);
        strcat(filename, ".txt");
        fid = fopen(filename,"w");
        if (fid == NULL){
            printf("ERROR: Could not open file E_field!");
        }
    }
    if(plotB){
        sprintf(filename, "B_field%d", index);
        strcat(filename, ".txt");
        fid2 = fopen(filename,"w");
        if (fid2 == NULL){
            printf("ERROR: Could not open file B_field!");
        }
    }
    
    else{
        int nx = Grid->numberOfGridPointsInX;
        int ny = Grid->numberOfGridPointsInY;
        int nz = Grid->numberOfGridPointsInZ;
        int k = 128;
        double Ex, Ey, Ez, Bx, By, Bz;
        
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                Bx = Grid->B[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 0];
                By = Grid->B[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1];
                Bz = Grid->B[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2];
                
                Ex = Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 0];
                Ey = Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1];
                Ez = Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2];
                if(plotE){
                    double Esq = Ex * Ex + Ey * Ey + Ez * Ez;
                    fprintf(fid, "%f\t", Esq);
                }
                if (plotB){
                    double Bsq = Bx * Bx + By * By + Bz * Bz;
                    fprintf(fid2, "%f\t", Bsq);
                }
            }
            if(plotE){
                fprintf(fid,"\n");
            }
            if(plotB){
                fprintf(fid2,"\n");
            }
        }
    }
    fclose(fid);
}

///@brief writes Grid parameters to file, in order for python to use them as plot parameters
///@throws ERROR: Could not open file gridParameters.txt
void writeGridParametersToFile(Grid *Grid){
    FILE *fid = fopen("gridParameters.txt", "w");
    if (fid == NULL){
        printf("ERROR: Could not open gridParameters.txt");
    }
    else{
        printf("writing grid parameters to file\n");
        fprintf(fid, "%f %f %f %d %d %d %f %f %f\n", Grid->dx, Grid->dy, Grid->dz, Grid->numberOfGridPointsForBoxInX, Grid->numberOfGridPointsForBoxInY, Grid->numberOfGridPointsForBoxInZ, Grid->lengthOfSimulationBoxInX, Grid->lengthOfSimulationBoxInY, Grid->lengthOfSimulationBoxInZ);
    }
    fclose(fid);
}
