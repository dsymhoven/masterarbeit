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

void allocateFieldsOnBoxBorders(Grid *Grid){
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    int numberOfGridPointsForBoxInX = Grid -> numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid -> numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid -> numberOfGridPointsForBoxInZ;
    
    Grid->Bz_im1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->By_im1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Bx_jm1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Bz_jm1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Bx_km1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->By_km1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Ey_ip1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Ez_ip1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Ez_jp1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Ex_jp1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Ex_kp1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Ey_kp1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    
    if (Grid->Bz_im1 == NULL || Grid->By_im1 == NULL || Grid->Bx_jm1 == NULL || Grid->Bz_jm1 == NULL
        || Grid->Bx_km1 == NULL || Grid->By_km1 == NULL || Grid->Ey_ip1 == NULL || Grid->Ez_ip1 == NULL
        || Grid->Ez_jp1 == NULL || Grid->Ex_jp1 == NULL || Grid->Ex_kp1 == NULL || Grid->Ey_kp1 == NULL)
        printf("cannot allocate memory for the fields on box borders!\n");
    
    for (int i = 0; i < numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ; i++)
    {
        if ((Grid->Bz_im1[i] = (double*) malloc(numberOfGridPointsForBoxInY * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->By_im1[i] = (double*) malloc(numberOfGridPointsForBoxInY * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Bx_jm1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Bz_jm1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Bx_km1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInY * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->By_km1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInY * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Ey_ip1[i] = (double*) malloc(numberOfGridPointsForBoxInY * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Ez_ip1[i] = (double*) malloc(numberOfGridPointsForBoxInY * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Ez_jp1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Ex_jp1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Ex_kp1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInY * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Ey_kp1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInY * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        
    }

}

///@brief method for releasing all previously allocated memory in struct Grid. Put all free() invokations in here
void freeMemoryOnGrid(Grid *Grid){
    printf("releasing allocated memory in Grid...\n");
    free(Grid->B);
    free(Grid->E);
    
    for (int i = 0; i < Grid->numberOfGridPointsForBoxInX * Grid->numberOfGridPointsForBoxInY * Grid->numberOfGridPointsForBoxInZ; i++)
    {
        
        free(Grid->Bz_im1[i]);
        Grid->Bz_im1[i] = NULL;
        free(Grid->By_im1[i]);
        Grid->By_im1[i] = NULL;
        free(Grid->Bx_jm1[i]);
        Grid->Bx_jm1[i] = NULL;
        free(Grid->Bz_jm1[i]);
        Grid->Bz_jm1[i] = NULL;
        free(Grid->Bx_km1[i]);
        Grid->Bx_km1[i] = NULL;
        free(Grid->By_km1[i]);
        Grid->By_km1[i] = NULL;
        free(Grid->Ey_ip1[i]);
        Grid->Ey_ip1[i] = NULL;
        free(Grid->Ez_ip1[i]);
        Grid->Ez_ip1[i] = NULL;
        free(Grid->Ez_jp1[i]);
        Grid->Ez_jp1[i] = NULL;
        free(Grid->Ex_jp1[i]);
        Grid->Ex_jp1[i] = NULL;
        free(Grid->Ex_kp1[i]);
        Grid->Ex_kp1[i] = NULL;
        free(Grid->Ey_kp1[i]);
        Grid->Ey_kp1[i] = NULL;
        
    }
    
    free(Grid->Bz_im1);
    Grid->Bz_im1 = NULL;
    free(Grid->By_im1);
    Grid->By_im1 = NULL;
    free(Grid->Bx_jm1);
    Grid->Bx_jm1 = NULL;
    free(Grid->Bz_jm1);
    Grid->Bz_jm1 = NULL;
    free(Grid->Bx_km1);
    Grid->Bx_km1 = NULL;
    free(Grid->By_km1);
    Grid->By_km1 = NULL;
    free(Grid->Ey_ip1);
    Grid->Ey_ip1 = NULL;
    free(Grid->Ez_ip1);
    Grid->Ez_ip1 = NULL;
    free(Grid->Ez_jp1);
    Grid->Ez_jp1 = NULL;
    free(Grid->Ex_jp1);
    Grid->Ex_jp1 = NULL;
    free(Grid->Ex_kp1);
    Grid->Ex_kp1 = NULL;
    free(Grid->Ey_kp1);
    Grid->Ey_kp1 = NULL;
}

///@brief maxwellPusher for E field. This method is only used for testing purposes. See pushEFieldInsideBoxes for the final version
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
void pushBFieldOnGrid(Grid *Grid, double dt){
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
void writeFieldsToFile(Grid *Grid, char *filename, int index, int planeForPlotting, bool plotE, bool plotB){
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
        int k = planeForPlotting;
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

///@brief Maxwell Pusher for EField inside Boxes. This method is used for the hybrid field method.
///@remark E field is calculated via negative curl. Therefore value of B on the left side of E on the grid is required. Thus start i,j,k with 1
void pushEFieldInsideBoxes(Grid *Grid, double dt){
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
                
                if (i % Grid->numberOfGridPointsForBoxInX == 0 || j % Grid->numberOfGridPointsForBoxInY == 0 || k % Grid->numberOfGridPointsForBoxInZ == 0){
                    continue;
                }
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
///@brief Maxwell Pusher for BField inside Boxes. This method is used for the hybrid field method.
///@remark B field is calculated via positive curl. Therefore value of E on the right side of B on the grid is required. Thus stop i,j,k with at n - 1 where n denotes the numberOfGridPoints
void pushBFieldInsideBoxes(Grid *Grid, double dt){
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
                if (i % Grid->numberOfGridPointsForBoxInX == 0 || j % Grid->numberOfGridPointsForBoxInY == 0 || k % Grid->numberOfGridPointsForBoxInZ == 0){
                    continue;
                }
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

void setBFieldOnBorders(Grid *Grid){
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    int numberOfGridPointsForBoxInX = Grid -> numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid -> numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid -> numberOfGridPointsForBoxInZ;
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    for (int ib = 0; ib < numberOfBoxesInX; ib ++){
        for (int jb = 0; jb < numberOfBoxesInY; jb ++){
            for (int kb = 0; kb < numberOfBoxesInZ; kb ++){
                
                if (ib != 0){
                    for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            Grid->Bz_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = Grid->B[3 * numberOfGridPointsInZ
                                                                                            * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX - 1) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd)
                                                                                            + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 2];
                            Grid->By_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = Grid->B[3 * numberOfGridPointsInZ
                                                                                            * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX - 1) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd)
                                                                                            + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 1];
                        }
                    }
                }
                else{
                    for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            Grid->Bz_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = 0;
                            Grid->By_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = 0;
                        }
                    }
                }
                
                if (jb != 0){
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            Grid->Bx_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = Grid->B[3 * numberOfGridPointsInZ
                                                                                            * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY - 1)
                                                                                            + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 0];
                            Grid->Bz_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = Grid->B[3 * numberOfGridPointsInZ
                                                                                            * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY - 1)
                                                                                            + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 2];
                        }
                    }
                }
                else{
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            Grid->Bx_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = 0;
                            Grid->Bz_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = 0;
                        }
                    }
                }
                
                if (kb != 0){
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                        
                            Grid->Bx_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = Grid->B[3 * numberOfGridPointsInZ
                                                                                            * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd)
                                                                                            + 3 * (kb * numberOfGridPointsForBoxInZ - 1) + 0];
                            Grid->By_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = Grid->B[3 * numberOfGridPointsInZ
                                                                                            * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd)
                                                                                            + 3 * (kb * numberOfGridPointsForBoxInZ - 1) + 1];
                        }
                    }
                }
                else{
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                            Grid->Bx_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = 0;
                            Grid->By_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = 0;
                        }
                    }
                }
            }

                
            }
        }
}

void SetEFieldOnBorders(Grid *Grid){
    
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    int numberOfBoxesInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfBoxesInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfBoxesInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int ib, jb, kb;
    
    for (ib = 0; ib < numberOfBoxesInX; ib++)
        for (jb = 0; jb < numberOfBoxesInY; jb++)
            for (kb = 0; kb < numberOfBoxesInZ; kb++)
            {
                
                if (ib != numberOfBoxesInX - 1)
                    for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++)
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++)
                        {
                            Grid->Ey_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = Grid->E[3 * numberOfGridPointsInZ
                                                                                            * numberOfGridPointsInY * ((ib + 1) * numberOfGridPointsForBoxInX) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd)
                                                                                            + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 1];
                            Grid->Ez_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = Grid->E[3 * numberOfGridPointsInZ
                                                                                            * numberOfGridPointsInY * ((ib + 1) * numberOfGridPointsForBoxInX) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd)
                                                                                            + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 2];
                        }
                else
                    for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++)
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++)
                        {
                            Grid->Ey_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = 0;
                            Grid->Ez_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = 0;
                        }
                
                if (jb != numberOfBoxesInY - 1)
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++)
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++)
                        {
                            Grid->Ez_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = Grid->E[3 * numberOfGridPointsInZ
                                                                                            * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * ((jb + 1) * numberOfGridPointsForBoxInY)
                                                                                            + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 2];
                            Grid->Ex_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = Grid->E[3 * numberOfGridPointsInZ
                                                                                            * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * ((jb + 1) * numberOfGridPointsForBoxInY)
                                                                                            + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 0];
                        }
                else
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++)
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++)
                        {
                            Grid->Ez_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = 0;
                            Grid->Ex_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = 0;
                        }
                
                if (kb != numberOfBoxesInZ - 1)
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++)
                        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++)
                        {
                            Grid->Ex_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = Grid->E[3 * numberOfGridPointsInZ
                                                                                            * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd)
                                                                                            + 3 * ((kb + 1) * numberOfGridPointsForBoxInZ) + 0];
                            Grid->Ey_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = Grid->E[3 * numberOfGridPointsInZ
                                                                                            * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd)
                                                                                            + 3 * ((kb + 1) * numberOfGridPointsForBoxInZ) + 1];
                        }
                else
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++)
                        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++)
                        {
                            Grid->Ex_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = 0;
                            Grid->Ey_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = 0;
                        }
                
            }
    
}



