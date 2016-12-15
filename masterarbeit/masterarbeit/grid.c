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
#include "calculations.h"


bool useUPML = true;


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
    
    Grid->D = (double *) malloc(arrayLength * sizeof(double));
    Grid->H = (double *) malloc(arrayLength * sizeof(double));
    
    if( Grid->B == NULL || Grid->E == NULL || Grid->D == NULL || Grid->H == NULL){
        printf("ERROR: allocation for E and B failed!\n");
    }
    else{
        for (int i = 0; i < arrayLength; i++){
            Grid->E[i] = 0;
            Grid->B[i] = 0;
            Grid->D[i] = 0;
            Grid->H[i] = 0;
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
    
    Grid->Hz_im1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Hy_im1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Hx_jm1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Hz_jm1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Hx_km1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Hy_km1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Ey_ip1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Ez_ip1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Ez_jp1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Ex_jp1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Ex_kp1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Grid->Ey_kp1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    
    if (Grid->Hz_im1 == NULL || Grid->Hy_im1 == NULL || Grid->Hx_jm1 == NULL || Grid->Hz_jm1 == NULL
        || Grid->Hx_km1 == NULL || Grid->Hy_km1 == NULL || Grid->Ey_ip1 == NULL || Grid->Ez_ip1 == NULL
        || Grid->Ez_jp1 == NULL || Grid->Ex_jp1 == NULL || Grid->Ex_kp1 == NULL || Grid->Ey_kp1 == NULL)
        printf("cannot allocate memory for the fields on box borders!\n");
    
    for (int i = 0; i < numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ; i++)
    {
        if ((Grid->Hz_im1[i] = (double*) malloc(numberOfGridPointsForBoxInY * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Hy_im1[i] = (double*) malloc(numberOfGridPointsForBoxInY * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Hx_jm1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Hz_jm1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Hx_km1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInY * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Grid->Hy_km1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInY * sizeof(double))) == NULL)
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

void allocateUPMLCoefficients(Grid *Grid){
    Grid->upml1E = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml2E = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml3E = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml4E = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml5E = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    Grid->upml6E = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    
    if(Grid->upml1E == NULL || Grid->upml2E == NULL || Grid->upml3E == NULL || Grid->upml4E == NULL || Grid->upml5E == NULL || Grid->upml6E == NULL){
        printf("ERROR: Could not allocate memory for UPML coefficients for E field!");
    }
    
    Grid->upml1H = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml2H = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml3H = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml4H = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml5H = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    Grid->upml6H = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    
    if(Grid->upml1H == NULL || Grid->upml2H == NULL || Grid->upml3H == NULL || Grid->upml4H == NULL || Grid->upml5H == NULL || Grid->upml6H == NULL){
        printf("ERROR: Could not allocate memory for UPML coefficients for H field!");
    }
}

///@brief calls allocateFieldsOnGrid() and allocateFieldsOnBoxBorders(). For details see their documentation
void allocateMemoryOnGrid(Grid *Grid){
    allocateFieldsOnGrid(Grid);
    allocateFieldsOnBoxBorders(Grid);
    allocateUPMLCoefficients(Grid);
    
}

///@brief method for releasing all previously allocated memory in struct Grid. Put all free() invokations in here
void freeMemoryOnGrid(Grid *Grid){
    printf("releasing allocated memory in Grid...\n");
    free(Grid->B);
    free(Grid->E);
    free(Grid->D);
    free(Grid->H);
    free(Grid->upml1E);
    free(Grid->upml2E);
    free(Grid->upml3E);
    free(Grid->upml4E);
    free(Grid->upml5E);
    free(Grid->upml6E);
    free(Grid->upml1H);
    free(Grid->upml2H);
    free(Grid->upml3H);
    free(Grid->upml4H);
    free(Grid->upml5H);
    free(Grid->upml6H);
    
    for (int i = 0; i < Grid->numberOfBoxesInX * Grid->numberOfBoxesInY * Grid->numberOfBoxesInZ; i++)
    {
        
        free(Grid->Hz_im1[i]);
        Grid->Hz_im1[i] = NULL;
        free(Grid->Hy_im1[i]);
        Grid->Hy_im1[i] = NULL;
        free(Grid->Hx_jm1[i]);
        Grid->Hx_jm1[i] = NULL;
        free(Grid->Hz_jm1[i]);
        Grid->Hz_jm1[i] = NULL;
        free(Grid->Hx_km1[i]);
        Grid->Hx_km1[i] = NULL;
        free(Grid->Hy_km1[i]);
        Grid->Hy_km1[i] = NULL;
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
    
    free(Grid->Hz_im1);
    Grid->Hz_im1 = NULL;
    free(Grid->Hy_im1);
    Grid->Hy_im1 = NULL;
    free(Grid->Hx_jm1);
    Grid->Hx_jm1 = NULL;
    free(Grid->Hz_jm1);
    Grid->Hz_jm1 = NULL;
    free(Grid->Hx_km1);
    Grid->Hx_km1 = NULL;
    free(Grid->Hy_km1);
    Grid->Hy_km1 = NULL;
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
    double Hx_ijk;
    double Hy_ijk;
    double Hz_ijk;
    double Hz_ijm1k;
    double Hy_ijkm1;
    double Hx_ijkm1;
    double Hz_im1jk;
    double Hy_im1jk;
    double Hx_ijm1k;
    
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
                
                Hx_ijk = Grid->H[3 * nz * ny * i + 3 * nz * j + 3 * k + 0];
                Hy_ijk = Grid->H[3 * nz * ny * i + 3 * nz * j + 3 * k + 1];
                Hz_ijk = Grid->H[3 * nz * ny * i + 3 * nz * j + 3 * k + 2];
                
                Hz_ijm1k = Grid->H[3 * nz * ny * i + 3 * nz * (j - 1) + 3 * k + 2];
                Hy_ijkm1 = Grid->H[3 * nz * ny * i + 3 * nz * j + 3 * (k - 1) + 1];
                Hx_ijkm1 = Grid->H[3 * nz * ny * i + 3 * nz * j + 3 * (k - 1) + 0];
                Hz_im1jk = Grid->H[3 * nz * ny * (i - 1) + 3 * nz * j + 3 * k + 2];
                Hy_im1jk = Grid->H[3 * nz * ny * (i - 1) + 3 * nz * (j) + 3 * (k) + 1];
                Hx_ijm1k = Grid->H[3 * nz * ny * (i) + 3 * nz * (j - 1) + 3 * (k) + 0];
                
                
                double val_x = cny * (Hz_ijk - Hz_ijm1k) - cnz * (Hy_ijk - Hy_ijkm1);
                double val_y = cnz * (Hx_ijk - Hx_ijkm1) - cnx * (Hz_ijk - Hz_im1jk);
                double val_z = cnx * (Hy_ijk - Hy_im1jk) - cny * (Hx_ijk - Hx_ijm1k);
                
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
                Grid->H[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1] =  exp(-(pow((i - mu) * Grid->dx, 2) + pow((j - mu) * Grid->dy, 2) + pow((k - mu) * Grid->dz, 2)));
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
                
                Grid->H[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 0] += add_x;
                Grid->H[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1] += add_y;
                Grid->H[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2] += add_z;
                
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
        double Ex, Ey, Ez, Hx, Hy, Hz;
        
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                Hx = Grid->H[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 0];
                Hy = Grid->H[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1];
                Hz = Grid->H[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2];
                
                Ex = Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 0];
                Ey = Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1];
                Ez = Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2];
                if(plotE){
                    double Esq = Ex * Ex + Ey * Ey + Ez * Ez;
                    fprintf(fid, "%f\t", Esq);
                }
                if (plotB){
                    double Bsq = Hx * Hx + Hy * Hy + Hz * Hz;
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


void pushEField(Grid *Grid, Particle *Particle, double t, double dt){
    pushEFieldInsideBoxes(Grid, dt);
    setHFieldOnBorders(Grid);
    adjustHFields(Grid, Particle, t);
    pushEFieldAtBorders(Grid, dt);
    
}

void pushHField(Grid *Grid, Particle *Particle, double t, double dt){
    pushHFieldInsideBoxes(Grid, dt);
    setEFieldOnBorders(Grid);
    adjustEFields(Grid, Particle, t);
    pushHFieldAtBorders(Grid, dt);
}

///@brief Maxwell Pusher for EField inside Boxes. This method is used for the hybrid field method.
///@remark E field is calculated via negative curl. Therefore value of H on the left side of E on the grid is required. Thus start i,j,k with 1. UPML is activated by default
void pushEFieldInsideBoxes(Grid *Grid, double dt){
    double Hx_ijk;
    double Hy_ijk;
    double Hz_ijk;
    double Hz_ijm1k;
    double Hy_ijkm1;
    double Hx_ijkm1;
    double Hz_im1jk;
    double Hy_im1jk;
    double Hx_ijm1k;
    
    double cnx = 0.5 * dt / Grid->dx;
    double cny = 0.5 * dt / Grid->dy;
    double cnz = 0.5 * dt / Grid->dz;
    
    double dOld;
    
    int numberOfGridPointsInX = Grid->numberOfGridPointsInX;
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    int i, j, k;
    for (i = 1; i < numberOfGridPointsInX - 1; i++){
        for (j = 1; j < numberOfGridPointsInY - 1; j++){
            for (k = 1; k < numberOfGridPointsInZ - 1; k++){
                
                if (i % Grid->numberOfGridPointsForBoxInX == 0 || j % Grid->numberOfGridPointsForBoxInY == 0 || k % Grid->numberOfGridPointsForBoxInZ == 0){
                    continue;
                }
                Hx_ijk = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                Hy_ijk = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                Hz_ijk = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                
                Hz_ijm1k = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * (j - 1) + 3 * k + 2];
                Hy_ijkm1 = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * (k - 1) + 1];
                Hx_ijkm1 = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * (k - 1) + 0];
                Hz_im1jk = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i - 1) + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                Hy_im1jk = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i - 1) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1];
                Hx_ijm1k = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j - 1) + 3 * (k) + 0];
                
                if (useUPML){
                    dOld = Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                    
                    Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Grid->upml1E[j] * Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Grid->upml2E[j] * ((Hz_ijk - Hz_ijm1k) / Grid->dy - (Hy_ijk - Hy_ijkm1) / Grid->dz);
                    Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Grid->upml3E[k] * Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Grid->upml4E[k] * (Grid->upml5E[i] * Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] - Grid->upml6E[i] * dOld);
                    
                    dOld = Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                    
                    Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Grid->upml1E[k] * Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Grid->upml2E[k] * ((Hx_ijk - Hx_ijkm1) / Grid->dz - (Hz_ijk - Hz_im1jk) / Grid->dx);
                    
                    Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Grid->upml3E[i] * Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Grid->upml4E[i] * (Grid->upml5E[j] * Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] - Grid->upml6E[j] * dOld);
                    
                    dOld = Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                    
                    Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Grid->upml1E[i] * Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Grid->upml2E[i] * ((Hy_ijk - Hy_im1jk) / Grid->dx - (Hx_ijk - Hx_ijm1k) / Grid->dy);
                    
                    
                    Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Grid->upml3E[j] * Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Grid->upml4E[j] * (Grid->upml5E[k] * Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] - Grid->upml6E[k] * dOld);

                }
                else{
                    double val_x = cny * (Hz_ijk - Hz_ijm1k) - cnz * (Hy_ijk - Hy_ijkm1);
                    double val_y = cnz * (Hx_ijk - Hx_ijkm1) - cnx * (Hz_ijk - Hz_im1jk);
                    double val_z = cnx * (Hy_ijk - Hy_im1jk) - cny * (Hx_ijk - Hx_ijm1k);
                    
                    Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 0] += val_x;
                    Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1] += val_y;
                    Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 2] += val_z;
                }
            }
        }
    }
    
}
///@brief Maxwell Pusher for HField inside Boxes. This method is used for the hybrid field method.
///@remark H field is calculated via positive curl. Therefore value of E on the right side of H on the grid is required. Thus stop i,j,k at n - 1 where n denotes the numberOfGridPoints. UPML is activated by default
void pushHFieldInsideBoxes(Grid *Grid, double dt){
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
    
    
    int numberOfGridPointsInX = Grid->numberOfGridPointsInX;
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    double bOld;
    
    int i, j, k;
    for (i = 1; i < numberOfGridPointsInX - 1; i++){
        for (j = 1; j < numberOfGridPointsInY - 1; j++){
            for (k = 1; k < numberOfGridPointsInZ - 1; k++){
                if ((i + 1) % Grid->numberOfGridPointsForBoxInX == 0 || (j + 1) % Grid->numberOfGridPointsForBoxInY == 0 || (k + 1) % Grid->numberOfGridPointsForBoxInZ == 0){
                    continue;
                }
                Ey_ijkp1 = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k + 1) + 1];
                Ey_ijk = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1];
                Ez_ijp1k = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j + 1) + 3 * (k) + 2];
                Ez_ijk = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 2];
                Ez_ip1jk = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i + 1) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 2];
                Ex_ijkp1 = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k + 1) + 0];
                Ex_ijk = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 0];
                Ex_ijp1k = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j + 1) + 3 * (k) + 0];
                Ey_ip1jk = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i + 1) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1];
                
                if(useUPML){
                    
                    bOld = Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                    
                    Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Grid->upml1H[j] * Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Grid->upml2H[j] * ((Ey_ijkp1 - Ey_ijk) / Grid->dz - (Ez_ijp1k - Ez_ijk) / Grid->dy);
                    
                    Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Grid->upml3H[k] * Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Grid->upml4H[k] * (Grid->upml5H[i] * Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] - Grid->upml6E[i] * bOld);
                    
                    bOld = Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                    
                    Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Grid->upml1H[k] * Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Grid->upml2H[k] * ((Ez_ip1jk - Ez_ijk) / Grid->dx - (Ex_ijkp1 - Ex_ijk) / Grid->dz);
                    
                    Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Grid->upml3H[i] * Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Grid->upml4H[i] * (Grid->upml5H[j] * Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] - Grid->upml6H[j] * bOld);
                    
                    bOld = Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                    
                    Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Grid->upml1H[i] * Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Grid->upml2H[i] * ((Ex_ijp1k - Ex_ijk) / Grid->dy - (Ey_ip1jk - Ey_ijk) / Grid->dx);
                    
                    Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Grid->upml3H[j] * Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Grid->upml4H[j] * (Grid->upml5H[k] * Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] - Grid->upml6H[k] * bOld);
                }
                else{
                    
                    double add_x = cnz * (Ey_ijkp1 - Ey_ijk) - cny * (Ez_ijp1k - Ez_ijk);
                    double add_y = cnx * (Ez_ip1jk - Ez_ijk) - cnz * (Ex_ijkp1 - Ex_ijk);
                    double add_z = cny * (Ex_ijp1k - Ex_ijk) - cnx * (Ey_ip1jk - Ey_ijk);
                    
                    Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 0] += add_x;
                    Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1] += add_y;
                    Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 2] += add_z;
                }
            }
        }
    }
}

///@brief this methods loops through all boxes and stores all values of H in the plane left, infront and below of the current box.  As can be seen from the Yee Scheme only Hy and Hz are needed from the left plane. Hz and Hx are needed from the plane infront and Hy and Hx are needed from the plane below. If the current box is a box where te left to the left does not exist (ib = 0) then the respective values for H are set to 0.
///@remark Hz_im1 and all others are matrices. The first index denotes the boxIndex and the second Index the gridPoints in the plane.
///@param Grid pointer to Grid struct
void setHFieldOnBorders(Grid *Grid){
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
                            Grid->Hz_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX - 1) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 2];
                            Grid->Hy_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX - 1) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 1];
                        }
                    }
                }
                else{
                    for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            Grid->Hz_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = 0;
                            Grid->Hy_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = 0;
                        }
                    }
                }
                
                if (jb != 0){
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            Grid->Hx_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY - 1) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 0];
                            Grid->Hz_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY - 1) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 2];
                        }
                    }
                }
                else{
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            Grid->Hx_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = 0;
                            Grid->Hz_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = 0;
                        }
                    }
                }
                
                if (kb != 0){
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                            
                            Grid->Hx_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * (kb * numberOfGridPointsForBoxInZ - 1) + 0];
                            Grid->Hy_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * (kb * numberOfGridPointsForBoxInZ - 1) + 1];
                        }
                    }
                }
                else{
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                            Grid->Hx_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = 0;
                            Grid->Hy_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = 0;
                        }
                    }
                }
            }
            
            
        }
    }
}

///@brief this methods loops through all boxes and stores all values of E in the plane right, behind and above of the current box.  As can be seen from the Yee Scheme only Ey and Ez are needed from the right plane. Ez and Ex are needed from the plane behind and Ey and Ex are needed from the plane above. If the current box is a box where the plane to the right does not exist (ib = numberOfBoxesInX) then the respective values for E are set to 0.
///@remark Ey_ip1 and all others are matrices. The first index denotes the boxIndex and the second Index the gridPoints in the plane.
///@param Grid pointer to Grid struct
void setEFieldOnBorders(Grid *Grid){
    
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int ib, jb, kb;
    
    for (ib = 0; ib < numberOfBoxesInX; ib++){
        for (jb = 0; jb < numberOfBoxesInY; jb++){
            for (kb = 0; kb < numberOfBoxesInZ; kb++){
                
                if (ib != numberOfBoxesInX - 1){
                    for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            
                            Grid->Ey_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * ((ib + 1) * numberOfGridPointsForBoxInX) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 1];
                            Grid->Ez_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * ((ib + 1) * numberOfGridPointsForBoxInX) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 2];
                        }
                    }
                }
                else{
                    for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            
                            Grid->Ey_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = 0;
                            Grid->Ez_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = 0;
                        }
                    }
                }
                
                if (jb != numberOfBoxesInY - 1){
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            
                            Grid->Ez_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * ((jb + 1) * numberOfGridPointsForBoxInY) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 2];
                            Grid->Ex_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * ((jb + 1) * numberOfGridPointsForBoxInY) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 0];
                        }
                    }
                }
                else{
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            
                            Grid->Ez_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = 0;
                            Grid->Ex_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = 0;
                        }
                    }
                }
                
                if (kb != numberOfBoxesInZ - 1){
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                            
                            Grid->Ex_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * ((kb + 1) * numberOfGridPointsForBoxInZ) + 0];
                            Grid->Ey_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * ((kb + 1) * numberOfGridPointsForBoxInZ) + 1];
                        }
                    }
                }
                else{
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                            
                            Grid->Ex_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = 0;
                            Grid->Ey_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = 0;
                        }
                    }
                }
            }
        }
    }
}

///@brief this method loops through all boxes and adjusts the H fields in the plane to the left, infront and below. The actual adjustemnt takes place in "adjustHyz_im1()", "adjustHxz_jm1()" and "adjustHxy_km1()" method.
void adjustHFields(Grid *Grid, Particle *Particle, const double t){
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    for (int ib = 0; ib < numberOfBoxesInX; ib++){
        for (int jb = 0; jb < numberOfBoxesInY; jb++){
            for (int kb = 0; kb < numberOfBoxesInZ; kb++){
                
                int boxIndex = ib * numberOfBoxesInY * numberOfBoxesInZ + jb * numberOfBoxesInZ + kb;
                
                if (ib != 0){
                    adjustHyz_im1(Grid, Particle, boxIndex, ib, jb, kb, t);
                }
                
                if (jb != 0){
                    adjustHxz_jm1(Grid, Particle, boxIndex, ib, jb, kb, t);
                }
                
                if (kb != 0){
                    adjustHxy_km1(Grid, Particle, boxIndex, ib, jb, kb, t);
                }
            }
        }
    }
}

///@brief this method adjusts the H values on the left and right side border of the near field.
///@param Grid pointer to Grid struct
///@param Particle pointer to Particle struct
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustHFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexIm1 is not. Then we are at the left border of the near field. Since we want to push a value for the Efield inside the near field we need (among others) the value to the left. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the right border of the near field, i.e. current box is the not in near field but boxIndexIm1 is, then we push EField values from the far field. For this we need values to the left, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the E field value in the far field correctly.
void adjustHyz_im1(Grid *Grid, Particle *Particle, const int boxIndex, const int ib, const int jb, const int kb, const double t){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int boxIndexIm1 = calcBoxIndexIm1(Grid, boxIndex);
    double xObserver[4];
    
    if (boxIsInNearFieldOfParticle(Grid, Particle, boxIndex)
        && !boxIsInNearFieldOfParticle(Grid, Particle, boxIndexIm1)){
        
        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
            for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX - 1);
                xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                
                subLWField(Grid, Particle, &Grid->Hy_im1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd], xObserver, 4);
                subLWField(Grid, Particle, &Grid->Hz_im1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd], xObserver, 5);
            }
        }
    }
    else if (!boxIsInNearFieldOfParticle(Grid, Particle, boxIndex)
             && boxIsInNearFieldOfParticle(Grid, Particle, boxIndexIm1)){
        
        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
            for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX - 1);
                xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                
                addLWField(Grid, Particle, &(Grid->Hy_im1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd]), xObserver, 4);
                addLWField(Grid, Particle, &(Grid->Hz_im1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd]), xObserver, 5);
            }
        }
        
    }
    
}

///@brief this method adjusts the H values infront and on the back side border of the near field.
///@param Grid pointer to Grid struct
///@param Particle pointer to Particle struct
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustHFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexJm1 is not. Then we are at the facing border of the near field. Since we want to push a value for the Efield inside the near field we need (among others) the value infront. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the back border of the near field, i.e. current box is the not in near field but boxIndexJm1 is, then we push EField values from the far field. For this we need values infront, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the E field value in the far field correctly.
void adjustHxz_jm1(Grid *Grid, Particle *Particle, const int boxIndex, const int ib, const int jb, const int kb, const double t){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int boxIndexJm1 = calcBoxIndexJm1(Grid, boxIndex);
    double xObserver[4];
    
    if (boxIsInNearFieldOfParticle(Grid, Particle, boxIndex)
        && !boxIsInNearFieldOfParticle(Grid, Particle, boxIndexJm1)){
        
        for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
            for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY - 1);
                xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                
                subLWField(Grid, Particle, &Grid->Hx_jm1[boxIndex][numberOfGridPointsForBoxInZ * id + kd], xObserver, 3);
                subLWField(Grid, Particle, &Grid->Hz_jm1[boxIndex][numberOfGridPointsForBoxInZ * id + kd], xObserver, 5);
            }
        }
    }
    else if (!boxIsInNearFieldOfParticle(Grid, Particle, boxIndex)
             && boxIsInNearFieldOfParticle(Grid, Particle, boxIndexJm1)){
        
        for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
            for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY - 1);
                xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                
                addLWField(Grid, Particle, &(Grid->Hx_jm1[boxIndex][numberOfGridPointsForBoxInZ * id + kd]), xObserver, 3);
                addLWField(Grid, Particle, &(Grid->Hz_jm1[boxIndex][numberOfGridPointsForBoxInZ * id + kd]), xObserver, 5);
            }
        }
        
    }
    
}

///@brief this method adjusts the H values on the top and bottom side border of the near field.
///@param Grid pointer to Grid struct
///@param Particle pointer to Particle struct
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustHFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexKm1 is not. Then we are at bottom border of the near field. Since we want to push a value for the Efield inside the near field we need (among others) the value below. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the top border of the near field, i.e. current box is the not in near field but boxIndexKm1 is, then we push EField values from the far field. For this we need values below, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the E field value in the far field correctly.
void adjustHxy_km1(Grid *Grid, Particle *Particle, const int boxIndex, const int ib, const int jb, const int kb, const double t){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int boxIndexKm1 = calcBoxIndexKm1(Grid, boxIndex);
    double xObserver[4];
    
    if (boxIsInNearFieldOfParticle(Grid, Particle, boxIndex)
        && !boxIsInNearFieldOfParticle(Grid, Particle, boxIndexKm1)){
        
        for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
            for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ - 1);
                
                subLWField(Grid, Particle, &Grid->Hx_km1[boxIndex][numberOfGridPointsForBoxInY * id + jd], xObserver, 3);
                subLWField(Grid, Particle, &Grid->Hy_km1[boxIndex][numberOfGridPointsForBoxInY * id + jd], xObserver, 4);
            }
        }
    }
    else if (!boxIsInNearFieldOfParticle(Grid, Particle, boxIndex)
             && boxIsInNearFieldOfParticle(Grid, Particle, boxIndexKm1)){
        
        for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
            for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ - 1);
                
                addLWField(Grid, Particle, &(Grid->Hx_km1[boxIndex][numberOfGridPointsForBoxInY * id + jd]), xObserver, 3);
                addLWField(Grid, Particle, &(Grid->Hy_km1[boxIndex][numberOfGridPointsForBoxInY * id + jd]), xObserver, 4);
            }
        }
        
    }
    
}

///@brief this method loops through all boxes and adjusts the E fields in the plane to the left, infront and below. The actual adjustemnt takes place in "adjustEyz_ip1()", "adjustExz_jp1()" and "adjustExy_kp1()" method.
void adjustEFields(Grid *Grid, Particle *Particle, const double t){
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    for (int ib = 0; ib < numberOfBoxesInX; ib++){
        for (int jb = 0; jb < numberOfBoxesInY; jb++){
            for (int kb = 0; kb < numberOfBoxesInZ; kb++){
                
                int boxIndex = ib * numberOfBoxesInY * numberOfBoxesInZ + jb * numberOfBoxesInZ + kb;
                
                if (ib != numberOfBoxesInX - 1){
                    adjustEyz_ip1(Grid, Particle, boxIndex, ib, jb, kb, t);
                }
                
                if (jb != numberOfBoxesInY - 1){
                    adjustExz_jp1(Grid, Particle, boxIndex, ib, jb, kb, t);
                }
                
                if (kb != numberOfBoxesInZ - 1){
                    adjustExy_kp1(Grid, Particle, boxIndex, ib, jb, kb, t);
                }
            }
        }
    }
}


///@brief this method adjusts the E values on the left and right side border of the near field.
///@param Grid pointer to Grid struct
///@param Particle pointer to Particle struct
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustEFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexIp1 is not. Then we are at the right border of the near field. Since we want to push a value for the Hfield inside the near field we need (among others) the E field value to the right. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the left border of the near field, i.e. current box is the not in near field but boxIndexIp1 is, then we push HField values from the far field. For this we need values to the rigth, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the H field value in the far field correctly.
void adjustEyz_ip1(Grid *Grid, Particle *Particle, const int boxIndex, const int ib, const int jb, const int kb, const double t){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int boxIndexIp1 = calcBoxIndexIp1(Grid, boxIndex);
    double xObserver[4];
    
    if (boxIsInNearFieldOfParticle(Grid, Particle, boxIndex)
        && !boxIsInNearFieldOfParticle(Grid, Particle, boxIndexIp1)){
        
        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
            for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * ((ib + 1) * numberOfGridPointsForBoxInX);
                xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                
                subLWField(Grid, Particle, &Grid->Ey_ip1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd], xObserver, 1);
                subLWField(Grid, Particle, &Grid->Ez_ip1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd], xObserver, 2);
            }
        }
    }
    else if (!boxIsInNearFieldOfParticle(Grid, Particle, boxIndex)
             && boxIsInNearFieldOfParticle(Grid, Particle, boxIndexIp1)){
        
        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
            for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * ((ib + 1) * numberOfGridPointsForBoxInX);
                xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                
                addLWField(Grid, Particle, &(Grid->Ey_ip1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd]), xObserver, 1);
                addLWField(Grid, Particle, &(Grid->Ez_ip1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd]), xObserver, 2);
            }
        }
        
    }
    
}

///@brief this method adjusts the E values infront and on the back side border of the near field.
///@param Grid pointer to Grid struct
///@param Particle pointer to Particle struct
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustEFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexJp1 is not. Then we are at the back side border of the near field. Since we want to push a value for the Hfield inside the near field we need (among others) the E field value behind. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the front side border of the near field, i.e. current box is the not in near field but boxIndexJp1 is, then we push HField values from the far field. For this we need values infront, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the H field value in the far field correctly.
void adjustExz_jp1(Grid *Grid, Particle *Particle, const int boxIndex, const int ib, const int jb, const int kb, const double t){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int boxIndexJp1 = calcBoxIndexJp1(Grid, boxIndex);
    double xObserver[4];
    
    if (boxIsInNearFieldOfParticle(Grid, Particle, boxIndex)
        && !boxIsInNearFieldOfParticle(Grid, Particle, boxIndexJp1)){
        
        for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
            for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                xObserver[2] = (Grid->dy) * ((jb + 1) * numberOfGridPointsForBoxInY);
                xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                
                subLWField(Grid, Particle, &Grid->Ex_jp1[boxIndex][numberOfGridPointsForBoxInZ * id + kd], xObserver, 0);
                subLWField(Grid, Particle, &Grid->Ez_jp1[boxIndex][numberOfGridPointsForBoxInZ * id + kd], xObserver, 2);
            }
        }
    }
    else if (!boxIsInNearFieldOfParticle(Grid, Particle, boxIndex)
             && boxIsInNearFieldOfParticle(Grid, Particle, boxIndexJp1)){
        
        for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
            for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                xObserver[2] = (Grid->dy) * ((jb + 1) * numberOfGridPointsForBoxInY);
                xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                
                addLWField(Grid, Particle, &(Grid->Ex_jp1[boxIndex][numberOfGridPointsForBoxInZ * id + kd]), xObserver, 0);
                addLWField(Grid, Particle, &(Grid->Ez_jp1[boxIndex][numberOfGridPointsForBoxInZ * id + kd]), xObserver, 2);
            }
        }
        
    }
    
}

///@brief this method adjusts the E values at top and bottom side border of the near field.
///@param Grid pointer to Grid struct
///@param Particle pointer to Particle struct
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustEFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexKp1 is not. Then we are at the top side border of the near field. Since we want to push a value for the Hfield inside the near field we need (among others) the E field value above. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the bottom side border of the near field, i.e. current box is the not in near field but boxIndexKp1 is, then we push HField values from the far field. For this we need values above, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the H field value in the far field correctly.
void adjustExy_kp1(Grid *Grid, Particle *Particle, const int boxIndex, const int ib, const int jb, const int kb, const double t){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int boxIndexKp1 = calcBoxIndexKp1(Grid, boxIndex);
    double xObserver[4];
    
    if (boxIsInNearFieldOfParticle(Grid, Particle, boxIndex)
        && !boxIsInNearFieldOfParticle(Grid, Particle, boxIndexKp1)){
        
        for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
            for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                xObserver[3] = (Grid->dz) * ((kb + 1) * numberOfGridPointsForBoxInZ);
                
                subLWField(Grid, Particle, &Grid->Ex_kp1[boxIndex][numberOfGridPointsForBoxInY * id + jd], xObserver, 0);
                subLWField(Grid, Particle, &Grid->Ey_kp1[boxIndex][numberOfGridPointsForBoxInY * id + jd], xObserver, 1);
            }
        }
    }
    else if (!boxIsInNearFieldOfParticle(Grid, Particle, boxIndex)
             && boxIsInNearFieldOfParticle(Grid, Particle, boxIndexKp1)){
        
        for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
            for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                xObserver[3] = (Grid->dz) * ((kb + 1) * numberOfGridPointsForBoxInZ);
                
                addLWField(Grid, Particle, &(Grid->Ex_kp1[boxIndex][numberOfGridPointsForBoxInY * id + jd]), xObserver, 0);
                addLWField(Grid, Particle, &(Grid->Ey_kp1[boxIndex][numberOfGridPointsForBoxInY * id + jd]), xObserver, 1);
            }
        }
        
    }
    
}

///@brief Maxwell Pusher for EField at box borders. This method is used for the hybrid field method.
///@remark E field is calculated via negative curl. Therefore value of H on the left side of E on the grid is required. Thus start i,j,k with 1. UPML is activated by default
void pushEFieldAtBorders(Grid *Grid, double dt){
    
    double Hx_ijk;
    double Hy_ijk;
    double Hz_ijk;
    double Hz_ijm1k;
    double Hy_ijkm1;
    double Hx_ijkm1;
    double Hz_im1jk;
    double Hy_im1jk;
    double Hx_ijm1k;
    
    double dOld;
    
    double cnx = 0.5 * dt / Grid->dx;
    double cny = 0.5 * dt / Grid->dy;
    double cnz = 0.5 * dt / Grid->dz;
    
    int numberOfGridPointsInX = Grid->numberOfGridPointsInX;
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    int i, j, k;
    
    for (i = 1; i < numberOfGridPointsInX - 1; i++){
        for (j = 1; j < numberOfGridPointsInY - 1; j++){
            for (k = 1; k < numberOfGridPointsInZ - 1; k++){
                
                if (i % numberOfGridPointsForBoxInX != 0 && j % numberOfGridPointsForBoxInY != 0 && k % numberOfGridPointsForBoxInZ != 0){
                    continue;
                }
                
                int ib = i / numberOfGridPointsForBoxInX;
                int jb = j / numberOfGridPointsForBoxInY;
                int kb = k / numberOfGridPointsForBoxInZ;
                
                Hx_ijk = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                Hy_ijk = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                Hz_ijk = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                
                
                /*still need to be adjusted...*/
                
                Hx_ijm1k = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * (j - 1) + 3 * k + 0];
                Hy_im1jk = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i - 1) + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                Hz_im1jk = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i - 1) + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                
                Hx_ijkm1 = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * (k - 1) + 0];
                Hy_ijkm1 = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * (k - 1) + 1];
                Hz_ijm1k = Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * (j - 1) + 3 * k + 2];
                
                
                /*adjust values*/
                if (i % numberOfGridPointsForBoxInX == 0){
                    Hz_im1jk = Grid->Hz_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (j % numberOfGridPointsForBoxInY) + (k % numberOfGridPointsForBoxInZ)];
                    Hy_im1jk = Grid->Hy_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (j % numberOfGridPointsForBoxInY) + (k % numberOfGridPointsForBoxInZ)];
                }
                
                if (j % numberOfGridPointsForBoxInY == 0){
                    Hx_ijm1k = Grid->Hx_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (i % numberOfGridPointsForBoxInX) + (k % numberOfGridPointsForBoxInZ)];
                    Hz_ijm1k = Grid->Hz_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (i % numberOfGridPointsForBoxInX) + (k % numberOfGridPointsForBoxInZ)];
                }
                
                if (k % numberOfGridPointsForBoxInZ == 0){
                    Hx_ijkm1 = Grid->Hx_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * (i % numberOfGridPointsForBoxInX) + (j % numberOfGridPointsForBoxInY)];
                    Hy_ijkm1 = Grid->Hy_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * (i % numberOfGridPointsForBoxInX) + (j % numberOfGridPointsForBoxInY)];
                }
                
                if (useUPML){
                    
                    
                    dOld = Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                    
                    Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Grid->upml1E[j] * Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Grid->upml2E[j] * ((Hz_ijk - Hz_ijm1k) / Grid->dy - (Hy_ijk - Hy_ijkm1) / Grid->dz);
                    
                    Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Grid->upml3E[k] * Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Grid->upml4E[k] * (Grid->upml5E[i] * Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] - Grid->upml6E[i] * dOld);
                    
                    dOld = Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                    
                    Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Grid->upml1E[k] * Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Grid->upml2E[k] * ((Hx_ijk - Hx_ijkm1) / Grid->dz - (Hz_ijk - Hz_im1jk) / Grid->dx);
                    
                    Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Grid->upml3E[i] * Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Grid->upml4E[i] * (Grid->upml5E[j] * Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] - Grid->upml6E[j] * dOld);
                    
                    dOld = Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                    
                    Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Grid->upml1E[i] * Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Grid->upml2E[i] * ((Hy_ijk - Hy_im1jk) / Grid->dx - (Hx_ijk - Hx_ijm1k) / Grid->dy);
                    
                    Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Grid->upml3E[j] * Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Grid->upml4E[j] * (Grid->upml5E[k] * Grid->D[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] - Grid->upml6E[k] * dOld);
                }
                else{
                    double val_x = cny * (Hz_ijk - Hz_ijm1k) - cnz * (Hy_ijk - Hy_ijkm1);
                    double val_y = cnz * (Hx_ijk - Hx_ijkm1) - cnx * (Hz_ijk - Hz_im1jk);
                    double val_z = cnx * (Hy_ijk - Hy_im1jk) - cny * (Hx_ijk - Hx_ijm1k);
                    
                    Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 0] += val_x;
                    Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1] += val_y;
                    Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 2] += val_z;
                }
            }
        }
    }
}


///@brief Maxwell Pusher for HField at box borders. This method is used for the hybrid field method.
///@remark H field is calculated via positive curl. Therefore value of E on the right side of H on the grid is required. Thus stop i,j,k at n - 1 where n denotes the numberOfGridPoints. UPML is activated by default
void pushHFieldAtBorders(Grid *Grid, double dt){
    
    double Ex_ijk;
    double Ey_ijk;
    double Ez_ijk;
    double Ey_ijkp1;
    double Ez_ijp1k;
    double Ez_ip1jk;
    double Ex_ijkp1;
    double Ex_ijp1k;
    double Ey_ip1jk;
    
    double bOld;
    
    double cnx = 0.5 * dt / Grid->dx;
    double cny = 0.5 * dt / Grid->dy;
    double cnz = 0.5 * dt / Grid->dz;
    
    int numberOfGridPointsInX = Grid->numberOfGridPointsInX;
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    int i, j, k;
    
    for (i = 1; i < numberOfGridPointsInX - 1; i++){
        for (j = 1; j < numberOfGridPointsInY - 1; j++){
            for (k = 1; k < numberOfGridPointsInZ - 1; k++){
                
                if ((i + 1) % numberOfGridPointsForBoxInX != 0 && (j + 1) % numberOfGridPointsForBoxInY != 0 && (k + 1) % numberOfGridPointsForBoxInZ != 0){
                    continue;
                }
                
                int ib = i / numberOfGridPointsForBoxInX;
                int jb = j / numberOfGridPointsForBoxInY;
                int kb = k / numberOfGridPointsForBoxInZ;
                
                Ex_ijk = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                Ey_ijk = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                Ez_ijk = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                
                
                /*still need to be adjusted...*/
                Ez_ip1jk = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i + 1) + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                Ey_ip1jk = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i + 1) + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                Ex_ijp1k = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * (j + 1) + 3 * k + 0];
                Ez_ijp1k = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * (j + 1) + 3 * k + 2];
                Ex_ijkp1 = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * (k + 1) + 0];
                Ey_ijkp1 = Grid->E[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * (k + 1) + 1];
                
                /*adjust values*/
                if ((i + 1) % numberOfGridPointsForBoxInX == 0){
                    Ez_ip1jk =
                    Grid->Ez_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (j % numberOfGridPointsForBoxInY) + (k % numberOfGridPointsForBoxInZ)];
                    Ey_ip1jk =
                    Grid->Ey_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (j % numberOfGridPointsForBoxInY) + (k % numberOfGridPointsForBoxInZ)];
                    
                }
                
                if ((j + 1) % numberOfGridPointsForBoxInY == 0){
                    Ex_ijp1k =
                    Grid->Ex_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (i % numberOfGridPointsForBoxInX) + (k % numberOfGridPointsForBoxInZ)];
                    Ez_ijp1k =
                    Grid->Ez_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (i % numberOfGridPointsForBoxInX) + (k % numberOfGridPointsForBoxInZ)];
                }
                
                if ((k + 1) % numberOfGridPointsForBoxInZ == 0){
                    Ex_ijkp1 =
                    Grid->Ex_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * (i % numberOfGridPointsForBoxInX) + (j % numberOfGridPointsForBoxInY)];
                    Ey_ijkp1 =
                    Grid->Ey_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * (i % numberOfGridPointsForBoxInX) + (j % numberOfGridPointsForBoxInY)];
                }
                
                if (useUPML){
                    
                    bOld = Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                    
                    Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Grid->upml1H[j] * Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Grid->upml2H[j] * ((Ey_ijkp1 - Ey_ijk) / Grid->dz - (Ez_ijp1k - Ez_ijk) / Grid->dy);
                    
                    Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Grid->upml3H[k] * Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Grid->upml4H[k] * (Grid->upml5H[i] * Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] - Grid->upml6E[i] * bOld);
                    
                    bOld = Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                    
                    Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Grid->upml1H[k] * Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Grid->upml2H[k] * ((Ez_ip1jk - Ez_ijk) / Grid->dx - (Ex_ijkp1 - Ex_ijk) / Grid->dz);
                    
                    Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Grid->upml3H[i] * Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Grid->upml4H[i] * (Grid->upml5H[j] * Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] - Grid->upml6H[j] * bOld);
                    
                    bOld = Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                    
                    Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Grid->upml1H[i] * Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Grid->upml2H[i] * ((Ex_ijp1k - Ex_ijk) / Grid->dy - (Ey_ip1jk - Ey_ijk) / Grid->dx);
                    
                    Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Grid->upml3H[j] * Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Grid->upml4H[j] * (Grid->upml5H[k] * Grid->B[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] - Grid->upml6H[k] * bOld);
                }
                else{
                    
                    double val_x = cnz * (Ey_ijkp1 - Ey_ijk) - cny * (Ez_ijp1k - Ez_ijk);
                    double val_y = cnx * (Ez_ip1jk - Ez_ijk) - cnz * (Ex_ijkp1 - Ex_ijk);
                    double val_z = cny * (Ex_ijp1k - Ex_ijk) - cnx * (Ey_ip1jk - Ey_ijk);
                    
                    Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 0] += val_x;
                    Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1] += val_y;
                    Grid->H[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 2] += val_z;
                }
            }
        }
    }
}


