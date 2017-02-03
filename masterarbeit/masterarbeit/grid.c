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
#include "calculations.h"



/// @brief initializes all properties of struct Grid.
void initGrid(Grid *Grid, double dx, double dy, double dz, int numberOfGridPointsForBoxInX, int numberOfGridPointsForBoxInY, int numberOfGridPointsForBoxInZ, int numberOfBoxesInX, int numberOfBoxesInY, int numberOfBoxesInZ){
    
    printf("initializing Grid ...\n");
    Grid->dx = dx;
    Grid->dy = dy;
    Grid->dz = dz;
    
    Grid->EMax = 0;
    Grid->HMax = 0;
    
    Grid->upmlLayerWidth = 10;
    
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
/// @param Grid pointer to Grid struct
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
/// @brief allocates memory for H and E fields on the planes at box borders
/// @param Grid pointer to Grid struct
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
/// @brief allocates memeory for UPML coefficients
/// @param Grid pointer to Grid struct
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
///@param Grid pointer to Grid struct
void allocateMemoryOnGrid(Grid *Grid){
    allocateFieldsOnGrid(Grid);
    allocateFieldsOnBoxBorders(Grid);
    allocateUPMLCoefficients(Grid);
    
}

///@brief method for releasing all previously allocated memory in struct Grid. Put all free() invokations in here
///@param Grid pointer to Grid struct
void freeMemoryOnGrid(Grid *Grid){
    printf("releasing allocated memory in Grid...\n");
    free(Grid->B);
    Grid->B = NULL;
    free(Grid->E);
    Grid->E = NULL;
    free(Grid->D);
    Grid->D = NULL;
    free(Grid->H);
    Grid->H = NULL;
    free(Grid->upml1E);
    Grid->upml1E = NULL;
    free(Grid->upml2E);
    Grid->upml2E = NULL;
    free(Grid->upml3E);
    Grid->upml3E = NULL;
    free(Grid->upml4E);
    Grid->upml4E = NULL;
    free(Grid->upml5E);
    Grid->upml5E = NULL;
    free(Grid->upml6E);
    Grid->upml6E = NULL;
    free(Grid->upml1H);
    Grid->upml1H = NULL;
    free(Grid->upml2H);
    Grid->upml2H = NULL;
    free(Grid->upml3H);
    Grid->upml3H = NULL;
    free(Grid->upml4H);
    Grid->upml4H = NULL;
    free(Grid->upml5H);
    Grid->upml5H = NULL;
    free(Grid->upml6H);
    Grid->upml6H = NULL;
    
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
///@param Grid pointer to Grid struct
///@param dt time increment
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

///@brief sets external time dependent E and H field as plane wave.
///@param x four vector containing actual position where you want to evaluate the external field
///@param Eextern vector containing external E field components
///@param Hextern vector containing external H field components
void externalPlaneWave(const double x[4], double tStart, double Eextern[3], double Hextern[3]){
    double E0 = 1;
    double H0 = 1;
    double factor = 1;
    
    Eextern[0] = 0;
    if(x[0] - tStart >= x[1] / factor){
        Eextern[1] = E0 * sin(x[0] - tStart -x[1]/factor);
    }
    else{
        Eextern[1] = 0;
    }
    Eextern[2] = 0;
    
    Hextern[0] = 0;
    Hextern[1] = 0;
    if(x[0] - tStart >= x[1] / factor){
        Hextern[2] = H0 * sin(x[0] - tStart -x[1]/factor);
    }
    else{
        Hextern[2] = 0;
    }
    
    
}

///@brief inits a sample E and B field onto the grid for testing purposes.
///@param Grid pointer to Grid struct
void initSamplePulseOnGrid(Grid *Grid){
    
    int ny = Grid->numberOfGridPointsInY;
    int nz = Grid->numberOfGridPointsInZ;
    int lengthInX = Grid->dx * Grid->numberOfGridPointsInX;
    
        for (int i = 16 * 4; i < 16 * 5; i++)
        {
            for (int j = 1; j < ny; j++)
            {
                for (int l = 1; l < nz; l++)
                {
                    Grid->E[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (l) + 2] = sin(i * 2 * M_PI / lengthInX)*sin(M_PI*j/ny);
    
                    Grid->H[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (l) + 1] = sin(i * 2 * M_PI / lengthInX)*sin(M_PI*j/ny);
                }
            }
    
        }
}

///@brief maxwellPusher for H field.
///@remark H field is calculated via positive curl. Therefore value of E on the right side of H on the grid is required. Thus stop i,j,k with at n - 1 where n denotes the numberOfGridPoints
///@param Grid pointer to Grid struct
///@param dt time increment
void pushHFieldOnGrid(Grid *Grid, double dt){
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


///@brief writes Grid parameters to file, in order for python to use them as plot parameters
///@throws ERROR: Could not open file gridParameters.txt
///@param Grid pointer to Grid struct
void writeGridParametersToFile(Grid *Grid){
    FILE *fid = fopen("gridParameters.txt", "w");
    if (fid == NULL){
        printf("ERROR: Could not open gridParameters.txt");
    }
    else{
        printf("writing grid parameters to file\n");
        fprintf(fid, "%f %f %f %d %d %d %f %f %f %f %f\n", Grid->dx, Grid->dy, Grid->dz, Grid->numberOfGridPointsForBoxInX, Grid->numberOfGridPointsForBoxInY, Grid->numberOfGridPointsForBoxInZ, Grid->lengthOfSimulationBoxInX, Grid->lengthOfSimulationBoxInY, Grid->lengthOfSimulationBoxInZ, Grid->EMax, Grid->HMax);
    }
    fclose(fid);
}

void clearFieldsFromGrid(Grid *Grid){
    printf("Clearing Fields from Grid ...\n");
    for(int i = 0; i < Grid->numberOfGridPointsInX * Grid->numberOfGridPointsInY * Grid->numberOfGridPointsInZ * 3; i++){
        Grid->E[i] = 0;
        Grid->H[i] = 0;
        Grid->D[i] = 0;
        Grid->B[i] = 0;
    }
}


///@brief loops through the specified grid plane and calcutes |B|^2 and |E|^2 of external fields. Values are stored in seperate files. File is structured similiar to the grid.
///@param filename pointer to a char. Gets modified inside the method
///@param index outer loop index. Is used to name the output file
///@param plotE set to true, if you want to write E-field to file
///@param plotB set to true, if you want to write B-field to file
///@throws ERROR: Could not open file for E or B field
void writeExternalFieldsToFile(Grid *Grid, double Eextern[3], double Bextern[3], double t, char *filename, int index, int planeForPlotting, bool plotE, bool plotB){
    printf("Writing external fields to file ...\n");
    FILE *fid = NULL;
    FILE *fid2 = NULL;
    
    if (plotE){
        sprintf(filename, "E_extern%d", index);
        strcat(filename, ".txt");
        fid = fopen(filename,"w");
        if (fid == NULL){
            printf("ERROR: Could not open file E_extern!");
        }
    }
    if(plotB){
        sprintf(filename, "BH_extern%d", index);
        strcat(filename, ".txt");
        fid2 = fopen(filename,"w");
        if (fid2 == NULL){
            printf("ERROR: Could not open file H_extern!");
        }
    }
    
    else{
        int nx = Grid->numberOfGridPointsInX;
        int ny = Grid->numberOfGridPointsInY;

        int k = planeForPlotting;
        
        double Ex, Ey, Ez, Hx, Hy, Hz, Esq, Hsq;
        double x[4]={0};
        x[0] = t;
        x[3] = k * Grid->dz;
        
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                x[1] = i * Grid->dx;
                x[2] = j * Grid->dy;

                externalPlaneWave(x, 0, Eextern, Bextern);
                Ex = Eextern[0];
                Ey = Eextern[1];
                Ez = Eextern[2];
                
                
                Hx = Bextern[0];
                Hy = Bextern[1];
                Hz = Bextern[2];
                
                if(plotE){
                    Esq = Ex * Ex + Ey * Ey + Ez * Ez;
                    fprintf(fid, "%f\t", Esq);
                }
                if (plotB){
                    Hsq = Hx * Hx + Hy * Hy + Hz * Hz;
                    fprintf(fid2, "%f\t", Hsq);
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
    fclose(fid2);
}

