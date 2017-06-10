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
    printf("allocating Fields on box borders... \n");
    
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    int numberOfGridPointsForBoxInX = Grid -> Box.numberOfGridPointsInX;
    int numberOfGridPointsForBoxInY = Grid -> Box.numberOfGridPointsInY;
    int numberOfGridPointsForBoxInZ = Grid -> Box.numberOfGridPointsInZ;
    
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
    
    printf("allocating UPML Coefficients ... \n");
    
    Grid->upml1Ex = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml2Ex = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml3Ex = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml4Ex = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml5Ex = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    Grid->upml6Ex = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    
    if(Grid->upml1Ex == NULL || Grid->upml2Ex == NULL || Grid->upml3Ex == NULL || Grid->upml4Ex == NULL || Grid->upml5Ex == NULL || Grid->upml6Ex == NULL){
        printf("ERROR: Could not allocate memory for UPML coefficients for Ex field!");
    }
    
    Grid->upml1Ey = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml2Ey = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml3Ey = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    Grid->upml4Ey = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    Grid->upml5Ey = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml6Ey = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    
    if(Grid->upml1Ey == NULL || Grid->upml2Ey == NULL || Grid->upml3Ey == NULL || Grid->upml4Ey == NULL || Grid->upml5Ey == NULL || Grid->upml6Ey == NULL){
        printf("ERROR: Could not allocate memory for UPML coefficients for Ey field!");
    }
    
    Grid->upml1Ez = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    Grid->upml2Ez = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    Grid->upml3Ez = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml4Ez = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml5Ez = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml6Ez = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    
    if(Grid->upml1Ez == NULL || Grid->upml2Ez == NULL || Grid->upml3Ez == NULL || Grid->upml4Ez == NULL || Grid->upml5Ez == NULL || Grid->upml6Ez == NULL){
        printf("ERROR: Could not allocate memory for UPML coefficients for Ez field!");
    }
    
    Grid->upml1Hx = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml2Hx = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml3Hx = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml4Hx = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml5Hx = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    Grid->upml6Hx = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    
    if(Grid->upml1Hx == NULL || Grid->upml2Hx == NULL || Grid->upml3Hx == NULL || Grid->upml4Hx == NULL || Grid->upml5Hx == NULL || Grid->upml6Hx == NULL){
        printf("ERROR: Could not allocate memory for UPML coefficients for Hx field!");
    }
    
    Grid->upml1Hy = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml2Hy = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml3Hy = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    Grid->upml4Hy = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    Grid->upml5Hy = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml6Hy = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    
    if(Grid->upml1Hy == NULL || Grid->upml2Hy == NULL || Grid->upml3Hy == NULL || Grid->upml4Hy == NULL || Grid->upml5Hy == NULL || Grid->upml6Hy == NULL){
        printf("ERROR: Could not allocate memory for UPML coefficients for Hy field!");
    }
    
    Grid->upml1Hz = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    Grid->upml2Hz = (double *) malloc(Grid->numberOfGridPointsInX * sizeof(double));
    Grid->upml3Hz = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml4Hz = (double *) malloc(Grid->numberOfGridPointsInY * sizeof(double));
    Grid->upml5Hz = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    Grid->upml6Hz = (double *) malloc(Grid->numberOfGridPointsInZ * sizeof(double));
    
    if(Grid->upml1Hz == NULL || Grid->upml2Hz == NULL || Grid->upml3Hz == NULL || Grid->upml4Hz == NULL || Grid->upml5Hz == NULL || Grid->upml6Hz == NULL){
        printf("ERROR: Could not allocate memory for UPML coefficients for Hz field!");
    }
}

///@brief calls allocateFieldsOnGrid() and allocateFieldsOnBoxBorders(). For details see their documentation
///@param Grid pointer to Grid struct
void allocateMemoryOnGrid(Grid *Grid){
    allocateFieldsOnGrid(Grid);
    allocateFieldsOnBoxBorders(Grid);
    
    if(Grid->useUPML){
        allocateUPMLCoefficients(Grid);
    }
    
}

///@brief initializes all properties of struct Grid and allocates memory for all arrays including UPML
///@param Grid pointer to Grid struct
///@param Resolution pointer to Resolution struct
///@param Box pointer to Box struct
///@param numberOfBoxesInX number of Boxes in x direction
///@param numberOfBoxesInY number of Boxes in y direction
///@param numberOfBoxesInZ number of Boxes in z direction
///@param useUPML set this value to true and UPML will be used
void initGrid(Grid *Grid, Resolution *Resolution, Box *Box, const int numberOfBoxesInX, const int numberOfBoxesInY, const int numberOfBoxesInZ, bool useUPML){
    
    printf("initializing Grid ...\n");

    Grid->EMax = 0;
    Grid->HMax = 0;
    
    Grid->Resolution = *Resolution;
    
    Grid->numberOfBoxesInX = numberOfBoxesInX;
    Grid->numberOfBoxesInY = numberOfBoxesInY;
    Grid->numberOfBoxesInZ = numberOfBoxesInZ;
    
    Grid->Box = *Box;
    
    Grid->numberOfGridPointsInX = Box->numberOfGridPointsInX * numberOfBoxesInX;
    Grid->numberOfGridPointsInY = Box->numberOfGridPointsInY * numberOfBoxesInY;
    Grid->numberOfGridPointsInZ = Box->numberOfGridPointsInZ * numberOfBoxesInZ;
    
    Grid->lengthOfSimulationBoxInX = Grid->numberOfGridPointsInX * Resolution->dx;
    Grid->lengthOfSimulationBoxInY = Grid->numberOfGridPointsInY * Resolution->dy;
    Grid->lengthOfSimulationBoxInZ = Grid->numberOfGridPointsInZ * Resolution->dz;
    
    Grid->useUPML = useUPML;
    Grid->upmlLayerWidth = 10;
    
    allocateMemoryOnGrid(Grid);
    
    if (useUPML){
        calcUPMLCoefficients(Grid);
    }
}

///@brief Releases all allocated memory for fields.
///@param Grid pointer to Grid struct
void freeFieldsOnGrid(Grid *Grid){
    printf("releasing allocated memory for fields...\n");
    free(Grid->B);
    Grid->B = NULL;
    free(Grid->E);
    Grid->E = NULL;
    free(Grid->D);
    Grid->D = NULL;
    free(Grid->H);
    Grid->H = NULL;
}

///@brief Releases all allocated memory for UPML Coefficients.
///@param Grid pointer to Grid struct
void freeUPMLCoefficients(Grid *Grid){
    printf("releasing allocated memory for UPML Coefficients...\n");
    
    free(Grid->upml1Ex);
    Grid->upml1Ex = NULL;
    free(Grid->upml2Ex);
    Grid->upml2Ex = NULL;
    free(Grid->upml3Ex);
    Grid->upml3Ex = NULL;
    free(Grid->upml4Ex);
    Grid->upml4Ex = NULL;
    free(Grid->upml5Ex);
    Grid->upml5Ex = NULL;
    free(Grid->upml6Ex);
    Grid->upml6Ex = NULL;
    free(Grid->upml1Hx);
    Grid->upml1Hx = NULL;
    free(Grid->upml2Hx);
    Grid->upml2Hx = NULL;
    free(Grid->upml3Hx);
    Grid->upml3Hx = NULL;
    free(Grid->upml4Hx);
    Grid->upml4Hx = NULL;
    free(Grid->upml5Hx);
    Grid->upml5Hx = NULL;
    free(Grid->upml6Hx);
    Grid->upml6Hx = NULL;
    
    free(Grid->upml1Ey);
    Grid->upml1Ey = NULL;
    free(Grid->upml2Ey);
    Grid->upml2Ey = NULL;
    free(Grid->upml3Ey);
    Grid->upml3Ey = NULL;
    free(Grid->upml4Ey);
    Grid->upml4Ey = NULL;
    free(Grid->upml5Ey);
    Grid->upml5Ey = NULL;
    free(Grid->upml6Ey);
    Grid->upml6Ey = NULL;
    free(Grid->upml1Hy);
    Grid->upml1Hy = NULL;
    free(Grid->upml2Hy);
    Grid->upml2Hy = NULL;
    free(Grid->upml3Hy);
    Grid->upml3Hy = NULL;
    free(Grid->upml4Hy);
    Grid->upml4Hy = NULL;
    free(Grid->upml5Hy);
    Grid->upml5Hy = NULL;
    free(Grid->upml6Hy);
    Grid->upml6Hy = NULL;
    
    free(Grid->upml1Ez);
    Grid->upml1Ez = NULL;
    free(Grid->upml2Ez);
    Grid->upml2Ez = NULL;
    free(Grid->upml3Ez);
    Grid->upml3Ez = NULL;
    free(Grid->upml4Ez);
    Grid->upml4Ez = NULL;
    free(Grid->upml5Ez);
    Grid->upml5Ez = NULL;
    free(Grid->upml6Ez);
    Grid->upml6Ez = NULL;
    free(Grid->upml1Hz);
    Grid->upml1Hz = NULL;
    free(Grid->upml2Hz);
    Grid->upml2Hz = NULL;
    free(Grid->upml3Hz);
    Grid->upml3Hz = NULL;
    free(Grid->upml4Hz);
    Grid->upml4Hz = NULL;
    free(Grid->upml5Hz);
    Grid->upml5Hz = NULL;
    free(Grid->upml6Hz);
    Grid->upml6Hz = NULL;
}

///@brief Releases all allocated memory for fields on box borders.
///@param Grid pointer to Grid struct
void freeFieldsOnBoxBorders(Grid *Grid){
    printf("releasing allocated memory for fields on box borders...\n");
    
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

///@brief method for releasing all previously allocated memory in struct Grid. Put all free() invokations in here
///@param Grid pointer to Grid struct
void freeMemoryOnGrid(Grid *Grid){
    freeFieldsOnGrid(Grid);
    freeFieldsOnBoxBorders(Grid);
    
    if (Grid->useUPML){
        freeUPMLCoefficients(Grid);
    }

}

///@brief maxwellPusher for E field. This method is only used for testing purposes. See pushEFieldInsideBoxes for the final version
///@remark E field is calculated via negative curl. Therefore value of B on the left side of E on the grid is required. Thus start i,j,k with 1
///@param Grid pointer to Grid struct
///@param dt time increment
void pushEFieldOnGrid(Grid *Grid, const double dt){
    double Hx_ijk;
    double Hy_ijk;
    double Hz_ijk;
    double Hz_ijm1k;
    double Hy_ijkm1;
    double Hx_ijkm1;
    double Hz_im1jk;
    double Hy_im1jk;
    double Hx_ijm1k;
    
    double cnx = 0.5 * dt / Grid->Resolution.dx;
    double cny = 0.5 * dt / Grid->Resolution.dy;
    double cnz = 0.5 * dt / Grid->Resolution.dz;
    
    
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
void externalPlaneWave(const double x[4], const double tStart, double Eextern[3], double Hextern[3]){
    double E0 = 0.01;
    double H0 = 0.01;
    double frequency = 2 * M_PI / 600;
    
    Eextern[0] = 0;
    if(x[0] - tStart >= x[1]){
        Eextern[1] = E0 * cos(frequency * (x[0] - tStart - x[1]));
    }
    else{
        Eextern[1] = 0;
    }
    Eextern[2] = 0;
    
    Hextern[0] = 0;
    Hextern[1] = 0;
    if(x[0] - tStart >= x[1]){
        Hextern[2] = H0 * cos(frequency * (x[0] - tStart - x[1]));
    }
    else{
        Hextern[2] = 0;
    }
    
    
}

///@brief sets external time dependent E and H field as plane wave.
///@param x four vector containing actual position where you want to evaluate the external field
///@param Eextern vector containing external E field components
///@param Hextern vector containing external H field components
void externalPulse(const double x[4], const double tStart, double Eextern[3], double Hextern[3], Grid *Grid){
    double E0 = 1.77 * pow(10, -8.0);
    double H0 = 1.77 * pow(10, -8.0);
    double frequency = 1.77 * pow(10, -8);
    
    double lengthOfSimulationBoxInX = Grid->lengthOfSimulationBoxInX;
    double lengthOfSimulationBoxInY = Grid->lengthOfSimulationBoxInY;
    double widthOfPulseInY = 0.5 * lengthOfSimulationBoxInY;
    double widthOfPulseInX = 0.5 * lengthOfSimulationBoxInX;
    
    Eextern[0] = 0;
    if(x[0] - tStart >= x[1] && x[1] > x[0] - tStart - widthOfPulseInX){
        Eextern[1] = E0 * cos(frequency * (x[0] - tStart - x[1])) * pow(sin(0.05 * frequency * (x[0] - tStart - x[1])),2) * exp(-pow(0.14 * frequency * (x[2] - widthOfPulseInY),2));
    }

    Eextern[2] = 0;
    
    Hextern[0] = 0;
    Hextern[1] = 0;
    if(x[0] - tStart >= x[1] && x[1] > x[0] - tStart - widthOfPulseInX){
        Hextern[2] = H0 * cos(frequency * (x[0] - tStart - x[1])) * pow(sin(0.05 * frequency * (x[0] - tStart - x[1])),2) * exp(-pow(0.14 * frequency * (x[2] - widthOfPulseInY),2));
    }
}

///@brief inits a sample E and B field onto the grid for testing purposes.
///@param Grid pointer to Grid struct
void initSamplePulseOnGrid(Grid *Grid){
    
    int ny = Grid->numberOfGridPointsInY;
    int nz = Grid->numberOfGridPointsInZ;
    int lengthInX = Grid->Resolution.dx * Grid->numberOfGridPointsInX;
    
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
void pushHFieldOnGrid(Grid *Grid, const double dt){
    double Ex_ijk;
    double Ey_ijk;
    double Ez_ijk;
    double Ey_ijkp1;
    double Ez_ijp1k;
    double Ez_ip1jk;
    double Ex_ijkp1;
    double Ex_ijp1k;
    double Ey_ip1jk;
    
    double cnx = 0.5 * dt / Grid->Resolution.dx;
    double cny = 0.5 * dt / Grid->Resolution.dy;
    double cnz = 0.5 * dt / Grid->Resolution.dz;
    
    
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

void clearFieldsFromGrid(Grid *Grid){
    printf("Clearing Fields from Grid ...\n");
    for(int i = 0; i < Grid->numberOfGridPointsInX * Grid->numberOfGridPointsInY * Grid->numberOfGridPointsInZ * 3; i++){
        Grid->E[i] = 0;
        Grid->H[i] = 0;
        Grid->D[i] = 0;
        Grid->B[i] = 0;
    }
}


