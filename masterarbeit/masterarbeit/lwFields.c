//
//  lwFields.c
//  masterarbeit
//
//  Created by David Symhoven on 16.12.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#include "lwFields.h"
#include "calculations.h"
#include "string.h"
#include "math.h"
#include "stdlib.h"

bool useUPML = true;

///@brief pushes E field inside Boxes. Sets H field at box borders. Adjusts H fields at box borders and finally pushes E field at box borders.
///@param Grid instance of Grid struct
///@param Particles struct containing all particles
///@param numberOfParticles number of particles
///@param t current simulation time
///@param dt time increment
void pushEField(Grid *Grid, Particle *Particles, int numberOfParticles, double t, double dt){
    pushEFieldInsideBoxes(Grid, dt);
    setHFieldOnBorders(Grid);
    adjustHFields(Grid, Particles, numberOfParticles, t);
    pushEFieldAtBorders(Grid, dt);
    
}
///@brief pushes H field inside Boxes. Sets E field at box borders. Adjusts E fields at box borders and finally pushes H field at box borders.
///@param Grid instance of Grid struct
///@param Particles struct containing all particles
///@param numberOfParticles number of particles
///@param t current simulation time
///@param dt time increment
void pushHField(Grid *Grid, Particle *Particles, int numberOfParticles, double t, double dt){
    pushHFieldInsideBoxes(Grid, dt);
    setEFieldOnBorders(Grid);
    adjustEFields(Grid, Particles, numberOfParticles, t);
    pushHFieldAtBorders(Grid, dt);
}

///@brief Maxwell Pusher for EField inside Boxes. This method is used for the hybrid field method.
///@remark E field is calculated via negative curl. Therefore value of H on the left side of E on the grid is required. Thus start i,j,k with 1. UPML is activated by default
///@param Grid instance of Grid struct
///@param dt time increment
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
///@param Grid instance of Grid struct
///@param dt time increment
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

///@brief this methods loops through all boxes and stores all values of H in the plane left, infront and below of the current box.  As can be seen from the Yee Scheme only Hy and Hz are needed from the left plane. Hz and Hx are needed from the plane infront and Hy and Hx are needed from the plane below. If the current box is a box where to the left does not exist (ib = 0) then the respective values for H are set to 0.
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
///@param Grid pointer to Grid struct
///@param Particles pointer to Particles struct array
///@param numberOfParticles number of particles
///@param t current simulation time
void adjustHFields(Grid *Grid, Particle *Particles, int numberOfParticles, const double t){
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    for (int ib = 0; ib < numberOfBoxesInX; ib++){
        for (int jb = 0; jb < numberOfBoxesInY; jb++){
            for (int kb = 0; kb < numberOfBoxesInZ; kb++){
                
                int boxIndex = ib * numberOfBoxesInY * numberOfBoxesInZ + jb * numberOfBoxesInZ + kb;
                
                if (ib != 0){
                    adjustHyz_im1(Grid, Particles, numberOfParticles, boxIndex, ib, jb, kb, t);
                }
                
                if (jb != 0){
                    adjustHxz_jm1(Grid, Particles, numberOfParticles, boxIndex, ib, jb, kb, t);
                }
                
                if (kb != 0){
                    adjustHxy_km1(Grid, Particles, numberOfParticles, boxIndex, ib, jb, kb, t);
                }
            }
        }
    }
}

///@brief this method adjusts the H values on the left and right side border of the near field.
///@param Grid pointer to Grid struct
///@param Particles pointer to Particles struct array
///@param numberOfParticles number of particles
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustHFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexIm1 is not. Then we are at the left border of the near field. Since we want to push a value for the Efield inside the near field we need (among others) the value to the left. Since this point is in the far field, the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the right border of the near field, i.e. current box is the not in near field but boxIndexIm1 is, then we push EField values from the far field. For this we need values to the left, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the E field value in the far field correctly.
void adjustHyz_im1(Grid *Grid, Particle *Particles, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int boxIndexIm1 = calcBoxIndexIm1(Grid, boxIndex);
    double xObserver[4];
    for (int p = 0; p < numberOfParticles; p++){
        if (boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndex)
            && !boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndexIm1)){
            
            for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX - 1);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    subLWField(Grid, &Particles[p], &Grid->Hy_im1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd], xObserver, 4);
                    subLWField(Grid, &Particles[p], &Grid->Hz_im1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd], xObserver, 5);
                }
            }
        }
        else if (!boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndex)
                 && boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndexIm1)){
            
            for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX - 1);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    addLWField(Grid, &Particles[p], &(Grid->Hy_im1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd]), xObserver, 4);
                    addLWField(Grid, &Particles[p], &(Grid->Hz_im1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd]), xObserver, 5);
                }
            }
            
        }
    }
    
}

///@brief this method adjusts the H values infront and on the back side border of the near field.
///@param Grid pointer to Grid struct
///@param Particles pointer to Particles struct array
///@param numberOfParticles number of particles
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustHFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexJm1 is not. Then we are at the facing border of the near field. Since we want to push a value for the Efield inside the near field we need (among others) the value infront. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the back border of the near field, i.e. current box is the not in near field but boxIndexJm1 is, then we push EField values from the far field. For this we need values infront, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the E field value in the far field correctly.
void adjustHxz_jm1(Grid *Grid, Particle *Particles, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int boxIndexJm1 = calcBoxIndexJm1(Grid, boxIndex);
    double xObserver[4];
    for (int p = 0; p < numberOfParticles; p++){
        if (boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndex)
            && !boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndexJm1)){
            
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY - 1);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    subLWField(Grid, &Particles[p], &Grid->Hx_jm1[boxIndex][numberOfGridPointsForBoxInZ * id + kd], xObserver, 3);
                    subLWField(Grid, &Particles[p], &Grid->Hz_jm1[boxIndex][numberOfGridPointsForBoxInZ * id + kd], xObserver, 5);
                }
            }
        }
        else if (!boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndex)
                 && boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndexJm1)){
            
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY - 1);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    addLWField(Grid, &Particles[p], &(Grid->Hx_jm1[boxIndex][numberOfGridPointsForBoxInZ * id + kd]), xObserver, 3);
                    addLWField(Grid, &Particles[p], &(Grid->Hz_jm1[boxIndex][numberOfGridPointsForBoxInZ * id + kd]), xObserver, 5);
                }
            }
            
        }
    }
    
}

///@brief this method adjusts the H values on the top and bottom side border of the near field.
///@param Grid pointer to Grid struct
///@param Particles pointer to Particles struct array
///@param numberOfParticles number of particles
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustHFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexKm1 is not. Then we are at bottom border of the near field. Since we want to push a value for the Efield inside the near field we need (among others) the value below. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the top border of the near field, i.e. current box is the not in near field but boxIndexKm1 is, then we push EField values from the far field. For this we need values below, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the E field value in the far field correctly.
void adjustHxy_km1(Grid *Grid, Particle *Particles, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int boxIndexKm1 = calcBoxIndexKm1(Grid, boxIndex);
    double xObserver[4];
    
    for (int p = 0; p < numberOfParticles; p++){
        if (boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndex)
            && !boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndexKm1)){
            
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ - 1);
                    
                    subLWField(Grid, &Particles[p], &Grid->Hx_km1[boxIndex][numberOfGridPointsForBoxInY * id + jd], xObserver, 3);
                    subLWField(Grid, &Particles[p], &Grid->Hy_km1[boxIndex][numberOfGridPointsForBoxInY * id + jd], xObserver, 4);
                }
            }
        }
        else if (!boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndex)
                 && boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndexKm1)){
            
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ - 1);
                    
                    addLWField(Grid, &Particles[p], &(Grid->Hx_km1[boxIndex][numberOfGridPointsForBoxInY * id + jd]), xObserver, 3);
                    addLWField(Grid, &Particles[p], &(Grid->Hy_km1[boxIndex][numberOfGridPointsForBoxInY * id + jd]), xObserver, 4);
                }
            }
            
        }
    }
    
}

///@brief this method loops through all boxes and adjusts the E fields in the plane to the left, infront and below. The actual adjustemnt takes place in "adjustEyz_ip1()", "adjustExz_jp1()" and "adjustExy_kp1()" method.
///@param Grid pointer to Grid struct
///@param Particles pointer to Particles struct array
///@param numberOfParticles number of particles
///@param t current simulation time
void adjustEFields(Grid *Grid, Particle *Particles, int numberOfParticles, const double t){
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    for (int ib = 0; ib < numberOfBoxesInX; ib++){
        for (int jb = 0; jb < numberOfBoxesInY; jb++){
            for (int kb = 0; kb < numberOfBoxesInZ; kb++){
                
                int boxIndex = ib * numberOfBoxesInY * numberOfBoxesInZ + jb * numberOfBoxesInZ + kb;
                
                if (ib != numberOfBoxesInX - 1){
                    adjustEyz_ip1(Grid, Particles, numberOfParticles, boxIndex, ib, jb, kb, t);
                }
                
                if (jb != numberOfBoxesInY - 1){
                    adjustExz_jp1(Grid, Particles, numberOfParticles, boxIndex, ib, jb, kb, t);
                }
                
                if (kb != numberOfBoxesInZ - 1){
                    adjustExy_kp1(Grid, Particles, numberOfParticles, boxIndex, ib, jb, kb, t);
                }
            }
        }
    }
}


///@brief this method adjusts the E values on the left and right side border of the near field.
///@param Grid pointer to Grid struct
///@param Particles pointer to Particles struct array
///@param numberOfParticles number of particles
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustEFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexIp1 is not. Then we are at the right border of the near field. Since we want to push a value for the Hfield inside the near field we need (among others) the E field value to the right. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the left border of the near field, i.e. current box is the not in near field but boxIndexIp1 is, then we push HField values from the far field. For this we need values to the rigth, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the H field value in the far field correctly.
void adjustEyz_ip1(Grid *Grid, Particle *Particles, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int boxIndexIp1 = calcBoxIndexIp1(Grid, boxIndex);
    double xObserver[4];
    for (int p = 0; p < numberOfParticles; p++){
        if (boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndex)
            && !boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndexIp1)){
            
            for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * ((ib + 1) * numberOfGridPointsForBoxInX);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    subLWField(Grid, &Particles[p], &Grid->Ey_ip1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd], xObserver, 1);
                    subLWField(Grid, &Particles[p], &Grid->Ez_ip1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd], xObserver, 2);
                }
            }
        }
        else if (!boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndex)
                 && boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndexIp1)){
            
            for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * ((ib + 1) * numberOfGridPointsForBoxInX);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    addLWField(Grid, &Particles[p], &(Grid->Ey_ip1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd]), xObserver, 1);
                    addLWField(Grid, &Particles[p], &(Grid->Ez_ip1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd]), xObserver, 2);
                }
            }
            
        }
    }
    
}

///@brief this method adjusts the E values infront and on the back side border of the near field.
///@param Grid pointer to Grid struct
///@param Particles pointer to Particles struct array
///@param numberOfParticles number of particles
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustEFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexJp1 is not. Then we are at the back side border of the near field. Since we want to push a value for the Hfield inside the near field we need (among others) the E field value behind. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the front side border of the near field, i.e. current box is the not in near field but boxIndexJp1 is, then we push HField values from the far field. For this we need values infront, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the H field value in the far field correctly.
void adjustExz_jp1(Grid *Grid, Particle *Particles, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int boxIndexJp1 = calcBoxIndexJp1(Grid, boxIndex);
    double xObserver[4];
    for (int p = 0; p < numberOfParticles; p++){
        if (boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndex)
            && !boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndexJp1)){
            
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * ((jb + 1) * numberOfGridPointsForBoxInY);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    subLWField(Grid, &Particles[p], &Grid->Ex_jp1[boxIndex][numberOfGridPointsForBoxInZ * id + kd], xObserver, 0);
                    subLWField(Grid, &Particles[p], &Grid->Ez_jp1[boxIndex][numberOfGridPointsForBoxInZ * id + kd], xObserver, 2);
                }
            }
        }
        else if (!boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndex)
                 && boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndexJp1)){
            
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * ((jb + 1) * numberOfGridPointsForBoxInY);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    addLWField(Grid, &Particles[p], &(Grid->Ex_jp1[boxIndex][numberOfGridPointsForBoxInZ * id + kd]), xObserver, 0);
                    addLWField(Grid, &Particles[p], &(Grid->Ez_jp1[boxIndex][numberOfGridPointsForBoxInZ * id + kd]), xObserver, 2);
                }
            }
            
        }
    }
    
}

///@brief this method adjusts the E values at top and bottom side border of the near field.
///@param Grid pointer to Grid struct
///@param Particles pointer to Particles struct array
///@param numberOfParticles number of particles
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustEFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexKp1 is not. Then we are at the top side border of the near field. Since we want to push a value for the Hfield inside the near field we need (among others) the E field value above. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the bottom side border of the near field, i.e. current box is the not in near field but boxIndexKp1 is, then we push HField values from the far field. For this we need values above, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the H field value in the far field correctly.
void adjustExy_kp1(Grid *Grid, Particle *Particles, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t){
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int boxIndexKp1 = calcBoxIndexKp1(Grid, boxIndex);
    double xObserver[4];
    for(int p = 0; p < numberOfParticles; p++){
        if (boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndex)
            && !boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndexKp1)){
            
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * ((kb + 1) * numberOfGridPointsForBoxInZ);
                    
                    subLWField(Grid, &Particles[p], &Grid->Ex_kp1[boxIndex][numberOfGridPointsForBoxInY * id + jd], xObserver, 0);
                    subLWField(Grid, &Particles[p], &Grid->Ey_kp1[boxIndex][numberOfGridPointsForBoxInY * id + jd], xObserver, 1);
                }
            }
        }
        else if (!boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndex)
                 && boxIsInNearFieldOfParticle(Grid, &Particles[p], boxIndexKp1)){
            
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * ((kb + 1) * numberOfGridPointsForBoxInZ);
                    
                    addLWField(Grid, &Particles[p], &(Grid->Ex_kp1[boxIndex][numberOfGridPointsForBoxInY * id + jd]), xObserver, 0);
                    addLWField(Grid, &Particles[p], &(Grid->Ey_kp1[boxIndex][numberOfGridPointsForBoxInY * id + jd]), xObserver, 1);
                }
            }
            
        }
    }
    
}

///@brief Maxwell Pusher for EField at box borders. This method is used for the hybrid field method.
///@remark E field is calculated via negative curl. Therefore value of H on the left side of E on the grid is required. Thus start i,j,k with 1. UPML is activated by default
///@param Grid pointer to Grid struct
///@param dt time increment
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
///@param Grid pointer to Grid struct
///@param dt time increment
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
        double Ex, Ey, Ez, Hx, Hy, Hz, Esq, Bsq;
        
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
                    Esq = Ex * Ex + Ey * Ey + Ez * Ez;
                    fprintf(fid, "%f\t", Esq);
                    if (Esq > Grid->EMax){
                        Grid->EMax = Esq;
                    }
                }
                if (plotB){
                    Bsq = Hx * Hx + Hy * Hy + Hz * Hz;
                    fprintf(fid2, "%f\t", Bsq);
                    if (Bsq > Grid->HMax){
                        Grid->HMax = Bsq;
                    }
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



///@brief loops through the entire grid and writes respective E and B values to E_initalFields.txt and H_initialFields.txt respectively.
///@param Grid instance of Grid struct
void writeFieldsFromCompleteGridToFile(Grid *Grid){
    printf("Writing initial fields to file ...\n");
    FILE *fid = fopen("E_initialField.txt", "w");
    if(fid == NULL){
        printf("ERROR: Could not open E_initialField.txt!\n");
    }
    FILE *fid2 = fopen("H_initialField.txt","w");
    if(fid2 == NULL){
        printf("ERROR: Could not open H_initialField.txt!\n");
    }

    int nx = Grid->numberOfGridPointsInX;
    int ny = Grid->numberOfGridPointsInY;
    int nz = Grid->numberOfGridPointsInZ;

    for(int i = 0; i < nx * ny * nz * 3; i++){
        fprintf(fid,"%f\n", Grid->E[i]);
        fprintf(fid2,"%f\n", Grid->H[i]);
    }
    fclose(fid);
    fclose(fid2);
    
}

///@brief adds LW fields in box specified by input parameter "boxIndex"
///@param Grid pointer to an instance of a Grid struct.
///@param Particle pointer to an instance of a Particle struct
///@param boxIndex specified box in which fields shall be calcualted
///@param t time at which the fields shall be calcualted. This is also the zeroth component of xObserver
///@remark translate box index into number of box in x,y and z direction first (ib, jb, kb). Then calculate the gridIndex of the lower left corner of the current box and start looping through all grid points in that box. Each gridpoint is the current observation point, where fields shall be calcualted. Also grid index is updated each time the current grid point changes in the loop. Then add Lw field at that point.
void addLWFieldsInBox(Grid *Grid, Particle *Particle, int boxIndex, double t){
    printf("adding LW fields in box %d\n", boxIndex);
    
    
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    double dx = Grid->dx;
    double dy = Grid->dy;
    double dz = Grid->dz;
    
    int adjustmentDueToUpmlInXRight = 0;
    int adjustmentDueToUpmlInXLeft = 0;
    int adjustmentDueToUpmlInYBack = 0;
    int adjustmentDueToUpmlInYFront = 0;
    int adjustmentDueToUpmlInZTop = 0;
    int adjustmentDueToUpmlInZBottom = 0;
    
    int ib = boxIndex / (numberOfBoxesInY * numberOfBoxesInZ);
    int jb = (boxIndex - ib * (numberOfBoxesInY * numberOfBoxesInZ)) / numberOfBoxesInY;
    int kb = (boxIndex - ib * (numberOfBoxesInZ * numberOfBoxesInZ)) - jb * numberOfBoxesInZ;
    
    if (useUPML){
        if(ib == 0) {
            adjustmentDueToUpmlInXLeft = Grid->upmlLayerWidth - 1;
        }
        else if(ib == Grid->numberOfBoxesInX - 1){
            adjustmentDueToUpmlInXRight = Grid->upmlLayerWidth - 1;
        }
        else{
            adjustmentDueToUpmlInXRight = 0;
            adjustmentDueToUpmlInXLeft = 0;
        }
        if(jb == 0){
            adjustmentDueToUpmlInYFront = Grid->upmlLayerWidth - 1;
        }
        else if(jb == Grid->numberOfBoxesInY - 1){
            adjustmentDueToUpmlInYBack = Grid->upmlLayerWidth - 1;
        }
        else{
            adjustmentDueToUpmlInYBack = 0;
            adjustmentDueToUpmlInYFront = 0;
        }
        if(kb == 0){
            adjustmentDueToUpmlInZBottom = Grid->upmlLayerWidth - 1;
        }
        else if(kb == Grid->numberOfBoxesInZ - 1){
            adjustmentDueToUpmlInZTop = Grid->upmlLayerWidth - 1;
        }
        else{
            adjustmentDueToUpmlInZBottom = 0;
            adjustmentDueToUpmlInZTop = 0;
        }
    }
    
    double xObserver[4] = {0};
    xObserver[0] = t;
    
    int lowerLeftGridIndexInBox = ib * numberOfGridPointsForBoxInX * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + jb * numberOfGridPointsForBoxInY * numberOfGridPointsInZ * 3 + kb * numberOfGridPointsForBoxInZ * 3;
    
    for (int id = 0 + adjustmentDueToUpmlInXLeft; id < numberOfGridPointsForBoxInX - adjustmentDueToUpmlInXRight; id++){
        for (int jd = 0 + adjustmentDueToUpmlInYFront; jd < numberOfGridPointsForBoxInY - adjustmentDueToUpmlInYBack; jd++){
            for (int kd = 0 + adjustmentDueToUpmlInZBottom; kd < numberOfGridPointsForBoxInZ - adjustmentDueToUpmlInZTop; kd++){
                
                int gridIndexInBox = lowerLeftGridIndexInBox + 3 * kd + 3 * jd * numberOfGridPointsInZ + 3 * id * numberOfGridPointsInZ * numberOfGridPointsInY;
                
                xObserver[1] = (ib * numberOfGridPointsForBoxInX + id) * dx;
                xObserver[2] = (jb * numberOfGridPointsForBoxInY + jd) * dy;
                xObserver[3] = (kb * numberOfGridPointsForBoxInZ + kd) * dz;
                
                addLWField(Grid, Particle, &Grid->H[gridIndexInBox], xObserver, 3);
                addLWField(Grid, Particle, &Grid->H[gridIndexInBox + 1], xObserver, 4);
                addLWField(Grid, Particle, &Grid->H[gridIndexInBox + 2], xObserver, 5);
                
                addLWField(Grid, Particle, &Grid->E[gridIndexInBox], xObserver, 0);
                addLWField(Grid, Particle, &Grid->E[gridIndexInBox + 1], xObserver, 1);
                addLWField(Grid, Particle, &Grid->E[gridIndexInBox + 2], xObserver, 2);
                
            }
        }
    }
}

///@brief subtracts LW fields in box specified by input parameter "boxIndex"
///@param Grid pointer to an instance of a Grid struct.
///@param Particle pointer to an instance of a Particle struct
///@param boxIndex specified box in which fields shall be calcualted
///@param t time at which the fields shall be calcualted. This is also the zeroth component of xObserver
///@remark translate box index into number of box in x,y and z direction first (ib, jb, kb). Then calculate the gridIndex of the lower left corner of the current box and start looping through all grid points in that box. Each gridpoint is the current observation point, where fields shall be calcualted. Also grid index is updated each time the current grid point changes in the loop. Then sub Lw field at that point.
void subLWFieldsInBox(Grid *Grid, Particle *Particle, int boxIndex, double t){
    printf("subtracting LW fields in box %d\n", boxIndex);
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsForBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsForBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsForBoxInZ;
    
    int numberOfGridPointsInY = Grid->numberOfGridPointsInY;
    int numberOfGridPointsInZ = Grid->numberOfGridPointsInZ;
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    double dx = Grid->dx;
    double dy = Grid->dy;
    double dz = Grid->dz;
    
    int adjustmentDueToUpmlInXRight = 0;
    int adjustmentDueToUpmlInXLeft = 0;
    int adjustmentDueToUpmlInYBack = 0;
    int adjustmentDueToUpmlInYFront = 0;
    int adjustmentDueToUpmlInZTop = 0;
    int adjustmentDueToUpmlInZBottom = 0;
    
    int ib = boxIndex / (numberOfBoxesInY * numberOfBoxesInZ);
    int jb = (boxIndex - ib * (numberOfBoxesInY * numberOfBoxesInZ)) / numberOfBoxesInY;
    int kb = (boxIndex - ib * (numberOfBoxesInZ * numberOfBoxesInZ)) - jb * numberOfBoxesInZ;
    
    if (useUPML){
        if(ib == 0) {
            adjustmentDueToUpmlInXLeft = Grid->upmlLayerWidth - 1;
        }
        else if(ib == Grid->numberOfBoxesInX){
            adjustmentDueToUpmlInXRight = Grid->upmlLayerWidth - 1;
        }
        else{
            adjustmentDueToUpmlInXRight = 0;
            adjustmentDueToUpmlInXLeft = 0;
        }
        if(jb == 0){
            adjustmentDueToUpmlInYFront = Grid->upmlLayerWidth - 1;
        }
        else if(jb == Grid->numberOfBoxesInY){
            adjustmentDueToUpmlInYBack = Grid->upmlLayerWidth - 1;
        }
        else{
            adjustmentDueToUpmlInYBack = 0;
            adjustmentDueToUpmlInYFront = 0;
        }
        if(kb == 0){
            adjustmentDueToUpmlInZBottom = Grid->upmlLayerWidth - 1;
        }
        else if(kb == Grid->numberOfBoxesInZ){
            adjustmentDueToUpmlInZTop = Grid->upmlLayerWidth - 1;
        }
        else{
            adjustmentDueToUpmlInZBottom = 0;
            adjustmentDueToUpmlInZTop = 0;
        }
    }
    
    double xObserver[4] = {0};
    xObserver[0] = t;
    
    int lowerLeftGridIndexInBox = ib * numberOfGridPointsForBoxInX * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + jb * numberOfGridPointsForBoxInY * numberOfGridPointsInZ * 3 + kb * numberOfGridPointsForBoxInZ * 3;
    
    for (int id = 0 + adjustmentDueToUpmlInXLeft; id < numberOfGridPointsForBoxInX - adjustmentDueToUpmlInXRight; id++){
        for (int jd = 0 + adjustmentDueToUpmlInYFront; jd < numberOfGridPointsForBoxInY - adjustmentDueToUpmlInYBack; jd++){
            for (int kd = 0 + adjustmentDueToUpmlInZBottom; kd < numberOfGridPointsForBoxInZ - adjustmentDueToUpmlInZTop; kd++){
                
                int gridIndexInBox = lowerLeftGridIndexInBox + 3 * kd + 3 * jd * numberOfGridPointsInZ + 3 * id * numberOfGridPointsInZ * numberOfGridPointsInY;
                
                xObserver[1] = (ib * numberOfGridPointsForBoxInX + id) * dx;
                xObserver[2] = (jb * numberOfGridPointsForBoxInY + jd) * dy;
                xObserver[3] = (kb * numberOfGridPointsForBoxInZ + kd) * dz;
                
                subLWField(Grid, Particle, &Grid->H[gridIndexInBox], xObserver, 3);
                subLWField(Grid, Particle, &Grid->H[gridIndexInBox + 1], xObserver, 4);
                subLWField(Grid, Particle, &Grid->H[gridIndexInBox + 2], xObserver, 5);
                
                subLWField(Grid, Particle, &Grid->E[gridIndexInBox], xObserver, 0);
                subLWField(Grid, Particle, &Grid->E[gridIndexInBox + 1], xObserver, 1);
                subLWField(Grid, Particle, &Grid->E[gridIndexInBox + 2], xObserver, 2);
                
            }
        }
    }
}

///@brief adds LW fields at grid point with index "gridIndexInBox". Each point on the grid is considered to be the observation point, where the zeroth component is the current simulation time and also the time at which we want to calculate the fields. Therefore loop through the particle xhistory vector up to the current simulation time and search for the pair of positions which fulfil the condition that the old position is inside and the new position is outside the backward lightcone of the observation point. Once that pair of positions is found the rest of the xHistory vector can be skipped because all following points will be outside as well.
///@param Grid pointer to an instance of a Grid struct.
///@param xObserver observation point at which the LW fields ahall be calculted
///@param component component in which the observation point shall be shifted.
///@param Particle pointer to an instance of a Particle struct
///@param destination pointer to where LW fields shall ne added.
///@remark the Yee-scheme requires that the E and B fields are shifted against each other. Maybe the switch case is not neccessary. But I didn't tested it yet. For the staggering a xObservationCopy vector is introduced, which is used for field calculations, so that the original xObservation vector is left untouched. After staggering took place the currentHistoryLength property of Particle struct is used as upper loop index. By doing this, this method (and all above the chain) can be used either within the outer simulation loop to calculate the fileds every time step (video) or afterwards to caluclate the fields just ones.
void addLWField(Grid *Grid, Particle *Particle, double *destination, double xObserver[4], int component){
    double xObserverCopy[4];
    memcpy(xObserverCopy, xObserver, 4 * sizeof(double));
    switch (component)
    {
        case 0:
            xObserverCopy[1] += 0.5 * (Grid->dx);
            break;
        case 1:
            xObserverCopy[2] += 0.5 * (Grid->dy);
            break;
        case 2:
            xObserverCopy[3] += 0.5 * (Grid->dz);
            break;
        case 3:
            xObserverCopy[2] += 0.5 * (Grid->dy);
            xObserverCopy[3] += 0.5 * (Grid->dz);
            break;
        case 4:
            xObserverCopy[3] += 0.5 * (Grid->dz);
            xObserverCopy[1] += 0.5 * (Grid->dx);
            break;
        case 5:
            xObserverCopy[1] += 0.5 * (Grid->dx);
            xObserverCopy[2] += 0.5 * (Grid->dy);
            break;
    }
    
    
    double beta[3] = {0};
    double intersectionPoint[4] = {0};
    double velocityAtIntersectionPoint[4] = {0};
    double gamma_sq;
    double R_sq;
    double R;
    double n[3] = {0};
    double betaDot[3] = {0};
    double dt = 0.5 * Grid->dx;
    double E[3] = {0};
    double B[3] = {0};
    int currentHistoryLength = Particle->currentHistoryLength;
    
    for (int index = 0; index < currentHistoryLength - 1; index ++){
        if(isInsideBackwardLightcone(Particle->xHistory[index], xObserverCopy) && !isInsideBackwardLightcone(Particle->xHistory[index+1], xObserverCopy)){
            calculateIntersectionPoint(Particle->xHistory[index], Particle->xHistory[index+1], Particle->uHistory[index], Particle->uHistory[index+1], xObserverCopy, intersectionPoint, velocityAtIntersectionPoint);
//            calculateBeta(Particle->xHistory[index], Particle->xHistory[index+1], beta);
            calculateLienardWiechertParameters(intersectionPoint, xObserverCopy, velocityAtIntersectionPoint, &gamma_sq, &R_sq, &R, n, beta);
            calculateBetaDot(Particle->uHistory[index], Particle->uHistory[index+1], dt, betaDot);
            calcuateLienardWiechertFields(gamma_sq, R_sq, R, n, beta, betaDot, Particle->charge, E, B);
            break;
        }
    }
    
    switch (component)
    {
        case 0:
            *destination += E[0];
            break;
        case 1:
            *destination += E[1];
            break;
        case 2:
            *destination += E[2];
            break;
        case 3:
            *destination += B[0];
            break;
        case 4:
            *destination += B[1];
            break;
        case 5:
            *destination += B[2];
            break;
    }
    
}

///@brief adds LW fields at grid point with index "gridIndexInBox". Each point on the grid is considered to be the observation point, where the zeroth component is the current simulation time and also the time at which we want to calculate the fields. Therefore loop through the particle xhistory vector up to the current simulation time and search for the pair of positions which fulfil the condition that the old position is inside and the new position is outside the backward lightcone of the observation point. Once that pair of positions is found the rest of the xHistory vector can be skipped because all following points will be outside as well.
///@param Grid pointer to an instance of a Grid struct.
///@param xObserver observation point at which the LW fields ahall be calculted
///@param component component in which the observation point shall be shifted.
///@param Particle pointer to an instance of a Particle struct
///@param destination pointer to where LW fields shall be subtracted.
///@remark the Yee-scheme requires that the E and B fields are shifted against each other. Maybe the switch case is not neccessary. But I didn't tested it yet. For the staggering a xObservationCopy vector is introduced, which is used for field calculations, so that the original xObservation vector is left untouched. After staggering took place the currentHistoryLength property of Particle struct is used as upper loop index. By doing this, this method (and all above the chain) can be used either within the outer simulation loop to calculate the fileds every time step (video) or afterwards to caluclate the fields just once.
void subLWField(Grid *Grid, Particle *Particle, double *destination, double xObserver[4], int component){
    double xObserverCopy[4];
    memcpy(xObserverCopy, xObserver, 4 * sizeof(double));
    switch (component)
    {
        case 0:
            xObserverCopy[1] += 0.5 * (Grid->dx);
            break;
        case 1:
            xObserverCopy[2] += 0.5 * (Grid->dy);
            break;
        case 2:
            xObserverCopy[3] += 0.5 * (Grid->dz);
            break;
        case 3:
            xObserverCopy[2] += 0.5 * (Grid->dy);
            xObserverCopy[3] += 0.5 * (Grid->dz);
            break;
        case 4:
            xObserverCopy[3] += 0.5 * (Grid->dz);
            xObserverCopy[1] += 0.5 * (Grid->dx);
            break;
        case 5:
            xObserverCopy[1] += 0.5 * (Grid->dx);
            xObserverCopy[2] += 0.5 * (Grid->dy);
            break;
    }
    
    
    double beta[3] = {0};
    double intersectionPoint[4] = {0};
    double velocityAtIntersectionPoint[4] = {0};
    double gamma_sq;
    double R_sq;
    double R;
    double n[3] = {0};
    double betaDot[3] = {0};
    double dt = 0.5 * Grid->dx;
    double E[3] = {0};
    double B[3] = {0};
    int currentHistoryLength = Particle->currentHistoryLength;
    
    
    for (int index = 0; index < currentHistoryLength - 1; index ++){
        if(isInsideBackwardLightcone(Particle->xHistory[index], xObserverCopy) && !isInsideBackwardLightcone(Particle->xHistory[index+1], xObserverCopy)){
            calculateIntersectionPoint(Particle->xHistory[index], Particle->xHistory[index+1], Particle->uHistory[index], Particle->uHistory[index+1], xObserverCopy, intersectionPoint, velocityAtIntersectionPoint);
            calculateBeta(Particle->xHistory[index], Particle->xHistory[index+1], beta);
            calculateLienardWiechertParameters(intersectionPoint, xObserverCopy, velocityAtIntersectionPoint, &gamma_sq, &R_sq, &R, n, beta);
            calculateBetaDot(Particle->uHistory[index], Particle->uHistory[index+1], dt, betaDot);
            calcuateLienardWiechertFields(gamma_sq, R_sq, R, n, beta, betaDot, Particle->charge, E, B);
            break;
        }
    }
    
    switch (component)
    {
        case 0:
            *destination -= E[0];
            break;
        case 1:
            *destination -= E[1];
            break;
        case 2:
            *destination -= E[2];
            break;
        case 3:
            *destination -= B[0];
            break;
        case 4:
            *destination -= B[1];
            break;
        case 5:
            *destination -= B[2];
            break;
    }
    
}

///@brief in order to push the current particle properly not just the external fields but also the field contributions from other particles propagated into the near field region need to be taken into account.
///@param Particle instance of Particle struct
///@param Grid instance of Grid struct
///@param Eextern vector containing external E field components
///@param Bextern vector containing external B field components
///@param E vector containg both external fields and those from other particles propagated into the near field region of the current particle
///@param B vector containg both external fields and those from other particles propagated into the near field region of the current particle
void updateFieldsForParticlePush(Particle *Particle, Grid *Grid, double Eextern[3], double Bextern[3], double E[3], double B[3]){
    
    interpolateFields(Grid, Particle, E, B);
    for (int i = 0; i < 3; i++){
        E[i] = Eextern[i] + E[i];
        B[i] = Bextern[i] + B[i];
    }
}




