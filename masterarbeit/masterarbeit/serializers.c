//
//  serializers.c
//  masterarbeit
//
//  Created by David Symhoven on 01.04.17.
//  Copyright © 2017 David Symhoven. All rights reserved.
//

#include "serializers.h"
#include "string.h"
#include "math.h"
#include "calculations.h"

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
    FILE *PoyntingX = fopen("Sx.txt", "w");
    FILE *PoyntingY = fopen("Sy.txt", "w");;
    FILE *PoyntingZ = fopen("Sz.txt", "w");;
    
    double Sx, Sy, Sz;
    
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
                
                Sx = Ey * Hz - Ez * Hy;
                Sy = Ez * Hx - Ex * Hz;
                Sz = Ex * Hy - Ey * Hx;
                
                fprintf(PoyntingX, "%.18f\t", Sx);
                fprintf(PoyntingY, "%.18f\t", Sy);
                fprintf(PoyntingZ, "%.18f\t", Sz);
                
                if(plotE){
                    Esq = Ex * Ex + Ey * Ey + Ez * Ez;
                    
                    fprintf(fid, "%.18f\t", Esq);
                    if (Esq > Grid->EMax){
                        Grid->EMax = Esq;
                    }
                }
                if (plotB){
                    Bsq = Hx * Hx + Hy * Hy + Hz * Hz;
                    fprintf(fid2, "%.18f\t", Bsq);
                    if (Bsq > Grid->HMax){
                        Grid->HMax = Bsq;
                    }
                }
            }
            fprintf(PoyntingX, "\n");
            fprintf(PoyntingY, "\n");
            fprintf(PoyntingZ, "\n");
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
    fclose(PoyntingX);
    fclose(PoyntingY);
    fclose(PoyntingZ);
}

///@brief loops through the entire E and B array and writes Bx,By,Bz and (or) Ex,Ey,Ez to seperate files. File is structured similiar to the grid.
///@param filename pointer to a char. Gets modified inside the method
///@param index outer loop index. Is used to name the output file
///@param plotE set to true, if you want to write E-field to file
///@param plotB set to true, if you want to write B-field to file
///@throws ERROR: Could not open file for E or B field
void writeFieldComponentsForFourierAnalysisToFile(Grid *Grid, char *filename, int index, int planeForPlotting, bool plotE, bool plotB){
    printf("Writing field components for Fourier Analysis to file ...\n");
    FILE *fid_Ex = NULL;
    FILE *fid_Ey = NULL;
    FILE *fid_Ez = NULL;
    FILE *fid_Hx = NULL;
    FILE *fid_Hy = NULL;
    FILE *fid_Hz = NULL;
    
    if (plotE){
        sprintf(filename, "E_field_x%d", index);
        strcat(filename, ".txt");
        fid_Ex = fopen(filename,"w");
        if (fid_Ex == NULL){
            printf("ERROR: Could not open file E_field_x!");
        }
        sprintf(filename, "E_field_y%d", index);
        strcat(filename, ".txt");
        fid_Ey = fopen(filename,"w");
        if (fid_Ey == NULL){
            printf("ERROR: Could not open file E_field_y!");
        }
        sprintf(filename, "E_field_z%d", index);
        strcat(filename, ".txt");
        fid_Ez = fopen(filename,"w");
        if (fid_Ez == NULL){
            printf("ERROR: Could not open file E_field_z!");
        }
    }
    if(plotB){
        sprintf(filename, "B_field_x%d", index);
        strcat(filename, ".txt");
        fid_Hx = fopen(filename,"w");
        if (fid_Hx == NULL){
            printf("ERROR: Could not open file B_field_x!");
        }
        sprintf(filename, "B_field_y%d", index);
        strcat(filename, ".txt");
        fid_Hy = fopen(filename,"w");
        if (fid_Hy == NULL){
            printf("ERROR: Could not open file B_field_y!");
        }
        sprintf(filename, "B_field_z%d", index);
        strcat(filename, ".txt");
        fid_Hz = fopen(filename,"w");
        if (fid_Hz == NULL){
            printf("ERROR: Could not open file B_field_z!");
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
                    fprintf(fid_Ex, "%.18f\t", Ex);
                    fprintf(fid_Ey, "%.18f\t", Ey);
                    fprintf(fid_Ez, "%.18f\t", Ez);
                }
                if (plotB){
                    fprintf(fid_Hx, "%.18f\t", Hx);
                    fprintf(fid_Hy, "%.18f\t", Hy);
                    fprintf(fid_Hz, "%.18f\t", Hz);
                }
            }
            if(plotE){
                fprintf(fid_Ex,"\n");
                fprintf(fid_Ey,"\n");
                fprintf(fid_Ez,"\n");
            }
            if(plotB){
                fprintf(fid_Hx,"\n");
                fprintf(fid_Hy,"\n");
                fprintf(fid_Hz,"\n");
            }
        }
    }
    fclose(fid_Ex);
    fclose(fid_Ey);
    fclose(fid_Ez);
    fclose(fid_Hx);
    fclose(fid_Hy);
    fclose(fid_Hz);
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
        fprintf(fid,"%.18f\n", Grid->E[i]);
        fprintf(fid2,"%.18f\n", Grid->H[i]);
    }
    fclose(fid);
    fclose(fid2);
    
}

///@brief writes current position, velocity vectors and near field info to file. First line is four vector x. Second line is four vector u. Third line is xMin up to yMax of near field box
///@param filename pointer to a char. Gets modified inside the method
///@param index outer loop index. Is used to name the output file
///@param Particles pointer to Particle struct containing all particles
///@param numberOfParticles number of particles
void writeParticlesToFile(Particle *Particles, int numberOfParticles, char *filename, int index){
    for(int p = 0; p < numberOfParticles; p++){
        
        printf("Writing Particle%d to file ...\n", p);
        FILE *fid = NULL;
        sprintf(filename, "Particle%d_%d", p, index);
        strcat(filename, ".txt");
        fid = fopen(filename,"w");
        if (fid == NULL){
            printf("ERROR: Could not open file Particle!\n");
        }
        
        for (int i = 0; i < 4; i++){
            fprintf(fid, "%f\t", Particles[p].x[i]);
        }
        fprintf(fid, "\n");
        
        for (int i = 0; i < 4; i++){
            fprintf(fid, "%f\t", Particles[p].u[i]);
        }
        fprintf(fid, "\n");
        for (int i = 0; i < 4; i++){
            fprintf(fid, "%f\t", Particles[p].edgesOfNearFieldBox[i]);
        }
        fprintf(fid, "\n");
        fclose(fid);
    }
    
}

void writeForcesToFile(Forces *Forces, double t){
    FILE *fid2 = fopen("dampingTermVSLorentzForce.txt","a");
    
    double absDampingTerm = sqrt(Forces->damping[0]*Forces->damping[0] + Forces->damping[1]*Forces->damping[1] + Forces->damping[2]*Forces->damping[2]);
    double absLorentzForce = sqrt(Forces->lorentz[0]*Forces->lorentz[0] + Forces->lorentz[1]*Forces->lorentz[1] + Forces->lorentz[2]*Forces->lorentz[2]);
    
    fprintf(fid2, "%f %.18f %.18f %.18f\n", t, absDampingTerm, absLorentzForce, absDampingTerm / absLorentzForce);
    fclose(fid2);
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
        fprintf(fid, "%f %f %f %d %d %d %f %f %f %.18f %.18f\n", Grid->Resolution.dx, Grid->Resolution.dy, Grid->Resolution.dz, Grid->Box.numberOfGridPointsInX, Grid->Box.numberOfGridPointsInY, Grid->Box.numberOfGridPointsInZ, Grid->lengthOfSimulationBoxInX, Grid->lengthOfSimulationBoxInY, Grid->lengthOfSimulationBoxInZ, Grid->EMax, Grid->HMax);
    }
    fclose(fid);
}

///@brief loops through the specified grid plane and calcutes |B|^2 and |E|^2 of external fields. Values are stored in seperate files. File is structured similiar to the grid.
///@param filename pointer to a char. Gets modified inside the method
///@param index outer loop index. Is used to name the output file
///@param plotE set to true, if you want to write E-field to file
///@param plotB set to true, if you want to write B-field to file
///@throws ERROR: Could not open file for E or B field
void writeExternalFieldsToFile(Grid *Grid, double Eextern[3], double Bextern[3], const double t, const double tStart, char *filename, const int index, const int planeForPlotting, bool plotE, bool plotB){
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
        x[3] = k * Grid->Resolution.dz;
        
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                x[1] = i * Grid->Resolution.dx;
                x[2] = j * Grid->Resolution.dy;
                
                externalPulse(x, tStart, Eextern, Bextern, Grid);
//                externalPlaneWave(x, tStart, Eextern, Bextern);
                Ex = Eextern[0];
                Ey = Eextern[1];
                Ez = Eextern[2];
                
                
                Hx = Bextern[0];
                Hy = Bextern[1];
                Hz = Bextern[2];
                
                if(plotE){
                    Esq = Ex * Ex + Ey * Ey + Ez * Ez;
                    fprintf(fid, "%.18f\t", Ey);
                }
                if (plotB){
                    Hsq = Hx * Hx + Hy * Hy + Hz * Hz;
                    fprintf(fid2, "%.18f\t", Hsq);
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


///@brief writes numberOfParticles and startTime of simulation to file. This is needed for python scripts, if simulation does not start at t = 0.
///@param numberOfParticles number of particles
///@param startTime start time of simulation
void writeSimulationInfoToFile(int numberOfParticles, int startTime){
    printf("Writing simulationInfo to file ...\n");
    
    FILE *fid = fopen("simulationInfo.txt","w");
    fprintf(fid, "%d %d", numberOfParticles, startTime);
    fclose(fid);
}

///@brief In order to have a unique characterization of each simulation all parameters are written into a file called "initialConditions.txt".
///@param Particles struct containing all particles
///@param Grid instance of a Grid struct
///@param numberOfParticles number of particles
///@param t current simulation time
///@param tEnd time where simulation should end
///@param Eextern vector containing external E-Field components
///@param Bextern vector containing external B-Field components
void writeInitialConditionsToFile(Grid *Grid, Particle *Particles, int numberOfParticles, double t, double tEnd, double Eextern[3], double Bextern[3]){
    FILE *fid = fopen("initialConditions.txt","w");
    fprintf(fid, "%f %f %f %d %d %d %d %d %d %d %f ", Grid->Resolution.dx, Grid->Resolution.dy, Grid->Resolution.dz, Grid->Box.numberOfGridPointsInX, Grid->Box.numberOfGridPointsInY, Grid->Box.numberOfGridPointsInZ, Grid->numberOfBoxesInX, Grid->numberOfBoxesInY, Grid->numberOfBoxesInZ, numberOfParticles, t);
    for(int p = 0; p < numberOfParticles; p++){
        fprintf(fid, "%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f ", Particles[p].x[0], Particles[p].x[1], Particles[p].x[2], Particles[p].x[3], Particles[p].u[0], Particles[p].u[1], Particles[p].u[2], Particles[p].u[3]);
    }
    fprintf(fid, "%f %f %f %f %f %f", Eextern[0], Eextern[1], Eextern[2], Bextern[0], Bextern[1], Bextern[2]);
    fclose(fid);
}
