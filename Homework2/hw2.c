
//
//  hw2.c
//  MSE760
//  This program is designed to simulate an Ar Crystal
//  Part 1: Rescale the volume by rescaling the lattice constant to 5.7
//  Part 2: Calculate the forces on the Argon atoms using a Lennard jones potential
//  Part 3: Assign random velocities to atoms so the temperature is 300K.
//  Part 4: Implement velocity Verlet algorithm
//
//  Created by Carter Francis on 10/10/18.
//  
//
//
//  main.c
//  MSE760
//
//  Created by Carter Francis on 9/6/18.
//  Copyright © 2018 Carter Francis. All rights reserved.
//
#include <math.h>
#include <stdio.h>
#include "hw2Extended.h"
#include <time.h>
#include <unistd.h>

int main() {
    double energyArray[48];

    int a=4;
    double sigma = 3.4;
    double uCell = 5.7/sigma;// in angstroms
    int nCells = a+1;
    double epsilon = 0.0104;
    int nAtoms = nCells*nCells*nCells*4;
    double kb=0.000086173;
    double X[nAtoms];
    double Y[nAtoms];
    double Z[nAtoms];
    double Xforce[nAtoms];
    double Yforce[nAtoms];
    double Zforce[nAtoms];
    double xVelocity[nAtoms];
    double yVelocity[nAtoms];
    double zVelocity[nAtoms];
    double kinE;
    double timestep =0.01;
    int timesteps =20/timestep;
    //Temp in reduced units
    //double initialTemp = 2.4857684495; // 300K
    double initialTemp = 4.14293; //500K
    //double initialTemp = 0; //0K
    double EStep[timesteps];
    double PEStep[timesteps];
    double KEStep[timesteps];
    double TempStep[timesteps];
    double T;

    double Length = uCell*nCells/2;
    createFCC(nCells, uCell, X, Y, Z);
    assignRandomVelocities(xVelocity, yVelocity, zVelocity, initialTemp,nAtoms);
    T = calculateTemperature(xVelocity, yVelocity, zVelocity, nAtoms);

    //getting the values from aPos
    // basis
    double en [nAtoms-1];
    calculate_Potential(X, Y, Z, en,Xforce,Yforce,Zforce, nAtoms,Length);
    printf("Temperature, TotalKE,AverageKE, step\n");
    for (int i=0; i<timesteps; i++)
    {
        double energy;
        static double PE;
        double KE;
        double T_during = calculateTemperature(xVelocity, yVelocity, zVelocity, nAtoms);
        printf("%g  ",T_during*epsilon/kb);
        TempStep[i]=T_during*epsilon/kb;
        energy = velocity_verlet(X, Y, Z, Xforce, Yforce, Zforce, xVelocity, yVelocity, zVelocity, en, nAtoms, timestep, Length, &PE,&KE);
        EStep[i] = (energy/nAtoms)*epsilon;
        printf(" %g ", (KE/nAtoms)*epsilon);
        KEStep[i]= (KE/nAtoms)*epsilon;
        PEStep[i]= (PE/nAtoms)*epsilon;
        printf("%i \n", i);


    }
    double T_f= calculateTemperature(xVelocity, yVelocity, zVelocity, nAtoms);
    FILE *out_energy_file = fopen("../Homework2/out/energies(,001).txt", "w");
    fprintf(out_energy_file, "Step, TotalE, KinE(ev/atom), PotE(ev/at), Temp(K)\n");
    for (int l=0; l<timesteps; l++)
    {
        double e =(2*KEStep[l])+PEStep[l];
        fprintf(out_energy_file, "%i, %g ,%g, %g, %g\n",l, e ,KEStep[l], PEStep[l], TempStep[l]);
    }
    FILE *out_file = fopen( "../Homework2/out/forces.txt", "w"); // write only
    fprintf(out_file, "The system temp is: %g \n", T);
    fprintf(out_file, "The final system temp is: %g \n", T_f);
    int i;
    fprintf(out_file, "Fx(N),            Fy(N),            Fz(N)\n");
    for (i=0; i<nAtoms;i++)
    {
        fprintf(out_file, "%g,    %g,    %g\n", ((Xforce[i]*epsilon*1.6*pow(10,-9))/sigma),((Yforce[i]*epsilon*1.6*pow(10,-9))/sigma),((Xforce[i]*epsilon*1.6*pow(10,-9))/sigma));//,xVelocity[i],yVelocity[i],zVelocity[i]);
    }
    int n_bins=1000;
    int pdf[1000-1]= {0};
    int total = pair_distribution_function(X,Y,Z,n_bins,nAtoms,Length,pdf);
    FILE *pairdf = fopen("../Homework2/out/pdf250K.txt", "w");
    fprintf(pairdf, "Distance, Counts \n");
    double bin_size=(Length/(n_bins-1));
    double volume = 4/3 * 3.1415* pow(Length,3);
    double num_density = total/volume;
    for (i=0; i<n_bins-1;i++)
    {
        double norm = 4/3 * 3.1415* (pow(((i+1) * bin_size),3)-pow((i * bin_size),3))*num_density;
        fprintf(pairdf, "%g, %g \n",(Length/n_bins)*i,(pdf[i]/norm));
    }
    return 1;
}