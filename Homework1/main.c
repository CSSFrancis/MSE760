//
//  main.c
//  MSE760
//
//  Created by Carter Francis on 9/6/18.
//  Copyright © 2018 Carter Francis. All rights reserved.
//
#include <math.h>
#include <stdio.h>

//functions

//Building a FFC Unit cell in reduced units
int createFCC (int numCells, double UnitCell, double X[], double Y[], double Z[])
{
    int count = 0;
    
    for (int i=0; i < numCells; i=i+1) {
        for (int j = 0; j < numCells;j=j+1){
            for (int k=0; k < numCells; k=k+1) {
                
                X[count] = i*1.0*UnitCell;
                printf("%g \n", X[count]);
                Y[count] = j*1.0*UnitCell;
                Z[count] = k*1.0*UnitCell;
                X[count+1] = i*1.0*UnitCell + (UnitCell/2.0);
                Y[count+1] = j*1.0*UnitCell + (UnitCell/2.0);
                Z[count+1] = k*1.0*UnitCell;
                X[count+2] = i*1.0*UnitCell + (UnitCell/2.0);
                Y[count+2] = j*1.0*UnitCell;
                Z[count+2] = k*1.0*UnitCell + (UnitCell/2.0);;
                X[count+3] = i*1.0*UnitCell;
                Y[count+3] = j*1.0*UnitCell + (UnitCell/2.0);
                Z[count+3] = k*1.0*UnitCell + (UnitCell/2.0);
                
                count =count +4;
            }
        }
    }
    return 0;
    }
//Calculating Potential

double calculate_Potential(double X[],double Y[], double Z[] , double energy[],int numAtoms, double L)
{
    //double potential;
    double dist=0;
    double energySum =0;
    //double cutoff = .2;
    for (int i=0; i <numAtoms-1; i=i+1){
        //determining atoms within some cutoff
        double energyHolder =0;
        for(int j=i+1; j <numAtoms; j=j+1){
            double difX =(X[i] - X[j]);
            double difY =(Y[i] - Y[j]);
            double difZ =(Z[i] - Z[j]);
            //For periodic boundry conditions
            /*
            if (difX > L){
                difX = difX+2.0*L;
            }
            if (difX <= -L){
                difX = difX+2.0*L;
            }
            if (difY > L){
                difY = 2.0*L+difY;
            }
            if (difY <= -L){
                difY = difY+2.0*L;
            } if (difZ > L){
                difZ = difZ+2.0*L;
            }
            if (difZ <= -L){
                difZ = difZ+2.0*L;
            }
             */
            dist = sqrt(pow(difX,2) + pow(difY,2) + pow(difZ,2));
            double redDist = dist/3.4;
            energyHolder = energyHolder + 4.0 * (pow((1/redDist),12) - pow((1/redDist),6));
            
                
        }
        energySum = energySum +energyHolder;
        energy[i] = energyHolder;
        
    }
    return energySum;
}

int main() {
    double energyArray[48];
    FILE *out_energy_file = fopen(getcwd(), "w");
    for (int a=0; a<15; a=a+1){
        double uCell = 5.26;// in angstroms
        int nCells = a+1;
        double epsilon = 0.0104;
        int nAtoms = nCells*nCells*nCells*4;
        double X[nAtoms];
        double Y[nAtoms];
        double Z[nAtoms];
        double sigma =3.4; // in angstroms
        double redUnitCell = uCell/sigma;
        double Length = uCell*nCells/2;
    
        createFCC(nCells, uCell, X, Y, Z);
        //getting the values from aPos
        // basis
        double en [nAtoms-1];
        double energy = calculate_Potential(X, Y, Z, en, nAtoms,Length);
        double energyPerAtom= (energy/nAtoms)*epsilon;
        char buf[0x100];
        snprintf(buf, sizeof(buf), "/Users/shaw/Shaw/MSE760/Out/%i.txt",a);
        FILE *out_file = fopen(buf, "w"); // write only
        for (int i=0; i<nAtoms;i++)
        {
            fprintf(out_file, "%g,%g,%g \n", X[i], Y[i], Z[i]);
            printf( "%g,%g,%g \n", X[i], Y[i], Z[i]);
        }
        energyArray[a] = energyPerAtom;
        fprintf(out_energy_file, "%i, %g \n",nCells, energyPerAtom);
    }
    
    
    return 0;
}
