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
int createFCC (int numCells, double reducedUnitCell, double X[], double Y[], double Z[])
{
    int count = 0;
    
    for (int i=0; i < numCells; i=i+1) {
        for (int j = 0; j < numCells;j=j+1){
            for (int k=0; k < numCells; k=k+1) {
                X[count] = i*reducedUnitCell;
                Y[count] = j*reducedUnitCell;
                Z[count] = k*reducedUnitCell;
                X[count+1] = i*reducedUnitCell + (reducedUnitCell/2);
                Y[count+1] = j*reducedUnitCell + (reducedUnitCell/2);
                Z[count+1] = k*reducedUnitCell;
                count =count +2;
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
            //returns a very large number for the energy because the distance is lower than one
            double difX =(X[i] - X[j]);
            double difY =(Y[i] - Y[j]);
            double difZ =(Z[i] - Z[j]);
            //For periodic boundry conditions
            if (difX > L){
                difX = difX+2*L;
            }
            if (difX <= -L){
                difX = difX+2*L;
            }
            if (difY > L){
                difY = 2*L+difY;
            }
            if (difY <= -L){
                difY = difY+2*L;
            } if (difZ > L){
                difZ = difZ+2*L;
            }
            if (difZ <= -L){
                difZ = difZ+2*L;
            }
            dist = sqrt(pow(difX,2) + pow(difY,2) + pow(difZ,2));
            energyHolder = energyHolder + 4 * (pow((1/dist),12) - pow((1/dist),6));
            
                
        }
        energySum = energySum +energyHolder;
        energy[i] = energyHolder;
        
    }
    return energySum;
}

int main() {
    double energyArray[48];
    for (int a=0; a<50; a=a+1){
        double uCell = 5.25;// in angstroms
        int nCells = a+1;
        double epsilon = 0.0104;
        int nAtoms = nCells*nCells*nCells*2;
        double X[nAtoms];
        double Y[nAtoms];
        double Z[nAtoms];
        double sigma =3.4; // in angstroms
        double redUnitCell = uCell/sigma;
        double Length = redUnitCell*nCells/2;
    
        createFCC(nCells, redUnitCell, X, Y, Z);
        //getting the values from aPos
        // basis
        double en [nAtoms-1];
        double energy = calculate_Potential(X, Y, Z, en, nAtoms,Length);
        double energyPerAtom= (energy/nAtoms)*epsilon;
        
        energyArray[a] = energyPerAtom;
        
    }
    
    return 0;
}
