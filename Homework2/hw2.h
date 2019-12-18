//
//  homework2.h
//  MSE760
//
//  Created by Carter Francis on 11/1/18.
//  Copyright © 2018 Carter Francis. All rights reserved.
//

#ifndef hw2_h
#define hw2_h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//Building a FFC Unit cell in reduced units
int createFCC (int numCells, double UnitCell, double X[], double Y[], double Z[])
{
    int count = 0;

    for (int i=0; i < numCells; i=i+1) {
        for (int j = 0; j < numCells;j=j+1){
            for (int k=0; k < numCells; k=k+1) {
                X[count] = i*1.0*UnitCell;
                //printf("%g \n", X[count]);
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

//Calculating the energy and the potential on some atom and moving that atom

double calculate_Potential(double X[],double Y[], double Z[], double energy[], int numAtoms, double L)
{
    //double potential;
    double dist;
    double energySum =0.0;
    //double cutoff = .2;
    for (int i=0; i <numAtoms-1; i=i+1){
        //determining atoms within some cutoff
        double energyHolder =0.0;
        for(int j=i+1; j <numAtoms; j=j+1){
            double difX =(X[i] - X[j]);
            double difY =(Y[i] - Y[j]);
            double difZ =(Z[i] - Z[j]);
            //For periodic boundry conditions
            if (difX > L){difX = difX-2.0*L;}
            if (difX <= -L){difX = difX+2.0*L;}
            if (difY > L){difY = difY-2.0*L;}
            if (difY <= -L){difY = difY+2.0*L;}
            if (difZ > L){difZ = difZ-2.0*L;}
            if (difZ <= -L){difZ = difZ+2.0*L;}
            //calculating the distance from atom i to j
            dist = sqrt(pow(difX,2.0) + pow(difY,2.0) + pow(difZ,2.0));
            //printf("%g \n", dist);
            // calculating the forces on the particles
            energyHolder = energyHolder + 4.0 * (pow((1/dist),12) - pow((1/dist),6));
        }
        energySum = energySum +energyHolder;
        energy[i] = energyHolder;

    }
    return energySum;
}

double calculateAtomEnergy(double oX[],double oY[], double oZ[],int rAtom, double X,double Y,double Z, int natom, double Len)
{
    double e = 0.0;
    double d;
    for (int i=0; (i<natom && i !=rAtom); i++)
    {
        double dX =(oX[i] - X);
        double dY =(oY[i] - Y);
        double dZ =(oZ[i] - Z);
        //For periodic boundry conditions
        if (dX > Len){dX = dX-2.0*Len;}
        if (dX <= -Len){dX = dX+2.0*Len;}
        if (dY > Len){dY = dY-2.0*Len;}
        if (dY <= -Len){dY = dY+2.0*Len;}
        if (dZ > Len){dZ = dZ-2.0*Len;}
        if (dZ <= -Len){dZ = dZ+2.0*Len;}
        d = sqrt(pow(dX,2.0) + pow(dY,2.0) + pow(dZ,2.0));
        e = e+ 4.0 * (pow((1/d),12) - pow((1/d),6));
    }
    return e;
}



#endif /* homework2_h */
