//  Auxilliary functions for Homework #3.  MSE 760.
//  hw3.h
//  MSE760
//
//  Created by Carter Francis on 11/1/18.
//  Copyright © 2018 Carter Francis. All rights reserved.
//

#ifndef hw3_h
#define hw3_h

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

double calculate_Potential(double X[],double Y[], double Z[] , double energy[],double xForces[], double yForces[],
        double zForces[],int numAtoms, double L)
{
    //Resetting the Forces back to zero to reassign
    for (int i = 0; i<numAtoms; i++)
    {
        xForces[i]=0.0; yForces[i]=0.0; zForces[i]=0.0;
    }
    //double potential;
    double dist=0.0;
    double energySum =0.0;
    //double cutoff = .2;
    for (int i=0; i <numAtoms-1; i=i+1){
        //determining atoms within some cutoff
        double energyHolder =0.0;
        double xForce =0.0;
        double yForce =0.0;
        double zForce =0.0;
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
            double redDist = dist;

            double a=(48.0/pow(redDist,2.0))*((pow((1.0/redDist),12.0) - 0.5*pow((1.0/redDist),6.0)));
            xForce= xForce + ((difX)*a);
            yForce= yForce + ((difY)*a);
            zForce= zForce + ((difZ)*a);
            xForces[j] = xForces[j] - (difX*a);
            yForces[j] = yForces[j] - (difY*a);
            zForces[j] = zForces[j] - (difZ*a);
            energyHolder = energyHolder + 4.0 * (pow((1/redDist),12) - pow((1/redDist),6));
        }
        energySum = energySum +energyHolder;
        xForces[i]=xForces[i] +xForce;
        yForces[i]=yForces[i] +yForce;
        zForces[i]=zForces[i] +zForce;
        energy[i] = energyHolder;

    }
    return energySum;
}


// Calculating the energy of a single atom. Uses some cutoff for calculating the energy
double calculateAtomEnergy(double X[],double Y[], double Z[], int rAtom, double Xrand, double Yrand,double Zrand,
        int natom, double Len)
{
    double e = 0.0; // Initializing the energu
    double d; // initializing the distance
    for (int i=0; (i<natom && i !=rAtom); i++) // Compare atom with every other atom.
    {
        // Distance from atom i to the a random atom
        double dX =(X[i] - Xrand); double dY =(Y[i] - Yrand); double dZ =(Z[i] - Zrand);
        //For periodic boundary conditions
        if (dX > Len){dX = dX-2.0*Len;}
        if (dX <= -Len){dX = dX+2.0*Len;}
        if (dY > Len){dY = dY-2.0*Len;}
        if (dY <= -Len){dY = dY+2.0*Len;}
        if (dZ > Len){dZ = dZ-2.0*Len;}
        if (dZ <= -Len){dZ = dZ+2.0*Len;}
        d = sqrt(pow(dX, 2.0) + pow(dY,2.0) + pow(dZ ,2.0)); // Calculating the distance
        e = e + 4.0 * (pow((1/d),12) - pow((1/d),6));
        //if (d < cutoff){e = e + 4.0 * (pow((1/d),12) - pow((1/d),6));} // Only calculate potential energy if inside some cutoff
    }
    return e;
}

// Calculating the energy of a every atom. Still needs a reasonable cutoff to speed up the calculations...
double calculate_all_atom_energy(double X[],double Y[], double Z[], double energy[], int numAtoms, double L)
{
    double dist;
    double energySum =0.0;
    // Need to do everything the long way just to make sure that the energies for each atom are calculated correctly.
    for (int i=0; i <numAtoms-1; i=i+1){
        //determining atoms within some cutoff
        double energyHolder =0.0;
        double xi= X[i];
        double yi= Y[i];
        double zi= Z[i];

        for(int j=0; j <numAtoms && j!=i; j=j+1){
            // Calculating the distance from atom i to atom j
            double difX =(xi - X[j]);
            double difY =(yi - Y[j]);
            double difZ =(zi - Z[j]);
            
            // For periodic boundary conditions
            if (difX > L){difX = difX-2.0*L;}
            if (difX <= -L){difX = difX+2.0*L;}
            if (difY > L){difY = difY-2.0*L;}
            if (difY <= -L){difY = difY+2.0*L;}
            if (difZ > L){difZ = difZ-2.0*L;}
            if (difZ <= -L){difZ = difZ+2.0*L;}
            //calculating the distance from atom i to j
            dist = sqrt(pow(difX,2.0) + pow(difY,2.0) + pow(difZ,2.0));
            //if (dist < cutoff) {energyHolder = energyHolder + 4.0 * (pow((1/dist),12) - pow((1/dist),6));}
            energyHolder = energyHolder + 4.0 * (pow((1/dist),12) - pow((1/dist),6));
        }
        energySum = energySum + energyHolder;
        energy[i] = energyHolder;

    }
    return energySum;
}

// One MonteCarlo move.  Takes a random atom and randomly displaces it. The uses calculateAtomEnergy and determines
// the probability the move will be made.  It then takes another random number and accepts or denys the move.

double monteCarloStep(double X[],double Y[], double Z[], double T, int natoms, double stepSize, double L,double cutoff)
{
    // Pick a random atom from all of the atoms
    int randomAtom = (rand() % natoms);
    // Old positions
    double rmax = RAND_MAX;
    // Delta change of x,y,z
    double dx = rand()/ rmax -0.5; double dy = rand()/rmax -0.5; double dz = rand()/rmax-0.5;
    // New positions
    double newX = X[randomAtom] + dx * stepSize;
    double newY = Y[randomAtom] + dy * stepSize;
    double newZ = Z[randomAtom] + dz * stepSize;
    // Periodic boundary conditions
    if (newX > L){newX = newX-2.0*L;}
    if (newX <= -L){newX = newX+2.0*L;}
    if (newY > L){newY = newY-2.0*L;}
    if (newY <= -L){newY = newY + 2.0 * L;}
    if (newZ > L){newZ = newZ - 2.0 * L;}
    if (newZ <= -L){newZ = newZ + 2.0 * L;}
    // New energy.  Just calculating the energy for the one atom that changed position.
    double oldE = calculateAtomEnergy(X,Y,Z, randomAtom,X[randomAtom],Y[randomAtom], Z[randomAtom], natoms, L);
    double newE = calculateAtomEnergy(X, Y, Z, randomAtom, newX, newY, newZ, natoms, L);
    double prob = exp(-2*(newE-oldE)/T);
    int accepted = (prob>1 || prob > rand()/rmax);
    printf("Prob: %g , OldE: %g NewE: %g, Accepted: %i \n",prob, oldE, newE, accepted);
    //if (prob>100){printf("Prob: %g , OldE: %g NewE: %g \n",prob, oldE, newE);}

    if (accepted==1)
    {
        X[randomAtom]=newX;Y[randomAtom]=newY;Z[randomAtom]=newZ;
    }
    return accepted;
}



#endif /* hw3.h */
