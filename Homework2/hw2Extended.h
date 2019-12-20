#ifndef MSE760_HW2EXTENDED_H
#define MSE760_HW2EXTENDED_H

//
//  Methods for Homework2. Still need some
//  MSE760
//
//  Created by Carter Francis on 10/10/18.
//  Revised 12/18/19
//


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


//Calculating forces and energy for some atom. It passes in the positions of the atoms

double calculate_Potential(double X[],double Y[], double Z[] , double energy[],double xForces[], double yForces[], double zForces[],int numAtoms, double L)
{
    //Resetting the Forces back to zero to reassign
    for (int i = 0; i<numAtoms; i++)
    {
        xForces[i]=0.0;
        yForces[i]=0.0;
        zForces[i]=0.0;

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


//assigning random velocities to the particles
int assignRandomVelocities(double xvel[], double yvel[], double zvel[], double temp,double atoms)
{
    srand (time(NULL));
    for (int i =0; i< atoms; i++)
    {
        int velocityNotAssigned =1;
        while (velocityNotAssigned)
        {

            double random1 = (double)rand()/RAND_MAX*2.0-1.0;//float in range -1 to 1
            double random2 = (double)rand()/RAND_MAX*2.0-1.0;
            double s2 = pow(random1,2)+ pow(random2,2);
            if (s2 <1)
            {
                xvel[i] = pow(temp*3.0, 0.5)*2.0*pow(1.0-s2, 0.5)*random1;
                yvel[i] = pow(temp*3.0, 0.5)*2.0*pow(1.0-s2, 0.5)*random2;
                zvel[i] = pow(temp*3.0, 0.5)*(1.0-2.0*s2);
                velocityNotAssigned =0;
            }
        }
    }
    return 0;
}

double calculateTemperature(double V_x[],double V_y[], double V_z[], int natoms)
{
    double temp=0.0;
    for (int i =0; i< natoms; i++)
    {
        temp += pow(V_x[i],2) + pow(V_y[i],2) +pow(V_z[i],2);
    }
    return temp*(1.0/(3.0*natoms));
}

double velocity_verlet(double X[], double Y[], double Z[], double xForce[], double yForce[], double zForce[], double xVel[],
                       double yVel[],double zVel[],double E[], int natoms, double timestep, double L, double* potentialE, double* kineticE)
{
    double Vhalfx[natoms];
    double Vhalfy[natoms];
    double Vhalfz[natoms];
    for (int i = 0; i< natoms; i++)
    {
        //updating positions
        Vhalfx[i] = xVel[i]+0.5*xForce[i]*timestep;
        Vhalfy[i] = yVel[i]+0.5*yForce[i]*timestep;
        Vhalfz[i] = zVel[i]+0.5*zForce[i]*timestep;
        double xholder = X[i] + (Vhalfx[i]*timestep);
        double yholder = Y[i] + (Vhalfy[i]*timestep);
        double zholder = Z[i] + (Vhalfz[i]*timestep);
        //Periodic boundry conditions
        if (xholder > L){xholder = xholder-2.0*L;}
        if (xholder <= -L){xholder = xholder+2.0*L;}
        if (yholder > L){yholder = yholder-2.0*L;}
        if (yholder <= -L){yholder= yholder+2.0*L;}
        if (zholder > L){zholder = zholder-2.0*L;}
        if (zholder <= -L){zholder = zholder+2.0*L;}
        X[i] = xholder;
        Y[i] = yholder;
        Z[i] = zholder;
    }
    //determine the acceleration by calculating the potential
    *potentialE = calculate_Potential(X, Y, Z, E, xForce, yForce, zForce, natoms,L);
    //calculating the velocity at time = t + timestep.
    *kineticE = 0;
    for (int i = 0; i< natoms; i++)
    {

        xVel[i]= Vhalfx[i]+0.5*xForce[i]*timestep;
        yVel[i]= Vhalfy[i]+0.5*yForce[i]*timestep;
        zVel[i]= Vhalfz[i]+0.5*zForce[i]*timestep;
        *kineticE = *kineticE + pow(pow(xVel[i],2.0)+pow(yVel[i],2.0)+pow(zVel[i],2.0),0.5);
    }
    printf("%g ", *kineticE);
    double energy = (2 * *kineticE)+*potentialE;


    return energy;
}

// Calculates the pair distribution function by determining the distance from
int pair_distribution_function(double X[],double Y[], double Z[], int num_bins, int numAtoms, double L, int pdf[])
{
    double bin_size = L/(num_bins-1);
    int total=0;
    for (int i = 0; i < numAtoms - 1; i = i + 1) {
        for (int j = i + 1; j < numAtoms; j = j + 1) {
            double difX = (X[i] - X[j]);
            double difY = (Y[i] - Y[j]);
            double difZ = (Z[i] - Z[j]);
            //For periodic boundry conditions
            if (difX > L) { difX = difX - 2.0 * L; }
            if (difX <= -L) { difX = difX + 2.0 * L; }
            if (difY > L) { difY = difY - 2.0 * L; }
            if (difY <= -L) { difY = difY + 2.0 * L; }
            if (difZ > L) { difZ = difZ - 2.0 * L; }
            if (difZ <= -L) { difZ = difZ + 2.0 * L; }
            //calculating the distance from atom i to j
            double dist = sqrt(pow(difX, 2.0) + pow(difY, 2.0) + pow(difZ, 2.0));
            if (dist<L) {
                int index = round(dist / bin_size);
                pdf[index]++;
                total++;
            }
        }
    }
    return total;
}
#endif