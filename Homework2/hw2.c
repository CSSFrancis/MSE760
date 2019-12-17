
//
//  hw2.c
//  MSE760
//
//  Created by Carter Francis on 10/10/18.
//  
//

#include <math.h>
#include <stdio.h>
#include "FunctionHeader.h"
/*
int main() {
    double energyArray[48];
    FILE *out_energy_file = fopen("/Users/shaw/Shaw/MSE760/Out/energiesNonPeriodic.txt", "w");
    for (int a=0; a<10; a=a+1){
        double sigma = 3.4;
        double uCell = 5.7/sigma;// in angstroms
        int nCells = a+1;
        double epsilon = 0.0104;
        int nAtoms = nCells*nCells*nCells*4;
        double X[nAtoms];
        double Y[nAtoms];
        double Z[nAtoms];
        double Xforce[nAtoms];
        double Yforce[nAtoms];
        double Zforce[nAtoms];
        double Length = uCell*nCells/2;
        
        createFCC(nCells, uCell, X, Y, Z);
        //getting the values from aPos
        // basis
        double en [nAtoms-1];
        double energy = calculate_Potential(X, Y, Z, en,Xforce,Yforce,Zforce, nAtoms,Length);
        double energyPerAtom= (energy/nAtoms)*epsilon;
        char buf[0x100];
        snprintf(buf, sizeof(buf), "/Users/shaw/Shaw/MSE760/Out/%i.txt",a);
        FILE *out_file = fopen(buf, "w"); // write only
        int i;
        for (i=0; i<nAtoms;i++)
        {
            fprintf(out_file, "%g,%g,%g \n", X[i], Y[i], Z[i]);
            printf( "%g,%g,%g \n", X[i], Y[i], Z[i]);
        }
        energyArray[a] = energyPerAtom;
        fprintf(out_energy_file, "%i, %g \n",nCells, energyPerAtom);
    }
    
    
    
    
    return 0;
    
}
*/
