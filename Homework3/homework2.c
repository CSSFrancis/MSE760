//
//  homework2.c
//  MSE760
//
//  Created by Carter Francis on 11/1/18.
//  Copyright © 2018 Carter Francis. All rights reserved.
//

#include "homework2.h"
int main() {
    //initalizing varibles
    int a=4;
    double sigma = 3.4;
    double uCell = 5.7/sigma;// in angstroms
    int nCells = a+1;
    double epsilon = 0.0104;
    int nAtoms = nCells*nCells*nCells*4;
    double X[nAtoms]; double Y[nAtoms]; double Z[nAtoms];
    int timesteps =200;
    double stepSize =0.05;
    double EStep[timesteps];
    double T =2.0;
    double Length = uCell*nCells/2;
    //initalizing simulation
    createFCC(nCells, uCell, X, Y, Z);
    double en [nAtoms-1];
    calculate_Potential(X, Y, Z, en, nAtoms,Length);
    FILE *out_energy_file = fopen("/Users/shaw/Shaw/MSE760/Out/energiesNonPeriodic.txt", "w");
    
    for (int l=0; l<timesteps; l++)
    {
        EStep[l]= monteCarloStep(X,Y,Z, T, nAtoms, stepSize,Length);
        fprintf(out_energy_file, "Step: %i, %g ev/atom \n",l, EStep[l]);
    }
    char buf[0x100];
    snprintf(buf, sizeof(buf), "/Users/shaw/Shaw/MSE760/Out/%i.txt",a);


    return 1;
}
