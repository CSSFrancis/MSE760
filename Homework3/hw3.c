//
//  hw3.c
//  MSE760
//
//  Creating a Monte Carlo code to equalibriate an Ar System
//  Copyright © 2018 Carter Francis. All rights reserved.
//

#include "hw3.h"
#include <unistd.h>

int main() {
    //initializing variables
    srand (time(NULL));
    double cutoff= 5.0;
    int a = 4; //num of unit cells
    double sigma = 3.4;
    double uCell = 6.0/sigma; // in angstroms, for a density of .84?
    int nCells = a+1;
    double epsilon = 0.0104;
    int nAtoms = nCells*nCells*nCells*4;
    double X[nAtoms]; double Y[nAtoms]; double Z[nAtoms]; double E[nAtoms];
    int timesteps = 100000;
    double stepSize = 0.1;
    double T = 2.0;
    double Length = uCell * nCells / 2;
    // Initializing simulation
    createFCC(nCells, uCell, X, Y, Z);
    double en [nAtoms - 1];
    // Calculate the energy of every single atom.
    calculate_all_atom_energy(X, Y, Z, en, nAtoms,Length, cutoff);
    FILE *out_energy_file = fopen("../Homework3/out/energies.txt", "w");
    double numChanged = 0;
    for (int l=0; l<timesteps; l++)
    {
        numChanged += monteCarloStep(X, Y, Z, T, nAtoms, stepSize,Length, en, cutoff);
        if (l%50==0){
            double energySum = calculate_all_atom_energy(X, Y, Z, en ,nAtoms,Length, cutoff);
            fprintf(out_energy_file, "%g\n", energySum);
        }
    }
    printf("The %g Percent of steps accepted",numChanged/timesteps*100);
    return 1;
}
