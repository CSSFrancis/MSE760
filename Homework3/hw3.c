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


    int a = 4; //num of unit cells
    double sigma = 3.4;
    double density= 0.84; // # of atoms in a volume sigma squared
    double uCell = pow((4.0 / density), .3333333333333); // 4 atoms per fcc. Volume per FCC ^1/3 for Length
    printf(" The unit cell dimensions are: %g\n", uCell);
    int nCells = a+1;
    double epsilon = 0.0104;
    int nAtoms = nCells*nCells*nCells*4;
    // Positions of the atoms and energy of each atom
    double X[nAtoms]; double Y[nAtoms]; double Z[nAtoms]; double E[nAtoms];
    int timesteps = 60000;
    double stepSize =.2;
    double T = 2.0;
    double Length = uCell * nCells / 2;
    double cutoff= 2.5; // 2.5 Sigma in real units
    // Initializing simulation
    createFCC(nCells, uCell, X, Y, Z);
    double en [nAtoms - 1];
    // Calculate the energy of every single atom.
    calculate_all_atom_energy(X, Y, Z, en, nAtoms,Length);
    FILE *out_energy_file = fopen("../Homework3/out/energies.txt", "w");
    double numChanged = 0;
    for (int l=0; l<timesteps; l++)
    {
        numChanged += monteCarloStep(X, Y, Z, T, nAtoms, stepSize,Length, cutoff);
        if (l%5==0){
            double energySum = calculate_all_atom_energy(X, Y, Z, en ,nAtoms,Length);
            printf( "Energy: %g\n", (energySum/nAtoms)*epsilon);
            fprintf(out_energy_file, "%g\n", (energySum/nAtoms)*epsilon);
        }
    }
    printf("The %g Percent of steps accepted",numChanged/timesteps*100);
    return 1;
}
