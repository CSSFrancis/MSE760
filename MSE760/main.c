//
//  main.c
//  MSE760
//
//  Created by Carter Francis on 9/6/18.
//  Copyright © 2018 Carter Francis. All rights reserved.
//
#include <math.h>
#include <stdio.h>

int main() {
    // insert code here...
    double uCell = .541;
    int nCells = 5;
    double aPos;
    aPos = createFCC(nCells, uCell);
    

    // basis
    
    return 0;
}

//functions
double createFCC (int numCells,double unitCell){
    int atomsPerCell =4;
    int numAtoms = (pow(numCells,3)*atomsPerCell);
    double atomPositions[3][500];
    static int count = 0;
    for (int i=0; i < numCells; i++) {
        for (int j = 0; j<numCells;j++){
            for (int k=0; k < numCells; k++) {
                atomPositions[0][count] = i*unitCell;
                atomPositions[1][count] = j*unitCell;
                atomPositions[2][count] = k*unitCell;
                atomPositions[0][count+1] = i*unitCell+ (unitCell/2);
                atomPositions[1][count+1] = j*unitCell+ (unitCell/2);
                atomPositions[2][count+1] = k*unitCell;
                count =count +2;
            }
        }
    }
    return **atomPositions;
    }
