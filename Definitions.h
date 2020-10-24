//
//  Definitions.h
//  Lid_Driven_Cavity
//
//  Created by Patryk K on 14/3/20.
//

#ifndef Definitions_h
#define Definitions_h


#endif /* Definitions_h */

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include "mpi.h"
#include "cblas.h"

using namespace std;

// MARK: "C" function definitions"

#define F77NAME(x) x##_

extern "C" {
    // Perform LU factorisation and store the result in A - banded
    void F77NAME(dgbtrf) (const int& m, const int& n, const int& kl,
                          const int& ku, double* A, const int& lda,
                          int* ipiv, int& info);
    
    
    // Solve a linear system using a previously factored matrix A - banded
    void F77NAME(dgbtrs) (const char& trans, const int& n, const int& kl,
                          const int &ku, const int& nrhs, const double* A,
                          const int& lda, const int* ipiv, double* b,
                          const int& ldb, int& info);
    
    // Solve a linear system storing result in y (BANDED MATRIX)
    void F77NAME(dgbsv)(const int& n, const int& kl, const int& ku, const int& nrhs, const double * AB,
                        const int& ldab, int * ipiv, double * B,
                        const int& ldb, int& info);
    
    
    // Blacs functions used to initialise the parallelised grid
    void blacs_get_(int*, int*, int*);
    void blacs_pinfo_(int*, int*);
    void blacs_gridinit_(int*, char*, int*, int*);
    void blacs_gridinfo_(int*, int*, int*, int*, int*);
    //    void descinit_(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
    //    void pdpotrf_(char*, int*, double*, int*, int*, int*, int*);
    void blacs_gridexit_(int*);
    int numroc_(int*, int*, int*, int*, int*);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Print a col-major matrix for verification
void PrintMatrix_col(int rows, int cols, double* A) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout.width(7);
            cout << A[j * rows + i];
        }
        cout << endl;
    }
}

// Print a row-major matrix
void PrintMatrix_row(int rows, int cols, double* A) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << A[i*cols+j] << "  ";
        }
        cout << endl;
    }
}


// Print a vector (double)
void PrintVector(int n, double* u) {
    for (int i = 0; i < n; ++i) {
        cout << u[i] << endl;
    }
}

// Print a vector (int)
void PrintVector(int n, int* u) {
    for (int i = 0; i < n; ++i) {
        cout << u[i] << endl;
    }
}
