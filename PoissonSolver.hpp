//
//  PoissonSolver.hpp
//  Lid_Driven_Cavity
//
//  Created by Patryk K on 14/3/20.
//  Copyright © 2020 Patryk K. All rights reserved.
//

#ifndef PoissonSolver_hpp
#define PoissonSolver_hpp

#endif /* PoissonSolver_hpp */

#include "Definitions.h"

/**
 * @class PoissonSolver
 * @brief Solves the 2D Poisson problem
 */

class PoissonSolver
{
public:
    /// Poisson constructor
    PoissonSolver(int nx, int ny, double* vorticity, double* psi, double Alpha, double Beta, double Gamma, double Nu, int MPA, int NQA, double* S_north, double* S_east, double* S_south, double* S_west, double* S_boundary_north, double* S_boundary_east, double* S_boundary_south, double* S_boundary_west, int Process_east, int Process_south, int Process_west, int Process_north, int Start, int Myrank_mpi, int Nprocs_mpi); // Should I have functions to set all these variables?
    /// Destructor
    ~PoissonSolver();
    /// Initialise the Serial Solver
    void InitialiseSerial();
    /// Run Jacobi Solver
    void Jacobi(int MAX_COUNTER);
    /// Initialise Conjugate Gradient Solver
    void InitialiseCG();
    /// Calculate Boundary Values for the Conjugate Gradient Solver
    void BoundaryCG();
    /// Run Conjugate Gradient Solver
    void ConjugateGradient();
    /// Run Serial Solver
    void Serial(double* P);
    
    
private:
    double* v = nullptr;    // vorticity array (defined in row major)
    double* s = nullptr;    // streamfunction array (defined in row major)
    
    int    Nx;  // Number of grid points in x-direction.
    int    Ny;  // Number of grid points in y-direction.
    
    // Coefficients of the stream function used in the poisson equation
    double alpha;
    double beta;
    double gamma;
    double nu;
    
    
    // MPI variables:
    
    int myrank_mpi;
    
    int mpA; // number of rows of the local matrix
    int nqA; // number of columns of the local matrix
    
    int counter;    // A counter of the Jacobi solver while loop
    int max_counter;
    int I;              // Global row of the current position in the local matrix
    int J;              // Global column of the current position in the local matrix
    int gpos = 0;       // Global position of the current local matrix position
    
    int start;  // Starting position of the local array relative to the global array
    
    
    
    // Define psi and w boundary arrays. First set of arrays will represent the point at the edges of local arrays, while the second set of arrays will the represent the points outside of a given local array that need to be received from another process
    
    double* s_north = nullptr;
    double* s_east = nullptr;
    double* s_south = nullptr;
    double* s_west = nullptr;
    //
    //
    double* s_boundary_north = nullptr;
    double* s_boundary_east = nullptr;
    double* s_boundary_south = nullptr;
    double* s_boundary_west = nullptr;
    
    
    // Define the rank of processes at the boundaries of the current local array. Process = -1 indicates that the given process does not exist.
    int process_east = -1;
    int process_south = -1;
    int process_west = -1;
    int process_north = -1;
    
    // Variables needed for Conjugate Gradient function:
    //-----------------------
    double* A = nullptr;
    double* temp_v = nullptr;
    double* temp_s = nullptr;
    
    int mx;
    int my;
    
    int nprocs_mpi;
    
    double* r = nullptr;
    double* temp = nullptr; //temp
    double* p = nullptr;
    double* p_big = nullptr;
    int k = 0;
    double rTr_local;
    double rTr_global;
    double rTr_prev;
    double pTAp_local;
    double pTAp_global;
    double coeff_1;
    double coeff_2;
    double eps;
    double tol = 0.00001;
    
    //-----------------------
    
    // Variables needed for Serial implementation:
    int* ipiv = nullptr;
    int info = 0;
    
    
};
