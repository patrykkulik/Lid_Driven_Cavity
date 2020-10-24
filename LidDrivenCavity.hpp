//
//  LidDrivenCavity.hpp
//  LidDrivenCavity
//
//  Created by Patryk K on 13/3/20.
//  Copyright © 2020 Patryk K. All rights reserved.
//

#ifndef LidDrivenCavity_hpp
#define LidDrivenCavity_hpp

#endif /* LidDrivenCavityhpp */

/**
 * @class LidDrivenCavity
 * @brief Initialises and Solves the Lid Driven Cavity problem.
 */

class LidDrivenCavity{
    
public:
    /// Class constructor
    LidDrivenCavity(int argc, char **argv);
    /// Class destructor
    ~LidDrivenCavity();
    
    /// Set domain size
    void SetDomainSize(double xlen, double ylen);
    /// Set grid size
    void SetGridSize(int nx, int ny);
    /// Set time step
    void SetTimeStep(double deltat);
    /// Set final time
    void SetFinalTime(double finalt);
    /// Set Reynolds number
    void SetReynoldsNumber(double re);
    /// Set grid split
    void SetGridSplit(int px, int py);
    /// Set boundary velocity
    void SetBoundaryVelocity(double u);
    
    /// Initialise the problem
    int  Initialise(int argc, char **argv);
    /// Compute vorticity in parallel
    void ComputeVorticity();
    /// Compute vorticity in serial
    void ComputeVorticitySerial();
    /// Update boundary vorticity
    void UpdateBoundaryVorticity(int i, int j);
    /// Communicate vorticity between processes
    void CommunicateVorticity();
    /// Integrate vorticity in parallel
    void Integrate();
    /// Integrate vorticity in serial
    void IntegrateSerial();
    /// Solve the problem
    void Solve();
    /// Gather the solution to the problem on the root process
    void Gather();
    /// Output the solution to file
    void Output();
    
    
private:
    double* v = nullptr;    // vorticity array (defined in row major)
    double* s = nullptr;    // streamfunction array (defined in row major)
    
    double dt;      // time step
    double T;       // Final time
    int    Nx;      // Number of grid points in x-direction.
    int    Ny;      // Number of grid points in y-direction.
    double Lx;      // Domain size in the x-direction
    double Ly;      // Domain size in the y-direction
    double Re;      // Reynolds number
    int Px;         // Number of partitions in the x-direction (parallel)
    int Py;         // Number of partitions in the x-direction (parallel)
    double U = 1;   // Top boundary velocity
    double dx;      // Grid spacing in x-direction
    double dy;      // Grid spacing in y-direction
    int nt;         // Number of time steps
    
    int mx;
    int my;
    
    // Coefficients of the stream function used in the final equation
    double alpha;
    double beta;
    double gamma;
    double nu;
    
    
    // MPI variables:
    
    int izero=0;
    int myrank_mpi; // rank of the current process in MPI
    int nprocs_mpi; // total number of processes in MPI
    
    int nb_r;       // (Global) Block size (rows). This will be used to define a grid of processes
    int nb_c;       // (Global) Block size (columns). Note both are integers.
    char layout='R';// Block cyclic, Row major processor mapping
    
    // BLACS variables:
    
    int iam;        // current blacs process (rank of the blacs process)
    int nprocs;     // total number of blacs processes
    int zero = 0;
    int ictxt;      // blacs context handle
    int myrow;      // row grid position relative to the total number of row processes (nprow)
    int mycol;      // column grid position relative to the tolal number of column processes (npcol)
    
    int mpA; // number of rows of the local matrix
    int nqA; // number of columns of the local matrix
    
    
    
    int* location_loc = nullptr;  // Reference to how many rows and columns each block has
    int* locations_0 = nullptr;   // Array that will define the sizes of all local matrices
    
    int* sizes = nullptr;   // An empty sizes array to store the full size of the local matrices
    
    int start=0;    // Starting position integer. This will indicate the starting global matrix position of each local array
    
    int* starts_0 = nullptr; // Array that contains all starting positions on the root process
    
    // Create temporary arrays for psi and w which will be populated within each process
    double* temp_s = nullptr;
    double* temp_s2 = nullptr;
    double* temp_v = nullptr;
    double* temp_v2 = nullptr;
    
    
    int counter = 0;    // A counter of the Jacobi solver while loop
    int I;              // Global row of the current position in the local matrix
    int J;              // Global column of the current position in the local matrix
    int gpos = 0;       // Global position of the current local matrix position
    
    int* displs = nullptr; // Offset of the local arrays in which the arrays should be comunicated to the root process
    
    
    
    // Define psi and w boundary arrays. First set of arrays will represent the point at the edges of local arrays, while the second set of arrays will the represent the points outside of a given local array that need to be received from another process
    
    double* s_north = nullptr;
    double* s_east = nullptr;
    double* s_south = nullptr;
    double* s_west = nullptr;
    
    double* v_north = nullptr;
    double* v_east = nullptr;
    double* v_south = nullptr;
    double* v_west = nullptr;
    //
    //
    double* s_boundary_north = nullptr;
    double* s_boundary_east = nullptr;
    double* s_boundary_south = nullptr;
    double* s_boundary_west = nullptr;
    
    double* v_boundary_north = nullptr;
    double* v_boundary_east = nullptr;
    double* v_boundary_south = nullptr;
    double* v_boundary_west = nullptr;
    
    
    // Define the rank of processes at the boundaries of the current local array. Process = -1 indicates that the given process does not exist.
    int process_east = -1;
    int process_south = -1;
    int process_west = -1;
    int process_north = -1;
    
    // Define vector that will be used to compute the serial solution.
    double* p = nullptr;
};
