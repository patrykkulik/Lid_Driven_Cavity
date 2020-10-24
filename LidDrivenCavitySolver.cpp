//
//  main.cpp
//  Lid_Driven_Cavity
//
//  Created by Patryk K on 9/3/20.
//

// run with make
// execute with mpiexec -np 25 ./LidDrivenCavitySolver --Lx 1 --Ly --Nx 161 --Ny 161 --Px 5 --Py 5 --dt 0.000975 --T 1 --Re 100

#include "cblas.h"
#include "LidDrivenCavity.cpp"


//--------------------------------------------------------------------------

/**
 * @brief Function that uses the LidDrivenCavity class to solve the Lid Driven Cavity problem posed in HPC assignment.
 */
int main(int argc, char **argv) {
    
    LidDrivenCavity problem (argc, argv); // Create an instance of the LidDrivenCavity class
    
    int check = 0;
    check = problem.Initialise(argc, argv); // Initialise the problem with inputs from the Command Line
    
    
    if (check != 0){ // Check if the problem was initialised correctly
        // Finalise MPI and exit if the problem was badly initialised
        MPI_Finalize();
        exit(0);
    }
    
    
    problem.Solve(); // Solve the problem (calls the PoissonSolver class within the Solve function)
    
    return 0;
}
