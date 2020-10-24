//
//  PoissonSolver.cpp
//  Lid_Driven_Cavity
//
//  Created by Patryk K on 14/3/20.
//  Copyright © 2020 Patryk K. All rights reserved.
//

#include "PoissonSolver.hpp"


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// MARK: Class function definitions:

PoissonSolver::PoissonSolver(int nx, int ny, double* vorticity, double* psi, double Alpha, double Beta, double Gamma, double Nu, int MPA, int NQA, double* S_north, double* S_east, double* S_south, double* S_west, double* S_boundary_north, double* S_boundary_east, double* S_boundary_south, double* S_boundary_west, int Process_east, int Process_south, int Process_west, int Process_north, int Start, int Myrank_mpi, int Nprocs_mpi){ // Constructor function
    
    // Variables below have the same definition as the equivalent variables in LidDrivenCavity class:
    Nx  = nx;
    Ny  = ny;
    v   = vorticity;
    s   = psi;
    //
    alpha   = Alpha;
    beta    = Beta;
    gamma   = Gamma;
    nu      = Nu;
    //
    mpA = MPA;
    nqA = NQA;
    //
    s_north = S_north;
    s_east  = S_east;
    s_south = S_south;
    s_west  = S_west;
    //
    s_boundary_north    = S_boundary_north;
    s_boundary_east     = S_boundary_east;
    s_boundary_south    = S_boundary_south;
    s_boundary_west     = S_boundary_west;
    //
    process_east    = Process_east;
    process_south   = Process_south;
    process_west    = Process_west;
    process_north   = Process_north;
    //
    start = Start;
    //
    myrank_mpi = Myrank_mpi;
    nprocs_mpi = Nprocs_mpi;
    
    // Create a temporary streamfunctionarray
    temp_s = new double[mpA*nqA];
    
    for (int i = 0; i<mpA;i++){
        for (int j =0; j<nqA; j++){
            temp_s[i*nqA + j] = 0;
        }
    }
    
    
}

PoissonSolver::~PoissonSolver(){ // Class destructor
    
    // Delete all pointers unique to this class
    delete A;
    delete temp_s;
    delete temp_v;
}

//****************************************************************************************************************
//********************************************SERIAL SOLVER*******************************************************
//****************************************************************************************************************

void PoissonSolver::InitialiseSerial(){
    // Define the inner matrix dimensions
    mx = Nx - 2;
    my = Ny - 2;
    
    // Define the banded matrix that describes the system of equations: A*s = v
    A = new double[mx*my*(1 + 2*my + my)];
    
    // Define the pivotes array and info
    ipiv = new int[mx*my];
    info = 0;
    
    for (int i = 0; i<mx*my; i++){
        ipiv[i] = 0;
    }
    
    // Populate the matrix A:
    for (int j = 0; j<mx*my; j++){
        for (int k = 0; k<my; k++){
            A[(1 + 2*my + my)*j + k] = 0;
        }
        
        A[(1 + 2*my + my)*j + my] = alpha;
        
        for (int k = my+1; k<2*my - 1; k++){
            A[(1 + 2*my + my)*j + k] = 0;
        }
        
        A[(1 + 2*my + my)*j + 2*my - 1] = gamma;
        A[(1 + 2*my + my)*j + 2*my] = beta+nu;
        A[(1 + 2*my + my)*j + 2*my + 1] = gamma;
        
        if ((j+1)%my == 0) {
            A[(1 + 2*my + my)*j + 2*my + 1] = 0;
        }
        else if ((j)%my == 0){
            A[(1 + 2*my + my)*j + 2*my - 1] = 0;
        }
        
        
        for (int k = 2*my+2; k<3*my; k++){
            A[(1 + 2*my + my)*j + k] = 0;
        }
        
        A[(1 + 2*my + my)*j + 3*my] = alpha;
        
    }
    
    
    // MARK: Perform LU factorisation of the banded matrix
    F77NAME(dgbtrf)(mx*my,mx*my, my, my, A, (1 + 2*my + my), ipiv, info);
    
    // Check if the factorisation was done correctly
    if (info != 0){
        cout << "A problem occured in dbtrf" << endl;
    }
}


void PoissonSolver::Serial(double* P){
    // define the vorticity array
    p = P;
    
    F77NAME(dgbtrs)('N', mx*my, my, my, 1, A, (2*my+my+1), ipiv, p, mx*my, info); // Initially, p holds the vorticity. After this line, it holds the inner points of phi
    
    // check for issues
    if (info != 0){
        cout << "A problem occured in dbtrs" << endl;
    }
    
    // Convert the p array back into streamfunction matrix
    for (int i = 1; i <= mx; i++){
        for (int j = 1; j <= my; j++){
            s[i*Ny + j] = p[(i-1)*my + j-1];
        }
    }
}


//****************************************************************************************************************
//********************************************JACOBI SOLVER*******************************************************
//****************************************************************************************************************


void PoissonSolver::Jacobi(int MAX_COUNTER){
    
    // Initialise the counter and epsilon
    counter = 0;
    max_counter = MAX_COUNTER;
    eps = 1000;
    
    // Perform the Jacobi computation. A reduction in max_counter may produce wrong results
    while (counter < max_counter){
        
        for (int i = 0; i < mpA; i++) { // local row
            for (int j = 0; j < nqA; j++) { // local col
                gpos = start + i*Ny + j;
                
                I = gpos / Ny;// global row
                J = gpos % Ny;// global col
                
                if (I == 0 || I== Nx - 1 || J == 0 || J == Ny - 1){ // Ignore the boundary points
                    continue;
                }
                
                // Perform the Jacobi iteration
                s[I*Ny + J] = s[I*Ny + J] + ((-alpha*s[(I+1)*Ny + J] - alpha*s[(I-1)*Ny + J] - beta*s[(I*Ny + J)])
                                             + (-gamma*s[(I)*Ny + J+1] - gamma*s[(I)*Ny + J-1] - nu*s[(I*Ny + J)]) + v[I*Ny + J])*1.5 / (beta + nu);
                
                
                // Update the boundary arrays
                if (i == 0){
                    s_boundary_north[j] = s[I*Ny + J];
                }
                if (i == mpA - 1){
                    s_boundary_south[j] = s[I*Ny + J];
                }
                if (j == 0){
                    s_boundary_west[i] = s[I*Ny + J];
                    //
                }
                if (j == nqA - 1){
                    s_boundary_east[i] = s[I*Ny + J];
                }
            }
        }
        
        
        // If the counter is just before being divisible by 50, store the streamfunction array in a temporary variable.
        if (counter % 50 == 49){
            for (int i =0 ; i<mpA; i++){
                for (int j=0; j<nqA;j++){
                    gpos = start + i*Ny + j;
                    
                    I = gpos / Ny;// global row
                    J = gpos % Ny;// global col
                    if (I == 0 || I== Nx - 1 || J == 0 || J == Ny - 1){ // Ignore the boundary points
                        continue;
                    }
                    
                    temp_s[i*nqA + j] = s[I*Ny + J];
                }
            }
        }
        else if (counter % 50 == 0){ // If the counter is divisible by 50, check the error between the current streamfunction and the previous streamfunction
            
            // Define error terms:
            eps = 0;
            tol = 0;
            for (int i =0 ; i<mpA; i++){
                for (int j=0; j<nqA;j++){
                    gpos = start + i*Ny + j; // global postion
                    
                    I = gpos / Ny;// global row
                    J = gpos % Ny;// global col
                    if (I == 0 || I== Nx - 1 || J == 0 || J == Ny - 1){ // Ignore the boundary points
                        continue;
                    }
                    
                    if (abs(s[I*Ny + J]) > 0.000000000001){ // If the streamfunction is big enough, check error
                        eps += abs((temp_s[i*nqA + j] - s[I*Ny + J])/s[I*Ny + J]);
                        
                    }
                }
            }
            
            eps = eps/(mpA*nqA);
            
            // Combine the error terms on all processes and take their average
            MPI_Allreduce(&eps, &tol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            tol = tol/nprocs_mpi;
            
            // Break the Jacobi loop if the error is sufficiently small
            if (tol < 0.0000002){
                break;
            }
        }
        
        
        
        // If a given neighbour exists, communicate the relevant boundary array
        
        if (process_east > -1){
            MPI_Send(s_boundary_east,mpA,MPI_DOUBLE,process_east,0,MPI_COMM_WORLD);
        }
        
        if (process_south > -1){
            MPI_Send(s_boundary_south,nqA,MPI_DOUBLE,process_south,0,MPI_COMM_WORLD);
        }
        
        if (process_west > -1){
            MPI_Send(s_boundary_west,mpA,MPI_DOUBLE,process_west,0,MPI_COMM_WORLD);
        }
        
        if (process_north > -1){
            MPI_Send(s_boundary_north,nqA,MPI_DOUBLE,process_north,0,MPI_COMM_WORLD);
        }
        
        
        // If a given neighbour exists, receive the boundary array and assign it to the full array
        if (process_west > -1){
            MPI_Recv(s_west, mpA, MPI_DOUBLE, process_west, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (int i = 0; i<mpA; i++){
                s[start - 1 + i*Ny] = s_west[i];
            }
            
        }
        
        if (process_north > -1){
            MPI_Recv(s_north, nqA, MPI_DOUBLE, process_north, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (int i = 0; i<nqA; i++){
                s[start - Ny + i] = s_north[i];
            }
        }
        
        if (process_east > -1){
            MPI_Recv(s_east, mpA, MPI_DOUBLE, process_east, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (int i = 0; i<mpA; i++){
                s[start + nqA + i*Ny] = s_east[i];
            }
            
        }
        
        if (process_south > -1){
            MPI_Recv(s_south, nqA, MPI_DOUBLE, process_south, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (int i = 0; i<nqA; i++){
                s[start + Ny*mpA + i] = s_south[i];
            }
        }
        
        
        // increment the counter
        counter++;
    }
    
}


//****************************************************************************************************************
//*****************************************CONJUGATE GRADIENT SOLVER**********************************************
//****************************************************************************************************************


void PoissonSolver::InitialiseCG(){// Function to initialise the Conjugate Gradient function
    
    // Create a banded matrix A for the conjugate gradient solver
    A = new double[mpA*nqA*(1+2*nqA)];
    
    for (int j = 0; j<(mpA)*(nqA); j++){
        
        A[j*(1+2*(nqA))] = alpha;
        for (int k = 1; k<(nqA)-1; k++){
            A[j*(1+2*(nqA)) + k] = 0;
        }
        
        A[j*(1+2*(nqA)) + (nqA)-1] = gamma;
        A[j*(1+2*(nqA)) + (nqA)] = beta + nu;
        A[j*(1+2*(nqA)) + (nqA)+1] = gamma;
        
        for (int k = 1; k<(nqA)-1; k++){
            A[j*(1+2*(nqA)) + (nqA)+1 + k] = 0;
        }
        
        A[j*(1+2*(nqA)) + 2*(nqA)] = alpha;
        
        if ((j+1)%(nqA) == 0 && j!=1) {
            A[j*(1+2*(nqA)) + (nqA)+1] = 0;
        }
        else if ((j)%(nqA) == 0){
            A[j*(1+2*(nqA)) + (nqA)-1] = 0;
        }
    }
    
    
    
    // Create a temporary arrays for vorticity and psi on the local processes
    temp_v = new double[(mpA)*(nqA)];
    //
    // Populate the temp_v and temp_s arrays
    for (int i = 0; i < (mpA); i++) { // local row
        for (int j = 0; j < (nqA); j++) { // local col
            gpos = start + i*Ny + j;
            //
            I = gpos / Ny;// global row
            J = gpos % Ny;// global col
            temp_v[i*(nqA) + j] = 0;
            
            temp_v[i*(nqA) + j] = v[I*Ny + J];
            
            temp_s[i*(nqA) + j] = 0.0001;
            
        }
    }
    
    
}


///////////////////////////////////////////////////////////////////////////////////////////

void PoissonSolver::BoundaryCG(){
    // Update the boundaries of the local arrays
    
    for (int i = 0; i < mpA; i++) { // local row
        for (int j = 0; j < nqA; j++) { // local col
            gpos = start + i*Ny + j;
            
            I = gpos / Ny;// global row
            J = gpos % Ny;// global col
            
            
            p_big[I*Ny + J] = p[i*nqA +j];
            
            if (i == 0){
                s_boundary_north[j] = p_big[I*Ny + J];
            }
            if (i == mpA - 1){
                s_boundary_south[j] = p_big[I*Ny + J];
            }
            if (j == 0){
                s_boundary_west[i]  = p_big[I*Ny + J];
                //
            }
            if (j == nqA - 1){
                s_boundary_east[i]  = p_big[I*Ny + J];
            }
            
        }
    }
    
    // If a given neighbour exists, communicate the relevant boundary array
    if (process_east > -1){
        MPI_Send(s_boundary_east,mpA,MPI_DOUBLE,process_east,0,MPI_COMM_WORLD);
    }
    
    if (process_south > -1){
        MPI_Send(s_boundary_south,nqA,MPI_DOUBLE,process_south,0,MPI_COMM_WORLD);
    }
    
    if (process_west > -1){
        MPI_Send(s_boundary_west,mpA,MPI_DOUBLE,process_west,0,MPI_COMM_WORLD);
    }
    
    if (process_north > -1){
        MPI_Send(s_boundary_north,nqA,MPI_DOUBLE,process_north,0,MPI_COMM_WORLD);
    }
    
    
    // If a given neighbour exists, receive the boundary array and assign it to the full array
    if (process_west > -1){
        MPI_Recv(s_west, mpA, MPI_DOUBLE, process_west, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        for (int i = 0; i<mpA; i++){
            p_big[start - 1 + i*Ny] = s_west[i];
        }
        
    }
    
    if (process_north > -1){
        MPI_Recv(s_north, nqA, MPI_DOUBLE, process_north, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        for (int i = 0; i<nqA; i++){
            p_big[start - Ny + i] = s_north[i];
        }
    }
    
    if (process_east > -1){
        MPI_Recv(s_east, mpA, MPI_DOUBLE, process_east, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        for (int i = 0; i<mpA; i++){
            p_big[start + nqA + i*Ny] = s_east[i];
        }
        
    }
    
    if (process_south > -1){
        MPI_Recv(s_south, nqA, MPI_DOUBLE, process_south, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        for (int i = 0; i<nqA; i++){
            p_big[start + Ny*mpA + i] = s_south[i];
        }
    }
    
}



void PoissonSolver::ConjugateGradient(){
    // Define variables needed for conjugate gradient
    r = new double[(mpA)*(nqA)];
    temp = new double[(mpA)*(nqA)];
    p = new double[(mpA)*(nqA)];
    p_big = new double[Nx*Ny];
    
    int k = 0;
    
    for (int i = 0; i < (mpA); i++) { // local row
        for (int j = 0; j < (nqA); j++) { // local col
            gpos = start + i*Ny + j;
            //
            I = gpos / Ny;// global row
            J = gpos % Ny;// global col
            temp_v[i*(nqA) + j] = 0;
            
            temp_v[i*(nqA) + j] = v[I*Ny + J];
            
            temp_s[i*(nqA) + j] = s[I*Ny + J];
            
        }
    }
    
    // Initialise the arrays
    for (int i = 0; i<mpA; i++){
        for (int j = 0; j<nqA; j++){
            r[i*(nqA) + j] = 0;
            p[i*(nqA) + j] = 0;
            temp[i*(nqA) + j] = 0;
        }
    }
    
    for (int i = 0; i<Nx; i++){
        for (int j = 0; j<Ny; j++){
            p_big[i*(Ny) + j] = 0;
        }
    }
    
    cblas_dcopy((mpA)*(nqA), temp_v, 1, r, 1);        // r_0 = temp_w
    
    cblas_dgbmv(CblasColMajor,CblasNoTrans,(mpA)*(nqA),(mpA)*(nqA),(nqA),(nqA),-1.0,A,2*(nqA)+1, temp_s,1,1.0,r,1); // r = r - A*temp_s
    
    cblas_dcopy((mpA)*(nqA), r, 1, p, 1);             // p_0 = r_0
    
    // Populate the p_big array
    for (int i = 0; i < mpA; i++) { // local row
        for (int j = 0; j < nqA; j++) { // local col
            gpos = start + i*Ny + j;
            
            I = gpos / Ny;// global row
            J = gpos % Ny;// global col
            
            p_big[I*Ny + J] = p[i*nqA +j];
            
        }
    }
    
    rTr_local = cblas_ddot((mpA)*(nqA), r, 1, r, 1);  // r_k^T r_k
    
    // Combine all local r dot products on all processes
    MPI_Allreduce(&rTr_local, &rTr_prev, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    
    while (k<50){      // Set a maximum number of iterations
        
        // Use iteration to compute the matrix vector product
        for (int i = 0; i < mpA; i++) { // local row
            for (int j = 0; j < nqA; j++) { // local col
                gpos = start + i*Ny + j;
                
                I = gpos / Ny;// global row
                J = gpos % Ny;// global col
                
                if (I == 0 || I== Nx - 1 || J == 0 || J == Ny - 1){ // Ignore the boundary points
                    continue;
                }
                
                
                temp[i*nqA +j] = (beta+nu)*p_big[I*Ny + J] + alpha*p_big[(I-1)*Ny + J] + alpha*p_big[((I+1)*Ny + J)]
                + gamma*p_big[I*Ny + J+1] + gamma*p_big[I*Ny + J - 1];
                
            }
        }
        
        
        pTAp_local = cblas_ddot((mpA)*(nqA), temp, 1, p, 1);  // p_k^T A p_k
        
        // Combine all local p_k^T A p_k products on all processes
        MPI_Allreduce(&pTAp_local, &pTAp_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        coeff_1 = rTr_prev / pTAp_global; // compute alpha_k
        
        cblas_daxpy((mpA)*(nqA), coeff_1, p, 1, temp_s, 1);  // x_{k+1} = x_k + alpha_k p_k
        cblas_daxpy((mpA)*(nqA), -coeff_1, temp, 1, r, 1); // r_{k+1} = r_k - alpha_k A p_k
        
        rTr_local = cblas_ddot((mpA)*(nqA), r, 1, r, 1);  // r_k^T r_k
        
        // Combine all local r dot products on all processes
        MPI_Allreduce(&rTr_local, &rTr_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // Break iteration if epsilon is small enough
        eps = rTr_global;
        if (eps < tol*tol) {
            cout << "Iteration broken" << endl;
            break;
        }
        coeff_2 = rTr_global/rTr_prev;
        
        rTr_prev = rTr_global;
        
        cblas_dcopy((mpA)*(nqA), r, 1, temp, 1);            // temp = r
        cblas_daxpy((mpA)*(nqA), coeff_2, p, 1, temp, 1);   // temp = r + beta * p
        cblas_dcopy((mpA)*(nqA), temp, 1, p, 1);            // p = temp
        
        k++;
        
        
        // Apply boundary conditions
        BoundaryCG();
        
        
    }
    
    // convert the temp array back to the general array
    for (int i = 0; i < mpA; i++) { // local row
        for (int j = 0; j < nqA; j++) { // local col
            gpos = start + i*Ny + j;
            
            I = gpos / Ny;// global row
            J = gpos % Ny;// global col
            
            s[I*Ny + J] = temp_s[i*nqA +j];
            
        }
    }
    
    // Clean up
    delete[] r;
    delete[] p;
    delete[] temp;
    delete[] p_big;
}
