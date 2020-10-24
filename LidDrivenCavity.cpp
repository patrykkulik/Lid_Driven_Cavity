//
//  LidDrivenCavity.cpp
//  Lid_Driven_Cavity
//
//  Created by Patryk K on 13/3/20.
//

#include "LidDrivenCavity.hpp"
#include "PoissonSolver.cpp"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// MARK: Class function definitions:

LidDrivenCavity::LidDrivenCavity(int argc, char **argv){ // Constructor function
    
    // Read the arguments provided on the command line and assign them to the appropriate variables:
    for (int i = 1; i<argc; i+=1){
        if (std::string(argv[i]) == "--Lx"){
            Lx = atof(argv[i+1]);   // Domain size in the x direction
        }
        else if (std::string(argv[i]) == "--Ly"){
            Ly = atof(argv[i+1]);   // Domain size in the y direction
        }
        else if (std::string(argv[i]) == "--Nx"){
            Nx = atof(argv[i+1]);   // Domain discretisation in the x direction
        }
        else if (std::string(argv[i]) == "--Ny"){
            Ny = atof(argv[i+1]);   // Domain discretisation in the y direction
        }
        else if (std::string(argv[i]) == "--Px"){
            Px = atof(argv[i+1]);   // Number of partitions in the x-direction (parallel)
        }
        else if (std::string(argv[i]) == "--Py"){
            Py = atof(argv[i+1]);   // Number of partitions in the y-direction (parallel)
        }
        else if (std::string(argv[i]) == "--dt"){
            dt = atof(argv[i+1]);   // Time step
        }
        else if (std::string(argv[i]) == "--T"){
            T = atof(argv[i+1]);    // Final time
        }
        else if (std::string(argv[i]) == "--Re"){
            Re = atof(argv[i+1]);   // Reynolds number
        }
    }
    
    
    // Call the class functions to override the private variables
    SetDomainSize(Lx, Ly);
    SetGridSize(Nx, Ny);
    SetGridSplit(Px, Py);
    SetTimeStep(dt);
    SetFinalTime(T);
    SetReynoldsNumber(Re);
    
}

LidDrivenCavity::~LidDrivenCavity(){ // Class destructor
    
    delete s_north;
    delete s_east;
    delete s_south;
    delete s_west;
    //
    delete v_north;
    delete v_east;
    delete v_south;
    delete v_west;
    //
    //
    delete s_boundary_north;
    delete s_boundary_east;
    delete s_boundary_south;
    delete s_boundary_west;
    //
    delete v_boundary_north;
    delete v_boundary_east;
    delete v_boundary_south;
    delete v_boundary_west;
    //
    delete displs;
    //
    delete location_loc;
    delete locations_0;
    //
    delete sizes;
    //
    delete starts_0;
    //
    delete temp_v;
    delete temp_v2;
    delete temp_s;
    delete temp_s2;
    //
    delete v;
    delete s;
    
    
    // Finalise blacs and MPI
    blacs_gridexit_(&ictxt);
    MPI_Finalize();
    
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen){ // Function to set the domain size variables
    Lx = xlen;
    Ly = ylen;
}

void LidDrivenCavity::SetGridSize(int nx, int ny){ // Function to set the grid size variables
    Nx = nx;
    Ny = ny;
}

void LidDrivenCavity::SetGridSplit(int px, int py){ // Function to set the parallel grid split variables
    Px = px;
    Py = py;
}

void LidDrivenCavity::SetTimeStep(double deltat){ // Function to set the time step variable
    dt = deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt){ // Function to set the final time variable
    T = finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re){ // Function to set the Reynolds number variable
    Re = re;
}

void LidDrivenCavity::SetBoundaryVelocity(double u){ // Function to set the boundary velocity variable
    U = u;
}


////////////////////////////////////////////////////////////////////////////////////////////////////


int LidDrivenCavity::Initialise(int argc, char **argv){ // Function to initialise all necessary problem variables and the MPI processes needed to solve the problem
    
    
    // MARK: Define variables
    dx = Lx / (Nx-1);       // Grid spacing in the x-direction
    dy = Ly / (Ny-1);       // Grid spacing in the y-direction
    nt = ceil(T/dt);        // Number of time steps
    
    
    // Coefficients of the stream function used in the final equation
    alpha = -1/(dx*dx);
    beta = 2/(dx*dx);
    gamma = -1/(dy*dy);
    nu = 2/(dy*dy);
    
    
    v = new double[Nx*Ny];  // Vorticity
    s = new double[Nx*Ny];  // Stream function
    
    // Initialise the arrays to 0 at all points:
    for (int i = 0; i<Nx*Ny; i++){
        v[i] = 0;
        s[i] = 0;
    }
    
    
    // MARK: Initialise MPI
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi); // myrank_mpi will now contain the rank of the current process
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi); // nprocs_mpi will now contain the total number of processes available
    
    
    // Check if the number of processes matches the grid discretisation:
    if (Px * Py != nprocs_mpi){
        if (myrank_mpi == 0){
            cout << "Error: The grid discretisation does not match the specified number of processes" << endl;
        }
        return 1;
    }
    
    nb_r = Nx/Px;    // (Global) Block size (rows). This will be used to define a grid of processes
    nb_c = Ny/Py;    // (Global) Block size (columns). Note both are integers.
    
    
    // MARK: Initialize BLACS grid
    
    blacs_pinfo_(&iam, &nprocs) ;       // initialise the BLACS world communicator
    blacs_get_(&zero, &zero, &ictxt );  // -> Create context - overrides ictxt
    
    blacs_gridinit_(&ictxt, &layout, &Px, &Py );          // Context -> Initialize the grid
    blacs_gridinfo_(&ictxt, &Px, &Py, &myrow, &mycol );   // Context -> Context grid info (# procs row/col, current procs row/col)
    
    
    // Compute the size of the local matrices
    mpA = numroc_( &Nx, &nb_r, &myrow, &izero, &Px ); // number of rows of the local matrix
    nqA = numroc_( &Ny, &nb_c, &mycol, &izero, &Py ); // number of columns of the local matrix
    
    
    ////////////////////////////////////////////
    // MARK: Compute additional variables needed for the problem definition
    
    location_loc = new int[2]; // Create a reference to how many rows and columns each block has
    location_loc[0] = mpA;
    location_loc[1] = nqA;
    
    locations_0 = new int[2*nprocs]; // Create an array that will contain all the local number of rows and columns of each process and will be available to each process
    
    MPI_Gather(location_loc,2,MPI_INT,locations_0,2,MPI_INT,0, MPI_COMM_WORLD); // Communicate the number of rows and column of each block to the root process and store them in locations_0
    MPI_Bcast(locations_0, 2*nprocs, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast the local sizes to all processes
    
    
    
    sizes = new int[nprocs]; // Define an empty sizes array to store the full size of the local matrices
    
    
    for (int i =0; i<nprocs; i++){
        sizes[i] = locations_0[i*2]*locations_0[i*2+1];     // Define the size based on the locations array. This will be the same for all processes.
    }
    
    // Populate start array to contain the starting global position of each local process. This is unique to each process:
    for (int i = 0; i<mycol; i++){
        start += locations_0[i*2 + 1];
    }
    
    for (int i = 0; i<myrow; i++){
        start += Ny*locations_0[i*2*Py];
    }
    
    
    starts_0 = new int[nprocs];
    
    MPI_Gather(&start,1,MPI_INT,starts_0,1,MPI_INT,0, MPI_COMM_WORLD); // Communicate the local starting positions to the root process
    
    
    ////////////////////////////////////////////
    
    // Create temporary arrays for vorticity and psi which will be populated within each process
    temp_s = new double[mpA*nqA];
    temp_s2 = new double[Nx*Ny];
    temp_v = new double[mpA*nqA];
    temp_v2 = new double[Nx*Ny];
    //
    //
    displs = new int[nprocs]; // Offset of the local arrays in which the arrays should be comunicated to the root process. Used in the GATHERV function
    
    
    // Initialise the arrays:
    for (int i = 0; i<mpA*nqA; i++){
        temp_s[i] = 0;
        temp_v[i] = 0;
    }
    
    for (int i = 0; i<Nx*Ny; i++){
        temp_s2[i] = 0;
        temp_v2[i] = 0;
    }
    
    for (int i =0 ; i<nprocs; i++){
        displs[i] = 0;
        for (int j = 0; j<i; j++){
            displs[i] += sizes[j];
        }
    }
    
    
    ////////////////////////////////////////////
    // Define psi and w boundary arrays. First set of arrays will represent the point at the edges of local arrays, while the second set of arrays will the represent the points outside of a given local array that need to be received from another process
    
    s_north = new double[nqA];
    s_east = new double[mpA];
    s_south = new double[nqA];
    s_west = new double[mpA];
    
    v_north = new double[nqA];
    v_east = new double[mpA];
    v_south = new double[nqA];
    v_west = new double[mpA];
    //
    //
    s_boundary_north = new double[nqA];
    s_boundary_east = new double[mpA];
    s_boundary_south = new double[nqA];
    s_boundary_west = new double[mpA];
    
    v_boundary_north = new double[nqA];
    v_boundary_east = new double[mpA];
    v_boundary_south = new double[nqA];
    v_boundary_west = new double[mpA];
    //
    //
    
    // Initialise the arrays
    for (int i = 0; i<nqA; i++){
        s_north[i] = 0;
        s_south[i] = 0;
        s_boundary_north[i] = 0;
        s_boundary_south[i] = 0;
        
        v_north[i] = 0;
        v_south[i] = 0;
        v_boundary_north[i] = 0;
        v_boundary_south[i] = 0;
    }
    
    
    for (int i = 0; i<mpA; i++){
        s_east[i] = 0;
        s_west[i] = 0;
        s_boundary_east[i] = 0;
        s_boundary_west[i] = 0;
        
        v_east[i] = 0;
        v_west[i] = 0;
        v_boundary_east[i] = 0;
        v_boundary_west[i] = 0;
    }
    
    
    // Define the location of the processes around each local process. The values for each cardinal direction are by set to -1 and overwritten if the given process exists
    
    if (mycol < Py - 1){
        process_east = myrank_mpi + 1;
    }
    //
    if (mycol > 0){
        process_west = myrank_mpi - 1;
    }
    //
    if (myrow < Px - 1){
        process_south = myrank_mpi + Py;
    }
    //
    if (myrow > 0){
        process_north = myrank_mpi - Py;
    }
    
    
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void LidDrivenCavity::UpdateBoundaryVorticity(int i, int j){ // Update Boundary vorticity (to be sent to other processes)
    
    if (i == 0){
        v_boundary_north[j] = v[I*Ny + J];
    }
    //
    if (i == mpA - 1){
        v_boundary_south[j] = v[I*Ny + J];
    }
    //
    if (j == 0){
        v_boundary_west[i] = v[I*Ny + J];
    }
    //
    if (j == nqA - 1){
        v_boundary_east[i] = v[I*Ny + J];
    }
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void LidDrivenCavity::CommunicateVorticity(){   // Use MPI to communicate the boundary vorticity values.
    
    // If a process at the boundary exists send the boundary information
    if (process_east > -1){
        MPI_Send(v_boundary_east,mpA,MPI_DOUBLE,process_east,0,MPI_COMM_WORLD);
    }
    //
    if (process_south > -1){
        MPI_Send(v_boundary_south,nqA,MPI_DOUBLE,process_south,0,MPI_COMM_WORLD);
    }
    //
    if (process_west > -1){
        MPI_Send(v_boundary_west,mpA,MPI_DOUBLE,process_west,0,MPI_COMM_WORLD);
    }
    //
    if (process_north > -1){
        MPI_Send(v_boundary_north,nqA,MPI_DOUBLE,process_north,0,MPI_COMM_WORLD);
    }
    
    
    // If a process at the boundary exists receive the boundary information and update the full voritcity array
    
    if (process_west > -1){
        MPI_Recv(v_west, mpA, MPI_DOUBLE, process_west, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        for (int i = 0; i<mpA; i++){
            v[start - 1 + i*Ny] = v_west[i];
        }
    }
    //
    //
    if (process_north > -1){
        MPI_Recv(v_north, nqA, MPI_DOUBLE, process_north, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        for (int i = 0; i<nqA; i++){
            v[start - Ny + i] = v_north[i];
        }
    }
    //
    //
    if (process_east > -1){
        MPI_Recv(v_east, mpA, MPI_DOUBLE, process_east, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        for (int i = 0; i<mpA; i++){
            v[start + nqA + i*Ny] = v_east[i];
        }
    }
    //
    //
    if (process_south > -1){
        MPI_Recv(v_south, nqA, MPI_DOUBLE, process_south, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        for (int i = 0; i<nqA; i++){
            v[start + Ny*mpA + i] = v_south[i];
        }
    }
    
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void LidDrivenCavity::ComputeVorticity(){
    
    // MARK: Calculate the boundary values of the vorticity vector
    
    for (int i = 0; i < mpA; i++) { // local row
        gpos = start + i*Ny;        // Global position
        
        I = gpos / Ny;              // Global row
        
        if (I == 0 || I == Nx -1){ // If the current global position is on the corner, continue to next iteration
            continue;
        }
        else{
            v[I*Ny + Ny-1]= (2 / (dy*dy)) * (0 - s[I*Ny + Ny-2]) - 2*U/dy; // Vorticity at the top wall
            v[I*Ny + 0] = (2 / (dy*dy)) * (0 - s[I*Ny + 1]);               // Vorticity at the bottom wall
        }
    }
    
    
    for (int j = 0; j < nqA; j++) { // local column
        gpos =  start + j;  // Global position
        J = gpos % Ny;      // Global column
        
        if (J == 0 || J == Ny - 1){ // If the current global position is on the corner, continue to next iteration
            continue;
        }
        else{
            v[J] = (2 / (dx*dx)) * (0 - s[1*Ny + J]);                  // Vorticity at the left wall
            v[(Nx-1)*Ny + J] = (2 / (dx*dx)) * (0 - s[(Nx-2)*Ny + J]); // Vorticity at the right wall
        }
    }
    
    
    
    // MARK: Compute the interior vorticity
    
    for (int i = 0; i < mpA; i++) { // local row
        
        
        for (int j = 0; j < nqA; j++) { // local col
            
            gpos = start + i*Ny + j;    // Global position
            
            I = gpos / Ny;              // Global row
            J = gpos % Ny;              // Global col
            
            
            if (I == 0 || I == Nx -1 || J == 0 || J == Ny-1){ // If the current global position is on the boundary, update boundary conditions and continue to the next iteration
                UpdateBoundaryVorticity(i,j);
                
                continue;
            }
            
            
            // Compute the inner vorticity
            v[I*Ny + J] = (2*s[I*Ny + J] - s[(I-1)*Ny + J] - s[(I+1)*Ny + J]) * (1/(dx*dx))
            + (2*s[I*Ny + J] - s[I*Ny + J - 1] - s[I*Ny + J + 1]) * (1/(dy*dy));
            
            
            // Update boundary conditions:
            UpdateBoundaryVorticity(i,j);
        }
    }
    
    // Communicate the Boundary Vorticity to neighbour processes
    CommunicateVorticity();
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void LidDrivenCavity::Integrate(){
    
    // MARK: Compute the interior vorticity at next time step
    
    for (int i = 0; i < mpA; i++) { // local row
        
        for (int j = 0; j < nqA; j++) { // local col
            
            gpos = start + i*Ny + j; //global position
            
            I = gpos / Ny;// global row
            J = gpos % Ny;// global col
            
            
            if (I == 0 || I == Nx - 1 || J == 0 || J == Ny - 1){ // If the current position is at the boundary, continue
                continue;
            }
            
            // Compute the vorticity at the next time step and assign it to temp_v array
            temp_v[i*nqA + j] = v[I*Ny + J] + dt*(
                                                  (1/Re)*((-2*v[I*Ny + J] + v[(I-1)*Ny + J] + v[(I+1)*Ny + J]) * (1/(dx*dx))
                                                          + (-2*v[I*Ny + J] + v[I*Ny + J - 1] + v[I*Ny + J + 1]) * (1/(dy*dy)))
                                                  - (s[I*Ny + J+1] - s[I*Ny + J-1])*(v[(I+1)*Ny + J] - v[(I-1)*Ny + J])*(1/(4*dx*dy))
                                                  + (s[(I+1)*Ny + J] - s[(I-1)*Ny + J])*(v[I*Ny + J+1] - v[I*Ny + J-1])*(1/(4*dx*dy))
                                                  );
        }
        
    }
    
    
    
    
    // Convert the temporary array back into the general array.
    for (int i = 0; i< mpA; i++){
        for (int j = 0; j<nqA; j++){
            gpos = start + i*Ny + j;
            
            I = gpos / Ny;// global row
            J = gpos % Ny;// global col
            
            
            if (I == 0 || I == Nx - 1 || J == 0 || J == Ny - 1){
                continue;
            }
            
            v[I*Ny + J] = temp_v[i*nqA + j];
            
        }
    }
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void LidDrivenCavity::Gather(){ // Function to gather the final streamfunction array onto the root process
    
    // Define local arrays of streamfunction and vorticity
    for (int i = 0; i < mpA; i++) { // local row
        for (int j = 0; j < nqA; j++) { // local col
            gpos = start + i*Ny + j;
            
            I = gpos / Ny;// global row
            J = gpos % Ny;// global col
            
            temp_s[i*nqA + j] = s[I*Ny + J];
            temp_v[i*nqA + j] = v[I*Ny + J];
        }
    }
    
    // Communivate the local vectors to the root process
    MPI_Gatherv(temp_s,mpA*nqA,MPI_DOUBLE,temp_s2,sizes,displs,MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Gatherv(temp_v,mpA*nqA,MPI_DOUBLE,temp_v2,sizes,displs,MPI_DOUBLE,0, MPI_COMM_WORLD);
    
    
    // On the root process, reassign the local arrays back to the global array, which will be output to file.
    if (myrank_mpi == 0){
        int loc = 0;
        int nr;
        int nc;
        for (int proc = 0; proc<nprocs; proc += 1){
            nr = locations_0[proc*2]; // number of rows of the current process
            nc = locations_0[proc*2 + 1]; // number of columns of the current process
            
            for (int i = 0; i<nr; i+=1){
                for (int j = 0; j<nc; j+=1){
                    s[starts_0[proc] + i*Ny + j] = temp_s2[loc + j + i*nc];
                    v[starts_0[proc] + i*Ny + j] = temp_v2[loc + j + i*nc];
                }
            }
            loc += nr*nc;
        }
        
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void LidDrivenCavity::Output(){
    
    // MARK: Output the result to file to allow for plotting
    
    if (myrank_mpi == 0){
        
        // Output streamfunction
        
        ofstream sOut("stream_function.csv", ios::out | ios::trunc);
        
        for (int i = 0; i<Ny*Nx-1; i++){
            sOut.precision(15);
            sOut << s[i] << endl;
        }
        
        sOut << s[Ny*Nx-1] << endl;
        
        sOut.close();
        
        /////////////////////////////////
        
        // Output vorticity
        
        ofstream vOut("vorticity.csv", ios::out | ios::trunc);
        
        for (int i = 0; i<Ny*Nx-1; i++){
            vOut.precision(15);
            vOut << v[i] << endl;
        }
        
        vOut << v[Ny*Nx-1] << endl;
        
        vOut.close();
        
        cout << "------------------" << endl;
    }
    
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// DEFINE A SERIAL CALCULATION METHOD, THIS WILL BE USED WHEN THE NUMBER OF PROCESSES IS EQUAL TO 1

void LidDrivenCavity::ComputeVorticitySerial(){
    
    mx = Nx - 2;
    my = Ny - 2;
    
    // MARK: Calculate the boundary values of the vorticity vector
    for (int i = 1; i <= mx; i++){
        v[i*Ny + Ny-1] = (2 / (dy*dy)) * (0 - s[i*Ny + Ny-2]) - 2*U/dy; // Vorticity at the top wall
        v[i*Ny + 0] = (2 / (dy*dy)) * (0 - s[i*Ny + 1]); // Vorticity at the bottom wall
    }
    
    for (int j = 1; j <= my; j++){
        v[j] = (2 / (dx*dx)) * (0 - s[1*Ny + j]); // Vorticity at the left wall
        v[(Nx-1)*Ny + j] = (2 / (dx*dx)) * (0 - s[(Nx-2)*Ny + j]); // Vorticity at the right wall
        
    }
    
    
    // MARK: Compute the interior vorticity
    
    for (int i = 1; i <= mx; i++){
        for (int j = 1; j <= my; j++){
            v[i*Ny + j] = (2*s[i*Ny + j] - s[(i-1)*Ny + j] - s[(i+1)*Ny + j]) * (1/(dx*dx))
            + (2*s[i*Ny + j] - s[i*Ny + j - 1] - s[i*Ny + j + 1]) * (1/(dy*dy));
        }
    }
    
}


void LidDrivenCavity::IntegrateSerial(){
    // MARK: Compute the interior vorticity at next time step
    for (int i = 1; i <= mx; i++){
        for (int j = 1; j <= my; j++){
            p[(i-1)*my + j-1] = v[i*Ny + j] + dt*( (1/Re)*((-2*v[i*Ny + j] + v[(i-1)*Ny + j] + v[(i+1)*Ny + j]) * (1/(dx*dx))
                                                           + (-2*v[i*Ny + j] + v[i*Ny + j - 1] + v[i*Ny + j + 1]) * (1/(dy*dy)))
                                                  - (s[i*Ny + j+1] - s[i*Ny + j-1])*(v[(i+1)*Ny + j] - v[(i-1)*Ny + j])*(1/(4*dx*dy))
                                                  + (s[(i+1)*Ny + j] - s[(i-1)*Ny + j])*(v[i*Ny + j+1] - v[i*Ny + j-1])*(1/(4*dx*dy))
                                                  );
        }
    }
    
    
}




////////////////////////////////////////////////////////////////////////////////////////////////////

void LidDrivenCavity::Solve(){
    
    // MARK: Create an instance of the Poisson Solver class
    PoissonSolver solver (Nx, Ny, v, s, alpha, beta, gamma, nu, mpA, nqA, s_north, s_east, s_south, s_west, s_boundary_north, s_boundary_east, s_boundary_south, s_boundary_west, process_east, process_south, process_west, process_north, start, myrank_mpi, nprocs_mpi);
    
    
    if (nprocs_mpi == 1){
        p = new double[(Nx-2)*(Ny-2)];
        
        solver.InitialiseSerial();
        
        // MARK: Iterate over the calculated number of time steps
        for (int t = 0; t < nt; t++){
            
            ComputeVorticitySerial(); // Compute vorticity at the current time step
            IntegrateSerial();        // Compute interior vorticity at the next time step
            
            solver.Serial(p); // Solve the Poisson equation using the serial method
            
        }
        delete p;
    }
    
    else{
        // solver.InitialiseCG();// Function to initialise the Conjugate Gradient function
        
        // MARK: Iterate over the calculated number of time steps
        for (int t = 0; t < nt; t++){
            
            ComputeVorticity(); // Compute vorticity at the current time step
            Integrate();        // Compute interior vorticity at the next time step
            
            solver.Jacobi(5000); // Solve the Poisson equation using the Jacobi method (500 iterations at each time step)
          
            // solver.ConjugateGradient();
        }
        
        // Gather all local vorticity and streamfunction arrays onto the root process
        Gather();
    }
    
    
    // Output the final results to file
    Output();
    
}

