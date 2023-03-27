/**
 * @file ShalloWaterMain.cpp
 * Solves Shallow Water PDE, all matrices are by default column major
 *  // steps: 
 *  // 1. Start with initial condition u,v,h(x,y,0)
 *  // 2. Propagate RK and SixOrd matrix in both x, y directions (2 matrix of SixOrd, 
 *  // Runge kutta is not matrix cuz its only forward time propagation)
 *  // 3. calculate nabla u, v, h with SixOrd spacial propagation (should be parallelised in 3 paths, further parallelised as smaller
 *  // sub matrices)
 *  // 4. calculate f=d(u,v,h)/dt using u,v,h and nabla, u, v, h of the current time step (duv/dx etc use product rule)
 *  // 5. use this to forward time propagate at each grid pt with runge kutta
 *     
 *  //CACHE LOCALITY
 *  //Think in vectors of columns of dimension Ny*Nx (column major), not matrices
 *  //
 *  // 1. All Matrices Are Designed To Be Column Majors
 *  //
 *  // 2. All Following Loops Has Been Designed To Encourage Cache Locality, Reduced Arithmetrics And Precomputation Of Values
 *  // Where Possible.
 *  //
 *  // 3. Column majors stores the values of a matrix column by column in the memory. In order to encourage cache locality,
 *  // we must design the for loops such that We operate along the rows inside each column before we move onto the next row
 *  // Therefore, the inner most loop with the fastest changing index must be such that 0<=j<Ny, 
 *  // where j is the fastest changing index and Ny is the number of rows inside a column
 *  //
 *  // 4. We always rearrange the loops according to 3, working along columns. Also, when it is possible to utilise cache locality,
 *  // we reorder the loop such that the fastest changing index is always at the inner most loop,
 *  // unless it significantly reduces the amount of possible precompuations.
 *  //
 *  // 5. Therefore, j, which is number of rows, which is also Ny in the y direction should be the fastest changing index
 *  // most of the times. An example of loop populating column major matrix M should look like:
 *  //
 *  //  for(int i=0;i<Nx;i++){
 *  //      for(int j=0;i<Ny;j++){
 *  //          M[i*Ny+j]=0;
 *  //      }
 *  //  }    
 */
 
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <stdexcept>

#include "ShallowWater.h"
#include <omp.h>
#include <boost/program_options.hpp>

using namespace std;

namespace po = boost::program_options;




/**
 * @brief Main function to solve the Shallow water PDEs using finite differencing of spacial grids and runge kutta forward time propagation
 * 
 */
int main(int argc, char* argv[]) {

    // Use boost to take in and parse coomandline arguments
    // total of 5 command line arguments
    po::options_description opts("Available options.");
    opts.add_options()
    ("Nx", po::value<int>()->default_value(100), "Integer number of grid points in x, default 100")
    ("Ny", po::value<int>()->default_value(100), "Integer number of grid points in y, default 100")
    ("ic", po::value<int>()->default_value(4), "Test case, ranging from 1-4, default 1")
    ("T", po::value<double>()->default_value(3), "Total time to propagate PDE, default 80.0")
    ("dt", po::value<double>()->default_value(0.1), "Step size for time, default 0.1")
    ("case1", po::value<int>()->default_value(0), "Case1: For loop method if 0, BLAS if 1, default 0")
    ("help", "Print help message.");
    

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    // print options instead if command line argument was help
    if (vm.count("help")) {
        std::cout << opts << std::endl;
    }

    // parsing of arguments
    int Nx = vm["Nx"].as<int>();
    int Ny = vm["Ny"].as<int>();
    int ic = vm["ic"].as<int>();
    double T = vm["T"].as<double>();
    double dt = vm["dt"].as<double>();
    int case1 = vm["case1"].as<int>();

    // Modify omp settings: allow nested parallelisation, change number of threads
    omp_set_nested(1);
    
    ShallowWater solution1(Nx, Ny, ic, T, dt);
    // set initial conditions
    solution1.SetInitialConditions();

    // time integrate (0 if for loop method) (1 if BLAS)
    solution1.TimeIntegrate(case1);

    // blocks to catch exceptions
    try {
    }
    catch (const std::logic_error& e) {
        std::cout << "An error occured: " << e.what() << std::endl;
    }
    
    cout<<"End of function."<<endl;
    

    return 0;
}
