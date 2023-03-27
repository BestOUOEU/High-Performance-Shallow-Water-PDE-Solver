#ifndef ShallowWater_h
#define ShallowWater_h
/**
* Header file for class ShallowWater
* @class ShallowWater
* @brief class ShallowWater to solve shallow water PDEs
* @param swater is the ShallowWater object identifier used for different member functions
* @returns will be the the matrix u, v or h depending on user input
*/
#include <iostream>


class ShallowWater{
    public:        

        /// 1. default constructor
        ShallowWater();

        /**
        * @brief 2. User defined constructor, Creates a ShallowWater with the user given parameters 
        * (Parameters given prior to execution, passed to main() , then main passes to constructor)
        * @param Nx number of grid points in x for the matrix
        * @param Ny number of grid points in y for the matrix
        * @param ic initial condition case for h picked
        * @param T number of grid points in y for the matrix
        * @param Ny number of grid points in y for the matrix
        */
        ShallowWater(int& Nx, int& Ny, int& ic, double& T, double& dt);

        /// 3. default destructor
        ~ShallowWater(); 

        /**
         * @brief 4. Sets initial conditions u(x,y,0)=0, v(x,y,0)=0 and depending on input, h(x,y,0)=g(x,y) on the PDE. 
         * Cache locality much more efficient when populating values in Ny direction, the number of rows in a column, 
         * compared to Nx, the number of columns
         * There are no arguments for calling this function, but it takes from Shallow Water object:
         * @param   ic      Choice of initial condition for h
         * @param   Nx      Matrix dimension in x, number of columns
         * @param   Ny      Matrix dimension in y, number of rows
         * @param   dx      Grid spacing in x
         * @param   dy      Grid spacing in y
         * @param   u       Pointer to u velocity matrix of size Nx*Ny
         * @param   v       Pointer to v velocity matrix of size Nx*Ny
         * @param   h       Pointer to h wave height matrix of size Nx*Ny
         */
        void SetInitialConditions();

        
        /**
        * @brief 5. Time integrate the shallow water PDE with the grids Nx Ny T and dt given by the user
        * when constructing the ShallowWater object.
        * Spacial integration with 6th order finite differencing while time integration with 4th order
        * Runge-Kutta explicit
        * @param   BLAS      Choice of for loop or BLAS method (0 if for loop, 1 if BLAS)
        */
        void TimeIntegrate(int BLAS);



        /**
         * @brief 6. File Output x, y, u, v, h for each grid point at final time
         * Takes from Shallow Water object:
         * @param   Nx       Number of Grid Points in x
         * @param   Ny       Number of Grid Points in y
         * @param   u        Pointer to the final time u velocity matrix of size Nx*Ny
         * @param   v        Pointer to the final time v velocity matrix of size Nx*Ny
         * @param   h        Pointer to the final time h wave height matrix of size Nx*Ny
         */
        void FileOutput();

    private:

        /// @brief  Privately held data to avoid misuse outside the class
        // If any needs changing, doesnt have to call setparam in main, just change here for testing.
        // Otherwise if taking stuff from command prompt, need a proper constructor
        double T = 5.0; // total time to propagate
        int ic = 1; // index to choose initial conditions for h

        // Steps between grid points
        double    dx   = 1.0;          // x spacial step 0.1
        double    dy   = 1.0;          // y spacial step 0.1
        double    dt   = 0.1;          // t time step, set as 0.1 (can be varied for numerical accuracy and speed)

        // Number of grid points
        int Nx = 100;                      
        int Ny = 100;
        int Nt = T/dt; // Number of forward time propagation points (integer), trailing decimals truncated
        int ldh = 7;         // Leading dimension of 6th order forward differencing scheme matrix
        

        // u,v,h matrices
        double* u = new double[Nx*Ny];     // x velocity of points - Dimension: N*Nfor all grid points
        double* v = new double[Nx*Ny];     // y velocity
        double* h = new double[Nx*Ny];     // h height of wave

        // extra set of u,v,h matrices for computing spacial derivatives
        // at different runge kutta points k2 k3 k4 without overwriting the original
        // the original u v and h matrices needs to be kept for += operation
        double* u2 = new double[Nx*Ny];     
        double* v2 = new double[Nx*Ny];     
        double* h2 = new double[Nx*Ny];     

        // Partial spacial derivatives matrices for different
        // runge kutta points, can simply overwrite as its is only an
        // used to bridge between u, v, h and d(u,v,h)/dt
        double* dudx = new double[Nx*Ny];    
        double* dvdx = new double[Nx*Ny];   
        double* dhdx = new double[Nx*Ny];    
        double* dudy = new double[Nx*Ny];    
        double* dvdy = new double[Nx*Ny];  
        double* dhdy = new double[Nx*Ny]; 

        // d/dt time derivative Matrices for Runge Kutta point k1
        double* dudt = new double[Nx*Ny];   
        double* dvdt = new double[Nx*Ny];    
        double* dhdt = new double[Nx*Ny]; 

        // d/dt time derivative Matrices for Runge Kutta points k2 k3 k4
        double* k2u = new double[Nx*Ny];   
        double* k2v = new double[Nx*Ny];    
        double* k2h = new double[Nx*Ny];
        double* k3u = new double[Nx*Ny];   
        double* k3v = new double[Nx*Ny];    
        double* k3h = new double[Nx*Ny];
        double* k4u = new double[Nx*Ny];   
        double* k4v = new double[Nx*Ny];    
        double* k4h = new double[Nx*Ny];


        // Banded matrices for sixth order finite differencing scheme
        double* SixOrd1 = new double[Nx*Nx];      // Sixth order banded matrix in x direction
        double* SixOrd2 = new double[Ny*Ny];      // Sixth order banded matrix in y direction
        
        // Private member functions

        /**
         * @brief 1. Calculates the values of u, v ,h at next time step using 4th order Runge-Kutta explicit scheme
         * There are no arguments for calling this function, but it takes from Shallow Water object:
         * @param   Nx       Number of Grid Points in x
         * @param   Ny       Number of Grid Points in y
         * @param   dt       Time step size Delta(t)
         * @param   u        Pointer to u velocity matrix of size Nx*Ny
         * @param   v        Pointer to v velocity matrix of size Nx*Ny
         * @param   h        Pointer to h wave height matrix of size Nx*Ny
         * @param   dudt     Pointer to the matrix of du/dt of size Nx*Ny
         * @param   dvdt     Pointer to the matrix of dv/dt of size Nx*Ny
         * @param   dhdt     Pointer to the matrix of dh/dt of size Nx*Ny
         * @param   k2u      Pointer to the matrix of runge kutta term k2u of size Nx*Ny
         * @param   k2v      Pointer to the matrix of runge kutta term k2v of size Nx*Ny
         * @param   k2h      Pointer to the matrix of runge kutta term k2h of size Nx*Ny
         * @param   k3u      Pointer to the matrix of runge kutta term k3u of size Nx*Ny
         * @param   k3v      Pointer to the matrix of runge kutta term k3v of size Nx*Ny
         * @param   k3h      Pointer to the matrix of runge kutta term k3h of size Nx*Ny
         * @param   k4u      Pointer to the matrix of runge kutta term k4u of size Nx*Ny
         * @param   k4v      Pointer to the matrix of runge kutta term k4v of size Nx*Ny
         * @param   k4h      Pointer to the matrix of runge kutta term k4h of size Nx*Ny
         */
        void FourthOrdRK();



        /**
         * @brief 2. Populate Gradients at grid points not affected by periodic boundary condition assumptions.
         * @param   dxx       Grid Spacing (could be x or y, although its named x)
         * @param   Nxx       Number of Grid Points in x
         * @param   Nyy       Number of Grid Points in y
         * @param   incy      The increment of index in the finite differencing scheme
         * (if populating derivatives in dx direction this shouldbe Ny, in y direction this should be 1)
         * @param   uu        Pointer to the matrix of function values of size Nx*Ny (could be u, v, h)
         * @param   dudxx     (Passed by reference, overwritten) Pointer to the matrix of function value derivatives 
         *  of size Nx*Ny (Numerator could be du, dv, dh, while denominator could be dx dy, 
         * matching that of dx parameter)
         * 
         */
        void GradElseWhere(double dxx, int Nxx, int Nyy, int incy, double* uu, double* &dudxx);



        /**
         * @brief 3. Calculate Gradients at boundary grid points affected by the periodic boundary conditions
         * @param   dxx       Grid Spacing (could be dx or dy, although its named x)
         * @param   Nxx       Number of Grid Points in x
         * @param   Nyy       Number of Grid Points in y
         * @param   incy      The increment of index in the finite differencing scheme
         * (if populating derivatives in dx direction this shouldbe Ny, in y direction this should be 1)
         * @param   uu        Pointer to the matrix of function values of size Nx*Ny (could be u, v, h)
         * @param   dudxx     (Passed by reference, overwritten) Pointer to the matrix of function value derivatives
         * of size Nx*Ny (Numerator could be du, dv, dh, while denominator could be dx dy, 
         * matching that of dx parameter)
         * 
         */
        void GradAtBound(double dxx, int Nxx, int Nyy, int incy, double* uu, double* &dudxx);
           


        /**
         * @brief 4. Calculates time derivatives from the spacial values and dervatives using Shallow Water PDE
         * There are no arguments for calling this function, but it takes from Shallow Water object:
         * @param   Nx       Number of Grid Points in x
         * @param   Ny       Number of Grid Points in y
         * @param   u        Pointer to u velocity matrix of size Nx*Ny
         * @param   v        Pointer to v velocity matrix of size Nx*Ny
         * @param   h        Pointer to h wave height matrix of size Nx*Ny
         * @param   dudx     Pointer to the matrix of du/dx of size Nx*Ny
         * @param   dudy     Pointer to the matrix of du/dy of size Nx*Ny
         * @param   dvdx     Pointer to the matrix of dv/dx of size Nx*Ny
         * @param   dvdy     Pointer to the matrix of dv/dy of size Nx*Ny
         * @param   dhdx     Pointer to the matrix of dh/dx of size Nx*Ny
         * @param   dhdy     Pointer to the matrix of dh/dy of size Nx*Ny
         * @param   dudt     (Passed by reference, overwritten) Pointer to the matrix of du/dt of size Nx*Ny (k1u)
         * @param   dvdt     (Passed by reference, overwritten) Pointer to the matrix of dv/dt of size Nx*Ny (k1v)
         * @param   dhdt     (Passed by reference, overwritten) Pointer to the matrix of dh/dt of size Nx*Ny (k1h)
         * @param   k2u      (Passed by reference, overwritten) Pointer to the matrix of runge kutta term k2u of size Nx*Ny
         * @param   k2v      (Passed by reference, overwritten) Pointer to the matrix of runge kutta term k2v of size Nx*Ny
         * @param   k2h      (Passed by reference, overwritten) Pointer to the matrix of runge kutta term k2h of size Nx*Ny
         * @param   k3u      (Passed by reference, overwritten) Pointer to the matrix of runge kutta term k3u of size Nx*Ny
         * @param   k3v      (Passed by reference, overwritten) Pointer to the matrix of runge kutta term k3v of size Nx*Ny
         * @param   k3h      (Passed by reference, overwritten) Pointer to the matrix of runge kutta term k3h of size Nx*Ny
         * @param   k4u      (Passed by reference, overwritten) Pointer to the matrix of runge kutta term k4u of size Nx*Ny
         * @param   k4v      (Passed by reference, overwritten) Pointer to the matrix of runge kutta term k4v of size Nx*Ny
         * @param   k4h      (Passed by reference, overwritten) Pointer to the matrix of runge kutta term k4h of size Nx*Ny
         */
        void Calcddt();



        /**
         * @brief 5. Populate Sixth order finite differencing general matrices SixOrd1 and SixOrd2
         * There are no arguments for calling this function, but it takes from Shallow Water object:
         * @param   SixOrd1       Pointer to the matrix storage of size Nx*Nx (for computing x derivatives)
         * @param   SixOrd2      Pointer to the matrix storage of size Ny*Ny (for computing y derivatives)
         * @param   dx           Grid spacing
         * @param   Nx           Other dimension of matrix, x
         * @param   Ny           Other dimension of matrix, y
         */
        void FillSixOrdMatrix();

        /**
         * @brief 6. Calculate Spacial gradients of u, v and h wrt. x and y using BLAS
         * @param   uu           u matrix used to calculate spacial gradient (could be u or u2)
         * @param   vv           v matrix used to calculate spacial gradient (could be v or v2)
         * @param   hh           h matrix used to calculate spacial gradient (could be h or h2)
         */
        void BLASGrad(double* uu, double* vv, double* hh);

        /**
         * @brief 7. Calculate the time derivatives k1=d/dt k2 k3 k4 using BLAS
         */
        void BLASCalcddt();

        /**
         * @brief 8. Runge Kutta forward time propagation of u,v,h using BLAS
         */
        void BLASRK();

};

#endif