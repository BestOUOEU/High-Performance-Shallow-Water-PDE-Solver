/* CPP file
 * ShallowWater class function implementations.
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <omp.h>



#include "ShallowWater.h"
using namespace std;

#define F77NAME(x) x##_
extern "C"{
    //BLAS Prototypes
    
    void F77NAME(daxpy) (const int &n, const double &alpha, const double *x, const int &incx, const double *y, const int &incy);
    
    void F77NAME(dcopy) (const int &n, const double *x, const int &incx, double *y, const int &incy);
    
    void F77NAME(dgemv) (const char &trans, const int &m, const int &n, const double &alpha, const double *a, 
    const int &lda, const double *x, const int &incx, const double &beta, double *y, const int &incy);
}


// Default Constructor
ShallowWater::ShallowWater(){
    // Nothing to do
}



// User Defined Constructor
ShallowWater::ShallowWater(int& Nx1, int& Ny1, int& ic1, double& T1, double& dt1){

    // Throw error if invalid values
    if (Nx1<=0 || Ny1<=0 || ic1<=0 || T1<=0 || dt1<=0 ) {
        throw std::logic_error("All parameters given prior to execution must be > 0.");
    }

    // Throw error if dt>T
    if (dt1>T1) {
        throw std::logic_error("Total time must be bigger than time step");
    }

    // Throw error if case doesnt exist
    if (ic1>4) {
        throw std::logic_error("Case doesnt exist");
    }

    // Assign the given values to variables in private, other parameters will be initialised as given by default
    Nx = Nx1;
    Ny = Ny1;
    ic = ic1;
    T = T1;
    dt = dt1;
    Nt = T1/dt1;
}



// Default Destructor
ShallowWater::~ShallowWater(){
    // Clean up the dynamically allocated arrays
    // Other parameters will be destroyed automatically when the object goes out of scope, no need to worry
    delete[] SixOrd1;
    delete[] SixOrd2;
    delete[] u;
    delete[] v;
    delete[] h;
    delete[] dudx;
    delete[] dvdx;
    delete[] dhdx;
    delete[] dudy;
    delete[] dvdy;
    delete[] dhdy;
    delete[] dudt;
    delete[] dvdt;
    delete[] dhdt;
    delete[] u2;
    delete[] v2;
    delete[] h2;
    delete[] k2u;
    delete[] k2v;
    delete[] k2h;
    delete[] k3u;
    delete[] k3v;
    delete[] k3h;
    delete[] k4u;
    delete[] k4v;
    delete[] k4h;
}



void ShallowWater::SetInitialConditions(){
    
    for(int i=0;i<Nx*Ny;i++){ 
        // populate zero initial values of u and v velocities in a single loop of N*N iterations
        u[i] = 0.0;
        v[i] = 0.0;
    }
    
    
    double by25 = 1.0/25.0;// precompute 1/25 for repeated use inside switch cases loops
    cout<<"ic:"<<ic<<endl;

    
    // First initialise useful indices outside of for loop
    int Nyti;
    // Next populate intial condition for h wave height, depending on the input chosen by user
    switch (ic) {
        case 1:

            for(int i=0;i<Nx;i++){
                Nyti = Ny*i; // Precompute Ny*i array indices
                double exp1 = (dx*i-50.0); // exp1=(x-50)
                for(int j=0;j<Ny;j++){
                    // Ny is the number of rows in a column, j is inner loop, 
                    // so populate all rows inside a column first, before populating the next column
                    h[Nyti+j] = 10.0+exp(-(exp1*exp1)*by25); // h=10+exp(-(x-50)^2/25)
                }
            }
                
            break;

        case 2:
        
        for(int i=0;i<Nx;i++){
                Nyti = Ny*i; // Precompute Ny*i array indices
                // Special case, could swap order of the loop such that we encourage recomputations, but
                // cache locality is more important here for populating values it seems
                for(int j=0;j<Ny;j++){
                    double exp1 = (dy*j-50.0); // exp1=(y-50)^2
                    h[Nyti+j] = 10.0+exp(-(exp1*exp1)*by25); // h=10+exp(-(y-50)^2/25)
                }
            }
            break;
            
        case 3:

            for(int i=0;i<Nx;i++){
                Nyti=Ny*i; // Precompute Ny*i array indices
                double exp1 = (dx*i-50.0); // exp1=(x-50)
                double exp2 = (exp1*exp1); // exp2=(x-50)^2

                for(int j=0;j<Ny;j++){
                    double exp3 = (dy*j-50.0); // exp3=(y-50)
                    h[Nyti+j] = 10.0+exp(-(exp2+exp3*exp3)*by25); // h=10+exp(-((x-50)^2+(y-50)^2)/25)
                }
            }
            break;

        case 4:
            //h0(x, y) = H(x, y) + exp(-((x - 25)2 + (y - 25)2)/25) + exp(-((x - 75)2 + (y - 75)2)/25)

            for(int i=0;i<Nx;i++){
                Nyti=Ny*i; // Precompute Ny*i array indices
                double exp1 = (dx*i-25.0); // exp1=(x-25)
                double exp2 = (dx*i-75.0); // exp2=(x-75)
                double exp3 = exp1*exp1; // exp3=(x-25)^2
                double exp4 = exp2*exp2; // exp4=(x-75)^2
            
                for(int j=0;j<Ny;j++){
                    double exp5 = (dy*j-25.0); // exp5=(y-25)
                    double exp6 = (dy*j-75.0); // exp6=(y-75)
                    h[Nyti+j] = 10.0+exp(-(exp3+exp5*exp5)*by25)+exp(-(exp4+exp6*exp6)*by25); //h0(x, y) = 10 + exp(-((x - 25)2 + (y - 25)2)/25) + exp(-((x - 75)2 + (y - 75)2)/25)
                }
            }
            break;
        
        
        default:
            cout<<"invalid choice of initial condition: no case for ic !"<<endl;
        
    }
    cout<<"Setinitcond ran"<<endl;
    
}



void ShallowWater::FileOutput(){

    // initialise a file for output
    ofstream vOut("output.txt", ios::out | ios::trunc);
    // set precision and width
    vOut.precision(5);
    vOut.width(7);

    // initialise useful indices used in for loop
    int Nytjpi;

    // Normally could output in a single for loop, however the format required is row major, 
    // the data is in column major
    for(int i=0;i<Ny;i++){
        for(int j=0;j<Nx;j++){
            // Precompute indices
            Nytjpi=i+j*Ny;

            // Output data in required format to file
            vOut << j << " " << i << " " << u[Nytjpi] << " " << v[Nytjpi] << " " << h[Nytjpi] <<endl;
        }
        vOut<<endl;
    }
    vOut.close();

}

        

void ShallowWater::TimeIntegrate(int BLAS){

    if(BLAS==0){ 
        cout<<"Nt:"<<Nt<<endl;
        // Repeat until final time

        // This for loop cannot be parallelised since it depends on previous timesteps, however, the functions being called have been paralleised
        for(int i=0;i<Nt;i++){
            
//          Printing statements commented out
//          cout<<"TimeIntegrate t: "<<dt*(i+1)<<endl;
//          cout<<"TimeIntegrate u: "<<u[5000]<<endl;
//          cout<<"TimeIntegrate h: "<<h[5000]<<endl;
            //time derivatives
            ShallowWater::Calcddt();
            //forward time propagation of u,v,h
            ShallowWater::FourthOrdRK();
//          Printing statements commented out
//
//          cout<<"TimeIntegrate dudx: "<<dudx[5000]<<endl;
//          cout<<"TimeIntegrate dudy: "<<dudy[5000]<<endl;
//          cout<<"TimeIntegrate dhdx: "<<dhdx[5000]<<endl;
//          cout<<"TimeIntegrate dhdy: "<<dhdy[5000]<<endl;            
//          cout<<"TimeIntegrate dudt: "<<dudt[5000]<<endl;
//          cout<<"TimeIntegrate dvdt: "<<dvdt[5000]<<endl;
//          cout<<"TimeIntegrate dhdt: "<<dhdt[5000]<<endl;
        }
    }
    else if(BLAS==1){
        //Populate sixth order finite differencing matrices
        ShallowWater::FillSixOrdMatrix();

        cout<<"Nt:"<<Nt<<endl;
        // Repeat until final time
        for(int i=0;i<Nt;i++){
//          Printing statements commented out
//          cout<<"TimeIntegrate t: "<<dt*(i+1)<<endl;
//          cout<<"TimeIntegrate h: "<<h[5000]<<endl;
            //time derivatives
            ShallowWater::BLASCalcddt();
            //forward time propagation of u,v,h
            ShallowWater::BLASRK();
//          Printing statements commented out
//          cout<<"TimeIntegrate dudx: "<<dudx[5000]<<endl;
//          cout<<"TimeIntegrate dudy: "<<dudy[5000]<<endl;
//          cout<<"TimeIntegrate dhdx: "<<dhdx[5000]<<endl;
//          cout<<"TimeIntegrate dhdy: "<<dhdy[5000]<<endl; 
//          cout<<"TimeIntegrate dudt: "<<dudt[5000]<<endl;
//          cout<<"TimeIntegrate dvdt: "<<dvdt[5000]<<endl;
//          cout<<"TimeIntegrate dhdt: "<<dhdt[5000]<<endl;
        }
    }
    else{
        throw std::logic_error("Invalid BLAS parameter for TimeIntegrate must be (0 for for loop method)(1 for BLAS method)"); // throw error if invalid method choice
    }



        // After finish forward time propagation, output results of u,v,h as .txt ile
        ShallowWater::FileOutput();
        cout<<"TimeIntegration complete, Check you folder for .txt file"<<endl;
}



// Private member functions below
void ShallowWater::FourthOrdRK(){
    // Precompute useful divisions
    double by6 = 1.0/6.0;
    // intialise useful indices outside of for loop
    int Nytipj;
    // initialise useful varaibles outside of loop

    // parallelised the for loop, letting omp to decide chunk, using collapse to parallelise both the outer and inner-loop
    int i = 0, j = 0; // for parallelisation
    #pragma omp parallel for\
    default(shared)         \
    schedule (static)       \
    private(i,j,Nytipj)     \
    collapse(2)
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            // Precompute index Ny*i+j
            Nytipj = Ny*i+j;

            // Overwrite Original values on the same matrix to avoid memory overuse and additional copying
            // The original entries are not used another time in the forward time propagation scheme
            u[Nytipj] += by6*(dudt[Nytipj]+k2u[Nytipj]*2.0+k3u[Nytipj]*2.0+k4u[Nytipj])*dt;
            v[Nytipj] += by6*(dvdt[Nytipj]+k2v[Nytipj]*2.0+k3v[Nytipj]*2.0+k4v[Nytipj])*dt;
            h[Nytipj] += by6*(dhdt[Nytipj]+k2h[Nytipj]*2.0+k3h[Nytipj]*2.0+k4h[Nytipj])*dt;
        }
    } 
}



void ShallowWater::GradElseWhere(double dxx, int Nxx, int Nyy, int incy, double* uu, double* &dudxx){

    // Precompute useful fractions used in du/dx for loop
    double bydx = 1.0/dxx; // 1/dx
    double bydxby60 = bydx/60.0; // pre-compute finite differencing coefficients
    double bydxby20t3 = bydx*3.0/20.0;
    double bydxby4t3 = bydx*3.0/4.0;
    
    int incy2;
    int incy3;
    int iterations1lower; // Lower bound of outer loop
    int iterations1upper; // upper bound of outer loop
    int iterations2lower; // lower bound of outer loop
    int iterations2upper; // upper bound of outer loop

    // Precompute Column Major Index locations for the differencing scheme
    // This depends on which direction the user chose to finite difference in, depends on incy.
    if (incy==Nyy){
        // case where finite differencing in x direction (at boundary), so increment in index for F-D scheme incy is Ny

        // Precompute all index increments for finite difference scheme
        // ui+1 ui+2 ui+3
        // incy is for ui+1
        incy2 = incy+incy;
        incy3 = incy2+incy;
        // ui-1 ui-2 ui-3
        // -incy is for ui-1
        // -incy2 is for ui-2
        // -incy3 is for ui-3

        // Precompute iteration properties
        iterations1lower = 3; // Lower bound of outer loop
        iterations1upper = Nxx-3; // upper bound of outer loop
        iterations2lower = 0; // lower bound of outer loop
        iterations2upper = Nyy; // upper bound of outer loop
    }
    else if(incy==1){
        // case where finite differencing in the y direction (at boundary), so increment in index for F-D scheme incy is 1

        // Precompute all index increments for finite difference scheme
        // ui+1 ui+2 ui+3
        // incy is for ui+1
        incy2 = incy+incy;
        incy3 = incy2+incy;
        // ui-1 ui-2 ui-3
        // -incy is for ui-1
        // -incy2 is for ui-2
        // -incy3 is for ui-3

        // Precompute iteration properties
        iterations1lower = 0; // Lower bound of outer loop
        iterations1upper = Nxx; // upper bound of outer loop
        iterations2lower = 3; // lower bound of outer loop
        iterations2upper = Nyy-3; // upper bound of outer loop
    }
    else{
        cout<<"incorrect value of incy, please make sure your finite differencing index increment is correct for boundary values! "<<endl;
    }
    
    // Finally, populate the gradients in a for loop (not at boundaries)
    // First initialise useful indices outside of loop
    int Nytipj;
    // Indexing and loops ordered to encourage cache locality (fastest changing index is at the inner most loop)

    // parallelised the for loop, letting omp to decide chunk, using collapse to parallelise both the outer and inner-loop
    int i = 0, j = 0; // for parallelisation
    #pragma omp parallel for\
    default(shared)         \
    schedule (static)       \
    private(i,j,Nytipj)     \
    shared(iterations1lower,iterations1upper,iterations2lower,iterations2upper,incy,incy2,incy3)\
    collapse(2)
    for(i = iterations1lower ; i < iterations1upper ; i++){
        for(j = iterations2lower; j < iterations2upper ; j++){ 
            // Precompute central point index to reduce arithmetics
            Nytipj=Ny*i+j;

            // Propulate derivatives with F-D scheme
            dudxx[Nytipj] = (+uu[Nytipj+incy]-uu[Nytipj-incy])*bydxby4t3
            +(-uu[Nytipj+incy2]+uu[Nytipj-incy2])*bydxby20t3
            +(+uu[Nytipj+incy3]-uu[Nytipj-incy3])*bydxby60;
            
        }   
    }
}
    


void ShallowWater::GradAtBound(double dxx, int Nxx, int Nyy, int incy, double* uu, double* &dudxx){
    // Function to calculate Gradients at boundary points

    // Precompute useful fractions used in for loop
    double bydx = 1.0/dxx; // 1/dx
    double bydxby60 = bydx/60.0; // pre-compute finite differencing coefficients
    double bydxby20t3 = bydx*3.0/20.0;
    double bydxby4t3 = bydx*3.0/4.0;
    
    int Nyt1; //is the column major index of the 1st row, 1st column entry
    int Nyt2; // is the column major index of the 1st row, 2nd column entry
    int Nyt3; // Column major index of the 1st row, 3rd column entry
    int Nyt4; // Column major index of the 1st row, 4th column entry
    int Nyt5; // Column major index of the 1st row, 5th column entry
    int Nyt6; // Column major index of the 1st row, 6th column entry
    // u-1 u-2 u-3 u-4 u-5 u-6
    int Nytn1; // Column major index of the 1st row, last column entry
    int Nytn2; // Column major index of the 1st row, 2nd last column entry
    int Nytn3; // Column major index of the 1st row, 3rd last column entry
    int Nytn4; // Column major index of the 1st row, 4th last column entry
    int Nytn5; // Column major index of the 1st row, 5th last column entry
    int Nytn6; // Column major index of the 1st row, 6th last column entry

    // Iteration properties
    int iterations1; // Number of iterations
    int itinc; // Iteration index increment of each iteration, if incy==Ny then this is 1, if incy==1 this is Ny
    
    // Precompute Column Major Index locations for the differencing scheme
    // This depends on which direction the user chose to finite difference in, depends on incy.
    if (incy==Nyy){
        // Case where finite differencing in x direction (at boundary), so increment in index for F-D scheme incy is Ny


        // Precompute index increments for the finite differnce scheme
        // u1 u2 u3 u4 u5 u6
        Nyt1 = 0; //is the column major index of the 1st row, 1st column entry
        Nyt2 = incy; // is the column major index of the 1st row, 2nd column entry
        Nyt3 = incy+incy; // Column major index of the 1st row, 3rd column entry
        Nyt4 = Nyt3+incy; // Column major index of the 1st row, 4th column entry
        Nyt5 = Nyt4+incy; // Column major index of the 1st row, 5th column entry
        Nyt6 = Nyt5+incy; // Column major index of the 1st row, 6th column entry
        // u-1 u-2 u-3 u-4 u-5 u-6
        Nytn1  = Nyy*(Nxx-1); // Column major index of the 1st row, last column entry
        Nytn2 = Nytn1-incy; // Column major index of the 1st row, 2nd last column entry
        Nytn3 = Nytn2-incy; // Column major index of the 1st row, 3rd last column entry
        Nytn4 = Nytn3-incy; // Column major index of the 1st row, 4th last column entry
        Nytn5 = Nytn4-incy; // Column major index of the 1st row, 5th last column entry
        Nytn6 = Nytn5-incy; // Column major index of the 1st row, 6th last column entry

        // Iteration properties
        iterations1=Nyy; // Number of iterations
        itinc=1; // Iteration index increment of each iteration, if incy==Ny then this is 1, if incy==1 this is Ny
    }
    else if(incy==1){
        //******************************************************************************************can write a function to accomodate this better
        // Case where finite differencing in the y direction (at boundary), so increment in index for F-D scheme incy is 1

        // Precompute index increments for the finite differnce scheme
        // u1 u2 u3 u4 u5 u6
        Nyt1 = 0; //is the column major index of the 1st row, 1st column entry
        Nyt2 = incy; // is the column major index of the 2nd row, 1st column entry
        Nyt3 = incy+incy; // Column major index of the 3rd row, 1st column entry
        Nyt4 = Nyt3+incy; // Column major index of the 4th row, 1st column entry
        Nyt5 = Nyt4+incy; // Column major index of the 5th row, 1st column entry
        Nyt6 = Nyt5+incy; // Column major index of the 6th row, 1st column entry
        // u-1 u-2 u-3 u-4 u-5 u-6
        Nytn1 = Nyy-1; // Column major index of the last row, 1st column entry
        Nytn2 = Nytn1-incy; // Column major index of the 2nd last row, 1st column entry
        Nytn3 = Nytn2-incy; // Column major index of the 3rd last row, 1st column entry
        Nytn4 = Nytn3-incy; // Column major index of the 4th last row, 1st column entry
        Nytn5 = Nytn4-incy; // Column major index of the 5th last row, 1st column entry
        Nytn6 = Nytn5-incy; // Column major index of the 6th last row, 1st column entry

        // Iteration properties
        iterations1=Nxx; // Number of iterations
        itinc=Nyy; // Iteration index increment of each interation, if incy==Ny then this is 1, if incy==1 this is Ny
    }
    else{
        cout<<"incorrect value of incy, please make sure your finite differencing index increment is correct for boundary values! "<<endl;
    }

    // Finally, populate the gradients at boundaries in a for loop
    // First iniitalise useful indices outside of loop
    // u1 u2 u3 u4 u5 u6
    int Nyt1i;
    int Nyt2i;
    int Nyt3i;
    int Nyt4i;
    int Nyt5i;
    int Nyt6i;
    // u-1 u-2 u-3 u-4 u-5 u-6
    int Nytn1i;
    int Nytn2i;
    int Nytn3i;
    int Nytn4i; 
    int Nytn5i; 
    int Nytn6i;
    
    int iitinc;

    // Can further split into 2 loops if to maximise cache, increasing the chance that cache hasnt bee overwritten yet 
    // parallelised the for loop, letting omp to decide chunk
    int i = 0; // for parallelisation
    #pragma omp parallel for    \
    default(shared)             \
    schedule(static)            \
    private(i,Nyt1i,Nyt2i,Nyt3i,Nyt4i,Nyt5i,Nyt6i,Nytn1i,Nytn2i,Nytn3i,Nytn4i,Nytn5i,Nytn6i,iitinc)\
    shared(Nyt1,Nyt2,Nyt3,Nyt4,Nyt5,Nyt6,Nytn1,Nytn2,Nytn3,Nytn4,Nytn5,Nytn6,itinc)
    for(i=0;i<iterations1;i++){
        // Further precompute index locations for the F-D scheme around the central point for the current iteration/position i
        iitinc=i*itinc; // Precompute common adder
        // u1 u2 u3 u4 u5 u6
        Nyt1i = Nyt1+iitinc;
        Nyt2i = Nyt2+iitinc;
        Nyt3i = Nyt3+iitinc;
        Nyt4i = Nyt4+iitinc;
        Nyt5i = Nyt5+iitinc; 
        Nyt6i = Nyt6+iitinc; 
        // u-1 u-2 u-3 u-4 u-5 u-6
        Nytn1i = Nytn1+iitinc;
        Nytn2i = Nytn2+iitinc; 
        Nytn3i = Nytn3+iitinc; 
        Nytn4i = Nytn4+iitinc; 
        Nytn5i = Nytn5+iitinc; 
        Nytn6i = Nytn6+iitinc;
    
        // The following overwrites dudx matrix for populating derivatives computed using finite differencing
        // First column/row du/dx
        dudxx[Nyt1i] = (uu[Nyt2i]-uu[Nytn1i])*bydxby4t3
        +(-uu[Nyt3i]+uu[Nytn2i])*bydxby20t3
        +(uu[Nyt4i]-uu[Nytn3i])*bydxby60;
        
        // Second column/row du/dx
        dudxx[Nyt2i] = (uu[Nyt3i]-uu[Nyt1i])*bydxby4t3
        +(-uu[Nyt4i]+uu[Nytn1i])*bydxby20t3
        +(+uu[Nyt5i]-uu[Nytn2i])*bydxby60;
        
        // Third column/row du/dx
        dudxx[Nyt3i] = (uu[Nyt4i]-uu[Nyt2i])*bydxby4t3
        +(-uu[Nyt5i]+uu[Nyt1i])*bydxby20t3
        +(uu[Nyt6i]-uu[Nytn1i])*bydxby60;
        
        // Last column/row du/dx
        dudxx[Nytn1i] = (uu[Nyt1i]-uu[Nytn2i])*bydxby4t3
        +(-uu[Nyt2i]+uu[Nytn3i])*bydxby20t3
        +(uu[Nyt3i]-uu[Nytn4i])*bydxby60;
        
        // Second last column/row du/dx
        dudxx[Nytn2i] = (uu[Nytn1i]-uu[Nytn3i])*bydxby4t3
        +(-uu[Nyt1i]+uu[Nytn4i])*bydxby20t3
        +(uu[Nyt2i]-uu[Nytn5i])*bydxby60;
        
        // Third last column/row du/dx
        dudxx[Nytn3i] = (uu[Nytn2i]-uu[Nytn4i])*bydxby4t3
        +(-uu[Nytn1i]+uu[Nytn5i])*bydxby20t3
        +(uu[Nyt1i]-uu[Nytn6i])*bydxby60;
    }
}



void ShallowWater::Calcddt(){
    // Precompute g=9.81 gravitational acceleration
    double g = 9.81;

    // First initialise useful indices outside of for loop
    int Nytipj;

    // Calculate spacial derivatives for calculating time derivatives for k1
    int incy=Ny; //for the computation of x derivatives
    ShallowWater::GradAtBound(dx, Nx, Ny, incy, u, dudx);
    ShallowWater::GradAtBound(dx, Nx, Ny, incy, v, dvdx);
    ShallowWater::GradAtBound(dx, Nx, Ny, incy, h, dhdx);
    ShallowWater::GradElseWhere(dx, Nx, Ny, incy, u, dudx);
    ShallowWater::GradElseWhere(dx, Nx, Ny, incy, v, dvdx);
    ShallowWater::GradElseWhere(dx, Nx, Ny, incy, h, dhdx);
    incy=1; //for the coomputation of y derivatives
    ShallowWater::GradAtBound(dy, Nx, Ny, incy, u, dudy);
    ShallowWater::GradAtBound(dy, Nx, Ny, incy, v, dvdy);
    ShallowWater::GradAtBound(dy, Nx, Ny, incy, h, dhdy);
    ShallowWater::GradElseWhere(dy, Nx, Ny, incy, u, dudy);
    ShallowWater::GradElseWhere(dy, Nx, Ny, incy, v, dvdy);
    ShallowWater::GradElseWhere(dy, Nx, Ny, incy, h, dhdy); 

    // loop to populate d/dt at k1
    // parallelised the for loop, letting omp to decide chunk, using collapse to parallelise both the outer and inner-loop
    int i = 0, j = 0; // for parallelisation
    #pragma omp parallel for    \
    default(shared)             \
    schedule (static)           \
    private(i,j,Nytipj)         \
    collapse(2)
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            // Precompute index Ny*i+j
            Nytipj = Ny*i+j;

            // Populate k1=f(y_n) time derivatives d(u,v and h)/dt  
            dudt[Nytipj] = -u[Nytipj]*dudx[Nytipj]-v[Nytipj]*dudy[Nytipj]-g*dhdx[Nytipj];
            dvdt[Nytipj] = -u[Nytipj]*dvdx[Nytipj]-v[Nytipj]*dvdy[Nytipj]-g*dhdy[Nytipj];
            dhdt[Nytipj] = -u[Nytipj]*dhdx[Nytipj]-h[Nytipj]*dudx[Nytipj]-
            v[Nytipj]*dhdy[Nytipj]-h[Nytipj]*dvdy[Nytipj];

            // compute y_n1=y_n+dt*k_1/2
            u2[Nytipj] = (u[Nytipj]+dt*dudt[Nytipj]*0.5); // 0.5 to avoid division
            v2[Nytipj] = (v[Nytipj]+dt*dvdt[Nytipj]*0.5);
            h2[Nytipj] = (h[Nytipj]+dt*dhdt[Nytipj]*0.5);

        }
    }

    // recalculate spacial derivatives for point k2
    incy=Ny; //for the computation of x derivatives
    ShallowWater::GradAtBound(dx, Nx, Ny, incy, u2, dudx);
    ShallowWater::GradAtBound(dx, Nx, Ny, incy, v2, dvdx);
    ShallowWater::GradAtBound(dx, Nx, Ny, incy, h2, dhdx);
    ShallowWater::GradElseWhere(dx, Nx, Ny, incy, u2, dudx);
    ShallowWater::GradElseWhere(dx, Nx, Ny, incy, v2, dvdx);
    ShallowWater::GradElseWhere(dx, Nx, Ny, incy, h2, dhdx);
    incy=1; //for the computation of y derivatives
    ShallowWater::GradAtBound(dy, Nx, Ny, incy, u2, dudy);
    ShallowWater::GradAtBound(dy, Nx, Ny, incy, v2, dvdy);
    ShallowWater::GradAtBound(dy, Nx, Ny, incy, h2, dhdy);
    ShallowWater::GradElseWhere(dy, Nx, Ny, incy, u2, dudy);
    ShallowWater::GradElseWhere(dy, Nx, Ny, incy, v2, dvdy);
    ShallowWater::GradElseWhere(dy, Nx, Ny, incy, h2, dhdy);  

    // another loop to populate d/dt at k2
    // parallelised the for loop, letting omp to decide chunk, using collapse to parallelise both the outer and inner-loop
    i = 0; // for parallelisation
    j = 0; // for parallelisation
    #pragma omp parallel for    \
    default(shared)             \
    schedule (static)           \
    private(i,j,Nytipj)         \
    collapse(2)
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            // Precompute index Ny*i+j
            Nytipj = Ny*i+j;


             // Populate k2=f(y_n+dt*k_1/2) time derivatives
            k2u[Nytipj] = -(u2[Nytipj])*dudx[Nytipj]-(v2[Nytipj])*dudy[Nytipj]-g*dhdx[Nytipj];
            k2v[Nytipj] = -(u2[Nytipj])*dvdx[Nytipj]-(v2[Nytipj])*dvdy[Nytipj]-g*dhdy[Nytipj];
            k2h[Nytipj] = -(u2[Nytipj])*dhdx[Nytipj]-(h2[Nytipj])*dudx[Nytipj]-
            (v2[Nytipj])*dhdy[Nytipj]-(h2[Nytipj])*dvdy[Nytipj];



            // compute y_n2=y_n+dt*k_2/2
            u2[Nytipj] = (u[Nytipj]+dt*k2u[Nytipj]*0.5); // avoid division
            v2[Nytipj] = (v[Nytipj]+dt*k2v[Nytipj]*0.5);
            h2[Nytipj] = (h[Nytipj]+dt*k2h[Nytipj]*0.5);

        }
    }

    // recalculate spacial derivatives for point k3
    incy=Ny; //for the computation of x derivatives
    ShallowWater::GradAtBound(dx, Nx, Ny, incy, u2, dudx);
    ShallowWater::GradAtBound(dx, Nx, Ny, incy, v2, dvdx);
    ShallowWater::GradAtBound(dx, Nx, Ny, incy, h2, dhdx);
    ShallowWater::GradElseWhere(dx, Nx, Ny, incy, u2, dudx);
    ShallowWater::GradElseWhere(dx, Nx, Ny, incy, v2, dvdx);
    ShallowWater::GradElseWhere(dx, Nx, Ny, incy, h2, dhdx);
    incy=1; //for the computation of y derivatives
    ShallowWater::GradAtBound(dy, Nx, Ny, incy, u2, dudy);
    ShallowWater::GradAtBound(dy, Nx, Ny, incy, v2, dvdy);
    ShallowWater::GradAtBound(dy, Nx, Ny, incy, h2, dhdy);
    ShallowWater::GradElseWhere(dy, Nx, Ny, incy, u2, dudy);
    ShallowWater::GradElseWhere(dy, Nx, Ny, incy, v2, dvdy);
    ShallowWater::GradElseWhere(dy, Nx, Ny, incy, h2, dhdy);  

    // another loop to populate d/dt at k3
    // parallelised the for loop, letting omp to decide chunk, using collapse to parallelise both the outer and inner-loop
    i = 0; // for parallelisation
    j = 0; // for parallelisation
    #pragma omp parallel for    \
    default(shared)             \
    schedule (static)           \
    private(i,j,Nytipj)         \
    collapse(2)
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            // Precompute index Ny*i+j
            Nytipj = Ny*i+j;


            // Populate k3=f(y_n+dt*k_2/2) time derivatives
            k3u[Nytipj] = -(u2[Nytipj])*dudx[Nytipj]-(v2[Nytipj])*dudy[Nytipj]-g*dhdx[Nytipj];
            k3v[Nytipj] = -(u2[Nytipj])*dvdx[Nytipj]-(v2[Nytipj])*dvdy[Nytipj]-g*dhdy[Nytipj];
            k3h[Nytipj] = -(u2[Nytipj])*dhdx[Nytipj]-(h2[Nytipj])*dudx[Nytipj]-
            (v2[Nytipj])*dhdy[Nytipj]-(h2[Nytipj])*dvdy[Nytipj];


            // precompute y_n3=y_n+dt*k_3
            u2[Nytipj] = (u[Nytipj]+dt*k3u[Nytipj]); // avoid division
            v2[Nytipj] = (v[Nytipj]+dt*k3v[Nytipj]);
            h2[Nytipj] = (h[Nytipj]+dt*k3h[Nytipj]);

        }
    }

    // recalculate spacial derivatives for point k3
    incy=Ny; //for the computation of x derivatives
    ShallowWater::GradAtBound(dx, Nx, Ny, incy, u2, dudx);
    ShallowWater::GradAtBound(dx, Nx, Ny, incy, v2, dvdx);
    ShallowWater::GradAtBound(dx, Nx, Ny, incy, h2, dhdx);
    ShallowWater::GradElseWhere(dx, Nx, Ny, incy, u2, dudx);
    ShallowWater::GradElseWhere(dx, Nx, Ny, incy, v2, dvdx);
    ShallowWater::GradElseWhere(dx, Nx, Ny, incy, h2, dhdx);
    incy=1; //for the computation of y derivatives
    ShallowWater::GradAtBound(dy, Nx, Ny, incy, u2, dudy);
    ShallowWater::GradAtBound(dy, Nx, Ny, incy, v2, dvdy);
    ShallowWater::GradAtBound(dy, Nx, Ny, incy, h2, dhdy);
    ShallowWater::GradElseWhere(dy, Nx, Ny, incy, u2, dudy);
    ShallowWater::GradElseWhere(dy, Nx, Ny, incy, v2, dvdy);
    ShallowWater::GradElseWhere(dy, Nx, Ny, incy, h2, dhdy);  

    // another loop to populate d/dt at k4
    // parallelised the for loop, letting omp to decide chunk, using collapse to parallelise both the outer and inner-loop
    i = 0; // for parallelisation
    j = 0; // for parallelisation
    #pragma omp parallel for\
    default(shared)         \
    schedule (static)       \
    private(i,j,Nytipj)     \
    collapse(2)
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            // Precompute index Ny*i+j
            Nytipj = Ny*i+j;


            // Populate k4=f(y_n+dt*k_3) time derivatives 
            k4u[Nytipj] = -(u2[Nytipj])*dudx[Nytipj]-(v2[Nytipj])*dudy[Nytipj]-g*dhdx[Nytipj];
            k4v[Nytipj] = -(u2[Nytipj])*dvdx[Nytipj]-(v2[Nytipj])*dvdy[Nytipj]-g*dhdy[Nytipj];
            k4h[Nytipj] = -(u2[Nytipj])*dhdx[Nytipj]-(h2[Nytipj])*dudx[Nytipj]-
            (v2[Nytipj])*dhdy[Nytipj]-(h2[Nytipj])*dvdy[Nytipj];

        }
    }

}



void ShallowWater::FillSixOrdMatrix() {
    // Populate SixOrd1 as column major
    // Initialse matrices entries with zeros
    for (int i = 0; i < Nx*Nx; ++i){ // fill banded diagonal entries as column major
        SixOrd1[i] = 0; //central diagonal
    }

    double bydx = 1.0/dx; // precompute divisions 1/dx
    double bydxby60 = bydx/60.0; // pre-compute finite differencing coefficients
    double bydxby20t3 = bydx*3.0/20.0; //
    double bydxby4t3 = bydx*3.0/4.0;

    for (int i = 0; i < Nx; ++i){ // fill banded diagonal entries as column major
        SixOrd1[i*Nx+i    ] = 0; //central diagonal
    }

    for (int i = 0; i < Nx-1; ++i){ // fill banded diagonal entries as column major
        SixOrd1[i*Nx+i+1    ] = -bydxby4t3; //1st lower diagonal
        SixOrd1[i*Nx+i+Nx   ] = bydxby4t3; //1st upper diagonal
    }

    for (int i = 0; i < Nx-2; ++i){ // fill banded diagonal entries as column major
        SixOrd1[i*Nx+i+2    ] = bydxby20t3; //2nd lower diagonal
        SixOrd1[i*Nx+i+Nx+Nx    ] = -bydxby20t3; //2nd upper diagonal
    }

    for (int i = 0; i < Nx-3; ++i){ // fill banded diagonal entries as column major
        SixOrd1[i*Nx+i+3    ] = -bydxby60; //3rd lower diagonal
        SixOrd1[i*Nx+i+Nx+Nx+Nx    ] = bydxby60; //3rd upper diagonal
    }

    //populate manually the repeated BC
    // at the top right corner
    SixOrd1[Nx*(Nx-1)] = -bydxby4t3;
    SixOrd1[Nx*(Nx-2)] = bydxby20t3;
    SixOrd1[Nx*(Nx-3)] = -bydxby60;
    SixOrd1[Nx*(Nx-1)+1] = +bydxby20t3;
    SixOrd1[Nx*(Nx-2)+1] = -bydxby60;
    SixOrd1[Nx*(Nx-1)+2] = -bydxby60;

    // at the bottom left corner
    SixOrd1[Nx-1] = bydxby4t3;
    SixOrd1[Nx-1+Nx] = -bydxby20t3;
    SixOrd1[Nx-1+Nx+Nx] = bydxby60;
    SixOrd1[Nx-2] = -bydxby20t3;
    SixOrd1[Nx-2+Nx] = +bydxby60;
    SixOrd1[Nx-3] = +bydxby60;



    // Now Populate sixOrd2 FD matrix
    // Initialse matrices entries with zeros
    for (int i = 0; i < Ny*Ny; ++i){ // fill banded diagonal entries as column major
        SixOrd2[i] = 0;
    }
    // Precompute expressions
    double bydy = 1.0/dy; // precompute divisions 
    double bydyby60 = bydy/60.0; // also precompute multiplications
    double bydyby20t3 = bydy*3.0/20.0;
    double bydyby4t3 = bydy*3.0/4.0;

    for (int i = 0; i < Ny; ++i){ // fill banded diagonal entries as column major
        SixOrd2[i*Ny+i    ] = 0; //central diagonal
    }

    for (int i = 0; i < Ny-1; ++i){ // fill banded diagonal entries as column major
        SixOrd2[i*Ny+i+1    ] = -bydyby4t3; //1st lower diagonal
        SixOrd2[i*Ny+i+Ny   ] = bydyby4t3; //1st upper diagonal
    }

    for (int i = 0; i < Ny-2; ++i){ // fill banded diagonal entries as column major
        SixOrd2[i*Ny+i+2    ] = bydyby20t3; //2nd lower diagonal
        SixOrd2[i*Ny+i+Ny+Ny    ] = -bydyby20t3; //2nd upper diagonal
    }

    for (int i = 0; i < Ny-3; ++i){ // fill banded diagonal entries as column major
        SixOrd2[i*Ny+i+3    ] = -bydyby60; //3rd lower diagonal
        SixOrd2[i*Ny+i+Ny+Ny+Ny    ] = bydyby60; //3rd upper diagonal
    }

    //populate manually the repeated BC
    // at the top right corner
    SixOrd2[Ny*(Ny-1)] = -bydyby4t3;
    SixOrd2[Ny*(Ny-2)] = bydyby20t3;
    SixOrd2[Ny*(Ny-3)] = -bydyby60;
    SixOrd2[Ny*(Ny-1)+1] = bydyby20t3;
    SixOrd2[Ny*(Ny-2)+1] = -bydyby60;
    SixOrd2[Ny*(Ny-1)+2] = -bydyby60;

    // at the bottom left corner
    SixOrd2[Ny-1] = bydyby4t3;
    SixOrd2[Ny-1+Ny] = -bydyby20t3;
    SixOrd2[Ny-1+Ny+Ny] = bydyby60;
    SixOrd2[Ny-2] = -bydyby20t3;
    SixOrd2[Ny-2+Ny] = +bydyby60;
    SixOrd2[Ny-3] = +bydyby60;
} 



void ShallowWater::BLASGrad(double* uu, double* vv, double* hh){

    int i = 0; // for parallelisation

    // no interdependancy to loop thru different rows of u for dudx in parallel // parallelise with 8 threads // make iteration number i private
    // static to avoid overhead, let omp to decide chunksize as it better distributes than us
    // each thread writes to different array indices, as they work on different rows, so race condition is not a matter, no need to atomic etc.
    // share everything except iteration number
    #pragma omp parallel sections\
    default(shared)
    { 
        // 2 parallel for loops are grouped in parallel sections, this increases the total amount of work and operations to be divided
        // thus providing possible benefits when originally the number of operations is not divisible by the number of threads, dividing 
        // the work more evenly between threads.
        // do not wait until other thread to complete, there is no need since there is no interdependancies

        #pragma omp section
        #pragma omp parallel for\
        default(shared)         \
        schedule (static)       \
        private(i)
        for(i=0 ; i<Ny ; i++){
            // Calculate the spacial gradients du/dx row by row using dgemv
            // increment set to Ny to jump to the next element in the same row for matrix uu and dudx 
            // used pointer arithmetics to move down an element after each iteration
            F77NAME(dgemv)('n', Nx, Nx, 1.0, SixOrd1, 
            Nx, uu+i, Ny, 0.0, dudx+i, Ny);

            F77NAME(dgemv)('n', Nx, Nx, 1.0, SixOrd1, 
            Nx, vv+i, Ny, 0.0, dvdx+i, Ny);       

            F77NAME(dgemv)('n', Nx, Nx, 1.0, SixOrd1, 
            Nx, hh+i, Ny, 0.0, dhdx+i, Ny); 

        }

        #pragma omp section
        #pragma omp parallel for\
        default(shared)         \
        schedule (static)       \
        private(i)
        for(i=0 ; i<Nx ; i++){

            //Populate spacial gradients in d/dy in for column by column using dgemv
            F77NAME(dgemv)('n', Ny, Ny, 1.0, SixOrd2, 
            Ny, uu+(i*Ny), 1, 0.0, dudy+(i*Ny), 1);

            F77NAME(dgemv)('n', Ny, Ny, 1.0, SixOrd2, 
            Ny, vv+(i*Ny), 1, 0.0, dvdy+(i*Ny), 1);       

            F77NAME(dgemv)('n', Ny, Ny, 1.0, SixOrd2, 
            Ny, hh+(i*Ny), 1, 0.0, dhdy+(i*Ny), 1);      

        }
    }
}



void ShallowWater::BLASCalcddt(){

    //First calculate spacial gradients d(u,v,h)/d(x,y) using BLASGrad with u,v and h
    ShallowWater::BLASGrad(u, v, h);

    //Cannot do elementwise multiplication using BLAS (can do but inefficient)
    //so we use the same method from for loop method to calculate d/dt
    //declare loop index variables 
    int Nytipj;
    //Pre assign g
    double g = 9.81;
    // loop to populate d/dt at k1
    // parallelised the for loop, letting omp to decide chunk, using collapse to parallelise both the outer and inner-loop
    int i = 0, j = 0; // for parallelisation
    #pragma omp parallel for\
    default(shared)         \
    schedule (static)       \
    private(i,j,Nytipj)     \
    collapse(2)
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            // Precompute index Ny*i+j
            Nytipj = Ny*i+j;


            // Populate k1 time derivatives = d(u,v and h)/dt
            dudt[Nytipj] = -u[Nytipj]*dudx[Nytipj]-v[Nytipj]*dudy[Nytipj]-g*dhdx[Nytipj];
            dvdt[Nytipj] = -u[Nytipj]*dvdx[Nytipj]-v[Nytipj]*dvdy[Nytipj]-g*dhdy[Nytipj];
            dhdt[Nytipj] = -u[Nytipj]*dhdx[Nytipj]-h[Nytipj]*dudx[Nytipj]-
            v[Nytipj]*dhdy[Nytipj]-h[Nytipj]*dvdy[Nytipj];
        }
    }

    // But we can use dcopy and daxpy to calculate the new u2=u+dt/2*du/dt
    // Precompute dt/2
    double dtby2=dt*0.5;
    F77NAME(dcopy) (Nx*Ny, u, 1, u2, 1); // u2=u;
    F77NAME(daxpy) (Nx*Ny, dtby2, dudt, 1, u2, 1); // u2=u+dt/2*dudt;
    // Repeat for v and h
    F77NAME(dcopy) (Nx*Ny, v, 1, v2, 1);
    F77NAME(daxpy) (Nx*Ny, dtby2, dvdt, 1, v2, 1);
    F77NAME(dcopy) (Nx*Ny, h, 1, h2, 1);
    F77NAME(daxpy) (Nx*Ny, dtby2, dhdt, 1, h2, 1);

    // After u2 v2 h2 has been computed we can start computing k2
    // Similarly, first populate the spacial gradient matrices, but this time using u2,v2,h2 at k2
    ShallowWater::BLASGrad(u2, v2, h2);

    // again cant use BLAS for elementwise multiplication so use for loop to populate d/dt at k2
    // parallelised the for loop, letting omp to decide chunk, using collapse to parallelise both the outer and inner-loop
    i = 0; // for parallelisation
    j = 0; // for parallelisation
    #pragma omp parallel for\
    default(shared)         \
    schedule (static)       \
    private(i,j,Nytipj)     \
    collapse(2)
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            // Precompute index Ny*i+j
            Nytipj = Ny*i+j;

            // Populate k2 time derivatives 
            k2u[Nytipj] = -(u2[Nytipj])*dudx[Nytipj]-(v2[Nytipj])*dudy[Nytipj]-g*dhdx[Nytipj];
            k2v[Nytipj] = -(u2[Nytipj])*dvdx[Nytipj]-(v2[Nytipj])*dvdy[Nytipj]-g*dhdy[Nytipj];
            k2h[Nytipj] = -(u2[Nytipj])*dhdx[Nytipj]-(h2[Nytipj])*dudx[Nytipj]-
            (v2[Nytipj])*dhdy[Nytipj]-(h2[Nytipj])*dvdy[Nytipj];
        }
    } 
    

    // we can use dcopy and daxpy to calculate and overwrite the new u2=u+dt/2*(k2)
    F77NAME(dcopy) (Nx*Ny, u, 1, u2, 1);
    F77NAME(daxpy) (Nx*Ny, dtby2, k2u, 1, u2, 1);
    // Repeat for v and h
    F77NAME(dcopy) (Nx*Ny, v, 1, v2, 1);
    F77NAME(daxpy) (Nx*Ny, dtby2, k2v, 1, v2, 1);
    F77NAME(dcopy) (Nx*Ny, h, 1, h2, 1);
    F77NAME(daxpy) (Nx*Ny, dtby2, k2h, 1, h2, 1);

    // After u3 v3 h3 has been computed we can start computing k3
    // Similarly, first populate the spacial gradient matrices, but this time using u2,v2,h2 at k3
    ShallowWater::BLASGrad(u2, v2, h2);

    // again cant use BLAS for elementwise multiplication so use for loop to populate d/dt at k3
    // parallelised the for loop, letting omp to decide chunk, using collapse to parallelise both the outer and inner-loop
    i = 0; // for parallelisation
    j = 0; // for parallelisation
    #pragma omp parallel for\
    default(shared)         \
    schedule (static)       \
    private(i,j,Nytipj)     \
    collapse(2)
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            // Precompute index Ny*i+j
            Nytipj = Ny*i+j;

            // Populate k3 time derivatives
            k3u[Nytipj] = -(u2[Nytipj])*dudx[Nytipj]-(v2[Nytipj])*dudy[Nytipj]-g*dhdx[Nytipj];
            k3v[Nytipj] = -(u2[Nytipj])*dvdx[Nytipj]-(v2[Nytipj])*dvdy[Nytipj]-g*dhdy[Nytipj];
            k3h[Nytipj] = -(u2[Nytipj])*dhdx[Nytipj]-(h2[Nytipj])*dudx[Nytipj]-
            (v2[Nytipj])*dhdy[Nytipj]-(h2[Nytipj])*dvdy[Nytipj];
        }
    }
    

    // we can use dcopy and daxpy to calculate and overwrite the new u2=u+dt*(k3)
    F77NAME(dcopy) (Nx*Ny, u, 1, u2, 1);
    F77NAME(daxpy) (Nx*Ny, dt, k3u, 1, u2, 1);
    // Repeat for v and h
    F77NAME(dcopy) (Nx*Ny, v, 1, v2, 1);
    F77NAME(daxpy) (Nx*Ny, dt, k3v, 1, v2, 1);
    F77NAME(dcopy) (Nx*Ny, h, 1, h2, 1);
    F77NAME(daxpy) (Nx*Ny, dt, k3h, 1, h2, 1);

    // After u3 v3 h3 has been computed we can start computing the final k4
    // Similarly, first populate the spacial gradient matrices, but this time using u2,v2,h2 at k4
    ShallowWater::BLASGrad(u2, v2, h2);

    // again cant use BLAS for elementwise multiplication so use for loop to populate d/dt at k4
    // parallelised the for loop, letting omp to decide chunk, using collapse to parallelise both the outer and inner-loop
    i = 0; // for parallelisation
    j = 0; // for parallelisation
    #pragma omp parallel for\
    default(shared)         \
    schedule (static)       \
    private(i,j,Nytipj)     \
    collapse(2)
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            // Precompute index Ny*i+j
            Nytipj = Ny*i+j;


            // Populate k4 time derivatives 
            k4u[Nytipj] = -(u2[Nytipj])*dudx[Nytipj]-(v2[Nytipj])*dudy[Nytipj]-g*dhdx[Nytipj];
            k4v[Nytipj] = -(u2[Nytipj])*dvdx[Nytipj]-(v2[Nytipj])*dvdy[Nytipj]-g*dhdy[Nytipj];
            k4h[Nytipj] = -(u2[Nytipj])*dhdx[Nytipj]-(h2[Nytipj])*dudx[Nytipj]-
            (v2[Nytipj])*dhdy[Nytipj]-(h2[Nytipj])*dvdy[Nytipj];
        }
    }
}


void ShallowWater::BLASRK(){
    //Precompute fractions
    double dtby6 = dt/6.0;
    double dtby6t2 = dtby6*2.0;
    //Forward time propagation by adding terms to the original u, v and h at (t_n) for u,v and h at (t_n+1)
    // u 
    F77NAME(daxpy) (Nx*Ny, dtby6, dudt, 1, u, 1);
    F77NAME(daxpy) (Nx*Ny, dtby6t2, k2u, 1, u, 1); 
    F77NAME(daxpy) (Nx*Ny, dtby6t2, k3u, 1, u, 1);
    F77NAME(daxpy) (Nx*Ny, dtby6, k4u, 1, u, 1); 
    // v
    F77NAME(daxpy) (Nx*Ny, dtby6, dvdt, 1, v, 1);
    F77NAME(daxpy) (Nx*Ny, dtby6t2, k2v, 1, v, 1); 
    F77NAME(daxpy) (Nx*Ny, dtby6t2, k3v, 1, v, 1);
    F77NAME(daxpy) (Nx*Ny, dtby6, k4v, 1, v, 1); 
    // h
    F77NAME(daxpy) (Nx*Ny, dtby6, dhdt, 1, h, 1);
    F77NAME(daxpy) (Nx*Ny, dtby6t2, k2h, 1, h, 1); 
    F77NAME(daxpy) (Nx*Ny, dtby6t2, k3h, 1, h, 1);
    F77NAME(daxpy) (Nx*Ny, dtby6, k4h, 1, h, 1);     
}
