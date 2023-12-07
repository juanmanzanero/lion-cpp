#ifndef MUMPS_INTERFACE_H__
#define MUMPS_INTERFACE_H__

#include <stdio.h>
#include <string.h>
#include <algorithm>
#include "coin-or/mumps/dmumps_c.h"
#include "lion/foundation/lion_exception.h"
#define USE_COMM_WORLD -987654

template<typename size_type>
inline std::vector<std::vector<double>> mumps_solve_linear_system(size_type n_cpp, size_type nnz_cpp, std::vector<size_type>& rows_cpp, std::vector<size_type>& cols_cpp, 
    std::vector<double>& lhs_cpp, const std::vector<std::vector<double>>& rhs_cpp, bool symmetric)
{
    assert(rows_cpp.size() == nnz_cpp);
    assert(cols_cpp.size() == nnz_cpp);
    assert(lhs_cpp.size()  == nnz_cpp);
    size_t n_rhs = rhs_cpp.size();

    assert(n_rhs > 0);
    for ( const auto& rhs_i : rhs_cpp ) {
      assert(rhs_i.size() == n_cpp);
    }
  
    DMUMPS_STRUC_C id;
    MUMPS_INT n = static_cast<MUMPS_INT>(n_cpp);
    MUMPS_INT8 nnz = static_cast<MUMPS_INT8>(nnz_cpp);
    std::vector<MUMPS_INT> irn(nnz), jcn(nnz);
    std::transform(rows_cpp.cbegin(), rows_cpp.cend(), irn.begin(), [](const auto& i) -> auto { return static_cast<MUMPS_INT>(i+1); });
    std::transform(cols_cpp.cbegin(), cols_cpp.cend(), jcn.begin(), [](const auto& i) -> auto { return static_cast<MUMPS_INT>(i+1); });
    auto x = rhs_cpp;
  
    int myid = 0;
  
    int error = 0;
  
    /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
    id.comm_fortran=USE_COMM_WORLD;
  
    // (1) Initialize mumps structure
  
    // Mumps settings:
      
    // Use the host to factorize the matrix
    id.par = 1;     
  
    // Solve general symmetric matrix
    id.sym = (symmetric ? 2 : 0);
  
    // Job -1 is the initialization
    id.job = -1;
  
    // (1.1) Call mumps to initialize
    dmumps_c(&id);
  
    /* Define the problem on the host */
    id.n   = n; 
    id.nnz = nnz; 
    id.irn = irn.data(); 
    id.jcn = jcn.data();
    id.a   = lhs_cpp.data(); 
  
  #define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
  
    // ICNTL(1) is the output stream for error messages. If it is negative or zero, these messages will be suppressed. Default value is 6.
    id.ICNTL(1)=-1; 
     
    // ICNTL(2) is the output stream for diagnostic printing, statistics, and warning messages. 
    // If it is negative or zero, these messages will be suppressed. Default value is 0.
    id.ICNTL(2)=-1; 
  
    // ICNTL(3) is the output stream for global information, collected on the host. 
    // If it is negative or zero, these messages will be suppressed. Default value is 6
    id.ICNTL(3)=-1; 
  
    // ICNTL(4) is the level of printing for error, warning, and diagnostic messages. Maximum value is 4 and default value is 2 (errors and warnings printed)
    id.ICNTL(4) = 0;
  
    // These options have been directly copied from Ipopt
    id.ICNTL(14) = 1000.0;
  
    id.ICNTL(6) = 7;
    id.ICNTL(7) = 7;
    id.ICNTL(8) = 77;
    id.ICNTL(10) = 0;   //no iterative refinement iterations
  
    id.ICNTL(13) = 1;   //avoid lapack bug, ensures proper inertia; mentioned to be very expensive in mumps manual
    id.cntl[0] = 1.0e-6;  // Set pivot tolerance
  
  
  
    // Call the MUMPS package (analyse, factorization)
    id.job=4;
    dmumps_c(&id);
    if (id.infog[0]<0) {
      printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
          myid, id.infog[0], id.infog[1]);
      error = 1;
      if ( error ) 
      {
       /* Terminate instance. */
       id.job = -2;
       dmumps_c(&id);
        throw lion_exception("[ERROR] Mumps has thrown an error");
      }
    }
  
    // (2) Solve the systems
    for (size_t i = 0; i < n_rhs; ++i)
    {
      error = 0;
      id.rhs = x[i].data();
      id.job=3;
      dmumps_c(&id);
      if (id.infog[0]<0) {
        printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
          myid, id.infog[0], id.infog[1]);
        error = 1;
  
      }
      if ( error ) 
      {
       /* Terminate instance. */
       id.job = -2;
       dmumps_c(&id);
        throw lion_exception("[ERROR] Mumps has thrown an error");
      }
    }
  
    /* Terminate instance. */
    id.job = -2;
    dmumps_c(&id);
  
    return x;
}

#endif
