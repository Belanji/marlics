#include "boundary_strong.h"
#include "geometry.h"
#include <petscts.h>

Strong_Boundary::Strong_Boundary(const Simulation_Parameters * sim_param, int boundary_id) : BOUNDARY ( sim_param, boundary_id )
{

  condition_name="strong";
  printf("\n");
  
  
} ;

double Strong_Boundary::surface_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const { return 0.0; }  ;
double Strong_Boundary::surface_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const { return 0.0; }  ;
double Strong_Boundary::surface_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const { return 0.0; }  ;
double Strong_Boundary::surface_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const { return 0.0; }  ;
double Strong_Boundary::surface_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const { return 0.0; }  ;

Strong_Boundary::~Strong_Boundary() {};

void Strong_Boundary::fill_jacobian_boundary(const PetscScalar * Qij,Mat Jac,Mat Jac_p, const PetscScalar * v, const int i, const int j, const int k) const 
{

  PetscInt idxm[1],idxn[1];


  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+0;
  MatSetValues(Jac,1,idxm,1,idxn,0,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+1;
  MatSetValues(Jac,1,idxm,1,idxn,0,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+2;
  MatSetValues(Jac,1,idxm,1,idxn,0,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+3;
  MatSetValues(Jac,1,idxm,1,idxn,0,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+4;
  MatSetValues(Jac,1,idxm,1,idxn,0,ADD_VALUES);

  
};


