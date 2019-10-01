#ifndef BOUNDARY_
#define BOUNDARY_

#include <petscts.h>
#include <string>

class BOUNDARY
{
 
  
 protected:

  void  assert_parameter_is_set(bool parameter, std::string parameter_name, bool has_standard_value =false);

  const double sigma;
  const double a;
  const double bb;
  const double cc;
  const double dx_1, dy_1, dz_1;
  const double L1;
  const double L2;
  const double L3;
  const double Lq;
  const double Ls;
  const double Lq_tilde;
  const double Lambda, Lambda_s;
  const double  S_eq;
  const double q0;
  const double * Wo;
  const int Nx;
  const int Ny;
  const int Nz;
  double theta_0, phi_0;
  const int boundary_id;
  std::string condition_name;
  
 public:
    
  BOUNDARY(const class Simulation_Parameters *, int);
  virtual ~BOUNDARY() {};



  virtual double surface_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  = 0;
  virtual double surface_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  = 0;
  virtual double surface_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  = 0;
  virtual double surface_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  = 0;
  virtual double surface_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  = 0;

  virtual void fill_jacobian_boundary(const PetscScalar *,Mat ,Mat , const PetscScalar *, const int , const int , const int ) const {};

  void fill_dFijQij(PetscScalar dFijQij[5][5],const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const {};

  void fill_dFijQijk(PetscScalar dFijdQij[5][5][3],const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const {};

  
  
};


#endif
