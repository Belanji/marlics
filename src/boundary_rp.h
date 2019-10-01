#ifndef BOUNDARY_RP_
#define BOUNDARY_RP_

#include <petscts.h>

class Boundary_Rp : public BOUNDARY
{
  
 protected:

  double Wo1;
  double Q0_00;
  double Q0_01;
  double Q0_02;
  double Q0_11;
  double Q0_12;
  double Lambda_s;
  
   public:

  virtual double surface_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double surface_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double surface_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double surface_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double surface_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual ~Boundary_Rp() {};


  virtual void fill_jacobian_boundary(const PetscScalar *Qij,Mat Jac,Mat Jac_pc, const PetscScalar * v, const int i, const int j, const int k) const;

  void fill_dFijQij(PetscScalar dFijdQij[5][5],const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const;

  void fill_dFijQijk(PetscScalar dFijdQijk[5][5][3],const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const;
  
  Boundary_Rp(const class Simulation_Parameters *, int);
};

 

#endif
