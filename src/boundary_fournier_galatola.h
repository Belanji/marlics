#ifndef BOUNDARY_Fg_
#define BOUNDARY_Fg_

class Boundary_Fg : public BOUNDARY
{
  
 protected:

  double Wo1;

   public:

  virtual double surface_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double surface_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double surface_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double surface_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double surface_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;

    virtual void fill_jacobian_boundary(const PetscScalar *Qij,Mat Jac,Mat Jac_pc, const PetscScalar * v, const int i, const int j, const int k) const {};

  
  virtual ~Boundary_Fg() {};
  Boundary_Fg(const class Simulation_Parameters *, int);
};

 

#endif
