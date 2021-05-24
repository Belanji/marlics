#ifndef BOUNDARY_STRONG_
#define BOUNDARY_STRONG_

#include "energy.h"
#include "boundary.h"
#include "geometry.h"

class Strong_Boundary : public BOUNDARY
{

  public:

  virtual double functional_derivative_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double functional_derivative_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double functional_derivative_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double functional_derivative_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double functional_derivative_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double energy_calculation      (const double  QN[5],const double  dQ[15],const double v[3]) const ;

  virtual double force_00(const double  QN[27*5],const double dQ[], const double v[3]) const ;
  virtual double force_01(const double  QN[27*5],const double dQ[], const double v[3]) const ;
  virtual double force_02(const double  QN[27*5],const double dQ[], const double v[3]) const ;
  virtual double force_11(const double  QN[27*5],const double dQ[], const double v[3]) const ;
  virtual double force_12(const double  QN[27*5],const double dQ[], const double v[3]) const ;

  virtual ~Strong_Boundary();

  Strong_Boundary(const Simulation_Parameters * , int );

};

#endif
