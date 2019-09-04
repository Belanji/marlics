#ifndef BOUNDARY_STRONG_
#define BOUNDARY_STRONG_

#include "boundary.h"
#include "geometry.h"

class Strong_Boundary : public BOUNDARY
{

  public:

  virtual double surface_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double surface_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double surface_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double surface_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double surface_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;

  virtual ~Strong_Boundary();

  Strong_Boundary(const Simulation_Parameters * , int );

};

#endif
