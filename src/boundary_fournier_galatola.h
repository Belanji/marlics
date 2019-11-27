#ifndef BOUNDARY_Fg_
#define BOUNDARY_Fg_

class Boundary_Fg : public BOUNDARY
{
  
 protected:

  double Wo1;

   public:

  virtual double functional_derivative_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double functional_derivative_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double functional_derivative_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double functional_derivative_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
  virtual double functional_derivative_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;

  virtual ~Boundary_Fg() {};
  Boundary_Fg(const class Simulation_Parameters *, int);
};

 

#endif
