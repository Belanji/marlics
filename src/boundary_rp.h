#ifndef BOUNDARY_RP_
#define BOUNDARY_RP_


class Boundary_Rp : public BOUNDARY
{
  
 protected:

  double Wo1;
  double Q0_00;
  double Q0_01;
  double Q0_02;
  double Q0_11;
  double Q0_12;
  double dx, dy, dz;
  
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
  virtual ~Boundary_Rp() {};

  Boundary_Rp(const Simulation_Parameters *, int);
};

 

#endif
