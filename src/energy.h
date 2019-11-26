#ifndef ENERGY__
#define ENERGY__

class ENERGY
{
  public:


  ENERGY(const struct Simulation_Parameters * lc);
  
  virtual  double functional_derivative_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const =0;
  virtual  double functional_derivative_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const =0;
  virtual  double functional_derivative_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const =0;
  virtual  double functional_derivative_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const =0;
  virtual  double functional_derivative_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const =0;

 

 protected:
  
  const double sigma;
  const double a;
  const double bb;
  const double cc;
  const double L1;
  const double L2;
  const double L3;
  const double Lq;
  const double Ls;
  const double Lq_tilde;
  const double Lambda;
  const double Lambda_s;
  const double p0;
  const double q0;
  const double  S_eq;      
  const double deltaepslon;
  const double elecfieldx;
  const double elecfieldy;
  const double elecfieldz;

  //void  assert_parameter_is_set(bool parameter, class std::string parameter_name, bool has_standard_value =false);
  
  double dt;

};


#endif
