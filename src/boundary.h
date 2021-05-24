#ifndef BOUNDARY_
#define BOUNDARY_

#include <string>

class BOUNDARY: public ENERGY
{
 
  
 protected:

  void  assert_parameter_is_set(bool parameter, std::string parameter_name, bool has_standard_value =false);
  

  const int boundary_id;
  const double Lambda_s;
  double * Wo;
  
  std::string condition_name;
  
 public:
    
  BOUNDARY(const class Simulation_Parameters *, int);
  virtual ~BOUNDARY() {};



  virtual double functional_derivative_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const=0  ;
  virtual double functional_derivative_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  =0;
  virtual double functional_derivative_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  =0;
  virtual double functional_derivative_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  =0;
  virtual double functional_derivative_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  =0;
  virtual double energy_calculation      (const double  QN[5],const double  dQ[15],const double v[3]) const  =0;
  virtual double force_00(const double  QN[5],const double dQ[], const double v[3]) const  =0;
  virtual double force_01(const double  QN[5],const double dQ[], const double v[3]) const  =0;
  virtual double force_02(const double  QN[5],const double dQ[], const double v[3]) const  =0;
  virtual double force_11(const double  QN[5],const double dQ[], const double v[3]) const  =0;
  virtual double force_12(const double  QN[5],const double dQ[], const double v[3]) const  =0;
    
};


#endif
