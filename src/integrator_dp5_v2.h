#ifndef DP5_V2_
#define DP5_V2_

#include "integrator.h"
#include <stdio.h>           
#include <stdlib.h>

class DP5_2 : public Integrator
{

 private:


  const double Atol;
  const double Rtol;
  const double prefac;
  const double facmin;
  const double facmax;    

 public:

  

  double *k_1, *k_2, *k_3, *k_4, *k_5, *k_6, *k_7;
  
  
  virtual void evolve( double *, double *, double);    
  virtual ~DP5_2(void) {};
  DP5_2( GEOMETRY  * , const struct Simulation_Parameters * );


};


#endif
