#ifndef DP5_NOISE_
#define DP5_NOISE_

#include "integrator.h"
#include <stdio.h>           
#include <stdlib.h>
#include <gsl/gsl_rng.h>

class NOISE_DP5 : public Integrator
{

 private:


  const double Atol;
  const double Rtol;
  const double prefac;
  const double facmin;
  const double facmax;    
  gsl_rng *rng;
  double noise_factor;

 public:

  

  double *k_1, *k_2, *k_3, *k_4, *k_5, *k_6, *k_7, *noise;
  
  
  virtual void evolve( double *, double *, double);    
  virtual ~NOISE_DP5(void) {};
  NOISE_DP5( GEOMETRY  * , const struct Simulation_Parameters * );


};


#endif
