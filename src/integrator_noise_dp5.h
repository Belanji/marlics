#ifndef NOISE_DP5_
#define NOISE_DP5_

#include "integrator.h"
#include <gsl/gsl_randist.h>  
#include <stdio.h>           
#include <stdlib.h>

class NoiseDP5 : public Integrator
{

 private:


  const double Atol;
  const double Rtol;
  const double prefac;
  const double facmin;
  const double facmax;    
  const double noise_factor;    
  const gsl_rng *w;
 public:

  

  double *k_1, *k_2, *k_3, *k_4, *k_5, *k_6, *k_7;
  
  
  virtual void evolve( double *, double *, double);    
  virtual ~NoiseDP5(void) {};
  NoiseDP5( GEOMETRY  * , const struct Simulation_Parameters * );


};


#endif
