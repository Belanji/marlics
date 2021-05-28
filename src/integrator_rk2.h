#ifndef RK2_
#define RK2_

#include "integrator.h"
#include <stdio.h>           
#include <stdlib.h>

class RK2 : public Integrator
{

 private:

 public:

  

  double *k_1, *k_2;
  
  
  virtual bool evolve( double *, double *, double);    
  virtual ~RK2(void) {};
  RK2( GEOMETRY  * , const struct Simulation_Parameters * );


};


#endif
