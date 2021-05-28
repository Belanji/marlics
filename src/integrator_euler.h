#ifndef EULER_
#define EULER_

#include "integrator.h"
#include <stdio.h>           
#include <stdlib.h>

class Euler : public Integrator
{

 private:

 public:

  

  double *k_1;
  
  
  virtual bool evolve( double *, double *, double);    
  virtual ~Euler(void) {};
  Euler( GEOMETRY  * , const struct Simulation_Parameters * );


};


#endif
