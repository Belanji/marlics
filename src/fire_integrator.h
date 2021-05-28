#ifndef FIRE_
#define FIRE_

#include "integrator.h"
#include <stdio.h>           
#include <stdlib.h>

class FIRE : public Integrator
{

 private:

 public:

  

  double *force, *vQij, v2, f2, *Qijt, *energy;
  double dtmin, potence, alphainit,scale, alphaMin, dt_init;
  double scaleAlpha, scaleP, dtMax, scaleM, alpha;
  int Np, Nmin, proceed;
  
  virtual bool evolve( double *, double *, double);    
  virtual ~FIRE(void) {};
  FIRE( GEOMETRY  * , const struct Simulation_Parameters * );


};


#endif
