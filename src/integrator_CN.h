#ifndef CN_
#define CN_

#include "integrator.h"
#include <stdio.h>           
#include <stdlib.h>
#include <petscts.h>
#include <petscsnes.h>


class CN : public Integrator
{

 private:


  Vec  Qsolution, Qtij, Rhs;
  Mat  Jac, Jac_Approx;
  TS   cranck_int;
  
  


  void getIcArray(double * Qij_input);
  void writeSolutionToArray(double * Qij_out);
  
  const double Atol;
  const double Rtol;
  const double prefac;
  const double facmin;
  const double facmax;    
 
 public:

  
  virtual void evolve( double *, double *, double);    
  virtual ~CN(void);
  CN( GEOMETRY  * , const struct Simulation_Parameters * );


};


#endif
