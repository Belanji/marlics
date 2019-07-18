#ifndef INTEGRATOR_
#define INTEGRATOR_

#include <cstdlib>          
#include "geometry.h"


class Integrator
{

 protected:

  const int Nx, Ny, Nz;
  double dt;
  
  GEOMETRY * sample_geometry;
  virtual ~Integrator() {};
 
  

 public:

  
  virtual void evolve(double *,  double *, double) = 0;
  
  double *Qtij;
  
 Integrator( GEOMETRY  * lc_pointer ) :  sample_geometry(lc_pointer),

    Nx(lc_pointer->Nx),
    Ny(lc_pointer->Ny),
    Nz(lc_pointer->Nz)
    {

      if((Qtij= (double *)calloc(5*sample_geometry->Nx*sample_geometry->Ny*sample_geometry->Nz, sizeof(double)))==NULL){ERROr};

    };
      

};
 
#endif
