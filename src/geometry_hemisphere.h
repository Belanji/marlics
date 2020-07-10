#ifndef HEMISPHERE_
#define HEMISPHERE_

#include "driver.h"
#include "geometry.h"

class Geometry_Hemisphere : public GEOMETRY
{
  

 private:

   double R_ex;
   double R_in;

  const double HNx;
  const double HNy;
  const double HNz;


  const double dx;
  const double dy;
  const double dz;

  
 public:

  
  virtual void fill_ki(double * k_i, const double * Qij) const ;
    

  Geometry_Hemisphere(const struct Simulation_Parameters * );
  ~Geometry_Hemisphere(void);

   private:

  virtual int * fill_point_type( void ) const ;  
  
  
};

#endif
