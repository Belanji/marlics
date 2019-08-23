#ifndef HOMOGENEOUS_SPHERE_
#define HOMOGENEOUS_SPHERE_

#include "driver.h"
#include "geometry.h"

class Geometry_Sphere : public GEOMETRY
{
  

 protected:

  const int R_out;
  const int R_in;

  const double HNx;
  const double HNy;
  const double HNz;

  
 public:

  
  virtual void fill_ki(double * k_i, const double * Qij,const int i,const int j,const int k) const ;
    

  Geometry_Sphere(const struct Simulation_Parameters * );
  ~Geometry_Sphere(void);

   private:

  virtual int * fill_point_type( void ) const ;  
  
  
};

#endif
