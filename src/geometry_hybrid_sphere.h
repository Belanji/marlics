#ifndef HYBRID_SPHERE_
#define HYBRID_SPHERE_

#include "driver.h"
#include "geometry.h"

class Geometry_Hybrid_Sphere : public GEOMETRY
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
  virtual void Energy_calc(double * k_i, const double * Qij) const ;
  virtual void compute_forces(double *, const double * ) const ;
     

  Geometry_Hybrid_Sphere(const struct Simulation_Parameters * );
  ~Geometry_Hybrid_Sphere(void);

   private:

  virtual int * fill_point_type( void ) const ;  
  
  
};

#endif
