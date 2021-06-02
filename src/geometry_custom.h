#ifndef GEOMETRY_CUSTOM_
#define GEOMETRY_CUSTOM_

#include "driver.h"
#include "geometry.h"

class Geometry_Custom : public GEOMETRY
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
  double **v;
    
 public:
 
  virtual void fill_ki(double * k_i, const double * Qij) const ;
  virtual void Energy_calc(double * k_i, const double * Qij) const ;
  virtual void compute_forces(double *, const double * ) const ;
  virtual void test_derivatives(const char integrator[]) const;
  
  Geometry_Custom(const struct Simulation_Parameters * );
  ~Geometry_Custom(void);

   private:

  virtual int * fill_point_type( void ) const ;  
  virtual int * fill_pt_and_normals(double **v, const char bound_file_name[], int &num_of_boundaries) const;
};

#endif
