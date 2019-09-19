#ifndef GEOMETRY_BULK_
#define GEOMETRY_BULK_

#include <petscts.h>
#include <petscsnes.h>


class Geometry_Bulk : public GEOMETRY
{
  
  
 public:


  virtual void fill_ki(double *, const double * ,const int ,const int ,const int) const ;
    
  Geometry_Bulk(const struct Simulation_Parameters *);
  ~Geometry_Bulk(void) {};

  PetscErrorCode RHsFunc(TS TS,PetscReal time,Vec Qij, Vec Rhs,void* someParams);

  
   private:

  virtual int * fill_point_type( void ) const ;  
  

  //void check_bulk_limits(int i, int j, int k) const ;
  //void check_surface_limits(int i, int j, int k) const ;
  
};

#endif
