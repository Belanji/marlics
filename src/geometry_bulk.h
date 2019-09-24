#ifndef GEOMETRY_BULK_
#define GEOMETRY_BULK_

#include <petscts.h>
#include <petscsnes.h>


class Geometry_Bulk : public GEOMETRY
{

  
 public:

  static PetscErrorCode RhsFunction(TS ts,PetscReal tiem,Vec Qij_in,Vec Rhs,void *sim_geometry);

  static PetscErrorCode RhsJacobian(TS ts,PetscReal time,Vec Qij_in,Mat Jac,Mat Jac_pc, void* sim_param);


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
