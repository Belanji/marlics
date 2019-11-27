#ifndef geometry_

#define geometry_
#include "energy.h"
#include "energy_ldg.h"

#include <vector>
#include <string>
#include <cstring>


class GEOMETRY
{

 public:
  
  const int Nx;
  const int Ny;
  const int Nz;
  const double dx_1;
  const double dy_1;
  const double dz_1;

  
  void read_check(int , int );
  
  virtual void fill_ki(double * k_i, const  double * Qij) const = 0;

  const int *point_type;
  int number_of_boundaries;  
  
  virtual void boundary_init( struct Simulation_Parameters * );
  virtual ~GEOMETRY() {};

 protected:


  std::string geometry_name;
  std::vector<class BOUNDARY *> bc_conditions;
  std::string boundary_needed_to_be_defined;


  const landau_de_gennes bulk_energy;
  
  virtual int * fill_point_type( void ) const = 0 ;  
  GEOMETRY * geometry_pointer;
  GEOMETRY(const struct Simulation_Parameters * );
  
};

#endif
