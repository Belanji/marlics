#include "driver.h"
#include "boundary.h"
#include <string>
#include <iostream>

BOUNDARY::BOUNDARY(const  Simulation_Parameters *  sim_par, int boundary_number):
  sigma(sim_par->T*sim_par->a), 
  a(sim_par->a),
  bb(sim_par->B),
  cc(sim_par->C),
  Nx(sim_par->Nx),
  Ny(sim_par->Nx),
  Nz(sim_par->Nx),
  dx_1(1/sim_par->dx),
  dy_1(1/sim_par->dy),
  dz_1(1/sim_par->dz),
  L1(sim_par->L1),
  L2(sim_par->L2),
  L3(sim_par->L3),
  Lq(sim_par->Lq),
  Ls(sim_par->Ls),
  Lq_tilde(2.0*sim_par->Lq*sim_par->q0),
  Lambda(1/sim_par->mu_1)    ,
  Lambda_s(1/sim_par->mu_1_s),
  S_eq(sim_par->S_eq),
  q0(sim_par->q0),
  boundary_id(boundary_number)
  { };


void BOUNDARY::assert_parameter_is_set(bool parameter, std::string parameter_name, bool has_standard_value )
{

  if(!parameter)
    {

      if(has_standard_value)
	{

	  std::cout <<"Parameter " << parameter_name << " is not set. Using standard value.\n";

	}
      else
	{

	  std::cout << "Paramter " << parameter_name << "not defined for the boundary condition number " << " number " << boundary_id << ".\n";
	  std::cout << "The boundary condition " << condition_name << " needs the aforementioned parameter defined.\n";
	  std::cout << "Please, set it in your input file, or check if it is mispelled.\n";
	  std::cout <<"Aborting the program.\n";
	  exit(0);
	    
	}
    }
  
}
