#include "driver.h"
#include "energy.h"
#include "boundary.h"
#include <string>
#include <iostream>
//base class to posterior boudary definition
BOUNDARY::BOUNDARY(const  Simulation_Parameters *  sim_par, int boundary_number): ENERGY(sim_par),
  Lambda_s(1/sim_par->mu_1_s),
  boundary_id(boundary_number)
  { };

/*Warn about the use of a standart value in the absence of a parameter value.
 Print a error menssage if the standart value doesn't exist.*/
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

          std::cout << "Parameter " << parameter_name << " not defined for the boundary condition " << "#" << boundary_id << ".\n";
          std::cout << "The boundary condition " << condition_name << " needs the aforementioned parameter defined.\n";
          std::cout << "Please, set it in your input file, or check if it is mispelled.\n";
          std::cout <<"Aborting the program.\n";
          exit(1);
            
        }
    }
  
}
