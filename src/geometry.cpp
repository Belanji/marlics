#include "geometry.h"
#include "driver.h"
#include "boundary.h"
#include "boundary_strong.h"
#include "boundary_rp.h"
#include "boundary_fournier_galatola.h"
#include <stdio.h> 
#include <vector>

double Pi=3.14159265359;
//Base class to latter geometry definition
GEOMETRY::GEOMETRY(const struct Simulation_Parameters * lc) :
    Nx(lc->Nx),
    Ny(lc->Ny),
    Nz(lc->Nz),
    dx_1(1/lc->dx),
    dy_1(1/lc->dy),
    dz_1(1/lc->dz)
{

  std::cout << "Nx=" << Nx << " \n";
  std::cout << "Ny=" << Ny << " \n";
  std::cout << "Nz=" << Nz << " \n\n";


};


  
void GEOMETRY::boundary_init( struct Simulation_Parameters * sim_param)
{

  std::cout <<"\nInitiating boundary conditions:\n\n";
  int ii=0;
  std::string anc_type;
  for (ii=0; ii< number_of_boundaries; ii++)
    {

      try
        {
          anc_type=  sim_param->anchoring_type.at(ii);
          std::cout << "Boundary " << ii << ":" << anc_type <<".\n";
      
        }
      catch(std::out_of_range dummy_var )
        {

          std::cout<< "You must define boundaries "<< boundary_needed_to_be_defined <<" in the " << geometry_name <<" geometry.\nPlease review your input file.\nAborting the program.\n\n";
                   
          exit(0);      
        }

      

      if( strcasecmp(anc_type.c_str(),"strong") == 0 || strcasecmp(anc_type.c_str(),"fixed") == 0  )

        {
    
    
          bc_conditions[ii]=new Strong_Boundary(sim_param, ii);


        }
      else if( strcasecmp(anc_type.c_str(),"rp") == 0 || strcasecmp(anc_type.c_str(),"Rapini-Papoular") == 0  )
        {
    
    
          bc_conditions[ii]=new Boundary_Rp(sim_param, ii);


        }
      else if( strcasecmp(anc_type.c_str(),"fg") == 0 || strcasecmp(anc_type.c_str(),"fournier-galatola") == 0  )
        {
    
    
          bc_conditions[ii]=new Boundary_Fg(sim_param, ii);


        }
      else
        {
          std::cout<<"There is no anchoring type named " << anc_type <<".\nPlease check your input file 'anchoring type' fields.\n\n Aborting the program.'\n";
          exit(0);
        } 
    }  
};
