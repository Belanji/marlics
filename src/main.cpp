#include "driver.h"
#include "geometry.h"
#include  "container.h"
#include "integrator.h"
#include <cstdio>           
#include <cstdlib>          
#include <cstring>
#include <cmath>            


int main(int argc, char **argv)
{

  //debuging options:
  //feenableexcept(FE_INVALID   | 
  //		 FE_DIVBYZERO |
  //		 FE_OVERFLOW  | 
  //		 FE_UNDERFLOW);

  double t,tf;
  driver Lc_driver=driver();
  
  printf("Welcome to MarLiCS software v0.1 \n");
  
  
  
  
  //Parsing line command options:
  
  Lc_driver.parse_input_file();
  Lc_driver.setup_LC();
  Lc_driver.setup_Simulation();

  fflush(stdout);

  t=Lc_driver.sim_param.ti;  
  tf=Lc_driver.sim_param.tf;
  


  //Check if the initial conditions were gotten from a file.
  if(Lc_driver.sim_param.ic_flag[1]==parameter_status::unset)
   {
      //Print a snaphot if the ic was not gotten from a file.
      Lc_driver.Data_Container->write_state(t, Lc_driver.Qij, Lc_driver.LcS_Geometry->point_type);
   }
  

  while(t<tf)
    {

      Lc_driver.LcS_Integrator->evolve(Lc_driver.Qij, & t , Lc_driver.sim_param.timeprint );//! timestep evolution
      Lc_driver.Data_Container->write_state(t, Lc_driver.Qij, Lc_driver.LcS_Geometry->point_type);
  
      Lc_driver.update_timeprint( ) ;
      
    }
 
  return 0;
  
}
