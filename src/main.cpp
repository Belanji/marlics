#include "driver.h"
#include "geometry.h"
//#include "integrator.h"
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
  
  printf("Welcome to MarLiCS software v0.04 \n");
  
  
  
  
  //Parsing line command options:
  
  Lc_driver.parse_input_file();
  Lc_driver.setup_LC();
  Lc_driver.setup_Simulation();


  t=Lc_driver.sim_param.ti;  
  tf=Lc_driver.sim_param.tf;
  



  //if(Lc_driver.lc_prop.ic_file_flag==0)
  //  {
  //    write_state(t, Lc_driver.Qij, Lc_driver.LcS_Geometry->point_type);
  //  }
  //else 
  //  {
  //    printf("%lf\n", t);
  //  }
  //
  //while(t<tf)
  //  {
  //
  //    Lc_driver.LcS_Integrator->evolve(Lc_driver.Qij, & t , Lc_driver.lc_prop.timeprint );//! timestep evolution
  //    write_state(t, Lc_driver.Qij, Lc_driver.LcS_Geometry->point_type);
  //
  //    Lc_driver.update_timeprint( ) ;
  //    
  //  }
  //
  //
  //
  //finish_container();
  //  //ncio_finish_container();

  return 0;
  
}


