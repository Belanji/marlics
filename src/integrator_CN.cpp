#include "driver.h"
#include "omp.h"
#include <cstdio>           
#include <cstdlib>          
#include <cstring>
#include <cmath>            
#include "geometry.h"
#include "integrator.h"
#include "integrator_CN.h"






CN::CN( GEOMETRY  * lc_pointer, const struct Simulation_Parameters *sim_param ) : Integrator( lc_pointer),
							 	   Atol(sim_param->Atol),
								   Rtol(sim_param->Rtol),
								   prefac(sim_param->prefac),
								   facmin(sim_param->facmin),
                                                                   facmax(sim_param->facmax)
{

  
  dt=sim_param->dt;
  //allocate the ith-stage array:
  
  

};



void CN::evolve( double * Qij, double *time, double tf )
{

  int i,j,k,ll,information_step=1;
  double local_error;
  double global_error; //, global_error_1=1.;
  double sc_i;
  double hfactor=1.0;
  double dt=this->dt;
  //const int chunk_size=0.06*(4*Ny*Nz)/omp_get_num_threads();
  //const int chunk_size=1;
  

  
  
  *time+=tf;

	      

}; 

CN::~CN(void)
{


  
}
