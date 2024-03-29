#include "driver.h"
#include "omp.h"
#include <cstdio>           
#include <cstdlib>          
#include <cstring>
#include <cmath>            
#include "geometry.h"
#include "integrator.h"
#include "integrator_rk2.h"
#define cl 2
#define new_chunk_size 4*Nx*Ny

#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)

//Dormand-Prince integrator
RK2::RK2( GEOMETRY  * lc_pointer, const struct Simulation_Parameters *sim_param ) : Integrator( lc_pointer)
{
  dt=sim_param->dt;
  //allocate the ith-stage array:
  if((k_1= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((k_2= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((energy = (double *)calloc(2*5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  
  std::cout << "dt=" << dt << " \n";

};


bool RK2::evolve( double * Qij, double *time, double tf )
{

  int ll,information_step=1;
  double Total_Energy;
  double dt=this->dt;
  //const int chunk_size=0.06*(4*Ny*Nz)/omp_get_num_threads();
  //const int chunk_size=1;
  
  #pragma omp parallel default(shared) private(ll)
    {
      
      while(*time<tf)
        {
        
        //1st Stage:
        sample_geometry->fill_ki(k_1,Qij);     
        
        //2nd Stage:
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size)
        for(ll=0;ll<5*Nx*Ny*Nz; ll++) Qtij[ll]=Qij[ll]+0.5*dt*k_1[ll]; 
        sample_geometry->fill_ki(k_2,Qtij);             

        #pragma omp for simd schedule(simd:dynamic,new_chunk_size)  nowait        
          for( ll=0; ll<5*Nx*Ny*Nz;ll++) Qij[ll]+=dt*k_2[ll]; 

        #pragma omp single 
            {
              *time+=dt;
              if( information_step%25==0 )
                {
                  std::cout << "time=" << *time << ", dt=" << dt << std::endl;
                  information_step=0;
                }
              information_step++;
              
            }
        }        
        Total_Energy=0;
        sample_geometry->Energy_calc(energy,Qij);
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size) reduction(+: Total_Energy)
            for( ll=0; ll<2*5*Nx*Ny*Nz;ll++) Total_Energy+=energy[ll];
            
        #pragma omp barrier
        #pragma omp single 
          std::cout << "time=" << *time << ", Energy=" << Total_Energy << std::endl;
    }
  return true;    
};
