#include "driver.h"
#include "omp.h"
#include <cstdio>           
#include <cstdlib>          
#include <cstring>
#include <cmath>            
#include "geometry.h"
#include "integrator.h"
#include "integrator_euler.h"
#define cl 2
#define cll 2
#define fixed_chunk_size 1
#define chunk_size 1
#define new_chunk_size 4*Nx*Ny

#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)

//Dormand-Prince integrator
Euler::Euler( GEOMETRY  * lc_pointer, const struct Simulation_Parameters *sim_param ) : Integrator( lc_pointer)
{


  
  
  dt=sim_param->dt;
  //allocate the ith-stage array:
  if((k_1= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  
  
  std::cout << "dt=" << dt << " \n";

};


void Euler::evolve( double * Qij, double *time, double tf )
{

  int i,j,k,ll,information_step=1;
  double dt=this->dt;
  //const int chunk_size=0.06*(4*Ny*Nz)/omp_get_num_threads();
  //const int chunk_size=1;
  
  #pragma omp parallel default(shared) private(i,j,k,ll)
    { 
      while(*time<tf)
        {
        
          #pragma omp for simd schedule(simd:dynamic,new_chunk_size)
          for(ll=0;ll<5*Nx*Ny*Nz; ll++)   Qtij[ll]=Qij[ll];
          
          #pragma omp for schedule(dynamic,fixed_chunk_size) collapse(cl)
          for( k= 0; k< Nz; k++)
            {
              for( j= 0; j< Ny; j++)
                {
                  for( i= 0; i< Nx; i++)
                    {         
                      sample_geometry->fill_ki(k_1,Qtij,i,j,k);           
                      
                    }       
                }
            }
    
          #pragma omp barrier
          
          #pragma omp for simd schedule(simd:dynamic,new_chunk_size)  nowait        
            for( ll=0; ll<5*Nx*Ny*Nz;ll++) Qij[ll]+=dt*k_1[ll];
            
          #pragma omp single nowait
            {
              *time+=dt;
             
              if( information_step%10==0 )
                {
                  std::cout << "time=" << *time << ", dt=" << dt << std::endl;
                  information_step=0;
                }
              information_step++;
              if( (tf-*time) < dt) dt=tf-*time;
            }
        } 
    }    
};
