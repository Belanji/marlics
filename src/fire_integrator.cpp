#include "driver.h"
#include "omp.h"
#include <cstdio>           
#include <cstdlib>          
#include <cstring>
#include <cmath>            
#include "geometry.h"
#include "integrator.h"
#include "fire_integrator.h"
#define cl 2
#define cll 2
#define fixed_chunk_size 1
#define chunk_size 1
#define new_chunk_size 4*Nx*Ny

#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)

//Dormand-Prince integrator
FIRE::FIRE( GEOMETRY  * lc_pointer, const struct Simulation_Parameters *sim_param ) : Integrator( lc_pointer)
{
  dt=sim_param->dt;
  //allocate the ith-stage array:
  if((force= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((vQij = (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((Qijt = (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((energy = (double *)calloc(2*5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  
  dtmin = dt*0.01;
  dtMax = dt * 100;
  alphainit = 0.3;
  scaleAlpha = 0.99;
  alphaMin = 0.1;
  scaleP =  1.1 ;
  scaleM =  0.9 ;
  Nmin = 10;
  min_force=sim_param->min_force;
  dt_init = dt;
  alpha = alphainit;
  printf("Initializing simulation with FIRE(Fast Inertial Relaxation Engine) method.\n");
  std::cout << "dt=" << dt << " \n";
};

bool FIRE::evolve( double * Qij, double *time, double tf )
{
  int ll,information_step=1, count=0, ccta=0, cc=0;
  double Total_Energy0, Total_Energy1, Qijtest, f2max;
  double dt=dt_init>dtmin?dt_init:dtmin;
  bool continue_condition=true;
  Np=0;
  std::cout << "dt=" << dt <<  ", dtmin=" << dtmin <<  ", dtmax=" << dtMax << " \n";
  #pragma omp parallel default(shared) private(ll)
  {
    #pragma omp barrier
    sample_geometry->compute_forces(force,Qij);
    #pragma omp barrier
    while(*time<tf && continue_condition)
    {
      potence=0; v2=0; f2=0;proceed=1;Qijtest=0;f2max=0;
        
      #pragma omp for simd schedule(simd:dynamic,new_chunk_size)          
        for( ll=0; ll<5*Nx*Ny*Nz;ll++) Qijt[ll]=Qij[ll]+dt*vQij[ll]+0.5*dt*dt*force[ll];
      
      #pragma omp for simd schedule(simd:dynamic,new_chunk_size)          
        for( ll=0; ll<5*Nx*Ny*Nz;ll++) vQij[ll]+=0.5*dt*force[ll];

      #pragma omp barrier
        sample_geometry->compute_forces(force,Qijt);
      #pragma omp barrier
            
      #pragma omp for simd schedule(simd:dynamic,new_chunk_size)          
        for( ll=0; ll<5*Nx*Ny*Nz;ll++) vQij[ll]+=0.5*dt*force[ll];
        
      #pragma omp for simd schedule(simd:dynamic,new_chunk_size) reduction(+: potence)
        for( ll=0; ll<5*Nx*Ny*Nz;ll++) potence+=(sample_geometry->point_type[ll/5]!=0)?vQij[ll]*force[ll]:0;
      
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size) reduction(+: v2)
          for( ll=0; ll<5*Nx*Ny*Nz;ll++) v2+=vQij[ll]*vQij[ll];
        
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size) reduction(+: f2)
          for( ll=0; ll<5*Nx*Ny*Nz;ll++) f2+=force[ll]*force[ll];
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size) reduction(+: Qijtest)
          for( ll=0; ll<5*Nx*Ny*Nz;ll++) Qijtest+=Qij[ll]*Qij[ll];
          
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size) reduction(max: f2max)
          for( ll=0; ll<5*Nx*Ny*Nz;ll++) f2max=MAX(f2max,force[ll]*force[ll]);
          
        if((1+v2+f2)==0||v2/(Nx*Ny*Nz)>0.5){//fprintf(stderr,"NaN found in the system.\nAborting Operation at t=%g %g!\n",*time, dt);
             proceed=0;}
        if (proceed==1)
        {
        scale=(f2>0)?sqrt(v2/f2):0;
      #pragma omp barrier
      if(potence>0){
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size)    nowait      
          for( ll=0; ll<5*Nx*Ny*Nz;ll++) Qij[ll]=Qijt[ll];
          
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size)    nowait      
          for( ll=0; ll<5*Nx*Ny*Nz;ll++) vQij[ll]=(1-alpha)*vQij[ll]+alpha*force[ll]*scale;
          
        if(sqrt(f2max)<min_force) continue_condition=false;
        #pragma omp single
          {
        *time+=dt;
            Np++;
            if(Np>Nmin)
            {
              alpha=MAX(alphaMin,scaleAlpha*alpha);
              dt=min(dt*scaleP,dtMax);
            }
            if( information_step%50==0 )
            {
              printf("time=%g, dt=%g, force=%g, velocity=%g, P=%g, Np=%d\n",*time, dt, sqrt(f2), sqrt(v2), potence, Np);fflush(stdout);
              information_step=0;
            }
            information_step++;
            if(!continue_condition)  printf("Stable condition achivied at time %g\n",*time);
          }
        
      }else{
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size)    nowait      
        for( ll=0; ll<5*Nx*Ny*Nz;ll++) vQij[ll]=0;
        
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size)    nowait      
        for( ll=0; ll<5*Nx*Ny*Nz;ll++) force[ll]=0;
        
        #pragma omp single
        { 
          dtMax*=0.95;
          dt=dtMax*0.95;
          alpha=alphainit;
          Np=0;
          if (dtMax<=dt_init) {printf("convergence probem, exiting the simmulation!"); exit(10);}
        }
      }
      if( (tf-*time) < dt) dt=tf-*time;
      }else{
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size)    nowait      
        for( ll=0; ll<5*Nx*Ny*Nz;ll++) vQij[ll]=0;
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size)    nowait      
        for( ll=0; ll<5*Nx*Ny*Nz;ll++) force[ll]=0;
        
        #pragma omp single
        {
          dt=MAX(0.0001, dt*scaleM);
          alpha=alphainit;
          Np=0;
        }
      }
    }
    Total_Energy0=0;
    sample_geometry->Energy_calc(energy,Qij);
    #pragma omp for simd schedule(simd:dynamic,new_chunk_size) reduction(+: Total_Energy0)
        for( ll=0; ll<2*5*Nx*Ny*Nz;ll++) Total_Energy0+=energy[ll];
        
    #pragma omp barrier
    #pragma omp single 
      std::cout << "time=" << *time << ", Energy=" << Total_Energy0 << std::endl;
            
  }
  return continue_condition;
};
