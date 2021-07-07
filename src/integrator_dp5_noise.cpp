#include "driver.h"
#include "omp.h"
#include <cstdio>           
#include <cstdlib>          
#include <cstring>
#include <cmath>            
#include <gsl/gsl_rng.h>            
#include <gsl/gsl_randist.h>            
#include "geometry.h"
#include "integrator.h"
#include "integrator_dp5_noise.h"
#define cl 2
#define cll 2
#define fixed_chunk_size 1
#define chunk_size 1
#define new_chunk_size 4*Nx*Ny

#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)

//Dormand-Prince integrator
NOISE_DP5::NOISE_DP5( GEOMETRY  * lc_pointer, const struct Simulation_Parameters *sim_param ) : Integrator( lc_pointer),
                                   Atol(sim_param->Atol),
                                   Rtol(sim_param->Rtol),
                                   prefac(sim_param->prefac),
                                   facmin(sim_param->facmin),
                                                                   facmax(sim_param->facmax)
{
  dt=sim_param->dt;
  //allocate the ith-stage array:
  if((k_1= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((k_2= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((k_3= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((k_4= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((k_5= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((k_6= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((k_7= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((noise= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((energy = (double *)calloc(2*5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  
  if (sim_param->ic_flag[4]==parameter_status::unset) 
    //~ gsl_rng_default_seed=time(NULL);
    gsl_rng_default_seed=1;
  else gsl_rng_default_seed= sim_param->rng_seed;
//  gsl_rng_default_seed=1570800053;
  
  std::cout << "seed = " <<gsl_rng_default_seed<<std::endl;;
  std::cout << "noise_factor = " <<noise_factor<<std::endl;
  noise_factor=sim_param->noise_factor;
  std::cout << "dt=" << dt << " \n";

  if(sim_param->integrator_parameters_flag<5)
    {
      std::cout << "You defined " <<  sim_param->integrator_parameters_flag << " of the 5 optional paramters of the DOrmand-Prince integrator.\nThe remaining values will be set to their standardvalues.\n";
    }
  
  std::cout << "Atol=" << Atol << " \n";
  std::cout << "Rtol=" << Rtol << " \n";
  std::cout << "facmin=" << facmin << "\n";
  std::cout << "facmax=" << facmax << "\n";
  std::cout << "prefac=" << prefac << "\n\n";

};


bool NOISE_DP5::evolve( double * Qij, double *time, double tf )
{
  int ll,information_step=1, counter=0;
  double local_error,n[3];
  double Total_Energy;
  double global_error; //, global_error_1=1.;
  double sc_i;
  double hfactor=1.0;
  double dt=this->dt;
  const int * point_type=sample_geometry->point_type;
  static const gsl_rng *w= gsl_rng_alloc(gsl_rng_taus);
  //const int chunk_size=0.06*(4*Ny*Nz)/omp_get_num_threads();
  //const int chunk_size=1;
  #pragma omp parallel default(shared) private(ll,local_error,sc_i)
    {
      //1st-Stage:
      #pragma omp for simd schedule(simd:dynamic,new_chunk_size)
        for(ll=0;ll<5*Nx*Ny*Nz; ll++)   Qtij[ll]=Qij[ll];
      sample_geometry->fill_ki(k_1,Qtij);             
      
      while(*time<tf)
        {
        //######### 2nd Stage:
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size)
          for(ll=0;ll<5*Nx*Ny*Nz; ll++) Qtij[ll]=Qij[ll]+0.222222222*dt*k_1[ll]; 
          
        sample_geometry->fill_ki(k_2,Qtij);

        //######### 3rd Stage:
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size)
          for(ll=0;ll<5*Nx*Ny*Nz; ll++)   Qtij[ll]=Qij[ll]+dt*(0.083333333*k_1[ll]+0.25*k_2[ll]); 
          
        sample_geometry->fill_ki(k_3,Qtij);           
          
        //######### 4th Stage:
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size)
          for( ll=0; ll< 5*Nx*Ny*Nz; ll++)  Qtij[ll]=Qij[ll]+dt*(0.169753086*k_1[ll]-0.231481481*k_2[ll]+0.617283951*k_3[ll]); 
          
        sample_geometry->fill_ki(k_4,Qtij);           
      
        //######### 5th Stage:
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size)                      
          for( ll=0; ll< 5*Nx*Ny*Nz; ll++) Qtij[ll]=Qij[ll]+dt*(0.251515152*k_1[ll]-0.590909091*k_2[ll]+0.924242424*k_3[ll]+0.081818182*k_4[ll]); 
          
        sample_geometry->fill_ki(k_5,Qtij);           
          
        //######### 6th stage:
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size) 
          for( ll=0; ll< 5*Nx*Ny*Nz; ll++)  Qtij[ll]=Qij[ll]+dt*(-0.678571429*k_1[ll]+2.25*k_2[ll]+0.142857143*k_3[ll]-3.857142857*k_4[ll]+3.142857143*k_5[ll]); 
          
        sample_geometry->fill_ki(k_6,Qtij);
        
        //######### 7th stage:
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size)
          for( ll=0; ll< 5*Nx*Ny*Nz; ll++) Qtij[ll]=Qij[ll]+dt*(0.095*k_1[ll]+0.6*k_3[ll]-0.6075*k_4[ll] +0.825*k_5[ll]+ 0.0875*k_6[ll]); 
          
        sample_geometry->fill_ki(k_7,Qtij);           
          
        //Restarting global error:
        global_error=0.;
        
        #pragma omp barrier
      
        //Calculating error:
        #pragma omp for simd schedule(simd:dynamic,new_chunk_size) reduction(max:global_error)
          for( ll=0; ll<5*Ny*Ny*Nz;ll++)
            {
          
              sc_i=Atol+fabs(Qtij[ll])*Rtol;            
              local_error=-0.0088*k_1[ll]
                          +0.066*k_3[ll]
                          -0.1782*k_4[ll]
                          +0.132*k_5[ll]
                          +0.009*k_6[ll]
                          -0.02*k_7[ll];                        
              global_error=MAX( fabs(local_error)/sc_i,global_error );
                      
            }
      
            #pragma omp single nowait
              {
              
                if(global_error<1.0)
                  {
                    *time+=dt;
              
                    if( information_step%10==0 )
                      {
                        std::cout << "time=" << *time << ", dt=" << dt << ", global_error=" << global_error <<
                        ", hfactor=" << hfactor << std::endl;
                        information_step=0;
                      }
                    information_step++;
              
                  };
              
                hfactor=min(facmax,MAX(facmin,prefac*pow(global_error,-0.2000000000)));
                dt=dt*hfactor;
                if(dt<1e-14)
                  {
                    printf("Convergence problem!!!\n Aborting execution!!!\n");
                    exit(5);
                  }      
                if( (tf-*time) < dt) dt=tf-*time;
                
              }
        
          if(global_error<1.0)
            {
              if (counter%50==0)
              {
                
              //~ #pragma omp for simd schedule(simd:dynamic,new_chunk_size)  nowait
  
                #pragma omp single
                {
                for(ll= 0; ll< Nx*Ny*Nz; ll++)
                {
                  if(point_type[ll] !=0 )
                    {
                      gsl_ran_dir_3d(w, &n[0], &n[1], &n[2]);
                      noise[5*ll+0]=(0.5*noise_factor*(3.0*n[0]*n[0]-1.0));
                      noise[5*ll+1]=(0.5*noise_factor*(3.0*n[0]*n[1]));
                      noise[5*ll+2]=(0.5*noise_factor*(3.0*n[0]*n[2]));
                      noise[5*ll+3]=(0.5*noise_factor*(3.0*n[1]*n[1]-1.0));
                      noise[5*ll+4]=(0.5*noise_factor*(3.0*n[1]*n[2]));
                      
                    }           
                  else
                    {                    
                      n[0]=0;
                      n[1]=0;
                      n[2]=1;
                      noise[5*ll+0]=0.0;
                      noise[5*ll+1]=0.0;
                      noise[5*ll+2]=0.0;
                      noise[5*ll+3]=0.0;
                      noise[5*ll+4]=0.0;  
                    
                    }     
                }
              }
                #pragma omp for simd schedule(simd:dynamic,new_chunk_size)          
                  for( ll=0; ll<5*Nx*Ny*Nz;ll++) 
                  {
                    Qij[ll]=Qtij[ll]+noise[ll];
                    Qtij[ll]=Qij[ll];
                  }
                #pragma omp barrier
                  sample_geometry->fill_ki(k_1,Qtij);
              }else{
                #pragma omp for simd schedule(simd:dynamic,new_chunk_size)  nowait      
                  for( ll=0; ll<5*Nx*Ny*Nz;ll++) Qij[ll]=Qtij[ll]; 
                
                #pragma omp for simd schedule(simd:dynamic,new_chunk_size)
                  for( ll=0; ll<5*Nx*Ny*Nz;ll++) k_1[ll]=k_7[ll]; 
              }
              #pragma omp single
              counter++;
                                
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
