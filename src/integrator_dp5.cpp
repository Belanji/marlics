#include "driver.h"
#include "omp.h"
#include <cstdio>           
#include <cstdlib>          
#include <cstring>
#include <cmath>            
#include "geometry.h"
#include "integrator.h"
#include "integrator_dp5.h"
#define cl 2
#define cll 2
#define fixed_chunk_size 1
#define chunk_size 1
#define new_chunk_size 4*Nx*Ny

#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)

//Dormand-Prince integrator
DP5::DP5( GEOMETRY  * lc_pointer, const struct Simulation_Parameters *sim_param ) : Integrator( lc_pointer),
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


void DP5::evolve( double * Qij, double *time, double tf )
{

  int i,j,k,ll,information_step=1;
  double local_error;
  double global_error; //, global_error_1=1.;
  double sc_i;
  double hfactor=1.0;
  double dt=this->dt;
  //const int chunk_size=0.06*(4*Ny*Nz)/omp_get_num_threads();
  //const int chunk_size=1;
  
  #pragma omp parallel default(shared) private(i,j,k,ll,local_error,sc_i)
    {
    
      #pragma omp for simd schedule(simd:dynamic,new_chunk_size)
        for(ll=0;ll<5*Nx*Ny*Nz; ll++) 	Qtij[ll]=Qij[ll];
    
    
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
    
      
      while(*time<tf)
        {
        
        //1st-Stage:
        
        //######### 2nd Stage:
          #pragma omp for simd schedule(simd:dynamic,new_chunk_size)
          for(ll=0;ll<5*Nx*Ny*Nz; ll++) Qtij[ll]=Qij[ll]+0.2*dt*k_1[ll]; 
        
            
          #pragma omp for schedule(dynamic,fixed_chunk_size) collapse(cl)
          for( k= 0; k< Nz; k++)
            {
              for( j= 0; j< Ny; j++)
                {
                  for( i= 0; i< Nx; i++)
                    {	    
            
                        sample_geometry->fill_ki(k_2,Qtij,i,j,k);		      
                    
                    }
                }
            }
        //######### 3rd Stage:
          #pragma omp for simd schedule(simd:dynamic,new_chunk_size)
          for(ll=0;ll<5*Nx*Ny*Nz; ll++)   Qtij[ll]=Qij[ll]+dt*(0.075*k_1[ll]+0.225*k_2[ll]); 
                  
        
        
          #pragma omp for schedule(dynamic,fixed_chunk_size) collapse(cl)
          for( k= 0; k< Nz; k++)
            {
              for( j= 0; j< Ny; j++)
                {
                  for( i= 0; i< Nx; i++)
                    {	    
                    
                      sample_geometry->fill_ki(k_3,Qtij,i,j,k);		      
                    
                    }		
                }
            }
        
        //######### 4th Stage:
          #pragma omp for simd schedule(simd:dynamic,new_chunk_size)
          for( ll=0; ll< 5*Nx*Ny*Nz; ll++)  Qtij[ll]=Qij[ll]+dt*(0.9777777777777779*k_1[ll]-3.7333333333333334*k_2[ll]+3.5555555555555554*k_3[ll]); 
        
          #pragma omp for schedule(dynamic,fixed_chunk_size) collapse(cl)
          for( k= 0; k< Nz; k++)
            {
              for( j= 0; j< Ny; j++)
                {
                  for( i= 0; i< Nx; i++)
                    {	    
            
                        sample_geometry->fill_ki(k_4,Qtij,i,j,k);		      
                    
                    }		
                }
            }
        
          //######### 5th Stage:
          #pragma omp for simd schedule(simd:dynamic,new_chunk_size)        			    
          for( ll=0; ll< 5*Nx*Ny*Nz; ll++) Qtij[ll]=Qij[ll]+dt*(2.9525986892242035*k_1[ll]-11.595793324188385*k_2[ll]+9.822892851699436*k_3[ll]-0.2908093278463649*k_4[ll]); 
            
              
          #pragma omp for schedule(dynamic,fixed_chunk_size) collapse(cl)
          for( k= 0; k< Nz; k++)
            {
              for( j= 0; j< Ny; j++)
                {
                  for( i= 0; i< Nx; i++)
                    {	    
              
                      sample_geometry->fill_ki(k_5,Qtij,i,j,k);		      
                      
                  
                    }		
                }
            }
        
          //6th stage:
          #pragma omp for simd schedule(simd:dynamic,new_chunk_size) 
          for( ll=0; ll< 5*Nx*Ny*Nz; ll++)  Qtij[ll]=Qij[ll]+dt*(2.8462752525252526*k_1[ll]-10.757575757575758*k_2[ll]+8.906422717743473*k_3[ll]+0.27840909090909094*k_4[ll]-0.2735313036020583*k_5[ll]); 
                    
        
          #pragma omp for schedule(dynamic,fixed_chunk_size) collapse(cl)
          for( k= 0; k< Nz; k++)
            {
              for( j= 0; j< Ny; j++)
                {
                  for( i= 0; i< Nx; i++)
                    {	    
              
                        sample_geometry->fill_ki(k_6,Qtij,i,j,k);
                  
                    }		
                }
            }
        
        //7th stage:
          #pragma omp for simd schedule(simd:dynamic,new_chunk_size)	
          for( ll=0; ll< 5*Nx*Ny*Nz; ll++) Qtij[ll]=Qij[ll]+dt*(0.09114583333333333*k_1[ll]+0.44923629829290207*k_3[ll]+0.6510416666666666*k_4[ll] -0.322376179245283*k_5[ll]+ 0.13095238095238093*k_6[ll]); 
          
              
          #pragma omp for schedule(dynamic,fixed_chunk_size) collapse(cl)
          for( k= 0; k< Nz; k++)
            {
              for( j= 0; j< Ny; j++)
                {
                  for( i= 0; i< Nx; i++)
                    {	    
              
                      sample_geometry->fill_ki(k_7,Qtij,i,j,k);		      
                    
                    }		
                }
            }
          
          //Restarting global error:
          global_error=0.;
        
          //Calculating error:
          #pragma omp for simd schedule(simd:dynamic,new_chunk_size) reduction(max:global_error)
          for( ll=0; ll<5*Ny*Ny*Nz;ll++)
            {
          
              sc_i=Atol+fabs(Qtij[ll])*Rtol;			
              local_error=0.0012326388888888873*k_1[ll]
                        -0.004252770290506136*k_3[ll]
                        +0.03697916666666656*k_4[ll]
                        -0.05086379716981132*k_5[ll]
                        +0.04190476190476189*k_6[ll]
                        -0.025*k_7[ll];						
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
            
              if( (tf-*time) < dt) dt=tf-*time;
              
            }
        
          if(global_error<1.0)
            {
              
              #pragma omp for simd schedule(simd:dynamic,new_chunk_size)  nowait	    
                for( ll=0; ll<5*Nx*Ny*Nz;ll++) Qij[ll]=Qtij[ll]; 
              
              #pragma omp for simd schedule(simd:dynamic,new_chunk_size)
                for( ll=0; ll<5*Nx*Ny*Nz;ll++) k_1[ll]=k_7[ll]; 
                                
            }
        } 
    }    
};
