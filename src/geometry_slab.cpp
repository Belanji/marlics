#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <ctime>
#include <cstring>

#include "driver.h"
#include "geometry.h"
#include "geometry_slab.h"
#include "boundary.h"
#include "energy_ldg.h"
#include <gsl/gsl_randist.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)


slab::slab(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param)
{
  point_type=fill_point_type( );
  geometry_pointer=&(*this);

  geometry_name="Slab";
  number_of_boundaries=2;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
  boundary_needed_to_be_defined="0 and 1";
  
  
};


void  slab::fill_ki(double * k_i,
                    const double * Qij)  const 
{

    
  int i ,j, k;
  int ip1,jp1,kp1, im1, jm1, km1, ll;
  double dQ[15];
  double ddQ[30];
  double QN[5];
  double v[3];

  k=0;
#pragma omp for schedule(dynamic,fixed_chunk_size) collapse(cl)
  for( j= 0; j< Ny; j++)
    {
      for( i= 0; i< Nx; i++)
	{
	  
	  ip1= (i+1)%Nx;
	  jp1= (j+1)%Ny;
	  kp1= (k+1);
	  im1= (i+Nx-1)%Nx;
	  jm1= (j+Ny-1)%Ny;
	  km1= k;
	  v[0]=0.0;
	  v[1]=0.0;
	  v[2]=-1.0;

	  for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
        

	  //Calcule first derivatives of Qij:
	  for(ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;

      
	  for(ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;

      
	  for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((k*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
      
	  k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions[0]->functional_derivative_00(QN,dQ,ddQ,v);
	  k_i[5*(Nx*(Ny*k+j)+i)+1]= bc_conditions[0]->functional_derivative_01(QN,dQ,ddQ,v);
	  k_i[5*(Nx*(Ny*k+j)+i)+2]= bc_conditions[0]->functional_derivative_02(QN,dQ,ddQ,v);
	  k_i[5*(Nx*(Ny*k+j)+i)+3]= bc_conditions[0]->functional_derivative_11(QN,dQ,ddQ,v);
	  k_i[5*(Nx*(Ny*k+j)+i)+4]= bc_conditions[0]->functional_derivative_12(QN,dQ,ddQ,v);
	

	}
    }
	  
#pragma omp for schedule(dynamic,fixed_chunk_size) collapse(cl)
  for( k= 1; k< Nz-1; k++)
    {
      for( j= 0; j< Ny; j++)
	{
	  for( i= 0; i< Nx; i++)
	    {	
	      

	      //check_bulk_limits( i,  j,  k);  

	      ip1= (i+1)%Nx;
	      jp1= (j+1)%Ny;
	      kp1= (k+1);
	      im1= (i+Nx-1)%Nx;
	      jm1= (j+Ny-1)%Ny;
	      km1= (k-1);

      
	      for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
        

	      //Calcule first derivatives of Qij:
	      for(ll=0; ll<=4;ll++) dQ[ll]= 0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;

      
	      for(ll=0; ll<=4;ll++) dQ[5+ll]= 0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;

      
	      for(ll=0; ll<=4;ll++) dQ[10+ll]= 0.5*(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;



	      //Calculate second derivatives of Qij:
	      for(ll=0; ll<=4;ll++) ddQ[ll]= (Qij[5*((k*Ny+j)*Nx+ip1)+ll]-2.0*QN[ll]+Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1*dx_1;
            
	      for(ll=0; ll<=4;ll++) ddQ[ll+5]= 0.25*(Qij[5*((k*Ny+jp1)*Nx+ip1)+ll]+Qij[5*((k*Ny+jm1)*Nx+im1)+ll]-Qij[5*((k*Ny+jp1)*Nx+im1)+ll]-Qij[5*((k*Ny+jm1)*Nx+ip1)+ll])*dx_1*dy_1;
      
	      for(ll=0; ll<=4;ll++) ddQ[10+ll]= 0.25*(Qij[5*((kp1*Ny+j)*Nx+ip1)+ll]+Qij[5*((km1*Ny+j)*Nx+im1)+ll]-Qij[5*((kp1*Ny+j)*Nx+im1)+ll]-Qij[5*((km1*Ny+j)*Nx+ip1)+ll])*dx_1*dz_1;
      

	      for(ll=0; ll<=4;ll++) ddQ[15+ll]= (Qij[5*((k*Ny+jp1)*Nx+i)+ll]-2.0*QN[ll]+Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1*dy_1;


	      for(ll=0; ll<=4;ll++) ddQ[20+ll]= 0.25*(Qij[5*((kp1*Ny+jp1)*Nx+i)+ll]+Qij[5*((km1*Ny+jm1)*Nx+i)+ll]-Qij[5*((kp1*Ny+jm1)*Nx+i)+ll]-Qij[5*((km1*Ny+jp1)*Nx+i)+ll])*dy_1*dz_1;
      
	      for(ll=0; ll<=4;ll++) ddQ[25+ll]= (Qij[5*((kp1*Ny+j)*Nx+i)+ll]-2.0*QN[ll]+Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1*dz_1;
      
	      k_i[5*(Nx*(Ny*k+j)+i)+0]=bulk_energy.functional_derivative_00(QN,dQ,ddQ,v); 
	      k_i[5*(Nx*(Ny*k+j)+i)+1]=bulk_energy.functional_derivative_01(QN,dQ,ddQ,v); 
	      k_i[5*(Nx*(Ny*k+j)+i)+2]=bulk_energy.functional_derivative_02(QN,dQ,ddQ,v); 
	      k_i[5*(Nx*(Ny*k+j)+i)+3]=bulk_energy.functional_derivative_11(QN,dQ,ddQ,v); 
	      k_i[5*(Nx*(Ny*k+j)+i)+4]=bulk_energy.functional_derivative_12(QN,dQ,ddQ,v);
	    }
	}
    }
      


  k=Nz-1;
#pragma omp for schedule(dynamic,fixed_chunk_size) collapse(cl)  
  for( j= 0; j< Ny; j++)
    {
      for( i= 0; i< Nx; i++)
	{		          
	  ip1= (i+1)%Nx;
	  jp1= (j+1)%Ny;
	  kp1= k;
	  im1= (i+Nx-1)%Nx;
	  jm1= (j+Ny-1)%Ny;
	  km1= (k-1);
	  v[0]=0.0;
	  v[1]=0.0;
	  v[2]=1.0;

      
	  for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
        

	  //Calcule first derivatives of Qij:
	  for(ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;

      
	  for(ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;

      
	  for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((k*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
      
	  k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions[1]->functional_derivative_00(QN,dQ,ddQ,v);
	  k_i[5*(Nx*(Ny*k+j)+i)+1]= bc_conditions[1]->functional_derivative_01(QN,dQ,ddQ,v);
	  k_i[5*(Nx*(Ny*k+j)+i)+2]= bc_conditions[1]->functional_derivative_02(QN,dQ,ddQ,v);
	  k_i[5*(Nx*(Ny*k+j)+i)+3]= bc_conditions[1]->functional_derivative_11(QN,dQ,ddQ,v);
	  k_i[5*(Nx*(Ny*k+j)+i)+4]= bc_conditions[1]->functional_derivative_12(QN,dQ,ddQ,v);
		 

	}
    }
	

	      
  //check_surface_limits(i,j,k);
}
      
      
int * slab::fill_point_type( void )  const 
{
  int i,j,k;
  int * point_kind;

  
  if((point_kind= (int *)calloc(Nx*Ny*Nz, sizeof(int)))==NULL){ERROr}

  for( i= 0; i< Nx; i++)
    {
      for( j= 0; j< Ny; j++)
        {

          k=0;
          point_kind[(k*Ny+j)*Nx+i]=2;

          for( k= 1; k< Nz-1; k++)
            {       
              
              point_kind[(k*Ny+j)*Nx+i]=1;
                                                              
            }
          
          k=Nz-1;
          point_kind[(k*Ny+j)*Nx+i]=3;
        }
    }     
    
  return point_kind;
}






