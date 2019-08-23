#include "driver.h"
#include "geometry.h"
#include "geometry_sphere.h"
#include "boundary.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <ctime>
#include <cstring>
#include <gsl/gsl_randist.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)


Geometry_Sphere::Geometry_Sphere(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param),
									  R_out ((Nx-Nx/2)+1.0),
									  R_in ((Nx-Nx/2)-0.1),
									  HNx ((Nx-1)/2),
									  HNy ((Ny-1)/2),
									  HNz ((Nz-1)/2)
{
  point_type=fill_point_type( );
  geometry_pointer=&(*this);
  geometry_name="Sphere";
  number_of_boundaries=1;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
  
};




Geometry_Sphere::~Geometry_Sphere()
{
  
  
}

 

void  Geometry_Sphere::fill_ki(double * k_i,
		      const double * Qij, 
		      const int i,
		      const int j,
		      const int k)  const 
{

    
  
  int ip1,jp1,kp1, im1, jm1, km1, ll;
  double dQ[15];
  double ddQ[30];
  double QN[5];

  if(point_type[(k*Ny+j)*Nx+i] == 1)
    {

      //check_bulk_limits( i,  j,  k);	
      
      
      
      ip1= (i+1);
      jp1= (j+1);
      kp1= (k+1);
      im1= (i-1);
      jm1= (j-1);
      km1= (k-1);


      for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
	

      //Calcule first derivatives of Qij:
      for(ll=0; ll<=4;ll++) dQ[ll]= 0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;

      
      for(ll=0; ll<=4;ll++) dQ[5+ll]= 0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dx_1;

      
      for(ll=0; ll<=4;ll++) dQ[10+ll]= 0.5*(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dx_1;



      //Calculate second derivatives of Qij:
      for(ll=0; ll<=4;ll++) ddQ[ll]= (Qij[5*((k*Ny+j)*Nx+ip1)+ll]-2.0*QN[ll]+Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1*dx_1;
            
      for(ll=0; ll<=4;ll++) ddQ[ll+5]= 0.25*(Qij[5*((k*Ny+jp1)*Nx+ip1)+ll]+Qij[5*((k*Ny+jm1)*Nx+im1)+ll]-Qij[5*((k*Ny+jp1)*Nx+im1)+ll]-Qij[5*((k*Ny+jm1)*Nx+ip1)+ll])*dx_1*dx_1;
      
      for(ll=0; ll<=4;ll++) ddQ[10+ll]= 0.25*(Qij[5*((kp1*Ny+j)*Nx+ip1)+ll]+Qij[5*((km1*Ny+j)*Nx+im1)+ll]-Qij[5*((kp1*Ny+j)*Nx+im1)+ll]-Qij[5*((km1*Ny+j)*Nx+ip1)+ll])*dx_1*dx_1;
      

      for(ll=0; ll<=4;ll++) ddQ[15+ll]= (Qij[5*((k*Ny+jp1)*Nx+i)+ll]-2.0*QN[ll]+Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dx_1*dx_1;


      for(ll=0; ll<=4;ll++) ddQ[20+ll]= 0.25*(Qij[5*((kp1*Ny+jp1)*Nx+i)+ll]+Qij[5*((km1*Ny+jm1)*Nx+i)+ll]-Qij[5*((kp1*Ny+jm1)*Nx+i)+ll]-Qij[5*((km1*Ny+jp1)*Nx+i)+ll])*dx_1*dx_1;
      
      for(ll=0; ll<=4;ll++) ddQ[25+ll]= (Qij[5*((kp1*Ny+j)*Nx+i)+ll]-2.0*QN[ll]+Qij[5*((km1*Ny+j)*Nx+i)+ll])*dx_1*dx_1;
      
  
      k_i[5*(Nx*(Ny*k+j)+i)+0]=bulk_00(QN,dQ,ddQ); 
      k_i[5*(Nx*(Ny*k+j)+i)+1]=bulk_01(QN,dQ,ddQ); 
      k_i[5*(Nx*(Ny*k+j)+i)+2]=bulk_02(QN,dQ,ddQ); 
      k_i[5*(Nx*(Ny*k+j)+i)+3]=bulk_11(QN,dQ,ddQ); 
      k_i[5*(Nx*(Ny*k+j)+i)+4]=bulk_12(QN,dQ,ddQ);
      
    }  
  else if( point_type[(k*Ny+j)*Nx+i] == 2 )
    {

      //check_surface_limits( i,  j,  k);
      double delta_x,delta_y,delta_z;
      double rr;
      double v[3];

      
      delta_x=i-HNx;
      delta_y=j-HNy;
      delta_z=k-HNz;
      rr=sqrt( delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
			
      v[0]=delta_x/rr;
      v[1]=delta_y/rr;
      v[2]=delta_z/rr;
      //Boundary condtions:

      ip1= i>=Nx/2 ? i   : i+1 ;
      jp1= j>=Ny/2 ? j   : j+1 ;
      kp1= k>=Nz/2 ? k   : k+1 ;
  
      im1= i>=Nx/2 ? i-1 : i ;
      jm1= j>=Ny/2 ? j-1 : j ;
      km1= k>=Nz/2 ? k-1 : k ;
  

      for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
	

      //Calcule first derivatives of Qij:
      for(ll=0; ll<=4;ll++) dQ[ll]=   (Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;

      
      for(ll=0; ll<=4;ll++) dQ[5+ll]= (Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dx_1;

      
      for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dx_1;
      
  
      k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions[0]->surface_00(QN,dQ,ddQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+1]= bc_conditions[0]->surface_01(QN,dQ,ddQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+2]= bc_conditions[0]->surface_02(QN,dQ,ddQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+3]= bc_conditions[0]->surface_11(QN,dQ,ddQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+4]= bc_conditions[0]->surface_12(QN,dQ,ddQ,v); 


      //check_surface_limits(i,j,k);
    }
  else
    {

      k_i[5*(Nx*(Ny*k+j)+i)+0]=0.; 
      k_i[5*(Nx*(Ny*k+j)+i)+1]=0.; 
      k_i[5*(Nx*(Ny*k+j)+i)+2]=0.; 
      k_i[5*(Nx*(Ny*k+j)+i)+3]=0.; 
      k_i[5*(Nx*(Ny*k+j)+i)+4]=0.;
    }
  
}
      
      
int * Geometry_Sphere::fill_point_type( void )  const 
{
  int i,j,k;
  double delta_x,delta_y,delta_z;
  double rr;
  int * point_kind;

  
  if((point_kind= (int *)calloc(Nx*Ny*Nz, sizeof(int)))==NULL){ERROr}

  for( i= 0; i< Nx; i++)
    {
      for( j= 0; j< Ny; j++)
	{
	  for( k= 0; k< Nz; k++)
	    {	    
	      delta_x=i-HNx;
	      delta_y=j-HNy;
	      delta_z=k-HNz;
	      rr=sqrt( delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
	      

	      if(rr < R_in)
		{		  		  

		  point_kind[(k*Ny+j)*Nx+i]=1;
			  
		}
	      else if(rr >=R_in && rr <=R_out)    
		{

		  point_kind[(k*Ny+j)*Nx+i]=2;
		  
		}
	      else
		{

		  point_kind[(k*Ny+j)*Nx+i]=0;
		      
		}
	
	    }
	      
	}	  
    }

  return point_kind;
}




  

