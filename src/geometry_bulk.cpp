#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <ctime>
#include <cstring>

#include "driver.h"
#include "geometry.h"
#include "geometry_bulk.h"
#include <gsl/gsl_randist.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)


Geometry_Bulk::Geometry_Bulk(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param)
{
  point_type=fill_point_type( );
  geometry_pointer=&(*this);

  geometry_name="Bulk";
  number_of_boundaries=0;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);

  
  
};


void  Geometry_Bulk::fill_ki(double * k_i,
                    const double * Qij, 
                    const int i,
                    const int j,
                    const int k)  const 
{

    
  
  int ip1,jp1,kp1, im1, jm1, km1, ll;
  double dQ[15];
  double ddQ[30];
  double QN[5];
  double v[3];
  
  
  ip1= (i+1)%Nx;
  jp1= (j+1)%Ny;
  kp1= (k+1)%Nz;
  im1= i-1+((Nx-1-i)/(Nx-1))*Nx;
  jm1= j-1+((Ny-1-i)/(Ny-1))*Ny;
  km1= k-1+((Nz-1-i)/(Nz-1))*Nz;

      
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
      
  
      k_i[5*(Nx*(Ny*k+j)+i)+0]=bulk_00(QN,dQ,ddQ); 
      k_i[5*(Nx*(Ny*k+j)+i)+1]=bulk_01(QN,dQ,ddQ); 
      k_i[5*(Nx*(Ny*k+j)+i)+2]=bulk_02(QN,dQ,ddQ); 
      k_i[5*(Nx*(Ny*k+j)+i)+3]=bulk_11(QN,dQ,ddQ); 
      k_i[5*(Nx*(Ny*k+j)+i)+4]=bulk_12(QN,dQ,ddQ);
  
}
      
      
int * Geometry_Bulk::fill_point_type( void )  const 
{
  int i,j,k;
  int * point_kind;

  
  if((point_kind= (int *)calloc(Nx*Ny*Nz, sizeof(int)))==NULL){ERROr}

  for( i= 0; i< Nx; i++)
    {
      for( j= 0; j< Ny; j++)
          {
            for( k= 0; k< Nz; k++)
              {     
                
                point_kind[(k*Ny+j)*Nx+i]=1;
                                                              
              }
          }    
    }     
    
  return point_kind;
}






