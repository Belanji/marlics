#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <ctime>
#include <cstring>

#include "driver.h"
#include "geometry.h"
#include "geometry_bulk.h"
#include "energy_ldg.h"
#include <gsl/gsl_randist.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)

#define chunk_size 4*Nx*Ny


Geometry_Bulk::Geometry_Bulk(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param)
{
  point_type=fill_point_type( );
  geometry_pointer=&(*this);

  geometry_name="Bulk";
  number_of_boundaries=0;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);

};


void  Geometry_Bulk::fill_ki(double * k_i,const double * Qij) const
{

#pragma omp for schedule(dynamic,1) collapse(2)
  for( int k= 0; k< Nz; k++)
    {
      for( int j= 0; j< Ny; j++)
        {
          for( int i= 0; i< Nx; i++)
            {	
              int ip1= (i+1)%Nx;
              int im1= (i+Nx-1)%Nx;

              int jp1= (j+1)%Ny;
              int jm1= (j+Ny-1)%Ny;

              int kp1= (k+1)%Nz;
              int km1= (k+Nz-1)%Nz;

              double dQ[15];
              double ddQ[30];
              double QN[5];
              double v[3];

              v[0]=0;
              v[1]=0;
              v[2]=0;
      
              for(int ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
      
              //Calcule first derivatives of Qij:
              for(int ll=0; ll<=4;ll++) dQ[ll]= 0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
      
              for(int ll=0; ll<=4;ll++) dQ[5+ll]= 0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
      
              for(int ll=0; ll<=4;ll++) dQ[10+ll]= 0.5*(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
      
              //Calculate second derivatives of Qij:
              for(int ll=0; ll<=4;ll++) ddQ[ll]= (Qij[5*((k*Ny+j)*Nx+ip1)+ll]-2.0*QN[ll]+Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1*dx_1;
                  
              for(int ll=0; ll<=4;ll++) ddQ[5+ll]= 0.25*(Qij[5*((k*Ny+jp1)*Nx+ip1)+ll]+Qij[5*((k*Ny+jm1)*Nx+im1)+ll]-Qij[5*((k*Ny+jp1)*Nx+im1)+ll]-Qij[5*((k*Ny+jm1)*Nx+ip1)+ll])*dx_1*dy_1;
            
              for(int ll=0; ll<=4;ll++) ddQ[10+ll]= 0.25*(Qij[5*((kp1*Ny+j)*Nx+ip1)+ll]+Qij[5*((km1*Ny+j)*Nx+im1)+ll]-Qij[5*((kp1*Ny+j)*Nx+im1)+ll]-Qij[5*((km1*Ny+j)*Nx+ip1)+ll])*dx_1*dz_1;
      
              for(int ll=0; ll<=4;ll++) ddQ[15+ll]= (Qij[5*((k*Ny+jp1)*Nx+i)+ll]-2.0*QN[ll]+Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1*dy_1;
      
              for(int ll=0; ll<=4;ll++) ddQ[20+ll]= 0.25*(Qij[5*((kp1*Ny+jp1)*Nx+i)+ll]+Qij[5*((km1*Ny+jm1)*Nx+i)+ll]-Qij[5*((kp1*Ny+jm1)*Nx+i)+ll]-Qij[5*((km1*Ny+jp1)*Nx+i)+ll])*dy_1*dz_1;
            
              for(int ll=0; ll<=4;ll++) ddQ[25+ll]= (Qij[5*((kp1*Ny+j)*Nx+i)+ll]-2.0*QN[ll]+Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1*dz_1;
        
              k_i[5*(Nx*(Ny*k+j)+i)+0]=bulk_energy.functional_derivative_00(QN,dQ,ddQ,v); 
              k_i[5*(Nx*(Ny*k+j)+i)+1]=bulk_energy.functional_derivative_01(QN,dQ,ddQ,v); 
              k_i[5*(Nx*(Ny*k+j)+i)+2]=bulk_energy.functional_derivative_02(QN,dQ,ddQ,v); 
              k_i[5*(Nx*(Ny*k+j)+i)+3]=bulk_energy.functional_derivative_11(QN,dQ,ddQ,v); 
              k_i[5*(Nx*(Ny*k+j)+i)+4]=bulk_energy.functional_derivative_12(QN,dQ,ddQ,v);

            }
        }
    }
}

void  Geometry_Bulk::compute_forces(double * k_i,const double * Qij) const
{
#pragma omp for schedule(dynamic,1) collapse(2)
  for( int k= 0; k< Nz; k++)
  {
    for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
	    {	
	      int ip1= (i+1)%Nx;
        int im1= (i+Nx-1)%Nx;
        int jp1= (j+1)%Ny;
        int jm1= (j+Ny-1)%Ny;
        int kp1= (k+1)%Nz;
        int km1= (k+Nz-1)%Nz;

        double dQ[105];
        double QN[35];
      
	      for(int ll=0; ll<=4;ll++) QN[ll+5*(0)]=Qij[5*(Nx*(Ny*km1+j  )+i  )+ll];
	      for(int ll=0; ll<=4;ll++) QN[ll+5*(1)]=Qij[5*(Nx*(Ny*k  +jm1)+i  )+ll];
	      for(int ll=0; ll<=4;ll++) QN[ll+5*(2)]=Qij[5*(Nx*(Ny*k  +j  )+im1)+ll];
	      for(int ll=0; ll<=4;ll++) QN[ll+5*(3)]=Qij[5*(Nx*(Ny*k  +j  )+i  )+ll];
	      for(int ll=0; ll<=4;ll++) QN[ll+5*(4)]=Qij[5*(Nx*(Ny*k  +j  )+ip1)+ll];
	      for(int ll=0; ll<=4;ll++) QN[ll+5*(5)]=Qij[5*(Nx*(Ny*k  +jp1)+i  )+ll];
	      for(int ll=0; ll<=4;ll++) QN[ll+5*(6)]=Qij[5*(Nx*(Ny*kp1+j  )+i  )+ll];
	    
        for(int ll=0; ll<=4;ll++) dQ[0+3*(ll+5*(0))]=dx_1*(Qij[5*(Nx*(Ny*km1+j  )+ip1)+ll]-Qij[5*(Nx*(Ny*km1+j  )+im1)+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[0+3*(ll+5*(1))]=dx_1*(Qij[5*(Nx*(Ny*k  +jm1)+ip1)+ll]-Qij[5*(Nx*(Ny*k  +jm1)+im1)+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[0+3*(ll+5*(2))]=dx_1*(Qij[5*(Nx*(Ny*k  +j  )+i  )+ll]-Qij[5*(Nx*(Ny*k  +j  )+im1)+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[0+3*(ll+5*(3))]=dx_1*(Qij[5*(Nx*(Ny*k  +j  )+ip1)+ll]-Qij[5*(Nx*(Ny*k  +j  )+im1)+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[0+3*(ll+5*(4))]=dx_1*(Qij[5*(Nx*(Ny*k  +j  )+ip1)+ll]-Qij[5*(Nx*(Ny*k  +j  )+i  )+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[0+3*(ll+5*(5))]=dx_1*(Qij[5*(Nx*(Ny*k  +jp1)+ip1)+ll]-Qij[5*(Nx*(Ny*k  +jp1)+im1)+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[0+3*(ll+5*(6))]=dx_1*(Qij[5*(Nx*(Ny*kp1+j  )+ip1)+ll]-Qij[5*(Nx*(Ny*kp1+j  )+im1)+ll]);
                                                                                          
        for(int ll=0; ll<=4;ll++) dQ[1+3*(ll+5*(0))]=dy_1*(Qij[5*(Nx*(Ny*km1+jp1)+i  )+ll]-Qij[5*(Nx*(Ny*km1+jm1)+i  )+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[1+3*(ll+5*(1))]=dy_1*(Qij[5*(Nx*(Ny*k  +j  )+i  )+ll]-Qij[5*(Nx*(Ny*k  +jm1)+i  )+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[1+3*(ll+5*(2))]=dy_1*(Qij[5*(Nx*(Ny*k  +jp1)+im1)+ll]-Qij[5*(Nx*(Ny*k  +jm1)+im1)+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[1+3*(ll+5*(3))]=dy_1*(Qij[5*(Nx*(Ny*k  +jp1)+i  )+ll]-Qij[5*(Nx*(Ny*k  +jm1)+i  )+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[1+3*(ll+5*(4))]=dy_1*(Qij[5*(Nx*(Ny*k  +jp1)+ip1)+ll]-Qij[5*(Nx*(Ny*k  +jm1)+ip1)+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[1+3*(ll+5*(5))]=dy_1*(Qij[5*(Nx*(Ny*k  +jp1)+i  )+ll]-Qij[5*(Nx*(Ny*k  +j  )+i  )+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[1+3*(ll+5*(6))]=dy_1*(Qij[5*(Nx*(Ny*kp1+jp1)+i  )+ll]-Qij[5*(Nx*(Ny*kp1+jm1)+i  )+ll]);
                                                                                          
        for(int ll=0; ll<=4;ll++) dQ[2+3*(ll+5*(0))]=dz_1*(Qij[5*(Nx*(Ny*k  +j  )+i  )+ll]-Qij[5*(Nx*(Ny*km1+j  )+i  )+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[2+3*(ll+5*(1))]=dz_1*(Qij[5*(Nx*(Ny*kp1+jm1)+i  )+ll]-Qij[5*(Nx*(Ny*km1+jm1)+i  )+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[2+3*(ll+5*(2))]=dz_1*(Qij[5*(Nx*(Ny*kp1+j  )+im1)+ll]-Qij[5*(Nx*(Ny*km1+j  )+im1)+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[2+3*(ll+5*(3))]=dz_1*(Qij[5*(Nx*(Ny*kp1+j  )+i  )+ll]-Qij[5*(Nx*(Ny*km1+j  )+i  )+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[2+3*(ll+5*(4))]=dz_1*(Qij[5*(Nx*(Ny*kp1+j  )+ip1)+ll]-Qij[5*(Nx*(Ny*km1+j  )+ip1)+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[2+3*(ll+5*(5))]=dz_1*(Qij[5*(Nx*(Ny*kp1+jp1)+i  )+ll]-Qij[5*(Nx*(Ny*km1+jp1)+i  )+ll]);
	      for(int ll=0; ll<=4;ll++) dQ[2+3*(ll+5*(6))]=dz_1*(Qij[5*(Nx*(Ny*kp1+j  )+i  )+ll]-Qij[5*(Nx*(Ny*k  +j  )+i  )+ll]);
	    
	      k_i[5*(Nx*(Ny*k+j)+i)+0]=-bulk_force.force_00(QN,dQ); 
	      k_i[5*(Nx*(Ny*k+j)+i)+1]=-bulk_force.force_01(QN,dQ); 
	      k_i[5*(Nx*(Ny*k+j)+i)+2]=-bulk_force.force_02(QN,dQ); 
	      k_i[5*(Nx*(Ny*k+j)+i)+3]=-bulk_force.force_11(QN,dQ); 
	      k_i[5*(Nx*(Ny*k+j)+i)+4]=-bulk_force.force_12(QN,dQ);

	    }
    }
  }
}

void  Geometry_Bulk::Energy_calc(double * k_i,const double * Qij) const
{
#pragma omp for schedule(dynamic,1) collapse(2)
  for( int k= 0; k< Nz; k++)
  {
    for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
	    {	
	      int ip1= (i+1)%Nx;
        int im1= (i+Nx-1)%Nx;

        int jp1= (j+1)%Ny;
        int jm1= (j+Ny-1)%Ny;

        int kp1= (k+1)%Nz;
        int km1= (k+Nz-1)%Nz;

        double dQ[15];
        double QN[5];
        double v[3];

        v[0]=0;
        v[1]=0;
        v[2]=0;

	      for(int ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];

	      //Calcule first derivatives of Qij:
	      for(int ll=0; ll<=4;ll++) dQ[ll]= 0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
	      for(int ll=0; ll<=4;ll++) dQ[5+ll]= 0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
	      for(int ll=0; ll<=4;ll++) dQ[10+ll]= 0.5*(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;

	      k_i[5*(Nx*(Ny*k+j)+i)]=bulk_energy.energy_calculation(QN,dQ,v); 
	    }
    }
  }
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






