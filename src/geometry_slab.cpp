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
  
  int k=0;
  #pragma omp for schedule(dynamic,1) collapse(2)
  for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
        {
          double dQ[15];
          double ddQ[30];
          double QN[5];
          double v[3];
          
          int ip1= (i+1)%Nx;
          int jp1= (j+1)%Ny;
          int kp1= (k+1);
          int im1= (i+Nx-1)%Nx;
          int jm1= (j+Ny-1)%Ny;
          int km1= k;
          v[0]=0.0;
          v[1]=0.0;
          v[2]=-1.0;
      
          for(int ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
      
          //Calcule first derivatives of Qij:
          for(int ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
            
          for(int ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
            
          for(int ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
            
          k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions[0]->functional_derivative_00(QN,dQ,ddQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+1]= bc_conditions[0]->functional_derivative_01(QN,dQ,ddQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+2]= bc_conditions[0]->functional_derivative_02(QN,dQ,ddQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+3]= bc_conditions[0]->functional_derivative_11(QN,dQ,ddQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+4]= bc_conditions[0]->functional_derivative_12(QN,dQ,ddQ,v);
      
        }
    }
	  
#pragma omp for schedule(dynamic,1) collapse(2) private(k)
  for( k= 1; k< Nz-1; k++)
    {
      for( int j= 0; j< Ny; j++)
        {
          for( int i= 0; i< Nx; i++)
            {	
              double dQ[15];
              double ddQ[30];
              double QN[5];
              double v[3];
      
              //check_bulk_limits( i,  j,  k);  
      
              int ip1= (i+1)%Nx;
              int jp1= (j+1)%Ny;
              int kp1= (k+1);
              int im1= (i+Nx-1)%Nx;
              int jm1= (j+Ny-1)%Ny;
              int km1= (k-1);
      
              for(int ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
      
              //Calcule first derivatives of Qij:
              for(int ll=0; ll<=4;ll++) dQ[ll]= 0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
            
              for(int ll=0; ll<=4;ll++) dQ[5+ll]= 0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
            
              for(int ll=0; ll<=4;ll++) dQ[10+ll]= 0.5*(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
      
              //Calculate second derivatives of Qij:
              for(int ll=0; ll<=4;ll++) ddQ[ll]= (Qij[5*((k*Ny+j)*Nx+ip1)+ll]-2.0*QN[ll]+Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1*dx_1;
                  
              for(int ll=0; ll<=4;ll++) ddQ[ll+5]= 0.25*(Qij[5*((k*Ny+jp1)*Nx+ip1)+ll]+Qij[5*((k*Ny+jm1)*Nx+im1)+ll]-Qij[5*((k*Ny+jp1)*Nx+im1)+ll]-Qij[5*((k*Ny+jm1)*Nx+ip1)+ll])*dx_1*dy_1;
            
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

  k=Nz-1;
#pragma omp for schedule(dynamic,1) collapse(2)  
  for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
        {
          double dQ[15];
          double ddQ[30];
          double QN[5];
          double v[3];
      
          int ip1= (i+1)%Nx;
          int jp1= (j+1)%Ny;
          int kp1= k;
          int im1= (i+Nx-1)%Nx;
          int jm1= (j+Ny-1)%Ny;
          int km1= (k-1);
          v[0]=0.0;
          v[1]=0.0;
          v[2]=1.0;
            
          for(int ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
      
          //Calcule first derivatives of Qij:
          for(int ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
            
          for(int ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
            
          for(int ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((k*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
            
          k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions[1]->functional_derivative_00(QN,dQ,ddQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+1]= bc_conditions[1]->functional_derivative_01(QN,dQ,ddQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+2]= bc_conditions[1]->functional_derivative_02(QN,dQ,ddQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+3]= bc_conditions[1]->functional_derivative_11(QN,dQ,ddQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+4]= bc_conditions[1]->functional_derivative_12(QN,dQ,ddQ,v);
		 
        }
    }
	
}

    
void  slab::compute_forces(double * k_i,const double * Qij) const
{
  int k=0;
  #pragma omp for schedule(dynamic,1) collapse(2)
  for( int j= 0; j< Ny; j++)
  {
    for( int i= 0; i< Nx; i++)
    {
	    int ip1= (i+1)%Nx;
      int im1= (i+Nx-1)%Nx;
      int jp1= (j+1)%Ny;
      int jm1= (j+Ny-1)%Ny;
      int kp1= (k+1)%Nz;
      int km1= k;//(k+Nz-1)%Nz;
    
      double dQ[105];
      double QN[7*5];
      double v[3];
    
      v[0]=0;
      v[1]=0;
      v[2]=-1;
      
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
	    
      k_i[5*(Nx*(Ny*k+j)+i)+0]= -bc_conditions[0]->force_00(QN,dQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+1]= -bc_conditions[0]->force_01(QN,dQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+2]= -bc_conditions[0]->force_02(QN,dQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+3]= -bc_conditions[0]->force_11(QN,dQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+4]= -bc_conditions[0]->force_12(QN,dQ,v);
      
    }
  }
	  
  #pragma omp for schedule(dynamic,1) collapse(2) private(k)
  for( k= 1; k< Nz-1; k++)
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
        double v[3];

        v[0]=0;
        v[1]=0;
        v[2]=0;

      
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
  
  k=Nz-1;
  #pragma omp for schedule(dynamic,1) collapse(2)  
  for( int j= 0; j< Ny; j++)
  {
    for( int i= 0; i< Nx; i++)
    {   
      
	      int ip1= (i+1)%Nx;
        int im1= (i+Nx-1)%Nx;
        int jp1= (j+1)%Ny;
        int jm1= (j+Ny-1)%Ny;
        int kp1= k;
        int km1= (k+Nz-1)%Nz;

        double dQ[105];
        double QN[7*5];
        double v[3];

        v[0]=0;
        v[1]=0;
        v[2]=1;
      
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
        //for(int ll=0; ll<=4;ll++) dQ[2+3*(ll+5*(0))]==0;
        //for(int ll=0; ll<=4;ll++) dQ[2+3*(ll+5*(6))]==0;
      
	      k_i[5*(Nx*(Ny*k+j)+i)+0]= -bc_conditions[1]->force_00(QN,dQ,v);
	      k_i[5*(Nx*(Ny*k+j)+i)+1]= -bc_conditions[1]->force_01(QN,dQ,v);
	      k_i[5*(Nx*(Ny*k+j)+i)+2]= -bc_conditions[1]->force_02(QN,dQ,v);
	      k_i[5*(Nx*(Ny*k+j)+i)+3]= -bc_conditions[1]->force_11(QN,dQ,v);
	      k_i[5*(Nx*(Ny*k+j)+i)+4]= -bc_conditions[1]->force_12(QN,dQ,v);
		 }
  }
}
      
void  slab::Energy_calc(double * k_i, const double * Qij)  const 
{
  int k=0;
  #pragma omp for schedule(dynamic,1) collapse(2)
  for( int j= 0; j< Ny; j++)
  {
    for( int i= 0; i< Nx; i++)
    {
      double dQ[15];
      double ddQ[30];
      double QN[5];
      double v[3];
      
      int ip1= (i+1)%Nx;
      int jp1= (j+1)%Ny;
      int kp1= (k+1);
      int im1= (i+Nx-1)%Nx;
      int jm1= (j+Ny-1)%Ny;
      int km1= k;
      v[0]=0.0;
      v[1]=0.0;
      v[2]=-1.0;

      for(int ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
      
      //Calcule first derivatives of Qij:
      for(int ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
      for(int ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
      for(int ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
        
      k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions[0]->energy_calculation(QN,dQ,v);
    }
  }
	  
  #pragma omp for schedule(dynamic,1) collapse(2) private(k)
  for( k= 1; k< Nz-1; k++)
  {
    for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
	    {	
        double dQ[15];
        double ddQ[30];
        double QN[5];
        double v[3];

	      //check_bulk_limits( i,  j,  k);  

	      int ip1= (i+1)%Nx;
	      int jp1= (j+1)%Ny;
	      int kp1= (k+1);
	      int im1= (i+Nx-1)%Nx;
	      int jm1= (j+Ny-1)%Ny;
	      int km1= (k-1);

	      for(int ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
        
        //Calcule first derivatives of Qij:
	      for(int ll=0; ll<=4;ll++) dQ[ll]= 0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
	      for(int ll=0; ll<=4;ll++) dQ[5+ll]= 0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
	      for(int ll=0; ll<=4;ll++) dQ[10+ll]= 0.5*(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;

	      k_i[5*(Nx*(Ny*k+j)+i)]=bulk_energy.energy_calculation(QN,dQ,v); 
	    }
    }
  }
  k=Nz-1;
  #pragma omp for schedule(dynamic,1) collapse(2)  
  for( int j= 0; j< Ny; j++)
  {
    for( int i= 0; i< Nx; i++)
    {
      double dQ[15];
      double ddQ[30];
      double QN[5];
      double v[3];
  
      int ip1= (i+1)%Nx;
      int jp1= (j+1)%Ny;
      int kp1= k;
      int im1= (i+Nx-1)%Nx;
      int jm1= (j+Ny-1)%Ny;
      int km1= (k-1);
      v[0]=0.0;
      v[1]=0.0;
      v[2]=1.0;
      
      for(int ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
  
      //Calcule first derivatives of Qij:
      for(int ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
      for(int ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
      for(int ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((k*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
        
      k_i[5*(Nx*(Ny*k+j)+i)]= bc_conditions[1]->energy_calculation(QN,dQ,v);
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






