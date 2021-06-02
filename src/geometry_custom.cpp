#include "driver.h"
#include "geometry.h"
#include "geometry_custom.h"
#include "boundary.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <ctime>
#include <cstring>
#include <algorithm>   

#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)


Geometry_Custom::Geometry_Custom(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param),                                                                                                                HNx ((Nx-1.)/2),
                                                                                   HNy ((Ny-1.)/2.),
                                                                                   HNz ((Nz-1.)/2.),
                                                                                   dx(sim_param->dx),
                                                                                   dy(sim_param->dy),
                                                                                   dz(sim_param->dz)
{

  geometry_pointer=&(*this);
  geometry_name="Custom";
  v=(double**)calloc(sim_param->Nx*sim_param->Ny*sim_param->Nz,sizeof(double*));
  for (int i=0;i< (Nx*Ny*Nz); i++)v[i]=(double*)calloc(3,sizeof(double));
  point_type=fill_pt_and_normals(v, sim_param->bound_fname,number_of_boundaries);
  test_derivatives(sim_param->integrator_type);
  
  boundary_needed_to_be_defined=std::to_string(number_of_boundaries);
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
};

Geometry_Custom::~Geometry_Custom()
{
  
}

void  Geometry_Custom::fill_ki(double * k_i,
			       const double * Qij)  const 
{

  #pragma omp for simd schedule(dynamic,1) collapse(2)  
  for( int k= 0; k< Nz; k++)
  {
    for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
      {
        
        int ip1,jp1,kp1, im1, jm1, km1, ll;
        double dQ[15];
        double ddQ[30];
        double QN[5];
        
        if(point_type[(k*Ny+j)*Nx+i] == 1)
        {

          //check_bulk_limits( i,  j,  k);  
          ip1= (i+1)%Nx;
          jp1= (j+1)%Ny;
          kp1= (k+1)%Nz;
          im1= (i-1+Nx)%Nx;
          jm1= (j-1+Ny)%Ny;
          km1= (k-1+Nz)%Nz;
          
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
          for(ll=0; ll<=4;ll++) ddQ[20+ll]= 0.25*(Qij[5*((kp1*Ny+jp1)*Nx+i)+ll]+Qij[5*((km1*Ny+jm1)*Nx+i)+ll]-Qij[5*((kp1*Ny+jm1)*Nx+i)+ll]-Qij[5*((km1*Ny+jp1)*Nx+i)+ll])*dz_1*dy_1;
          for(ll=0; ll<=4;ll++) ddQ[25+ll]= (Qij[5*((kp1*Ny+j)*Nx+i)+ll]-2.0*QN[ll]+Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1*dz_1;
    
          k_i[5*(Nx*(Ny*k+j)+i)+0]=bulk_energy.functional_derivative_00(QN,dQ,ddQ,v[(k*Ny+j)*Nx+i]); 
          k_i[5*(Nx*(Ny*k+j)+i)+1]=bulk_energy.functional_derivative_01(QN,dQ,ddQ,v[(k*Ny+j)*Nx+i]); 
          k_i[5*(Nx*(Ny*k+j)+i)+2]=bulk_energy.functional_derivative_02(QN,dQ,ddQ,v[(k*Ny+j)*Nx+i]); 
          k_i[5*(Nx*(Ny*k+j)+i)+3]=bulk_energy.functional_derivative_11(QN,dQ,ddQ,v[(k*Ny+j)*Nx+i]); 
          k_i[5*(Nx*(Ny*k+j)+i)+4]=bulk_energy.functional_derivative_12(QN,dQ,ddQ,v[(k*Ny+j)*Nx+i]);
        }          
        else if (point_type[(k*Ny+j)*Nx+i] == 0)
        {
          k_i[5*(Nx*(Ny*k+j)+i)+0]=0.; 
          k_i[5*(Nx*(Ny*k+j)+i)+1]=0.; 
          k_i[5*(Nx*(Ny*k+j)+i)+2]=0.; 
          k_i[5*(Nx*(Ny*k+j)+i)+3]=0.; 
          k_i[5*(Nx*(Ny*k+j)+i)+4]=0.;
        }
        else 
        {
          //check_surface_limits( i,  j,  k);
          double delta_x,delta_y,delta_z;
          double rr;
          int bound=point_type[(k*Ny+j)*Nx+i]-2;
          
          //Boundary condtions:
          
          ip1= (v[(k*Ny+j)*Nx+i][0] >=0) ? i : (i+1)%Nx ;
          jp1= (v[(k*Ny+j)*Nx+i][1] >=0) ? j : (j+1)%Ny ;
          kp1= (v[(k*Ny+j)*Nx+i][2] >=0) ? k : (k+1)%Nz ;
          
          im1= (v[(k*Ny+j)*Nx+i][0] >=0) ? (i-1+Nx)%Nx : i ;
          jm1= (v[(k*Ny+j)*Nx+i][1] >=0) ? (j-1+Ny)%Ny : j ;
          km1= (v[(k*Ny+j)*Nx+i][2] >=0) ? (k-1+Nz)%Nz : k ;
      
    
          for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
          //Calcule first derivatives of Qij:
          for(ll=0; ll<=4;ll++) dQ[ll]=   (Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
          for(ll=0; ll<=4;ll++) dQ[5+ll]= (Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
          for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
    
          k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions[bound]->functional_derivative_00(QN,dQ,ddQ,v[(k*Ny+j)*Nx+i]);
          k_i[5*(Nx*(Ny*k+j)+i)+1]= bc_conditions[bound]->functional_derivative_01(QN,dQ,ddQ,v[(k*Ny+j)*Nx+i]);
          k_i[5*(Nx*(Ny*k+j)+i)+2]= bc_conditions[bound]->functional_derivative_02(QN,dQ,ddQ,v[(k*Ny+j)*Nx+i]);
          k_i[5*(Nx*(Ny*k+j)+i)+3]= bc_conditions[bound]->functional_derivative_11(QN,dQ,ddQ,v[(k*Ny+j)*Nx+i]);
          k_i[5*(Nx*(Ny*k+j)+i)+4]= bc_conditions[bound]->functional_derivative_12(QN,dQ,ddQ,v[(k*Ny+j)*Nx+i]);
        }
	    }
    }
  }
}

void  Geometry_Custom::compute_forces(double * k_i, const double * Qij)  const 
{
  #pragma omp for simd schedule(dynamic,1) collapse(2)  
  for( int k= 0; k< Nz; k++)
  {
    for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
	    {	
        int ip1,jp1,kp1, im1, jm1, km1, ll;
        double dQ[15];
        double ddQ[30];
        double QN[27*5];
              
	      if(point_type[(k*Ny+j)*Nx+i] == 1)
        {
          //check_bulk_limits( i,  j,  k);
          ip1= (i+1)%Nx;
          jp1= (j+1)%Ny;
          kp1= (k+1)%Nz;
          im1= (i-1+Nx)%Nx;
          jm1= (j-1+Ny)%Ny;
          km1= (k-1+Nz)%Nz;
          
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
	      else if( point_type[(k*Ny+j)*Nx+i] == 0 )
        {
          k_i[5*(Nx*(Ny*k+j)+i)+0]=0.; 
          k_i[5*(Nx*(Ny*k+j)+i)+1]=0.; 
          k_i[5*(Nx*(Ny*k+j)+i)+2]=0.; 
          k_i[5*(Nx*(Ny*k+j)+i)+3]=0.; 
          k_i[5*(Nx*(Ny*k+j)+i)+4]=0.;
        }
        else 
        {
          
          int bound=point_type[(k*Ny+j)*Nx+i]-2;
          
          //Boundary condtions:
          
          ip1= (v[(k*Ny+j)*Nx+i][0] >=0) ? i : (i+1)%Nx ;
          jp1= (v[(k*Ny+j)*Nx+i][1] >=0) ? j : (j+1)%Ny ;
          kp1= (v[(k*Ny+j)*Nx+i][2] >=0) ? k : (k+1)%Nz ;
          
          im1= (v[(k*Ny+j)*Nx+i][0] >=0) ? (i-1+Nx)%Nx : i ;
          jm1= (v[(k*Ny+j)*Nx+i][1] >=0) ? (j-1+Ny)%Ny : j ;
          km1= (v[(k*Ny+j)*Nx+i][2] >=0) ? (k-1+Nz)%Nz : k ;
	    
          ip1= (i+1)%Nx;
          im1= (i+Nx-1)%Nx;
          jp1= (j+1)%Ny;
          jm1= (j+Ny-1)%Ny;
          
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
	    
          k_i[5*(Nx*(Ny*k+j)+i)+0]= -bc_conditions[bound]->force_00(QN,dQ,v[Nx*(Ny*k+j)+i]);
          k_i[5*(Nx*(Ny*k+j)+i)+1]= -bc_conditions[bound]->force_01(QN,dQ,v[Nx*(Ny*k+j)+i]);
          k_i[5*(Nx*(Ny*k+j)+i)+2]= -bc_conditions[bound]->force_02(QN,dQ,v[Nx*(Ny*k+j)+i]);
          k_i[5*(Nx*(Ny*k+j)+i)+3]= -bc_conditions[bound]->force_11(QN,dQ,v[Nx*(Ny*k+j)+i]);
          k_i[5*(Nx*(Ny*k+j)+i)+4]= -bc_conditions[bound]->force_12(QN,dQ,v[Nx*(Ny*k+j)+i]);
       
        }
	    }
    }
  }
}

void  Geometry_Custom::test_derivatives( const char integrator_type[])  const 
{
  int counter=0;

//  #pragma omp master
  for( int k= 0; k< Nz; k++)
  {
    for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
      {
        
        int ip,jp,kp, im, jm, km, ll;
        double dQ[15];
        double ddQ[30];
        double QN[5];
        
        if(point_type[(k*Ny+j)*Nx+i] == 1)
        {

          //check_bulk_limits( i,  j,  k);  
          ip= (i+1)%Nx;
          jp= (j+1)%Ny;
          kp= (k+1)%Nz;
          im= (i-1+Nx)%Nx;
          jm= (j-1+Ny)%Ny;
          km= (k-1+Nz)%Nz;
          
          if(point_type[(k*Ny+j)*Nx+ip]*point_type[(k*Ny+j)*Nx+im]*
             point_type[(k*Ny+jp)*Nx+i]*point_type[(k*Ny+jm)*Nx+i]*
             point_type[(kp*Ny+j)*Nx+i]*point_type[(km*Ny+j)*Nx+i]*
             point_type[(k*Ny+jp)*Nx+ip]*point_type[(k*Ny+jp)*Nx+im]*
             point_type[(k*Ny+jm)*Nx+ip]*point_type[(k*Ny+jm)*Nx+im]*
             point_type[(kp*Ny+j)*Nx+ip]*point_type[(kp*Ny+j)*Nx+im]*
             point_type[(km*Ny+j)*Nx+ip]*point_type[(km*Ny+j)*Nx+im]*
             point_type[(kp*Ny+jp)*Nx+i]*point_type[(kp*Ny+jm)*Nx+i]*
             point_type[(km*Ny+jp)*Nx+i]*point_type[(km*Ny+jm)*Nx+i]==0)
             {
               printf("Warning!!! The bulk cell (%d,%d,%d) is trying to interact with a empty cell (type 0)!\n",i,j,k);
               counter++;
             }
             
        }          
        else if (point_type[(k*Ny+j)*Nx+i] == 0)
        {
          
        }
        else 
        {
          //check_surface_limits( i,  j,  k);
          double delta_x,delta_y,delta_z;
          double rr;
          int bound=point_type[(k*Ny+j)*Nx+i]-2;
          
          //Boundary condtions:
          
          ip= (v[(k*Ny+j)*Nx+i][0] >=0) ? i : (i+1)%Nx ;
          jp= (v[(k*Ny+j)*Nx+i][1] >=0) ? j : (j+1)%Ny ;
          kp= (v[(k*Ny+j)*Nx+i][2] >=0) ? k : (k+1)%Nz ;
          
          im= (v[(k*Ny+j)*Nx+i][0] >=0) ? (i-1+Nx)%Nx : i ;
          jm= (v[(k*Ny+j)*Nx+i][1] >=0) ? (j-1+Ny)%Ny : j ;
          km= (v[(k*Ny+j)*Nx+i][2] >=0) ? (k-1+Nz)%Nz : k ;
          if(strcasecmp(integrator_type,"fire"))
          {
            if(point_type[(k*Ny+j)*Nx+ip]*point_type[(k*Ny+j)*Nx+im]*
               point_type[(k*Ny+jp)*Nx+i]*point_type[(k*Ny+jm)*Nx+i]*
               point_type[(kp*Ny+j)*Nx+i]*point_type[(km*Ny+j)*Nx+i]*
               point_type[(k*Ny+jp)*Nx+ip]*point_type[(k*Ny+jp)*Nx+im]*
               point_type[(k*Ny+jm)*Nx+ip]*point_type[(k*Ny+jm)*Nx+im]*
               point_type[(kp*Ny+j)*Nx+ip]*point_type[(kp*Ny+j)*Nx+im]*
               point_type[(km*Ny+j)*Nx+ip]*point_type[(km*Ny+j)*Nx+im]*
               point_type[(kp*Ny+jp)*Nx+i]*point_type[(kp*Ny+jm)*Nx+i]*
               point_type[(km*Ny+jp)*Nx+i]*point_type[(km*Ny+jm)*Nx+i]==0)
              {
                printf("Warning!!! The surface cell (%d,%d,%d) is trying to interact with a empty cell (type 0)!\n",i,j,k);
                counter++;
              }
          }else{
            if(point_type[(k*Ny+j)*Nx+ip]*point_type[(k*Ny+j)*Nx+im]*
               point_type[(k*Ny+jp)*Nx+i]*point_type[(k*Ny+jm)*Nx+i]*
               point_type[(kp*Ny+j)*Nx+i]*point_type[(km*Ny+j)*Nx+i]==0)
              {
                printf("Warning!!! The surface cell (%d,%d,%d) is trying to interact with a empty cell (type 0)!\n",i,j,k);
                counter++;
              }
          }
        }
	    }
    }
  }
  #pragma omp master
  {
    if(counter!=0)
    {
    printf("%d cells are trying to interact with an empty cell!!\n Please, check your boundary file.\n",counter);
    exit(10);
    }
  }
}

void  Geometry_Custom::Energy_calc(double * k_i, const double * Qij)  const 
{
  #pragma omp for simd schedule(dynamic,1) collapse(2)  
  for( int k= 0; k< Nz; k++)
  {
    for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
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
          for(ll=0; ll<=4;ll++) dQ[5+ll]= 0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
          for(ll=0; ll<=4;ll++) dQ[10+ll]= 0.5*(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
    
          k_i[2*(Nx*(Ny*k+j)+i)]=bulk_energy.energy_calculation(QN,dQ,v[(k*Ny+j)*Nx+i]); 
        }  
	      else if( point_type[(k*Ny+j)*Nx+i] == 0 )
        {
          k_i[5*(Nx*(Ny*k+j)+i)+0]=0.; 
          k_i[5*(Nx*(Ny*k+j)+i)+1]=0.; 
          k_i[5*(Nx*(Ny*k+j)+i)+2]=0.; 
          k_i[5*(Nx*(Ny*k+j)+i)+3]=0.; 
          k_i[5*(Nx*(Ny*k+j)+i)+4]=0.;
        }
        else 
        {
          //check_surface_limits( i,  j,  k);
          double delta_x,delta_y,delta_z;
          double rr;
          int bound=point_type[(k*Ny+j)*Nx+i]-1;
          
          //Boundary condtions:
          
          ip1= (v[(k*Ny+j)*Nx+i][0] >=0) ? i : (i+1)%Nx ;
          jp1= (v[(k*Ny+j)*Nx+i][1] >=0) ? j : (j+1)%Ny ;
          kp1= (v[(k*Ny+j)*Nx+i][2] >=0) ? k : (k+1)%Nz ;
          
          im1= (v[(k*Ny+j)*Nx+i][0] >=0) ? (i-1+Nx)%Nx : i ;
          jm1= (v[(k*Ny+j)*Nx+i][1] >=0) ? (j-1+Ny)%Ny : j ;
          km1= (v[(k*Ny+j)*Nx+i][2] >=0) ? (k-1+Nz)%Nz : k ;
    
          for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
    
          //Calcule first derivatives of Qij:
          for(ll=0; ll<=4;ll++) dQ[ll]=   (Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
          for(ll=0; ll<=4;ll++) dQ[5+ll]= (Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
          for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
    
          k_i[2*(Nx*(Ny*k+j)+i)+0]=bulk_energy.energy_calculation(QN,dQ,v[(k*Ny+j)*Nx+i]); 
          k_i[2*(Nx*(Ny*k+j)+i)+1]=bc_conditions[0]->energy_calculation(QN,dQ,v[(k*Ny+j)*Nx+i]);
//~           printf("%g %g %d %d %d\n",k_i[2*(Nx*(Ny*k+j)+i)+1],k_i[2*(Nx*(Ny*k+j)+i)],i,j,k);
          //check_surface_limits(i,j,k);
        }
	    }
    }
  }
}
       
int * Geometry_Custom::fill_pt_and_normals( double **v , const char bound_file_name[], int &number_of_boundaries)  const 
{
  int i, j, k, ijk, nn;
  double v_temp[3];
  int * point_kind, max_point_kind, line_number, read_tester=1;
  int NoV,pos_nx, pos_ny, pos_nz, pos_pt;
  FILE *bound_input = fopen(bound_file_name,"r");
  if(bound_input==0){perror(bound_file_name);exit(2);}
  printf("Reading boundaries from %s\n",bound_file_name);
  char line[500];
  std::string testline;
  fgets(line,500,bound_input);
  testline=line;
  NoV=1+std::count(testline.begin(),testline.end(),',');
  double var[NoV];
  char strvar[500];
  rewind(bound_input);
  for (i=0; i<NoV;i++) 
  {
    char names[500];
    fscanf(bound_input,"%[^,\n ],",names);
    if(strcasecmp(names,"nx")==0)pos_nx=i;
    if(strcasecmp(names,"ny")==0)pos_ny=i;
    if(strcasecmp(names,"nz")==0)pos_nz=i;
    if(strcasecmp(names,"pt")==0)pos_pt=i;
  }
  printf("Using the positions of nx %d ny %d nz %d and pt %d.\n ",pos_nx,pos_ny,pos_nz,pos_pt);
  if((point_kind= (int *)calloc(Nx*Ny*Nz, sizeof(int)))==NULL){ERROr}
  line_number=0;
  fgets(line,500,bound_input);
  for (i=0; i<Nx*Ny*Nz; i++){
    point_kind[i]=1;
  }  
  while (read_tester!=EOF){
    for (nn=0;nn<NoV;nn++) 
    {
      read_tester=fscanf(bound_input,"%[^,\n ],",strvar);
      if(read_tester==1)var[nn]=atof(strvar);
      else continue;
    }
    int ii=(int)var[0];
    int jj=(int)var[1];
    int kk=(int)var[2];
    ijk=(Ny*kk+jj)*Nx+ii;
    line_number++;
    if (ii>=Nx|| jj>=Ny||kk>=Nz){fprintf(stderr,"Point out of the simulation box at line %d\n",line_number);exit(10);}
    if ((int)var[pos_pt]<0){fprintf(stderr,"Negative Point Type founde at line %d\n\
Please, use 0 as empty space, 1 as bulk and 2+ as surface point type.\n",line_number);exit(10);}
    fgets(line,500,bound_input);
             v[(Ny*kk+jj)*Nx+ii][0]=var[pos_nx];
             v[(Ny*kk+jj)*Nx+ii][1]=var[pos_ny];
             v[(Ny*kk+jj)*Nx+ii][2]=var[pos_nz];
    point_kind[(Ny*kk+jj)*Nx+ii]=(int)var[pos_pt];
    max_point_kind=MAX(max_point_kind,point_kind[(Ny*kk+jj)*Nx+ii]);
  }
  number_of_boundaries=max_point_kind-1;
  printf("%d boundary(boundaries) found in %s!!\n\n",number_of_boundaries, bound_file_name);
  fflush(stdout);
  return point_kind;
}

int * Geometry_Custom::fill_point_type( void )  const  { return 0;}
