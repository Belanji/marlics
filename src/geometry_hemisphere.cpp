#include "driver.h"
#include "geometry.h"
#include "geometry_hemisphere.h"
#include "boundary.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <ctime>
#include <cstring>

#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)


Geometry_Hemisphere::Geometry_Hemisphere(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param),
                                                                                   HNx ((Nx-1.)/2),
                                                                                   HNy ((Ny-1.)/2.),
                                                                                   HNz (0.0),
                                                                                   dx(sim_param->dx),
                                                                                   dy(sim_param->dy),
                                                                                   dz(sim_param->dz)
{
  geometry_pointer=&(*this);
  geometry_name="Hemisphere";
  number_of_boundaries=2;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
  boundary_needed_to_be_defined="0 and 1";
  
  switch( sim_param->radius_flag[0] )
    {
     case parameter_status::unset:
         
       std::cout <<"Sphere internal radius not seted.\n"
                 <<"Using standard value.\n\n";

       R_in=( (Nx-1.)/2. - 0.4)*dx;
     break;

     case parameter_status::set:

       R_in=sim_param->R_in;
       break;
    }

  
    switch (sim_param->radius_flag[1])
     {
       case parameter_status::unset:
   
         std::cout <<"Sphere external radius not seted.\n"
                   <<"Using standard value.\n\n";
         R_ex=( (Nx-1)/2.+0.95)*dx;

         break;
   
       case parameter_status::set:
   
         R_ex=sim_param->R_ex;
         break;
     }

   std::cout <<"R_in=" << R_in <<std::endl
             <<"R_ex=" << R_ex <<std::endl;

   point_type=fill_point_type( );
}

Geometry_Hemisphere::~Geometry_Hemisphere()
{
  
  
}

void  Geometry_Hemisphere::fill_ki(double * k_i,
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
        double v[3];
      
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
    
          //Calculate second derivatives of Qij:
          for(ll=0; ll<=4;ll++) ddQ[ll]= (Qij[5*((k*Ny+j)*Nx+ip1)+ll]-2.0*QN[ll]+Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1*dx_1;
          for(ll=0; ll<=4;ll++) ddQ[ll+5]= 0.25*(Qij[5*((k*Ny+jp1)*Nx+ip1)+ll]+Qij[5*((k*Ny+jm1)*Nx+im1)+ll]-Qij[5*((k*Ny+jp1)*Nx+im1)+ll]-Qij[5*((k*Ny+jm1)*Nx+ip1)+ll])*dx_1*dy_1;
          for(ll=0; ll<=4;ll++) ddQ[10+ll]= 0.25*(Qij[5*((kp1*Ny+j)*Nx+ip1)+ll]+Qij[5*((km1*Ny+j)*Nx+im1)+ll]-Qij[5*((kp1*Ny+j)*Nx+im1)+ll]-Qij[5*((km1*Ny+j)*Nx+ip1)+ll])*dx_1*dz_1;
          for(ll=0; ll<=4;ll++) ddQ[15+ll]= (Qij[5*((k*Ny+jp1)*Nx+i)+ll]-2.0*QN[ll]+Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1*dy_1;
          for(ll=0; ll<=4;ll++) ddQ[20+ll]= 0.25*(Qij[5*((kp1*Ny+jp1)*Nx+i)+ll]+Qij[5*((km1*Ny+jm1)*Nx+i)+ll]-Qij[5*((kp1*Ny+jm1)*Nx+i)+ll]-Qij[5*((km1*Ny+jp1)*Nx+i)+ll])*dz_1*dy_1;
          for(ll=0; ll<=4;ll++) ddQ[25+ll]= (Qij[5*((kp1*Ny+j)*Nx+i)+ll]-2.0*QN[ll]+Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1*dz_1;
    
          k_i[5*(Nx*(Ny*k+j)+i)+0]=bulk_energy.functional_derivative_00(QN,dQ,ddQ,v); 
          k_i[5*(Nx*(Ny*k+j)+i)+1]=bulk_energy.functional_derivative_01(QN,dQ,ddQ,v); 
          k_i[5*(Nx*(Ny*k+j)+i)+2]=bulk_energy.functional_derivative_02(QN,dQ,ddQ,v); 
          k_i[5*(Nx*(Ny*k+j)+i)+3]=bulk_energy.functional_derivative_11(QN,dQ,ddQ,v); 
          k_i[5*(Nx*(Ny*k+j)+i)+4]=bulk_energy.functional_derivative_12(QN,dQ,ddQ,v);
          
        }  
	      else if( point_type[(k*Ny+j)*Nx+i] == 2 )
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
        //check_surface_limits(i,j,k);
        } 
        else if( point_type[(k*Ny+j)*Nx+i] == 3 )
        {
  
          //check_surface_limits( i,  j,  k);
          double delta_x,delta_y,delta_z;
          double rr;
          double v[3];
          
          delta_x=(i-HNx)*dx;
          delta_y=(j-HNy)*dy;
          delta_z=(k-HNz)*dz;
          rr=sqrt( delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
                            
          v[0]=delta_x/rr;
          v[1]=delta_y/rr;
          v[2]=delta_z/rr;
          //Boundary condtions:
    
          ip1= i>=HNx ? i   : i+1 ;
          jp1= j>=HNy ? j   : j+1 ;
          kp1= k>=HNz ? k   : k+1 ;
      
          im1= i>=HNx ? i-1 : i ;
          jm1= j>=HNy ? j-1 : j ;
          km1= k>=HNz ? k-1 : k ;
      
    
          for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
            
    
          //Calcule first derivatives of Qij:
          for(ll=0; ll<=4;ll++) dQ[ll]=   (Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
          for(ll=0; ll<=4;ll++) dQ[5+ll]= (Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
          for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
    
          k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions[1]->functional_derivative_00(QN,dQ,ddQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+1]= bc_conditions[1]->functional_derivative_01(QN,dQ,ddQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+2]= bc_conditions[1]->functional_derivative_02(QN,dQ,ddQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+3]= bc_conditions[1]->functional_derivative_11(QN,dQ,ddQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+4]= bc_conditions[1]->functional_derivative_12(QN,dQ,ddQ,v);
    
    
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
    }
  }
  
}

void  Geometry_Hemisphere::compute_forces(double * k_i, const double * Qij)  const 
{
  #pragma omp for simd schedule(dynamic,1) collapse(2)  
  for( int k= 0; k< Nz; k++)
  {
    for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
	    {	
        int ip1,jp1,kp1, im1, jm1, km1, ll;
      double dQ[105];
      double QN[7*5];
      double v[3];
              
	      if(point_type[(k*Ny+j)*Nx+i] == 1)
        {
          //check_bulk_limits( i,  j,  k);                    
          ip1= (i+1);
          jp1= (j+1);
          kp1= (k+1);
          im1= (i-1);
          jm1= (j-1);
          km1= (k-1);
          if(i*j*k==0||i==Nx-1 ||j==Ny-1 ||k==Nz-1) printf("Point type violation at (%d,%d,%d)\n",i,j,k);
  
 
  
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
        else if( point_type[(k*Ny+j)*Nx+i] == 2 )
        {
          
          int ip1= (i+1)%Nx;
          int jp1= (j+1)%Ny;
          int kp1= (k+1);
          int im1= (i+Nx-1)%Nx;
          int jm1= (j+Ny-1)%Ny;
          int km1= k;
          
          v[0]=0.0;
          v[1]=0.0;
          v[2]=-1.0;
          
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
        
        
          //check_surface_limits(i,j,k);
        }
        else if( point_type[(k*Ny+j)*Nx+i] == 3 )
        {
  
          //check_surface_limits( i,  j,  k);
          double delta_x,delta_y,delta_z;
          double rr;
          double v[3];
          
          delta_x=(i-HNx)*dx;
          delta_y=(j-HNy)*dy;
          delta_z=(k-HNz)*dz;
          rr=sqrt( delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
                            
          v[0]=delta_x/rr;
          v[1]=delta_y/rr;
          v[2]=delta_z/rr;
          //Boundary condtions:
    
          ip1= i>=HNx ? i   : i+1 ;
          jp1= j>=HNy ? j   : j+1 ;
          kp1= k>=HNz ? k   : k+1 ;
      
          im1= i>=HNx ? i-1 : i ;
          jm1= j>=HNy ? j-1 : j ;
          km1= k>=HNz ? k-1 : k ;
          
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
	    
          k_i[5*(Nx*(Ny*k+j)+i)+0]= -bc_conditions[1]->force_00(QN,dQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+1]= -bc_conditions[1]->force_01(QN,dQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+2]= -bc_conditions[1]->force_02(QN,dQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+3]= -bc_conditions[1]->force_11(QN,dQ,v);
          k_i[5*(Nx*(Ny*k+j)+i)+4]= -bc_conditions[1]->force_12(QN,dQ,v);
        
        
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
    }
  }
}

void  Geometry_Hemisphere::Energy_calc(double * k_i,
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
        double v[3];
      
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
    
          k_i[5*(Nx*(Ny*k+j)+i)]=bulk_energy.energy_calculation(QN,dQ,v); 
        }  
	      else if( point_type[(k*Ny+j)*Nx+i] == 2 )
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
            
          k_i[5*(Nx*(Ny*k+j)+i)]= bc_conditions[0]->energy_calculation(QN,dQ,v);
        //check_surface_limits(i,j,k);
        } 
        else if( point_type[(k*Ny+j)*Nx+i] == 3 )
        {
  
          //check_surface_limits( i,  j,  k);
          double delta_x,delta_y,delta_z;
          double rr;
          double v[3];
          
          delta_x=(i-HNx)*dx;
          delta_y=(j-HNy)*dy;
          delta_z=(k-HNz)*dz;
          rr=sqrt( delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
                            
          v[0]=delta_x/rr;
          v[1]=delta_y/rr;
          v[2]=delta_z/rr;
          //Boundary condtions:
    
          ip1= i>=HNx ? i   : i+1 ;
          jp1= j>=HNy ? j   : j+1 ;
          kp1= k>=HNz ? k   : k+1 ;
      
          im1= i>=HNx ? i-1 : i ;
          jm1= j>=HNy ? j-1 : j ;
          km1= k>=HNz ? k-1 : k ;
    
          for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
          //Calcule first derivatives of Qij:
          for(ll=0; ll<=4;ll++) dQ[ll]=   (Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
          for(ll=0; ll<=4;ll++) dQ[5+ll]= (Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
          for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
    
          k_i[5*(Nx*(Ny*k+j)+i)]= bc_conditions[1]->energy_calculation(QN,dQ,v);
          //check_surface_limits(i,j,k);
        }
          else
        {
          k_i[5*(Nx*(Ny*k+j)+i)]=0.; 
        }

      }
    }
  }
  
}
   
int * Geometry_Hemisphere::fill_point_type( void )  const 
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
      {       point_kind[(k*Ny+j)*Nx+i]=2;//Setting surface points
            
        delta_x=(i-HNx)*dx;
        delta_y=(j-HNy)*dy;
        delta_z=(k-HNz)*dz;
        rr=sqrt( delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
        
    
        if(k ==0 && rr <=R_ex)
        {
          point_kind[(k*Ny+j)*Nx+i]=2;//Setting bottom points
        }
        else if(rr < R_in)
        {                                 
          point_kind[(k*Ny+j)*Nx+i]=1;//Setting bulk points
        }
        else if(rr >=R_in && rr <=R_ex)    
        {
          point_kind[(k*Ny+j)*Nx+i]=3;//Setting surface points
        }
        else
        {
          point_kind[(k*Ny+j)*Nx+i]=0;//Setting outside points
        }
      }  
    }         
  }

  return point_kind;
}




  

