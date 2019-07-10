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
//#include "boundary_fournier_galatola.h"
//#include "boundary_strong.h"
#include <gsl/gsl_randist.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)


slab::slab(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param)
{
  point_type=fill_point_type( );
  geometry_pointer=&(*this);

  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
  
};

void slab::boundary_init( struct Simulation_Parameters * sim_param)
{
  

//  printf("Boundary type:  %s\n", sim_param->upper_anc_type);
//
// 
//  //Upper boundary:
//
//
//  if( strcmp(sim_param->upper_anc_type,"strong") == 0 || strcmp(sim_param->upper_anc_type,"fixed") == 0 || strcmp(sim_param->upper_anc_type,"rapini_papoular") == 0 )
//    {
//
//      bc_conditions=new Strong_Boundary(geometry_pointer);
//      printf("\n");
//    }
//  else
//    {
//      printf("There is no anchoring type named %s.\nPlease check your input file 'upper anchoring field'.\n\n Aborting the program.'\n",sim_param->upper_anc_type);
//      exit(0);
//    }

 
  
};


void  slab::fill_ki(double * k_i,
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
  
  if(point_type[(k*Ny+j)*Nx+i] == 1)
    {

      //check_bulk_limits( i,  j,  k);	
                 
      ip1= i+1;
      jp1= j+1;
      kp1= k+1;
      im1= i-1;
      jm1= j-1;
      km1= k-1;


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
      ip1= (i+1)%(Nx);
      jp1= (j+1)%(Ny);
      kp1= (k+1);
      im1= i-1+((Nx-1-i)/(Nx-1))*Nx;
      jm1= j-1+((Ny-1-i)/(Ny-1))*Ny;
      km1= (k-1);
      v[0]=0.0;
      v[1]=0.0;
      v[2]=-1.0;
      
      for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
	

      //Calcule first derivatives of Qij:
      for(ll=0; ll<=4;ll++) dQ[ll]=   (Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;

      
      for(ll=0; ll<=4;ll++) dQ[5+ll]= (Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dx_1;

      
      for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dx_1;
      
  
      //k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions->surface_00(QN,dQ,ddQ,v);
      //k_i[5*(Nx*(Ny*k+j)+i)+1]= bc_conditions->surface_01(QN,dQ,ddQ,v);
      //k_i[5*(Nx*(Ny*k+j)+i)+2]= bc_conditions->surface_02(QN,dQ,ddQ,v);
      //k_i[5*(Nx*(Ny*k+j)+i)+3]= bc_conditions->surface_11(QN,dQ,ddQ,v);
      //k_i[5*(Nx*(Ny*k+j)+i)+4]= bc_conditions->surface_12(QN,dQ,ddQ,v); 


      
    }
  else if( point_type[(k*Ny+j)*Nx+i] == 3 )
    {

      
      ip1= (i+1)%(Nx);
      jp1= (j+1)%(Ny);
      kp1= (k+1);
      im1= i-1+((Nx-1-i)/(Nx-1))*Nx;
      jm1= j-1+((Ny-1-i)/(Ny-1))*Ny;
      km1= (k-1);
      v[0]=0.0;
      v[1]=0.0;
      v[2]=1.0;

      
      for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
	

      //Calcule first derivatives of Qij:
      for(ll=0; ll<=4;ll++) dQ[ll]=(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;

      
      for(ll=0; ll<=4;ll++) dQ[5+ll]= (Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dx_1;

      
      for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dx_1;
      
  
      //k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions->surface_00(QN,dQ,ddQ,v);
      //k_i[5*(Nx*(Ny*k+j)+i)+1]= bc_conditions->surface_01(QN,dQ,ddQ,v);
      //k_i[5*(Nx*(Ny*k+j)+i)+2]= bc_conditions->surface_02(QN,dQ,ddQ,v);
      //k_i[5*(Nx*(Ny*k+j)+i)+3]= bc_conditions->surface_11(QN,dQ,ddQ,v);
      //k_i[5*(Nx*(Ny*k+j)+i)+4]= bc_conditions->surface_12(QN,dQ,ddQ,v); 


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



void slab::ic(struct Simulation_Parameters * sim_param,double * Qij)
{


  switch(sim_param->ic_flag[0])
    {
    case parameter_status::set:

      printf("\nInitial Conditions: %s.\n",sim_param->initial_conditions);
      
      if(strcasecmp(sim_param->initial_conditions,"random") == 0 || strcmp(sim_param->initial_conditions,"Random") == 0 )
	{

	  random_ic( sim_param, Qij );
	  

	}
      if(strcasecmp(sim_param->initial_conditions,"homogeneous") == 0 )
	{

	  homogeneous_ic( sim_param, Qij );
	  

	}
      else if(strcasecmp(sim_param->initial_conditions,"read_from_file") == 0 )	     
	{

	  if( sim_param->ic_flag[1] == parameter_status::unset )
	    {
	      printf("Missing the \"initial_conditions_file\" in your in put file.\n Aborting the program.\n\n");
	      exit(0);
	    }
	  
	  read_from_file_ic( sim_param, Qij );

	}
      else
	{
      
	  printf("\n The program did not recognize the initial condition option \"%s\".\nPlease review your input file.\n\nAborting the program.\n",sim_param->initial_conditions);

	  exit(0);
	}
      break;
      
    case parameter_status::unset:

      printf("Parameter \"initial_conditions\" not set in your in put file.\nPlease, set the parameter for one of the available initial consitions in this geometry:\n");
      printf("random,read_from_file\n\n");
      printf("Aborting the program.\n\n");
      exit(0);
    }

}

