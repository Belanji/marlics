#include "geometry.h"
#include "driver.h"
#include "boundary.h"
#include "boundary_strong.h"
#include "boundary_rp.h"
#include "boundary_fournier_galatola.h"
#include <stdio.h> 
#include <gsl/gsl_randist.h>
#include <ctime>
#include <math.h>
#include <vector>
#include <string>
#include <cstring>
#include <iostream>

double Pi=3.14159265359;

GEOMETRY::GEOMETRY(const struct Simulation_Parameters * lc) :
    Lambda(1/lc->mu_1),
    Lambda_s(1/(lc->mu_1_s)),
    a(lc->a),
    sigma(lc->a*lc->T),
    bb(lc->B),
    cc(lc->C),
    L1(lc->L1),
    L2(lc->L2),
    L3(lc->L3),
    Lq(lc->Lq),
    Ls(lc->Ls),
    Lq_tilde(lc->Lq*2.0*lc->q0),
    S_eq(lc->S_eq),
    Nx(lc->Nx),
    Ny(lc->Ny),
    Nz(lc->Nz),
    p0(lc->p0),
    q0(lc->q0),
    dx_1(1/lc->dx),
    dy_1(1/lc->dy),
    dz_1(1/lc->dz),
    Nx_1(1./lc->Nx)
{

  printf("Nx=%i \n",Nx);
  printf("Ny=%i \n",Ny);
  printf("Nz=%i \n\n",Nz);
	
	
  printf("L1=%lf \n",L1);
  printf("L2=%lf \n",L2);
  printf("L3=%lf \n",L3);
  printf("Lq=%lf \n",Lq);
  printf("Ls=%lf \n\n",Ls);
  printf("p0=%lf \n",p0);
  printf("dx=%lf \n",lc->dx);
  printf("dy=%lf \n",lc->dy);
  printf("dz=%lf \n",lc->dz);
  printf("S_eq=%lf \n\n",S_eq);
	

	
  printf("a=%lf \n",a);
  printf("b=%lf \n",bb);
  printf("c=%lf \n",cc);
  printf("T=%lf \n\n",lc->T);

  printf("Lambda=%lf \n",1/Lambda);
  printf("Lambda_s=%lf \n\n",1/Lambda_s);
  

  RhsPtr=NULL;
  JacobianPtr=NULL;

  
};


#define QN00 QN[0] 
#define QN01 QN[1]
#define QN02 QN[2]
#define QN11 QN[3]
#define QN12 QN[4]

#define  Q_00_0 dQ[0]
#define  Q_01_0 dQ[1]
#define  Q_02_0 dQ[2]
#define  Q_11_0 dQ[3]
#define  Q_12_0 dQ[4]

#define  Q_00_1 dQ[5]
#define  Q_01_1 dQ[6]
#define  Q_02_1 dQ[7]
#define  Q_11_1 dQ[8]
#define  Q_12_1 dQ[9]

#define  Q_00_2 dQ[10]
#define  Q_01_2 dQ[11]
#define  Q_02_2 dQ[12]
#define  Q_11_2 dQ[13]
#define  Q_12_2 dQ[14]


#define  Q_00_00 ddQ[0]
#define  Q_01_00 ddQ[1]
#define  Q_02_00 ddQ[2]
#define  Q_11_00 ddQ[3]
#define  Q_12_00 ddQ[4]
		
#define  Q_00_01 ddQ[5]
#define  Q_01_01 ddQ[6]
#define  Q_02_01 ddQ[7]
#define  Q_11_01 ddQ[8]
#define  Q_12_01 ddQ[9]
		
#define  Q_00_02 ddQ[10]
#define  Q_01_02 ddQ[11]
#define  Q_02_02 ddQ[12]
#define  Q_11_02 ddQ[13]
#define  Q_12_02 ddQ[14]

#define  Q_00_11 ddQ[15]
#define  Q_01_11 ddQ[16]
#define  Q_02_11 ddQ[17]
#define  Q_11_11 ddQ[18]
#define  Q_12_11 ddQ[19]

#define  Q_00_12 ddQ[20]
#define  Q_01_12 ddQ[21]
#define  Q_02_12 ddQ[22]
#define  Q_11_12 ddQ[23]
#define  Q_12_12 ddQ[24]

#define  Q_00_22 ddQ[25]
#define  Q_01_22 ddQ[26]
#define  Q_02_22 ddQ[27]
#define  Q_11_22 ddQ[28]
#define  Q_12_22 ddQ[29]


double GEOMETRY::bulk_00(const double  QN[5],const double  dQ[15],const double  ddQ[30]) const
  { 
    return Lambda*((-3.*sigma*QN00 - 6.*cc*QN00*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12)) - bb*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) - 2*QN00*QN11 - 2.*((QN11*QN11) + (QN12*QN12))) + 6.*Lq_tilde*(Q_01_2 - Q_02_1) + (3.*L1 + 2*(L2 + Ls))*Q_00_00 + 3.*L1*Q_00_11 + 3.*L1*Q_00_22 + L2*Q_00_22 + Ls*Q_00_22 + L3*((Q_00_0*Q_00_0) + (Q_00_1*Q_00_1) - 2.*(Q_00_2*Q_00_2) - 2.*(Q_01_0*Q_01_0) + (Q_01_1*Q_01_1) + (Q_01_2*Q_01_2) - 2.*(Q_02_0*Q_02_0) + (Q_02_1*Q_02_1) + (Q_02_2*Q_02_2) + Q_00_0*(3.*(Q_01_1 + Q_02_2) - 2.*Q_11_0) - 2.*(Q_11_0*Q_11_0) + (Q_11_1*Q_11_1) + (Q_11_2*Q_11_2) - 2.*(Q_12_0*Q_12_0) + (Q_12_1*Q_12_1) + Q_00_2*(3.*Q_02_0 - 2.*Q_11_2 + 3.*Q_12_1) + (Q_12_2*Q_12_2) + Q_00_1*(3.*Q_01_0 + 4.*Q_11_1 + 3.*Q_12_2) + 3.*QN00*Q_00_00 + 6.*QN01*Q_00_01 + 6.*QN02*Q_00_02 + 3.*QN11*Q_00_11 + 6.*QN12*Q_00_12 - 3.*(QN00 + QN11)*Q_00_22) + L2*Q_01_01 + Ls*Q_01_01 + L2*Q_02_02 + Ls*Q_02_02 - L2*Q_11_11 - Ls*Q_11_11 + L2*Q_11_22 + Ls*Q_11_22 - 2.*(L2 + Ls)*Q_12_12)/3.);
  }

double GEOMETRY::bulk_01(const double  QN[5],const double  dQ[15],const double  ddQ[30]) const
  { 
    return Lambda*((-2*Q_00_2*(Lq_tilde + L3*Q_01_2) - L3*Q_00_0*(2*Q_00_1 - 2*Q_01_0 + Q_11_1) + 2*Lq_tilde*(Q_02_0 + Q_11_2 - Q_12_1) - 2*(sigma*QN01 + bb*QN01*(QN00 + QN11) + bb*QN02*QN12 + 2*cc*QN01*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12)) - L1*(Q_01_00 + Q_01_11 + Q_01_22)) + L3*(-2*Q_02_0*Q_02_1 + 2*Q_01_0*(Q_01_1 + Q_02_2) - Q_11_0*(Q_00_1 + 2*Q_11_1) + 2*Q_01_2*(Q_02_0 - Q_11_2 + Q_12_1) + 2*(-(Q_12_0*Q_12_1) + Q_01_1*(Q_11_1 + Q_12_2) + QN00*Q_01_00 + 2*QN01*Q_01_01 + 2*QN02*Q_01_02 + QN11*Q_01_11 + 2*QN12*Q_01_12 - (QN00 + QN11)*Q_01_22)) + (L2 + Ls)*(Q_00_01 + Q_01_00 + Q_01_11 + Q_02_12 + Q_11_01 + Q_12_02))/2.);
  }

double GEOMETRY::bulk_02(const double  QN[5],const double  dQ[15],const double  ddQ[30]) const
  { 
    return Lambda*((4*Lq_tilde*Q_00_1 + 2*Lq_tilde*(-Q_01_0 + Q_11_1 + Q_12_2) - 2*(sigma*QN02 - bb*QN02*QN11 + bb*QN01*QN12 + 2*cc*QN02*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12)) - L1*(Q_02_00 + Q_02_11 + Q_02_22)) + L3*(-(Q_00_2*(2*Q_02_2 + Q_11_0)) - Q_00_0*(2*Q_00_2 - 2*Q_02_0 + Q_11_2) + 2*(Q_01_1*Q_02_0 + Q_02_0*Q_02_2 - (Q_02_2 + Q_11_0)*Q_11_2 + Q_02_2*Q_12_1) + 2*(Q_01_0*(-Q_01_2 + Q_02_1) - Q_12_0*Q_12_2 + Q_02_1*(Q_11_1 + Q_12_2) + QN00*Q_02_00 + 2*QN01*Q_02_01 + 2*QN02*Q_02_02 + QN11*Q_02_11 + 2*QN12*Q_02_12 - (QN00 + QN11)*Q_02_22)) + (L2 + Ls)*(Q_01_12 + Q_02_00 + Q_02_22 - Q_11_02 + Q_12_01))/2.);
  }

double GEOMETRY::bulk_11(const double  QN[5],const double  dQ[15],const double  ddQ[30]) const
  { 
    return Lambda*((bb*(2*(QN00*QN00) - (QN01*QN01) + 2*(QN02*QN02) + 2*QN00*QN11 - (QN11*QN11) - (QN12*QN12)) - 3*QN11*(sigma + 2*cc*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12))) + 6*Lq_tilde*(-Q_01_2 + Q_12_0) - (L2 + Ls)*Q_00_00 + 3*L1*(Q_11_00 + Q_11_11 + Q_11_22) + L3*((Q_00_0*Q_00_0) - 2*(Q_00_1*Q_00_1) + (Q_00_2*Q_00_2) + (Q_01_0*Q_01_0) - 2*(Q_01_1*Q_01_1) + (Q_01_2*Q_01_2) + (Q_02_0*Q_02_0) - 2*(Q_02_1*Q_02_1) + (Q_02_2*Q_02_2) + (4*Q_00_0 + 3*(Q_01_1 + Q_02_2))*Q_11_0 + (Q_11_0*Q_11_0) - 2*Q_00_1*Q_11_1 + 3*Q_01_0*Q_11_1 + (Q_11_1*Q_11_1) - 2*Q_00_2*Q_11_2 + 3*Q_02_0*Q_11_2 - 2*(Q_11_2*Q_11_2) + (Q_12_0*Q_12_0) + 3*Q_11_2*Q_12_1 - 2*(Q_12_1*Q_12_1) + 3*Q_11_1*Q_12_2 + (Q_12_2*Q_12_2) + 3*QN00*Q_11_00 + 6*QN01*Q_11_01 + 6*QN02*Q_11_02 + 3*QN11*Q_11_11 + 6*QN12*Q_11_12 - 3*(QN00 + QN11)*Q_11_22) + (L2 + Ls)*(Q_00_22 + Q_01_01 - 2*Q_02_02 + 2*Q_11_11 + Q_11_22 + Q_12_12))/3.);
  }

double GEOMETRY::bulk_12(const double  QN[5],const double  dQ[15],const double  ddQ[30]) const
  { 
    return Lambda*((2*Lq_tilde*(Q_01_1 - Q_02_2 - 2.*Q_11_0) - 2.*Q_00_0*(Lq_tilde - L3*Q_12_0) - (L2 + Ls)*Q_00_12 + (L2 + Ls)*(Q_01_02 + Q_02_01 + Q_12_11 + Q_12_22) - 2*(bb*QN01*QN02 + (sigma - bb*QN00 + 2*cc*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11)))*QN12 + 2.*cc*(QN12*QN12*QN12) - L1*(Q_12_00 + Q_12_11 + Q_12_22)) + L3*(-(Q_00_1*(2.*Q_00_2 + Q_11_2)) + 2.*Q_02_0*Q_12_2 - 2.*(Q_11_2 - Q_12_1)*(Q_11_1 + Q_12_2) - Q_00_2*(Q_11_1 + 2.*Q_12_2) + 2.*(Q_01_1*(-Q_01_2 + Q_12_0) + Q_02_2*(-Q_02_1 + Q_12_0) + Q_01_0*Q_12_1 + QN00*Q_12_00 + 2*QN01*Q_12_01 + 2.*QN02*Q_12_02 + QN11*Q_12_11 + 2.*QN12*Q_12_12 - (QN00 + QN11)*Q_12_22)))/2.);
  }



void GEOMETRY::random_ic( struct Simulation_Parameters * sim_param,double * Qij )
{
  int i,j,k;
  double n[3];	
  
  
  gsl_rng_default_seed=time(NULL);
  gsl_rng *w= gsl_rng_alloc(gsl_rng_taus);
  printf("seed = %lu\n", gsl_rng_default_seed);


  for(i= 0; i< Nx; i++)
	{
      

      
	  for(j= 0; j< Ny; j++)
	    {
	  

	  
	      for(k= 0; k< Nz; k++)
		{
	    
		  
		  if(point_type[(k*Ny+j)*Nx+i] !=0 )
		    {

	
		      gsl_ran_dir_3d(w, &n[0], &n[1], &n[2]);
		      Qij[5*(Nx*(Ny*k+j)+i)+0]=0.15*(0.5*S_eq*(3.0*n[0]*n[0]-1.0));
		      Qij[5*(Nx*(Ny*k+j)+i)+1]=0.15*(0.5*S_eq*(3.0*n[0]*n[1]));
		      Qij[5*(Nx*(Ny*k+j)+i)+2]=0.15*(0.5*S_eq*(3.0*n[0]*n[2]));
		      Qij[5*(Nx*(Ny*k+j)+i)+3]=0.15*(0.5*S_eq*(3.0*n[1]*n[1]-1.0));
		      Qij[5*(Nx*(Ny*k+j)+i)+4]=0.15*(0.5*S_eq*(3.0*n[1]*n[2]));     		  


			  
		    }		
		  else
		    {

		      n[0]=0.0;
		      n[1]=0.0;
		      n[2]=1.0;
				      
		      Qij[5*(Nx*(Ny*k+j)+i)+0]=(0.5*(3.0*n[0]*n[0]-1.0));
		      Qij[5*(Nx*(Ny*k+j)+i)+1]=(0.5*(3.0*n[0]*n[1]));
		      Qij[5*(Nx*(Ny*k+j)+i)+2]=(0.5*(3.0*n[0]*n[2]));
		      Qij[5*(Nx*(Ny*k+j)+i)+3]=(0.5*(3.0*n[1]*n[1]-1.0));
		      Qij[5*(Nx*(Ny*k+j)+i)+4]=(0.5*(3.0*n[1]*n[2]));     		  
		    }	  
		}		       
	    }
	}
  gsl_rng_free(w);
}


void GEOMETRY::homogeneous_ic( struct Simulation_Parameters * sim_param,double * Qij )
{
  int i,j,k;
  double n[3];	

  printf("Initiating homogeneous initial conditions:\n\n");
  switch(sim_param->ic_flag[2])
    {
    case parameter_status::set:
      
      theta_i=sim_param->theta_i;
      printf("theta_i=%lf\n",theta_i);
      break;

    case parameter_status::unset:

      printf("Parameter \"theta_i\" not set.\n");
      printf("The initial condition named \"homogeneous\" needs the paramters \"theta_i\" and \"phi_i\" set for use.\n");
      printf("Please set them in your in your input file.");
      printf("Aborting the program.\n");
      exit(0);
      break;
    }

  
  switch(sim_param->ic_flag[3])
    {
    case parameter_status::set:
      
      phi_i=sim_param->phi_i;
      printf("phi_i=%lf\n",phi_i);
      break;

    case parameter_status::unset:

      printf("Parameter \"phi_i\" not set.\n");
      printf("The initial condition named \"homogeneous\" needs the paramters \"theta_i\" and \"phi_i\" set for use.\n");
      printf("Please set them in your in your input file.");
      printf("Aborting the program.\n");
      exit(0);
      break;
    }

      
  theta_i=theta_i*Pi/180;
  phi_i=phi_i*Pi/180;
  

  for(i= 0; i< Nx; i++)
	{
      

      
	  for(j= 0; j< Ny; j++)
	    {
	  

	  
	      for(k= 0; k< Nz; k++)
		{
	    
		  
		  if(point_type[(k*Ny+j)*Nx+i] !=0 )
		    {

		      n[0]=cos(phi_i)*sin(theta_i);
		      n[1]=sin(phi_i)*sin(theta_i);
		      n[2]=cos(theta_i);


		      Qij[5*(Nx*(Ny*k+j)+i)+0]=(0.5*S_eq*(3.0*n[0]*n[0]-1.0));
		      Qij[5*(Nx*(Ny*k+j)+i)+1]=(0.5*S_eq*(3.0*n[0]*n[1]));
		      Qij[5*(Nx*(Ny*k+j)+i)+2]=(0.5*S_eq*(3.0*n[0]*n[2]));
		      Qij[5*(Nx*(Ny*k+j)+i)+3]=(0.5*S_eq*(3.0*n[1]*n[1]-1.0));
		      Qij[5*(Nx*(Ny*k+j)+i)+4]=(0.5*S_eq*(3.0*n[1]*n[2]));     		  
			  
		    }		
		  else
		    {

		      n[0]=0.0;
		      n[1]=0.0;
		      n[2]=1.0;
				      
		      Qij[5*(Nx*(Ny*k+j)+i)+0]=(0.5*(3.0*n[0]*n[0]-1.0));
		      Qij[5*(Nx*(Ny*k+j)+i)+1]=(0.5*(3.0*n[0]*n[1]));
		      Qij[5*(Nx*(Ny*k+j)+i)+2]=(0.5*(3.0*n[0]*n[2]));
		      Qij[5*(Nx*(Ny*k+j)+i)+3]=(0.5*(3.0*n[1]*n[1]-1.0));
		      Qij[5*(Nx*(Ny*k+j)+i)+4]=(0.5*(3.0*n[1]*n[2]));     		  
		    }	  
		}		       
	    }
	}

}



void GEOMETRY::homogeneous_easy_axis_ic( struct Simulation_Parameters * sim_param,double * Qij )
{

  std::cout<< "Starting homogeneous easy axis at the boundaries.\n";

  int i,j,k;
  double n[3];	
  std::vector<double> theta_0(number_of_boundaries);      
  std::vector<double> phi_0(number_of_boundaries);    

  for(int ii=0; ii<number_of_boundaries; ii++)
    {
  
      try
	{
      
	  theta_0[ii]=sim_param->theta_0.at(ii);
	  phi_0[ii]=sim_param->phi_0.at(ii);
     
	}
      catch(std::out_of_range dummy_var)
	{

	  printf("\n Easy axis angle (theta_0 or phi_0) number %i not defined.\n",ii);
	  printf("This initial condition needs both of them defined for use.\n");
	  printf("Please define them in your in your input file.");
	  printf("Aborting the program.\n");
      
	}
    }

  

  for(i= 0; i< Nx; i++)
	{
      

      
	  for(j= 0; j< Ny; j++)
	    {
	  

	  
	      for(k= 0; k< Nz; k++)
		{
	    
		  
		  if(point_type[(k*Ny+j)*Nx+i] >1 )
		    {
		      int ptype=point_type[(k*Ny+j)*Nx+i]-2;
		      
		      n[0]=cos(phi_0[ptype]*Pi/180)*sin(theta_0[ptype]*Pi/180);
		      n[1]=sin(phi_0[ptype]*Pi/180)*sin(theta_0[ptype]*Pi/180);
		      n[2]=cos(theta_0[ptype]*Pi/180);


		      Qij[5*(Nx*(Ny*k+j)+i)+0]=(0.5*S_eq*(3.0*n[0]*n[0]-1.0));
		      Qij[5*(Nx*(Ny*k+j)+i)+1]=(0.5*S_eq*(3.0*n[0]*n[1]));
		      Qij[5*(Nx*(Ny*k+j)+i)+2]=(0.5*S_eq*(3.0*n[0]*n[2]));
		      Qij[5*(Nx*(Ny*k+j)+i)+3]=(0.5*S_eq*(3.0*n[1]*n[1]-1.0));
		      Qij[5*(Nx*(Ny*k+j)+i)+4]=(0.5*S_eq*(3.0*n[1]*n[2]));     		  
			  
		    }		
		  
		}		       
	    }
	}

}

  
  void GEOMETRY::random_bulk_homogeneous_easy_axis_ic( struct Simulation_Parameters * sim_param,double * Qij )
  {

    std::cout <<"Starting composite initial conditions: Random bulk, homogeneous on boundary.\n\n";
    random_ic( sim_param, Qij );
    homogeneous_easy_axis_ic( sim_param, Qij );
    
  }

void GEOMETRY::read_from_file_ic( struct Simulation_Parameters * sim_param, double * Qij )
{
  int i,j,k,ii,jj,kk;
  double n[3],l[3],m[3],S,P;	
  FILE * ic_file;
  char string_placeholder[400];
  int read_status;
  int reading_line=1;
  
  ic_file=fopen(sim_param->ic_file_name,"r");
  if (ic_file== NULL)
    {
      printf("Unable to find the file \"%s\".\nPlease check your initial condition file name.\n\nAborting the program.\n\n",sim_param->ic_file_name);
      exit(0);
    }

  //get the file header:
  
  printf("\nReading initial conditions from \"%s\".\n",sim_param->ic_file_name);
  
  
  fgets(string_placeholder,400,ic_file);
  reading_line++;


  //Let's work:
  
  for(k= 0; k< Nz; k++)
    {
      for(j= 0; j< Ny; j++)
	{
	  for(i= 0; i< Nx; i++)
	    {

	  
	      fgets(string_placeholder,400,ic_file);
	      read_status=sscanf(string_placeholder,"%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",&ii,&jj,&kk,&n[0],&n[1],&n[2],&l[0],&l[1],&l[2],&S,&P);
	      read_check(read_status,reading_line);


	      m[0]=n[1]*l[2]-n[2]*l[1];
	      m[1]=n[2]*l[0]-n[0]*l[2];
	      m[2]=n[0]*l[1]-n[1]*l[0];
	      
	      Qij[5*(Nx*(Ny*kk+jj)+ii)+0]=(0.5*S*(3.0*n[0]*n[0]-1.0))+(0.5*P*(l[0]*l[0]-m[0]*m[0]));
	      Qij[5*(Nx*(Ny*kk+jj)+ii)+1]=(0.5*S*(3.0*n[0]*n[1]))+(0.5*P*(l[0]*l[1]-m[0]*m[1]));
	      Qij[5*(Nx*(Ny*kk+jj)+ii)+2]=(0.5*S*(3.0*n[0]*n[2]))+(0.5*P*(l[0]*l[2]-m[0]*m[2]));
	      Qij[5*(Nx*(Ny*kk+jj)+ii)+3]=(0.5*S*(3.0*n[1]*n[1]-1.0))+(0.5*P*(l[1]*l[1]-m[1]*m[1]));
	      Qij[5*(Nx*(Ny*kk+jj)+ii)+4]=(0.5*S*(3.0*n[1]*n[2]))+(0.5*P*(l[1]*l[2]-m[1]*m[2]));

	      
	      reading_line++;
	    }
	}

    }

  
  if (reading_line-1 != Nx*Ny*Nz+1)
    {
      int expected_lines=Nx*Ny*Nz+1;
      printf("Warning: The initial condition reader read %d lines, while the it was expected to read %d lines.\n This could imply inconsistence in the size of your actyual system and size of the initial condition system.\n\n",reading_line,expected_lines);
    }
  
}


//Auxiliary routines:
void GEOMETRY::read_check(int read_status, int line)
{

  if(read_status != 11)
    {
      printf("Failed to read the inicitial condtion file line number %d between field numbers %d-%d.\n Aborting the program.\n\n",line,read_status,read_status+1);
      exit(0);
    }

}

void GEOMETRY::boundary_init( struct Simulation_Parameters * sim_param)
{

  std::cout <<"\nInitiating boundary conditions:\n\n";
  int ii=0;
  std::string anc_type;
  for (ii=0; ii< number_of_boundaries; ii++)
    {

      try
	{
	  anc_type=  sim_param->anchoring_type.at(ii);
	  std::cout << "Boundary " << ii << ":" << anc_type <<".\n";
      
	}
      catch(std::out_of_range dummy_var )
	{

	  std::cout<< "You must define boundaries "<< boundary_needed_to_be_defined <<" in the " << geometry_name <<" geometry.\nPlease review your input file.\nAborting the program.\n\n";
		   
	  exit(0);      
	}

      

      if( strcasecmp(anc_type.c_str(),"strong") == 0 || strcasecmp(anc_type.c_str(),"fixed") == 0  )

	{
    
    
	  bc_conditions[ii]=new Strong_Boundary(sim_param, ii);


	}
      else if( strcasecmp(anc_type.c_str(),"rp") == 0 || strcasecmp(anc_type.c_str(),"Rapini-Papoular") == 0  )
	{
    
    
	  bc_conditions[ii]=new Boundary_Rp(sim_param, ii);


	}
      else if( strcasecmp(anc_type.c_str(),"fg") == 0 || strcasecmp(anc_type.c_str(),"fournier-galatola") == 0  )
	{
    
    
	  bc_conditions[ii]=new Boundary_Fg(sim_param, ii);


	}
      else
	{
	  std::cout<<"There is no anchoring type named " << anc_type <<".\nPlease check your input file 'anchoring type' fields.\n\n Aborting the program.'\n";
	  exit(0);
	}

      
    }

};


void GEOMETRY::ic(struct Simulation_Parameters * sim_param,double * Qij)
{


  switch(sim_param->ic_flag[0])
    {
    case parameter_status::set:

      printf("\nInitial Conditions: %s.\n",sim_param->initial_conditions);
       if(strcasecmp(sim_param->initial_conditions,"read_from_file") == 0 )	     
	{

	  if( sim_param->ic_flag[1] == parameter_status::unset )
	    {
	      printf("Missing the \"initial_conditions_file\" in your in put file.\n Aborting the program.\n\n");
	      exit(0);
	    }
	  
	  read_from_file_ic( sim_param, Qij );

	}
      else if(strcasecmp(sim_param->initial_conditions,"homogeneous") == 0 )
	{

	  homogeneous_ic( sim_param, Qij );
	  
	}
      else if(strcasecmp(sim_param->initial_conditions,"random_bulk_homogeneous_easy_axis") == 0 )
	{

	  random_bulk_homogeneous_easy_axis_ic( sim_param, Qij );
	  
	}
      else
      if(strcasecmp(sim_param->initial_conditions,"random") == 0 || strcmp(sim_param->initial_conditions,"Random") == 0 )
	{

	  random_ic( sim_param, Qij );
	  
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



void GEOMETRY::fill_jac_bulk(const PetscScalar *Qij,Mat Jac,Mat Jac_pc, int i, int j, int k)
{


  int ip1,jp1,kp1, im1, jm1, km1, ll;
  PetscInt idxm[1],idxn[1];
  PetscScalar Jac_values[1];
  PetscErrorCode ierr;

  
  ip1= (i+1)%Nx;
  jp1= (j+1)%Ny;
  kp1= (k+1)%Nz;
  im1= i-1+((Nx-1-i)/(Nx-1))*Nx;
  jm1= j-1+((Ny-1-j)/(Ny-1))*Ny;
  km1= k-1+((Nz-1-k)/(Nz-1))*Nz;

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+0;
  Jac_values[0]=Lambda*(-2*((dx_1*dx_1) + (dy_1*dy_1) + (dz_1*dz_1))*L1 - sigma - (2*(9*cc*(Qij[5*(Nx*(Ny*k+j)+i)+0]*Qij[5*(Nx*(Ny*k+j)+i)+0]) - bb*Qij[5*(Nx*(Ny*k+j)+i)+3] + Qij[5*(Nx*(Ny*k+j)+i)+0]*(bb + 6*cc*Qij[5*(Nx*(Ny*k+j)+i)+3]) + 3*cc*((Qij[5*(Nx*(Ny*k+j)+i)+1]*Qij[5*(Nx*(Ny*k+j)+i)+1]) + (Qij[5*(Nx*(Ny*k+j)+i)+2]*Qij[5*(Nx*(Ny*k+j)+i)+2]) + (Qij[5*(Nx*(Ny*k+j)+i)+3]*Qij[5*(Nx*(Ny*k+j)+i)+3]) + (Qij[5*(Nx*(Ny*k+j)+i)+4]*Qij[5*(Nx*(Ny*k+j)+i)+4]))))/3.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*kp1+j)+i)+0;
  Jac_values[0]=Lambda*((dz_1*dz_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*km1+j)+i)+0;
  Jac_values[0]=Lambda*((dz_1*dz_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*k+jp1)+i)+0;
  Jac_values[0]=Lambda*((dy_1*dy_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*k+jm1)+i)+0;
  Jac_values[0]=Lambda*((dy_1*dy_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*k+j)+ip1)+0;
  Jac_values[0]=Lambda*((dx_1*dx_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*k+j)+im1)+0;
  Jac_values[0]=Lambda*((dx_1*dx_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+1;
  Jac_values[0]=Lambda*((-2*(bb + 6*cc*Qij[5*(Nx*(Ny*k+j)+i)+0])*Qij[5*(Nx*(Ny*k+j)+i)+1])/3.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*kp1+j)+i)+1;
  Jac_values[0]=Lambda*(dz_1*Lq_tilde);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*km1+j)+i)+1;
  Jac_values[0]=Lambda*(-(dz_1*Lq_tilde));
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+2;
  Jac_values[0]=Lambda*((-2*(bb + 6*cc*Qij[5*(Nx*(Ny*k+j)+i)+0])*Qij[5*(Nx*(Ny*k+j)+i)+2])/3.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*k+jp1)+i)+2;
  Jac_values[0]=Lambda*(-(dy_1*Lq_tilde));
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*k+jm1)+i)+2;
  Jac_values[0]=Lambda*(dy_1*Lq_tilde);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+3;
  Jac_values[0]=Lambda*((2*(bb - 3*cc*Qij[5*(Nx*(Ny*k+j)+i)+0])*(Qij[5*(Nx*(Ny*k+j)+i)+0] + 2*Qij[5*(Nx*(Ny*k+j)+i)+3]))/3.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+0;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+4;
  Jac_values[0]=Lambda*((4*(bb - 3*cc*Qij[5*(Nx*(Ny*k+j)+i)+0])*Qij[5*(Nx*(Ny*k+j)+i)+4])/3.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+0;
  Jac_values[0]=Lambda*(-(Qij[5*(Nx*(Ny*k+j)+i)+1]*(bb + 4*cc*Qij[5*(Nx*(Ny*k+j)+i)+0] + 2*cc*Qij[5*(Nx*(Ny*k+j)+i)+3])));
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*kp1+j)+i)+0;
  Jac_values[0]=Lambda*(-(dz_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*km1+j)+i)+0;
  Jac_values[0]=Lambda*((dz_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+1;
  Jac_values[0]=Lambda*(-2*((dx_1*dx_1) + (dy_1*dy_1) + (dz_1*dz_1))*L1 - sigma - 2*cc*(Qij[5*(Nx*(Ny*k+j)+i)+0]*Qij[5*(Nx*(Ny*k+j)+i)+0]) - bb*Qij[5*(Nx*(Ny*k+j)+i)+3] - Qij[5*(Nx*(Ny*k+j)+i)+0]*(bb + 2*cc*Qij[5*(Nx*(Ny*k+j)+i)+3]) - 2*cc*(3*(Qij[5*(Nx*(Ny*k+j)+i)+1]*Qij[5*(Nx*(Ny*k+j)+i)+1]) + (Qij[5*(Nx*(Ny*k+j)+i)+2]*Qij[5*(Nx*(Ny*k+j)+i)+2]) + (Qij[5*(Nx*(Ny*k+j)+i)+3]*Qij[5*(Nx*(Ny*k+j)+i)+3]) + (Qij[5*(Nx*(Ny*k+j)+i)+4]*Qij[5*(Nx*(Ny*k+j)+i)+4])));
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*kp1+j)+i)+1;
  Jac_values[0]=Lambda*((dz_1*dz_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*km1+j)+i)+1;
  Jac_values[0]=Lambda*((dz_1*dz_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+jp1)+i)+1;
  Jac_values[0]=Lambda*((dy_1*dy_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+jm1)+i)+1;
  Jac_values[0]=Lambda*((dy_1*dy_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+j)+ip1)+1;
  Jac_values[0]=Lambda*((dx_1*dx_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+j)+im1)+1;
  Jac_values[0]=Lambda*((dx_1*dx_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+2;
  Jac_values[0]=Lambda*(-4*cc*Qij[5*(Nx*(Ny*k+j)+i)+1]*Qij[5*(Nx*(Ny*k+j)+i)+2] - bb*Qij[5*(Nx*(Ny*k+j)+i)+4]);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+j)+ip1)+2;
  Jac_values[0]=Lambda*((dx_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+j)+im1)+2;
  Jac_values[0]=Lambda*(-(dx_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+3;
  Jac_values[0]=Lambda*(-(Qij[5*(Nx*(Ny*k+j)+i)+1]*(bb + 2*cc*Qij[5*(Nx*(Ny*k+j)+i)+0] + 4*cc*Qij[5*(Nx*(Ny*k+j)+i)+3])));
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*kp1+j)+i)+3;
  Jac_values[0]=Lambda*((dz_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*km1+j)+i)+3;
  Jac_values[0]=Lambda*(-(dz_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+4;
  Jac_values[0]=Lambda*(-(bb*Qij[5*(Nx*(Ny*k+j)+i)+2]) - 4*cc*Qij[5*(Nx*(Ny*k+j)+i)+1]*Qij[5*(Nx*(Ny*k+j)+i)+4]);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+jp1)+i)+4;
  Jac_values[0]=Lambda*(-(dy_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+1;
  idxn[0]=5*(Nx*(Ny*k+jm1)+i)+4;
  Jac_values[0]=Lambda*((dy_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+0;
  Jac_values[0]=Lambda*(-2*cc*Qij[5*(Nx*(Ny*k+j)+i)+2]*(2*Qij[5*(Nx*(Ny*k+j)+i)+0] + Qij[5*(Nx*(Ny*k+j)+i)+3]));
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+jp1)+i)+0;
  Jac_values[0]=Lambda*(dy_1*Lq_tilde);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+jm1)+i)+0;
  Jac_values[0]=Lambda*(-(dy_1*Lq_tilde));
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+1;
  Jac_values[0]=Lambda*(-4*cc*Qij[5*(Nx*(Ny*k+j)+i)+1]*Qij[5*(Nx*(Ny*k+j)+i)+2] - bb*Qij[5*(Nx*(Ny*k+j)+i)+4]);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+j)+ip1)+1;
  Jac_values[0]=Lambda*(-(dx_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+j)+im1)+1;
  Jac_values[0]=Lambda*((dx_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+2;
  Jac_values[0]=Lambda*(-2*((dx_1*dx_1) + (dy_1*dy_1) + (dz_1*dz_1))*L1 - sigma + bb*Qij[5*(Nx*(Ny*k+j)+i)+3] - 2*cc*((Qij[5*(Nx*(Ny*k+j)+i)+0]*Qij[5*(Nx*(Ny*k+j)+i)+0]) + (Qij[5*(Nx*(Ny*k+j)+i)+1]*Qij[5*(Nx*(Ny*k+j)+i)+1]) + 3*(Qij[5*(Nx*(Ny*k+j)+i)+2]*Qij[5*(Nx*(Ny*k+j)+i)+2]) + Qij[5*(Nx*(Ny*k+j)+i)+0]*Qij[5*(Nx*(Ny*k+j)+i)+3] + (Qij[5*(Nx*(Ny*k+j)+i)+3]*Qij[5*(Nx*(Ny*k+j)+i)+3]) + (Qij[5*(Nx*(Ny*k+j)+i)+4]*Qij[5*(Nx*(Ny*k+j)+i)+4])));
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*kp1+j)+i)+2;
  Jac_values[0]=Lambda*((dz_1*dz_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*km1+j)+i)+2;
  Jac_values[0]=Lambda*((dz_1*dz_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+jp1)+i)+2;
  Jac_values[0]=Lambda*((dy_1*dy_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+jm1)+i)+2;
  Jac_values[0]=Lambda*((dy_1*dy_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+j)+ip1)+2;
  Jac_values[0]=Lambda*((dx_1*dx_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+j)+im1)+2;
  Jac_values[0]=Lambda*((dx_1*dx_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+3;
  Jac_values[0]=Lambda*(Qij[5*(Nx*(Ny*k+j)+i)+2]*(bb - 2*cc*(Qij[5*(Nx*(Ny*k+j)+i)+0] + 2*Qij[5*(Nx*(Ny*k+j)+i)+3])));
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+jp1)+i)+3;
  Jac_values[0]=Lambda*((dy_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+jm1)+i)+3;
  Jac_values[0]=Lambda*(-(dy_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+4;
  Jac_values[0]=Lambda*(-(bb*Qij[5*(Nx*(Ny*k+j)+i)+1]) - 4*cc*Qij[5*(Nx*(Ny*k+j)+i)+2]*Qij[5*(Nx*(Ny*k+j)+i)+4]);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*kp1+j)+i)+4;
  Jac_values[0]=Lambda*((dz_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+2;
  idxn[0]=5*(Nx*(Ny*km1+j)+i)+4;
  Jac_values[0]=Lambda*(-(dz_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+0;
  Jac_values[0]=Lambda*((2*(2*Qij[5*(Nx*(Ny*k+j)+i)+0] + Qij[5*(Nx*(Ny*k+j)+i)+3])*(bb - 3*cc*Qij[5*(Nx*(Ny*k+j)+i)+3]))/3.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+1;
  Jac_values[0]=Lambda*((-2*Qij[5*(Nx*(Ny*k+j)+i)+1]*(bb + 6*cc*Qij[5*(Nx*(Ny*k+j)+i)+3]))/3.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*kp1+j)+i)+1;
  Jac_values[0]=Lambda*(-(dz_1*Lq_tilde));
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*km1+j)+i)+1;
  Jac_values[0]=Lambda*(dz_1*Lq_tilde);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+2;
  Jac_values[0]=Lambda*((4*Qij[5*(Nx*(Ny*k+j)+i)+2]*(bb - 3*cc*Qij[5*(Nx*(Ny*k+j)+i)+3]))/3.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+3;
  Jac_values[0]=Lambda*(-2*((dx_1*dx_1) + (dy_1*dy_1) + (dz_1*dz_1))*L1 - sigma - (2*(3*cc*(Qij[5*(Nx*(Ny*k+j)+i)+0]*Qij[5*(Nx*(Ny*k+j)+i)+0]) + bb*Qij[5*(Nx*(Ny*k+j)+i)+3] - Qij[5*(Nx*(Ny*k+j)+i)+0]*(bb - 6*cc*Qij[5*(Nx*(Ny*k+j)+i)+3]) + 3*cc*((Qij[5*(Nx*(Ny*k+j)+i)+1]*Qij[5*(Nx*(Ny*k+j)+i)+1]) + (Qij[5*(Nx*(Ny*k+j)+i)+2]*Qij[5*(Nx*(Ny*k+j)+i)+2]) + 3*(Qij[5*(Nx*(Ny*k+j)+i)+3]*Qij[5*(Nx*(Ny*k+j)+i)+3]) + (Qij[5*(Nx*(Ny*k+j)+i)+4]*Qij[5*(Nx*(Ny*k+j)+i)+4]))))/3.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*kp1+j)+i)+3;
  Jac_values[0]=Lambda*((dz_1*dz_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*km1+j)+i)+3;
  Jac_values[0]=Lambda*((dz_1*dz_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*k+jp1)+i)+3;
  Jac_values[0]=Lambda*((dy_1*dy_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*k+jm1)+i)+3;
  Jac_values[0]=Lambda*((dy_1*dy_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*k+j)+ip1)+3;
  Jac_values[0]=Lambda*((dx_1*dx_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*k+j)+im1)+3;
  Jac_values[0]=Lambda*((dx_1*dx_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+4;
  Jac_values[0]=Lambda*((-2*(bb + 6*cc*Qij[5*(Nx*(Ny*k+j)+i)+3])*Qij[5*(Nx*(Ny*k+j)+i)+4])/3.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*k+j)+ip1)+4;
  Jac_values[0]=Lambda*(dx_1*Lq_tilde);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+3;
  idxn[0]=5*(Nx*(Ny*k+j)+im1)+4;
  Jac_values[0]=Lambda*(-(dx_1*Lq_tilde));
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+0;
  Jac_values[0]=Lambda*((bb - 2*cc*(2*Qij[5*(Nx*(Ny*k+j)+i)+0] + Qij[5*(Nx*(Ny*k+j)+i)+3]))*Qij[5*(Nx*(Ny*k+j)+i)+4]);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+j)+ip1)+0;
  Jac_values[0]=Lambda*(-(dx_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+j)+im1)+0;
  Jac_values[0]=Lambda*((dx_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+1;
  Jac_values[0]=Lambda*(-(bb*Qij[5*(Nx*(Ny*k+j)+i)+2]) - 4*cc*Qij[5*(Nx*(Ny*k+j)+i)+1]*Qij[5*(Nx*(Ny*k+j)+i)+4]);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+jp1)+i)+1;
  Jac_values[0]=Lambda*((dy_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+jm1)+i)+1;
  Jac_values[0]=Lambda*(-(dy_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+2;
  Jac_values[0]=Lambda*(-(bb*Qij[5*(Nx*(Ny*k+j)+i)+1]) - 4*cc*Qij[5*(Nx*(Ny*k+j)+i)+2]*Qij[5*(Nx*(Ny*k+j)+i)+4]);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*kp1+j)+i)+2;
  Jac_values[0]=Lambda*(-(dz_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*km1+j)+i)+2;
  Jac_values[0]=Lambda*((dz_1*Lq_tilde)/2.);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+3;
  Jac_values[0]=Lambda*(-2*cc*(Qij[5*(Nx*(Ny*k+j)+i)+0] + 2*Qij[5*(Nx*(Ny*k+j)+i)+3])*Qij[5*(Nx*(Ny*k+j)+i)+4]);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+j)+ip1)+3;
  Jac_values[0]=Lambda*(-(dx_1*Lq_tilde));
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+j)+im1)+3;
  Jac_values[0]=Lambda*(dx_1*Lq_tilde);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+j)+i)+4;
  Jac_values[0]=Lambda*(-2*((dx_1*dx_1) + (dy_1*dy_1) + (dz_1*dz_1))*L1 - sigma - 2*cc*(Qij[5*(Nx*(Ny*k+j)+i)+0]*Qij[5*(Nx*(Ny*k+j)+i)+0]) + Qij[5*(Nx*(Ny*k+j)+i)+0]*(bb - 2*cc*Qij[5*(Nx*(Ny*k+j)+i)+3]) - 2*cc*((Qij[5*(Nx*(Ny*k+j)+i)+1]*Qij[5*(Nx*(Ny*k+j)+i)+1]) + (Qij[5*(Nx*(Ny*k+j)+i)+2]*Qij[5*(Nx*(Ny*k+j)+i)+2]) + (Qij[5*(Nx*(Ny*k+j)+i)+3]*Qij[5*(Nx*(Ny*k+j)+i)+3]) + 3*(Qij[5*(Nx*(Ny*k+j)+i)+4]*Qij[5*(Nx*(Ny*k+j)+i)+4])));
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*kp1+j)+i)+4;
  Jac_values[0]=Lambda*((dz_1*dz_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*km1+j)+i)+4;
  Jac_values[0]=Lambda*((dz_1*dz_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+jp1)+i)+4;
  Jac_values[0]=Lambda*((dy_1*dy_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+jm1)+i)+4;
  Jac_values[0]=Lambda*((dy_1*dy_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+j)+ip1)+4;
  Jac_values[0]=Lambda*((dx_1*dx_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);

  idxm[0]=5*(Nx*(Ny*k+j)+i)+4;
  idxn[0]=5*(Nx*(Ny*k+j)+im1)+4;
  Jac_values[0]=Lambda*((dx_1*dx_1)*L1);
  MatSetValues(Jac,1,idxm,1,idxn,Jac_values,ADD_VALUES);


}


