#include "geometry.h"
#include "driver.h"
#include "boundary.h"
#include <stdio.h> 
#include <gsl/gsl_randist.h>
#include <ctime>
#include <math.h>
#include <vector>

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

