
#include "driver.h"
#include "omp.h"
#include <iostream>           
#include <cstdio>           
#include <cstdlib>          
#include <cstring>
#include <cmath>            
#include "geometry.h"
#include "integrator.h"
#include "integrator_CN.h"

#define NUMBEROFCLOSEPOINTS 27*5


static char help[] ="Solves one dimensional Burger's equation compares with exact solution\n\n";

PetscErrorCode RhsFunction(TS ts,PetscReal tiem,Vec Qij_in,Vec Rhs,void *sim_geometry);



CN::CN( GEOMETRY  * lc_pointer, const struct Simulation_Parameters *sim_param ) : Integrator( lc_pointer),
							 	   Atol(sim_param->Atol),
								   Rtol(sim_param->Rtol),
								   prefac(sim_param->prefac),
								   facmin(sim_param->facmin),
                                                                   facmax(sim_param->facmax)
{

  PetscErrorCode ierr;


  int argc = sim_param->argc;
  char ** argv=sim_param->argv;
  ierr = PetscInitialize(&argc,& argv,(char*)0,help);
  
  dt=sim_param->dt;
  

  // Creation of the solution vextors
  VecCreate( PETSC_COMM_SELF , &Qsolution  );
  VecSetSizes(Qsolution, PETSC_DECIDE, 5*Nx*Ny*Nz );
  VecSetType(Qsolution, VECSEQ  );
  
  VecDuplicate( Qsolution , &Qtij  );
  VecDuplicate( Qsolution , &Rhs  );


  //Creating Matrixes:

  //MatCreateSeqAIJ(PETSC_COMM_SELF, 5*Nx*Ny*Nz, 5*Nx*Ny*Nz, NUMBEROFCLOSEPOINTS , NULL, &Jac );

  TSAdapt * adaptor_ref;
  
  //Creating the Integrator:
  TSCreate(PETSC_COMM_SELF, &cranck_int);
  TSSetProblemType(cranck_int , TS_NONLINEAR);
  TSSetType( cranck_int , TSRK );
  TSRKSetType(cranck_int , TSRK5BS); 

    
  TSGetAdapt(cranck_int, adaptor_ref);
  TSAdaptSetType( *adaptor_ref, TSADAPTBASIC);
  TSAdaptSetClip( *adaptor_ref, 0.2, 5);
    
  TSSetTolerances(cranck_int , Atol,NULL,Rtol,NULL);

    

    


    
    
};


void CN::getIcArray(double * Qij_input)
{
  PetscScalar * Qij;
  int i,j, k, ll;

  VecGetArray(Qsolution,&Qij);
     
    for( k= 0; k< Nz; k++)
      {
	for( j= 0; j< Ny; j++)
	  {
               
	    for( i= 0; i< Nx; i++)
	      {
		for (ll=0; ll<5; ll++)
		  {
		   
                    Qij[5*((k*Ny+j)*Nx+i)+ll]=Qij_input[5*((k*Ny+j)*Nx+i)+ll];
                    
                  }
		   
	      }
            
	  }

      };
     

    VecRestoreArray(Qsolution,&Qij);

};


void CN::writeSolutionToArray(double * Qij_out)
{
  PetscScalar * Qij;
  int i,j, k, ll;

  VecGetArray(Qsolution,&Qij);
     
    for( k= 0; k< Nz; k++)
      {
	for( j= 0; j< Ny; j++)
	  {
               
	    for( i= 0; i< Nx; i++)
	      {
		for (ll=0; ll<5; ll++)
		  {
		   
                    Qij_out[5*((k*Ny+j)*Nx+i)+ll]=Qij[5*((k*Ny+j)*Nx+i)+ll];
                    
                  }
		   
	      }
            
	  }

      };
     

    VecRestoreArray(Qsolution,&Qij);

};


void CN::evolve( double * Qij, double *time, double tf )
{

  int i,j,k,ll,information_step=1;
  double local_error;
  double global_error; //, global_error_1=1.;
  double sc_i;
  double hfactor=1.0;
  double dt=this->dt;
  //const int chunk_size=0.06*(4*Ny*Nz)/omp_get_num_threads();
  //const int chunk_size=1;
  

  getIcArray(Qij);

  
  TSSetSolution(cranck_int, Qsolution);

  
  TSSetTime(cranck_int , *time);
  TSGetTimeStep(cranck_int , &dt);
  TSSetMaxTime(cranck_int , tf);
  TSSetExactFinalTime(cranck_int , TS_EXACTFINALTIME_MATCHSTEP);

  
  TSSetRHSFunction(cranck_int,Rhs, RhsFunction,  sample_geometry);
  TSSolve(cranck_int , NULL);

  *time=tf;
  


  writeSolutionToArray(Qij);


	      

}; 

CN::~CN(void)
{

PetscFinalize();
  
}



PetscErrorCode RhsFunction(TS ts,PetscReal time,Vec Qij_in,Vec Rhs, void *sim_geometry)
{

  GEOMETRY * sample_geometry=(GEOMETRY *) sim_geometry;
  int i, j, k;  
  int ip1,jp1,kp1, im1, jm1, km1, ll;
  double dQ[15];
  double ddQ[30];
  double QN[5];
  double v[3];
  PetscScalar *Qij, *k_i;
  
  const int Nx=sample_geometry->Nx;
  const int Ny=sample_geometry->Ny;
  const int Nz=sample_geometry->Nz;  

  const PetscScalar dx_1= sample_geometry->dx_1;
  const PetscScalar dy_1= sample_geometry->dy_1;
  const PetscScalar dz_1= sample_geometry->dz_1;  

  
  PetscErrorCode ierr;

  VecGetArray(Qij_in,&Qij);
  VecGetArray(Rhs, & k_i);

  for( i= 0; i< Nx; i++)
    {
      for( j= 0; j< Ny; j++)
	{
	  for( k= 0; k< Nz; k++)
	    {


	      ip1= (i+1)%Nx;
	      jp1= (j+1)%Ny;
	      kp1= (k+1)%Nz;
	      im1= i-1+((Nx-1-i)/(Nx-1))*Nx;
	      jm1= j-1+((Ny-1-j)/(Ny-1))*Ny;
	      km1= k-1+((Nz-1-k)/(Nz-1))*Nz;

	      //Remember to correct this bug in the production version. (i,j,k) in wrong places.
      
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
      
  
	      k_i[5*(Nx*(Ny*k+j)+i)+0]=sample_geometry->bulk_00(QN,dQ,ddQ); 
	      k_i[5*(Nx*(Ny*k+j)+i)+1]=sample_geometry->bulk_01(QN,dQ,ddQ); 
	      k_i[5*(Nx*(Ny*k+j)+i)+2]=sample_geometry->bulk_02(QN,dQ,ddQ); 
	      k_i[5*(Nx*(Ny*k+j)+i)+3]=sample_geometry->bulk_11(QN,dQ,ddQ); 
	      k_i[5*(Nx*(Ny*k+j)+i)+4]=sample_geometry->bulk_12(QN,dQ,ddQ);


	      
	    }

	}

    }
  
  VecRestoreArray(Rhs, & k_i);
  VecRestoreArray(Qij_in,&Qij);
  return ierr;
}
