#include "driver.h"
#include "omp.h"
#include <iostream>           
#include <cstdio>           
#include <cstdlib>          
#include <cstring>
#include <cmath>            
#include "geometry.h"
#include "geometry_slab.h"

#include "integrator.h"
#include "integrator_CN.h"
#include <petscsnes.h>
#include "geometry_bulk.h"

#define NUMBEROFCLOSEPOINTS 75


static char help[] ="Nothing for Now";


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

  MatCreateSeqAIJ(PETSC_COMM_SELF, 5*Nx*Ny*Nz, 5*Nx*Ny*Nz, NUMBEROFCLOSEPOINTS , NULL, &Jac );
  MatSetUp(Jac);

  
  //TSAdapt * adaptor_ref;
  
  //Creating the Integrator:
  TSCreate(PETSC_COMM_SELF, &cranck_int);
  TSSetProblemType(cranck_int , TS_NONLINEAR);
  TSSetType( cranck_int , TSCN );



  //Setuping the non-linear solver:
  
  SNES  snes;
  SNESLineSearch linesearch;
  
  TSGetSNES(cranck_int , &snes);
  SNESSetTolerances(snes, Atol, Rtol, (Atol+Rtol) , 10,10);
  SNESSetType(snes, SNESNEWTONLS);
  
  KSP ksp;
  SNESGetKSP( snes, &ksp);
  KSPSetType(ksp, KSPGMRES);
  KSPGMRESSetRestart(ksp, 100);
  KSPSetPCSide( ksp,  PC_LEFT );

  
  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc,  PCSOR);
  PCSORSetIterations(pc, 2,1);



  TSSetExactFinalTime(cranck_int , TS_EXACTFINALTIME_MATCHSTEP);
  TSSetRHSFunction(cranck_int,Rhs, sample_geometry->RhsPtr,  sample_geometry);
  TSSetRHSJacobian(cranck_int,Jac,Jac,sample_geometry->JacobianPtr,sample_geometry);
    
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
  
  
  TSSolve(cranck_int, NULL );
  *time=tf;
  


  writeSolutionToArray(Qij);


	      

}; 

CN::~CN(void)
{

PetscFinalize();
  
}



