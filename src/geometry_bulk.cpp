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
#include <petscts.h>

Geometry_Bulk::Geometry_Bulk(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param)
{
  point_type=fill_point_type( );
  geometry_pointer=&(*this);

  geometry_name="Bulk";
  number_of_boundaries=0;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);

  RhsPtr=RhsFunction;
  JacobianPtr=RhsJacobian;

  
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






PetscErrorCode Geometry_Bulk::RhsFunction(TS ts,PetscReal time,Vec Qij_in,Vec Rhs, void *sim_geometry)
{

  GEOMETRY * sample_geometry=(GEOMETRY *) sim_geometry;
  int i, j, k;  
  int ip1,jp1,kp1, im1, jm1, km1, ll;
  double dQ[15];
  double ddQ[30];
  double QN[5];
  double v[3];
  const PetscScalar *Qij;
  PetscScalar * k_i;
  
  const int Nx=sample_geometry->Nx;
  const int Ny=sample_geometry->Ny;
  const int Nz=sample_geometry->Nz;  

  const PetscScalar dx_1= sample_geometry->dx_1;
  const PetscScalar dy_1= sample_geometry->dy_1;
  const PetscScalar dz_1= sample_geometry->dz_1;  

  
  PetscErrorCode ierr;

  VecGetArrayRead(Qij_in,&Qij);
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
  VecRestoreArrayRead(Qij_in,&Qij);
  return ierr;
}

PetscErrorCode Geometry_Bulk::RhsJacobian(TS ts,PetscReal time,Vec Qij_in,Mat Jac,Mat Jac_pc, void* sim_geometry)
{



  GEOMETRY * sample_geometry=(GEOMETRY *) sim_geometry;
  int i, j, k;  
  int ip1,jp1,kp1, im1, jm1, km1, ll;
  double dQ[15];
  double ddQ[30];
  double QN[5];
  double v[3];
  const PetscScalar *Qij;

  PetscInt idxm[1],idxn[1];
  PetscScalar Jac_values[1];

  const PetscScalar sigma=sample_geometry->sigma;
  const PetscScalar bb=sample_geometry->bb;
  const PetscScalar cc=sample_geometry->cc;
  const PetscScalar Lambda=sample_geometry->Lambda;
  const PetscScalar L1=sample_geometry->L1;
  const PetscScalar Lq_tilde=sample_geometry->Lq_tilde;
  
  const int Nx=sample_geometry->Nx;
  const int Ny=sample_geometry->Ny;
  const int Nz=sample_geometry->Nz;  

  const PetscScalar dx_1= sample_geometry->dx_1;
  const PetscScalar dy_1= sample_geometry->dy_1;
  const PetscScalar dz_1= sample_geometry->dz_1;  

  
  PetscErrorCode ierr;

  VecGetArrayRead(Qij_in,&Qij);
 

  for( i= 0; i< Nx; i++)
    {
      for( j= 0; j< Ny; j++)
	{
	  for( k= 0; k< Nz; k++)
	    {

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


	}


    }


  MatAssemblyBegin(Jac, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Jac, MAT_FINAL_ASSEMBLY);
  VecRestoreArrayRead(Qij_in,&Qij);
  return ierr;

}
