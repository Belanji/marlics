#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <ctime>
#include <cstring>
#include <iostream>

#include "driver.h"
#include "geometry.h"
#include "geometry_slab.h"
#include "boundary.h"
#include <gsl/gsl_randist.h>
#include <petscts.h>

#define MAXX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)


slab::slab(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param)
{
  point_type=fill_point_type( );
  geometry_pointer=&(*this);

  geometry_name="Slab";
  number_of_boundaries=2;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
  boundary_needed_to_be_defined="0 and 1";



  RhsPtr=slab::RhsFunction;
  JacobianPtr=slab::Jacobian;  


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

      ip1= (i+1)%Nx;
      jp1= (j+1)%Ny;
      kp1= (k+1);
      im1= i-1+((Nx-1-i)/(Nx-1))*Nx;
      jm1= j-1+((Ny-1-i)/(Ny-1))*Ny;
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


      for(ll=0; ll<=4;ll++) ddQ[20+ll]= 0.25*(Qij[5*((kp1*Ny+jp1)*Nx+i)+ll]+Qij[5*((km1*Ny+jm1)*Nx+i)+ll]-Qij[5*((kp1*Ny+jm1)*Nx+i)+ll]-Qij[5*((km1*Ny+jp1)*Nx+i)+ll])*dy_1*dz_1;
      
      for(ll=0; ll<=4;ll++) ddQ[25+ll]= (Qij[5*((kp1*Ny+j)*Nx+i)+ll]-2.0*QN[ll]+Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1*dz_1;
      
  
      k_i[5*(Nx*(Ny*k+j)+i)+0]=bulk_00(QN,dQ,ddQ); 
      k_i[5*(Nx*(Ny*k+j)+i)+1]=bulk_01(QN,dQ,ddQ); 
      k_i[5*(Nx*(Ny*k+j)+i)+2]=bulk_02(QN,dQ,ddQ); 
      k_i[5*(Nx*(Ny*k+j)+i)+3]=bulk_11(QN,dQ,ddQ); 
      k_i[5*(Nx*(Ny*k+j)+i)+4]=bulk_12(QN,dQ,ddQ);
      
    }  
  else if( point_type[(k*Ny+j)*Nx+i] == 2 )
    {

      //check_surface_limits( i,  j,  k);
      ip1= (i+1)%Nx;
      jp1= (j+1)%Ny;
      kp1= (k+1);
      im1= i-1+((Nx-1-i)/(Nx-1))*Nx;
      jm1= j-1+((Ny-1-i)/(Ny-1))*Ny;
      km1= k;
      v[0]=0.0;
      v[1]=0.0;
      v[2]=-1.0;
      
      for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
	

      //Calcule first derivatives of Qij:
      for(ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;

      
      for(ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;

      
      for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((k*Ny+j)*Nx+i)+ll])*dz_1;
      
  
      k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions[0]->surface_00(QN,dQ,ddQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+1]= bc_conditions[0]->surface_01(QN,dQ,ddQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+2]= bc_conditions[0]->surface_02(QN,dQ,ddQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+3]= bc_conditions[0]->surface_11(QN,dQ,ddQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+4]= bc_conditions[0]->surface_12(QN,dQ,ddQ,v); 


      
    }
  else if( point_type[(k*Ny+j)*Nx+i] == 3 )
    {

      
      ip1= (i+1)%Nx;
      jp1= (j+1)%Ny;
      kp1= k;
      im1= i-1+((Nx-1-i)/(Nx-1))*Nx;
      jm1= j-1+((Ny-1-i)/(Ny-1))*Ny;
      km1= (k-1);
      v[0]=0.0;
      v[1]=0.0;
      v[2]=1.0;

      
      for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
	

      //Calcule first derivatives of Qij:
      for(ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;

      
      for(ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dx_1;

      
      for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((k*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dx_1;
      
  
      k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions[1]->surface_00(QN,dQ,ddQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+1]= bc_conditions[1]->surface_01(QN,dQ,ddQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+2]= bc_conditions[1]->surface_02(QN,dQ,ddQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+3]= bc_conditions[1]->surface_11(QN,dQ,ddQ,v);
      k_i[5*(Nx*(Ny*k+j)+i)+4]= bc_conditions[1]->surface_12(QN,dQ,ddQ,v); 


      //check_surface_limits(i,j,k);
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





PetscErrorCode slab::Jacobian(TS ts,PetscReal time,Vec Qij_in,Mat Jac,Mat Jac_pc, void* sim_geometry)
{
  GEOMETRY * sample_geometry=(GEOMETRY *) sim_geometry;

  int i, j, k;
  int ip1,jp1,kp1, im1, jm1, km1, ll;
  PetscScalar dQ[15];
  PetscScalar ddQ[30];
  PetscScalar QN[5];
  PetscScalar v[3];
  PetscScalar dFijdQij[5][5];
  PetscScalar dFijdQijk[5][5][3];
  
  const PetscScalar * Qij;
  PetscErrorCode ierr;

  const int *point_type= sample_geometry->point_type;;


  const int Nx=sample_geometry->Nx;
  const int Ny=sample_geometry->Ny;
  const int Nz=sample_geometry->Nz;

  const PetscScalar dx_1=sample_geometry->dx_1;
  const PetscScalar dy_1=sample_geometry->dy_1;
  const PetscScalar dz_1=sample_geometry->dz_1;


  
  VecGetArrayRead(Qij_in, &Qij);

  for( i= 0; i< Nx; i++)
    {
      for( j= 0; j< Ny; j++)
	{
	  for( k= 0; k< Nz; k++)
	    {

	      if(point_type[(k*Ny+j)*Nx+i] == 1)
                {

                  sample_geometry->fill_jac_bulk( Qij, Jac, Jac_pc,  i,  j,  k);      
      
                }
//              else if(point_type[(k*Ny+j)*Nx+i] == 2)
//                {
//
//                  PetscInt idxm[5];
//		  PetscInt idxn[5];
//		  PetscScalar jacobian_values[25];
//		  PetscInt ii, jj, kk, ll;
//		  
//                  ip1= (i+1)%Nx;
//                  jp1= (j+1)%Ny;
//                  kp1= (k+1);
//                  im1= i-1+((Nx-1-i)/(Nx-1))*Nx;
//                  jm1= j-1+((Ny-1-i)/(Ny-1))*Ny;
//                  km1= k;
//                  v[0]=0.0;
//                  v[1]=0.0;
//                  v[2]=-1.0;
//
//		  PetscInt iis[2]={ip1,im1};
//		  PetscScalar dQijxdQij[2]={0.5*dx_1,-0.5*dx_1};
//
//		  PetscInt jjs[2]={jp1,jm1};
//		  PetscScalar dQijydQij[2]={0.5*dy_1,-0.5*dy_1};
//
//		  PetscInt kks[2]={kp1,km1};
//		  PetscScalar dQijzdQij[2]={dz_1,-dz_1};
//
//		  
//		  //Clculating rows indexes:
//		  for(ll=0; ll<=4;ll++) idxm[ll]=5*(Nx*(Ny*k+j)+i)+ll;
//
//		  
//                  for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];	
//
//                  //Calcule first derivatives of Qij:
//                  for(ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
//                  for(ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
//                  for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((k*Ny+j)*Nx+i)+ll])*dz_1;
//
//		  //calculate de derivatves of Fij respect to Qij,k:
//                  sample_geometry->bc_conditions[0]->fill_dFijQij(dFijdQij,QN, dQ,ddQ,v);
//                  sample_geometry->bc_conditions[0]->fill_dFijQijk(dFijdQijk, QN, dQ, ddQ, v);
//
//
//		  //Filling the dFij/dQij:
//                  for(ll=0; ll<=4;ll++) idxn[ll]=5*(Nx*(Ny*k+j)+i)+ll;
//
//		  for(ii=0; ii<=4; ii++)
//		    {
//		      for(jj=0; jj<=4; jj++)
//			{
//
//			  jacobian_values[5*ii+jj]=dFijdQij[ii][jj];
//			  
//			}
//
//		    }
//
//		  MatSetValues(Jac,5,idxm,5,idxn,jacobian_values,ADD_VALUES);
//
//
//
//
//		  //Filling the dFij/dQij_x:
//		  for(int pm=0; pm<2; pm++)
//		    {
//		      for(ll=0; ll<=4;ll++) idxn[ll]=5*(Nx*(Ny*k+j)+iis[pm])+ll;
//
//		      for(ii=0; ii<=4; ii++)
//			{
//			  for(jj=0; jj<=4; jj++)
//			    {
//
//			      jacobian_values[5*ii+jj]=dFijdQijk[ii][jj][0]*dQijxdQij[pm];;
//			  
//			    }
//
//			}
//
//		      MatSetValues(Jac,5,idxm,5,idxn,jacobian_values,ADD_VALUES);
//		    }
//
//
//
//		  
//		  //Filling the dFij/dQij_y:
//		  for(int pm=0; pm<2; pm++)
//		    {
//		      for(ll=0; ll<=4;ll++) idxn[ll]=5*(Nx*(Ny*k+jjs[pm])+i)+ll;
//
//		      for(ii=0; ii<=4; ii++)
//			{
//			  for(jj=0; jj<=4; jj++)
//			    {
//
//			      jacobian_values[5*ii+jj]=dFijdQijk[ii][jj][1]*dQijydQij[pm];;
//			  
//			    }
//
//			}
//
//		      MatSetValues(Jac,5,idxm,5,idxn,jacobian_values,ADD_VALUES);
//		    }
//
//
//		  //Filling the dFij/dQij_z:
//		  for(int pm=0; pm<2; pm++)
//		    {
//		      for(ll=0; ll<=4;ll++) idxn[ll]=5*(Nx*(Ny*kks[pm]+j)+i)+ll;
//
//		      for(ii=0; ii<=4; ii++)
//			{
//			  for(jj=0; jj<=4; jj++)
//			    {
//
//			      jacobian_values[5*ii+jj]=dFijdQijk[ii][jj][2]*dQijzdQij[pm];;
//			  
//			    }
//
//			}
//
//		      MatSetValues(Jac,5,idxm,5,idxn,jacobian_values,ADD_VALUES);
//		    }
//
//
//
//                }
//	      else if(point_type[(k*Ny+j)*Nx+i] == 3)
//                {
//
//                  PetscInt idxm[5];
//		  PetscInt idxn[5];
//		  PetscScalar jacobian_values[25];
//		  PetscInt ii, jj, kk, ll;
//		  
//                  ip1= (i+1)%Nx;
//                  jp1= (j+1)%Ny;
//                  kp1= k;
//                  im1= i-1+((Nx-1-i)/(Nx-1))*Nx;
//                  jm1= j-1+((Ny-1-i)/(Ny-1))*Ny;
//                  km1= k-1;
//                  v[0]=0.0;
//                  v[1]=0.0;
//                  v[2]=1.0;
//
//		  PetscInt iis[2]={ip1,im1};
//		  PetscScalar dQijxdQij[2]={0.5*dx_1,-0.5*dx_1};
//
//		  PetscInt jjs[2]={jp1,jm1};
//		  PetscScalar dQijydQij[2]={0.5*dy_1,-0.5*dy_1};
//
//		  PetscInt kks[2]={kp1,km1};
//		  PetscScalar dQijzdQij[2]={dz_1,-dz_1};
//
//		  
//		  //Clculating rows indexes:
//		  for(ll=0; ll<=4;ll++) idxm[ll]=5*(Nx*(Ny*k+j)+i)+ll;
//
//		  
//                  for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];	
//
//                  //Calcule first derivatives of Qij:
//                  for(ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
//                  for(ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
//                  for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((k*Ny+j)*Nx+i)+ll])*dz_1;
//
//		  //calculate de derivatves of Fij respect to Qij,k:
//                  sample_geometry->bc_conditions[1]->fill_dFijQij(dFijdQij,QN, dQ,ddQ,v);
//                  sample_geometry->bc_conditions[1]->fill_dFijQijk(dFijdQijk, QN, dQ, ddQ, v);
//
//
//		  //Filling the dFij/dQij:
//                  for(ll=0; ll<=4;ll++) idxn[ll]=5*(Nx*(Ny*k+j)+i)+ll;
//
//		  for(ii=0; ii<=4; ii++)
//		    {
//		      for(jj=0; jj<=4; jj++)
//			{
//
//			  jacobian_values[5*ii+jj]=dFijdQij[ii][jj];
//			  
//			}
//
//		    }
//
//		  MatSetValues(Jac,5,idxm,5,idxn,jacobian_values,ADD_VALUES);
//
//
//
//
//		  //Filling the dFij/dQij_x:
//		  for(int pm=0; pm<2; pm++)
//		    {
//		      for(ll=0; ll<=4;ll++) idxn[ll]=5*(Nx*(Ny*k+j)+iis[pm])+ll;
//
//		      for(ii=0; ii<=4; ii++)
//			{
//			  for(jj=0; jj<=4; jj++)
//			    {
//
//			      jacobian_values[5*ii+jj]=dFijdQijk[ii][jj][0]*dQijxdQij[pm];;
//			  
//			    }
//
//			}
//
//		      MatSetValues(Jac,5,idxm,5,idxn,jacobian_values,ADD_VALUES);
//		    }
//
//
//
//		  
//		  //Filling the dFij/dQij_y:
//		  for(int pm=0; pm<2; pm++)
//		    {
//		      for(ll=0; ll<=4;ll++) idxn[ll]=5*(Nx*(Ny*k+jjs[pm])+i)+ll;
//
//		      for(ii=0; ii<=4; ii++)
//			{
//			  for(jj=0; jj<=4; jj++)
//			    {
//
//			      jacobian_values[5*ii+jj]=dFijdQijk[ii][jj][1]*dQijydQij[pm];;
//			  
//			    }
//
//			}
//
//		      MatSetValues(Jac,5,idxm,5,idxn,jacobian_values,ADD_VALUES);
//		    }
//
//
//		  //Filling the dFij/dQij_z:
//		  for(int pm=0; pm<2; pm++)
//		    {
//		      for(ll=0; ll<=4;ll++) idxn[ll]=5*(Nx*(Ny*kks[pm]+j)+i)+ll;
//
//		      for(ii=0; ii<=4; ii++)
//			{
//			  for(jj=0; jj<=4; jj++)
//			    {
//
//			      jacobian_values[5*ii+jj]=dFijdQijk[ii][jj][2]*dQijzdQij[pm];;
//			  
//			    }
//
//			}
//
//		      MatSetValues(Jac,5,idxm,5,idxn,jacobian_values,ADD_VALUES);
//		    }
//
//
//
//                }

	    }

	}
    }

  MatAssemblyBegin(Jac, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Jac, MAT_FINAL_ASSEMBLY);
  VecRestoreArrayRead(Qij_in,&Qij);
  return ierr;

  
}

PetscErrorCode slab::RhsFunction(TS ts,PetscReal time,Vec Qij_in,Vec Rhs,void *sim_geometry)
{

  GEOMETRY * sample_geometry=(GEOMETRY *) sim_geometry;
  PetscInt i, j, k;  
  PetscInt ip1,jp1,kp1, im1, jm1, km1, ll;
  double dQ[15];
  double ddQ[30];
  double QN[5];
  double v[3];
  const PetscScalar *Qij;
  PetscScalar * k_i;

  const int *point_type= sample_geometry->point_type;
  
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
  
	      if(point_type[(k*Ny+j)*Nx+i] == 1)
		{

		  //check_bulk_limits( i,  j,  k);	

		  ip1= (i+1)%Nx;
		  jp1= (j+1)%Ny;
		  kp1= (k+1);
		  im1= i-1+((Nx-1-i)/(Nx-1))*Nx;
		  jm1= j-1+((Ny-1-i)/(Ny-1))*Ny;
		  km1= (k-1);

      
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
//	      else if( point_type[(k*Ny+j)*Nx+i] == 2 )
//		{
//
//		  //check_surface_limits( i,  j,  k);
//		  ip1= (i+1)%Nx;
//		  jp1= (j+1)%Ny;
//		  kp1= (k+1);
//		  im1= i-1+((Nx-1-i)/(Nx-1))*Nx;
//		  jm1= j-1+((Ny-1-i)/(Ny-1))*Ny;
//		  km1= k;
//		  v[0]=0.0;
//		  v[1]=0.0;
//		  v[2]=-1.0;
//      
//		  for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
//	
//
//		  //Calcule first derivatives of Qij:
//		  for(ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
//
//      
//		  for(ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
//
//      
//		  for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((k*Ny+j)*Nx+i)+ll])*dz_1;
//      
//  
//		  k_i[5*(Nx*(Ny*k+j)+i)+0]= sample_geometry->bc_conditions[0]->surface_00(QN,dQ,ddQ,v);
//		  k_i[5*(Nx*(Ny*k+j)+i)+1]= sample_geometry->bc_conditions[0]->surface_01(QN,dQ,ddQ,v);
//		  k_i[5*(Nx*(Ny*k+j)+i)+2]= sample_geometry->bc_conditions[0]->surface_02(QN,dQ,ddQ,v);
//		  k_i[5*(Nx*(Ny*k+j)+i)+3]= sample_geometry->bc_conditions[0]->surface_11(QN,dQ,ddQ,v);
//		  k_i[5*(Nx*(Ny*k+j)+i)+4]= sample_geometry->bc_conditions[0]->surface_12(QN,dQ,ddQ,v); 
//
//		}
//	      else if( point_type[(k*Ny+j)*Nx+i] == 3 )
//		{
//
//      
//		  ip1= (i+1)%Nx;
//		  jp1= (j+1)%Ny;
//		  kp1= k;
//		  im1= i-1+((Nx-1-i)/(Nx-1))*Nx;
//		  jm1= j-1+((Ny-1-i)/(Ny-1))*Ny;
//		  km1= (k-1);
//		  v[0]=0.0;
//		  v[1]=0.0;
//		  v[2]=1.0;
//
//      
//		  for(ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
//	
//
//		  //Calcule first derivatives of Qij:
//		  for(ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
//
//      
//		  for(ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dx_1;
//
//      
//		  for(ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((k*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dx_1;
//      
//  
//		  k_i[5*(Nx*(Ny*k+j)+i)+0]= sample_geometry->bc_conditions[1]->surface_00(QN,dQ,ddQ,v);
//		  k_i[5*(Nx*(Ny*k+j)+i)+1]= sample_geometry->bc_conditions[1]->surface_01(QN,dQ,ddQ,v);
//		  k_i[5*(Nx*(Ny*k+j)+i)+2]= sample_geometry->bc_conditions[1]->surface_02(QN,dQ,ddQ,v);
//		  k_i[5*(Nx*(Ny*k+j)+i)+3]= sample_geometry->bc_conditions[1]->surface_11(QN,dQ,ddQ,v);
//		  k_i[5*(Nx*(Ny*k+j)+i)+4]= sample_geometry->bc_conditions[1]->surface_12(QN,dQ,ddQ,v); 
//
//
//		}
	    }
	}
    }
  
  VecRestoreArrayRead(Qij_in,&Qij);
  VecRestoreArray(Rhs, & k_i);
  return ierr;
  
}
  
