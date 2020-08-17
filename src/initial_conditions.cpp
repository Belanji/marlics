#include "driver.h"
#include "geometry.h"
#include "initial_conditions.h"
#include <gsl/gsl_randist.h>
#include <ctime>
#include <math.h>
#include <string>
#include <cstring>
#include <iostream>

const double Pi=M_PI;

void apply_initial_conditions(struct Simulation_Parameters * sim_param,double * Qij, const GEOMETRY * geometry)
{


  switch(sim_param->ic_flag[0])
    {
      case parameter_status::set:
      
        std::cout << "\nInitial Conditions: " << sim_param->initial_conditions << ".\n";
         if(strcasecmp(sim_param->initial_conditions,"read_from_file") == 0 )          
          {
      
            if( sim_param->ic_flag[1] == parameter_status::unset )
              {
                std::cout<<"Missing the \"initial_conditions_file\" in your in put file.\n Aborting the program.\n\n";
                exit(0);
              }
            
            read_from_file_ic( sim_param, Qij, geometry );
      
          }
        else if(strcasecmp(sim_param->initial_conditions,"homogeneous_easy_axis") == 0 )
          {
            homogeneous_easy_axis_ic( sim_param, Qij, geometry );
            
          }
        else if(strcasecmp(sim_param->initial_conditions,"homogeneous") == 0 )
          {
      
            homogeneous_ic( sim_param, Qij, geometry );
            
          }
        else if(strcasecmp(sim_param->initial_conditions,"random_bulk_homogeneous_easy_axis") == 0 )
          {
      
            random_bulk_homogeneous_easy_axis_ic( sim_param, Qij, geometry );
            
          }
        else if(strcasecmp(sim_param->initial_conditions,"random_hemis") == 0 )
          {
      
            random_hemis_ic( sim_param, Qij, geometry );
            
          }
        else if(strcasecmp(sim_param->initial_conditions,"random") == 0 )
          {
      
            random_ic( sim_param, Qij, geometry );
            
          }
        else if(strcasecmp(sim_param->initial_conditions,"cholesteric") == 0 )
          {
      
            cholesteric_ic( sim_param, Qij, geometry );
            
          }
        else
          {
        
            std::cout << "\n The program did not recognize the initial condition option \"" << sim_param->initial_conditions << "\".\nPlease review your input file.\n\nAborting the program.\n";
      
            exit(2);
          }
        break;
        
      case parameter_status::unset:
       std::cout<< "Parameter \"initial_conditions\" not set in your in put file.\n"<<
                   "Please, set the parameter for one of the available initial consitions in this geometry:\n"<<
                   "random,read_from_file\n\nAborting the program.\n\n";
        exit(1);
    }

}




void random_ic( struct Simulation_Parameters * sim_param,double * Qij, const GEOMETRY * geometry )
{
  int i,j,k;
  double n[3];  

  const double S_eq=sim_param->S_eq;
  const int Nx=geometry->Nx;
  const int Ny=geometry->Ny;
  const int Nz=geometry->Nz;
  const int * point_type=geometry->point_type;
  
  
  if (sim_param->ic_flag[4]==parameter_status::unset) 
    gsl_rng_default_seed=time(NULL);
  else gsl_rng_default_seed= sim_param->rng_seed;
//  gsl_rng_default_seed=1570800053;
  
  gsl_rng *w= gsl_rng_alloc(gsl_rng_taus);
  std::cout << "seed = " <<gsl_rng_default_seed;


  for(i= 0; i< Nx; i++)
    {
      for(j= 0; j< Ny; j++)
        {
          for(k= 0; k< Nz; k++)
            {
              if(point_type[(k*Ny+j)*Nx+i] !=0 )
                {

    
                  gsl_ran_dir_3d(w, &n[0], &n[1], &n[2]);
                  Qij[5*(Nx*(Ny*k+j)+i)+0]=0.1*(0.5*S_eq*(3.0*n[0]*n[0]-1.0));
                  Qij[5*(Nx*(Ny*k+j)+i)+1]=0.1*(0.5*S_eq*(3.0*n[0]*n[1]));
                  Qij[5*(Nx*(Ny*k+j)+i)+2]=0.1*(0.5*S_eq*(3.0*n[0]*n[2]));
                  Qij[5*(Nx*(Ny*k+j)+i)+3]=0.1*(0.5*S_eq*(3.0*n[1]*n[1]-1.0));
                  Qij[5*(Nx*(Ny*k+j)+i)+4]=0.1*(0.5*S_eq*(3.0*n[1]*n[2]));                


                      
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
void random_hemis_ic( struct Simulation_Parameters * sim_param,double * Qij, const GEOMETRY * geometry )
{
  int i,j,k;
  double n[3];  

  const double S_eq=sim_param->S_eq;
  const int Nx=geometry->Nx;
  const int Ny=geometry->Ny;
  const int Nz=geometry->Nz;
  const int * point_type=geometry->point_type;
  
  
  if (sim_param->ic_flag[4]==parameter_status::unset) 
    gsl_rng_default_seed=time(NULL);
  else gsl_rng_default_seed= sim_param->rng_seed;
//  gsl_rng_default_seed=1570800053;
  
  gsl_rng *w= gsl_rng_alloc(gsl_rng_taus);
  std::cout << "seed = " <<gsl_rng_default_seed;


  for(i= 0; i< Nx; i++)
    {
      for(j= 0; j< Ny; j++)
        {
          for(k= 0; k< Nz; k++)
            {
              if(point_type[(k*Ny+j)*Nx+i] ==1 )
                {

    
                  gsl_ran_dir_3d(w, &n[0], &n[1], &n[2]);
                  Qij[5*(Nx*(Ny*k+j)+i)+0]=0.1*(0.5*S_eq*(3.0*n[0]*n[0]-1.0));
                  Qij[5*(Nx*(Ny*k+j)+i)+1]=0.1*(0.5*S_eq*(3.0*n[0]*n[1]));
                  Qij[5*(Nx*(Ny*k+j)+i)+2]=0.1*(0.5*S_eq*(3.0*n[0]*n[2]));
                  Qij[5*(Nx*(Ny*k+j)+i)+3]=0.1*(0.5*S_eq*(3.0*n[1]*n[1]-1.0));
                  Qij[5*(Nx*(Ny*k+j)+i)+4]=0.1*(0.5*S_eq*(3.0*n[1]*n[2]));                


                      
                }
              else if(point_type[(k*Ny+j)*Nx+i] ==2 )
                {

    
                  gsl_ran_dir_2d(w, &n[0], &n[1]); n[2]=0;
                  Qij[5*(Nx*(Ny*k+j)+i)+0]=0.1*(0.5*S_eq*(3.0*n[0]*n[0]-1.0));
                  Qij[5*(Nx*(Ny*k+j)+i)+1]=0.1*(0.5*S_eq*(3.0*n[0]*n[1]));
                  Qij[5*(Nx*(Ny*k+j)+i)+2]=0.1*(0.5*S_eq*(3.0*n[0]*n[2]));
                  Qij[5*(Nx*(Ny*k+j)+i)+3]=0.1*(0.5*S_eq*(3.0*n[1]*n[1]-1.0));
                  Qij[5*(Nx*(Ny*k+j)+i)+4]=0.1*(0.5*S_eq*(3.0*n[1]*n[2]));                
                      
                }
              else if(point_type[(k*Ny+j)*Nx+i] ==3 )
                {
                  
                  double delta_x=(i-(Nx-1)/2);
                  double delta_y=(j-(Ny-1)/2);
                  double delta_z=(k);
                  double rr=sqrt( delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
                  
                  n[0]=delta_x/rr;
                  n[1]=delta_y/rr;
                  n[2]=delta_z/rr;
                  
                  Qij[5*(Nx*(Ny*k+j)+i)+0]=0.1*(0.5*S_eq*(3.0*n[0]*n[0]-1.0));
                  Qij[5*(Nx*(Ny*k+j)+i)+1]=0.1*(0.5*S_eq*(3.0*n[0]*n[1]));
                  Qij[5*(Nx*(Ny*k+j)+i)+2]=0.1*(0.5*S_eq*(3.0*n[0]*n[2]));
                  Qij[5*(Nx*(Ny*k+j)+i)+3]=0.1*(0.5*S_eq*(3.0*n[1]*n[1]-1.0));
                  Qij[5*(Nx*(Ny*k+j)+i)+4]=0.1*(0.5*S_eq*(3.0*n[1]*n[2]));                
                      
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
////build a homogeneous director field
void homogeneous_ic( struct Simulation_Parameters * sim_param,double * Qij, const GEOMETRY * geometry )
{
  int i,j,k;
  double n[3];  

  const double S_eq=sim_param->S_eq;
  const int Nx=geometry->Nx;
  const int Ny=geometry->Ny;
  const int Nz=geometry->Nz;
  const int * point_type=geometry->point_type;

  const double theta_i=sim_param->theta_i*M_PI/180;
  const double phi_i=sim_param->phi_i*M_PI/180;
  
  
  std::cout<<"Initiating homogeneous initial conditions:\n\n";
  switch(sim_param->ic_flag[2])
    {
    case parameter_status::set:
      

      std::cout << "theta_i=" << theta_i << "\n";
      break;

    case parameter_status::unset:
      std::cout<<"Parameter \"theta_i\" not set.\nThe initial condition named \"homogeneous\" needs the paramters \"theta_i\" and \"phi_i\" set for use.\n"<<
      "Please set them in your in your input file.\n \nAborting the program.\n";
     
      exit(1);
      break;
    }

  
  switch(sim_param->ic_flag[3])
    {
    case parameter_status::set:
      
      std::cout << "phi_i=" << phi_i << "\n";
      break;

    case parameter_status::unset:

      std::cout<<"Parameter \"phi_i\" not set.\nThe initial condition named \"homogeneous\" needs the paramters \"theta_i\" and \"phi_i\" set for use.\n"<<
      "Please set them in your in your input file.\n \nAborting the program.\n";
      exit(1);
      break;
    }
        

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


void cholesteric_ic( struct Simulation_Parameters * sim_param,double * Qij, const GEOMETRY * geometry )
{
  //Initiatiate a homogeneously twisted cholesteric along the z sxis
  int i,j,k;
  double n[3];  

  const double S_eq=sim_param->S_eq;
  const int Nx=geometry->Nx;
  const int Ny=geometry->Ny;
  const int Nz=geometry->Nz;
  const int * point_type=geometry->point_type;

  const double theta_i=sim_param->theta_i*M_PI/180;
  const double phi_i=sim_param->phi_i*M_PI/180;
  const double p0_i=sim_param->p0_i=p0_i/(M_PI*sim_param->dz);;
  
  std::cout<<"Initiating cholesteric initial conditions:\n\n";
  switch(sim_param->ic_flag[2])
    {
    case parameter_status::set:
      
  
      std::cout << "theta_i=" << theta_i << "\n";
      break;

    case parameter_status::unset:
    
      std::cout<<"Parameter \"theta_i\" not set.\nThe initial condition named \"cholesteric\" needs the paramters \"theta_i\", \"phi_i\" and \"p0_i\" set for use.\n"<<
      "Please set them in your in your input file.\n \nAborting the program.\n";
      exit(1);
      break;
    }

  
  switch(sim_param->ic_flag[3])
    {
    case parameter_status::set:
      
  
      std::cout << "phi_i=" << phi_i << "\n";
      break;

    case parameter_status::unset:

      std::cout<<"Parameter \"phi_i\" not set.\nThe initial condition named \"cholesteric\" needs the paramters \"theta_i\", \"phi_i\" and \"p0_i\" set for use.\n"<<
      "Please set them in your in your input file.\n \nAborting the program.\n";
      exit(1);
      break;
    }
    
  switch(sim_param->ic_flag[5])
    {
    case parameter_status::set:
      
  
      std::cout << "p0_i=" << p0_i*sim_param->dz << "\n";
      break;

    case parameter_status::unset:

      std::cout<<"Parameter \"p0_i\" not set.\nThe initial condition named \"cholesteric\" needs the paramters \"theta_i\", \"phi_i\" and \"p0_i\" set for use.\n"<<
      "Please set them in your in your input file.\n \nAborting the program.\n";
      exit(1);
      break;
    }

      

  
  for(i= 0; i< Nx; i++)
    {
      for(j= 0; j< Ny; j++)
        {      
          for(k= 0; k< Nz; k++)
            {
             
              if(point_type[(k*Ny+j)*Nx+i] !=0 )
                {
                  n[0]=sin(theta_i)*cos(phi_i+k/p0_i);
                  n[1]=sin(theta_i)*sin(phi_i+k/p0_i);
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

void homogeneous_easy_axis_ic( struct Simulation_Parameters * sim_param,double * Qij, const GEOMETRY * geometry )
{
  homogeneous_ic(sim_param,Qij, geometry);
  homogeneous_boundary(sim_param,Qij, geometry);
  
}
//Initiate a homogeneous_easy_axis at strong boundaries

void homogeneous_boundary( struct Simulation_Parameters * sim_param, double * Qij, const GEOMETRY * geometry )
{
  
  std::cout<< "Starting homogeneous easy axis at the boundaries.\n";
  
  int i,j,k;
  double n[3];
  const int number_of_boundaries=geometry->number_of_boundaries;
  std::vector<double> theta_0(number_of_boundaries);      
  std::vector<double> phi_0(number_of_boundaries);    


  const double S_eq=sim_param->S_eq;
  const int Nx=geometry->Nx;
  const int Ny=geometry->Ny;
  const int Nz=geometry->Nz;
  const int * point_type=geometry->point_type;

  
  for(int ii=0; ii<number_of_boundaries; ii++)
    {
  
      try
        {
      
          theta_0[ii]=sim_param->theta_0.at(ii);
          phi_0[ii]=sim_param->phi_0.at(ii);
     
        }
      catch(std::out_of_range dummy_var)
        {

          std::cout << "\n Easy axis angle (theta_0 or phi_0) number " << ii <<
          " not defined.\n This initial condition needs both of them defined for use.\n" <<
      "Please define them in your in your input file.\n \n Aborting the program.\n";
      
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

//Build a random director field with homogeneous border
void random_bulk_homogeneous_easy_axis_ic( struct Simulation_Parameters * sim_param,double * Qij, const GEOMETRY * geometry )
{

  std::cout <<"Starting composite initial conditions: Random bulk, homogeneous on boundary.\n\n";
  random_ic( sim_param, Qij, geometry );
  homogeneous_boundary( sim_param, Qij, geometry );
  
}

void read_from_file_ic( struct Simulation_Parameters * sim_param, double * Qij, const GEOMETRY * geometry )
{
  int i,j,k,ii,jj,kk;
  double n[3],l[3],m[3],S,P;    
  FILE * ic_file;
  char string_placeholder[400];
  int read_status;
  int reading_line=1;

  const int Nx=geometry->Nx;
  const int Ny=geometry->Ny;
  const int Nz=geometry->Nz;
  const int * point_type=geometry->point_type;



  
  ic_file=fopen(sim_param->ic_file_name,"r");
  if (ic_file== NULL)
    {
      perror(sim_param->ic_file_name);
      exit(EXIT_FAILURE);
    }

  //get the file header:
  
  std::cout << "\nReading initial conditions from \"" << sim_param->ic_file_name << "\".\n";
  
  
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
 
              if(ii>=Nx|jj>=Ny|kk>=Nz)
                {
                  printf("ERROR: in line %d of %s.\n",reading_line+1,sim_param->ic_file_name);
                  printf("You are trying to assing the lattice (%d,%d,%d) in a cell of dimensions (%d,%d,%d).\n",ii,jj,kk,Nx,Ny,Nz);
                  printf("Note that the cell count starts at lattice (0,0,0)\n");
                  printf("Please, check your input file and your ic file.\n");
                  exit(4);
                }
          

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
      std::cout <<"Warning: The initial condition reader read "<< reading_line<< " lines, while the it was expected to read "<<
      expected_lines<<" lines.\n This could imply inconsistence in the size of your actyual system and size of the initial condition system.\n\n";
    }
    fclose(ic_file);
  
}


//Auxiliary routines:
void read_check(int read_status, int line)
{

  if(read_status != 11)
    {
      std::cout <<"Failed to read the initial condtion file line number "<<line<<" between field numbers "<<
      read_status<<".\n Aborting the program.\n\n";
      exit(2);
    }

}



