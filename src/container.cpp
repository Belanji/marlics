#include "container.h"
#include "driver.h"
#include <cstdio>
#include <cmath>
#include <cstdlib> 
#include <ccomplex>
#include <string>
#include <sys/stat.h>
#include <gsl/gsl_eigen.h>
//#include "geometry.h"
#define FILE_NUMBER_DIGITS 5

using namespace std;

//Manage the program output
container::container(struct Simulation_Parameters * sim_param)
{
  struct stat buffer;  
  sprintf(output_folder,"%s", sim_param->output_folder);
  output_fname.assign(sim_param->output_fname);

  if (stat(output_folder,&buffer)!=0)
  {
      std::cout<<output_folder<<": directory not found"<<endl;
      perror(output_folder);
      exit(4);
  }  
    switch(sim_param->output_name_status[0])
      {
      case   parameter_status::unset:
        
        std::cout << "You didn't define the \"output_folder\". Using the default output folder -> \".\" (here)."<< std::endl;
        break;  
      }
    switch(sim_param->output_name_status[1])
      {
      case   parameter_status::unset:
        
        std::cout << "You didn't define the \"output_name\". Using the standart output name -> director_field."<< std::endl;
        break; 
      }
  std::cout << "Using \""<< output_folder<<"\\"<<output_fname<<" as output file pattern.\n";
  switch(sim_param->timeprint_status[3])
    {
      case parameter_status::set:
        {
          char first_name[100];
          char fn[3];
          sprintf(fn,"%d",sim_param->first_output_file_number);
          file_number=sim_param->first_output_file_number;
          std::string output_name=output_fname;
          std::size_t found = output_name.rfind("$$");
          if (found!=std::string::npos)
            output_name.replace (found,2,fn);
            
          sprintf(first_name,"%s/%s", output_folder, output_name.c_str());
          cout << "First snapshot will written at \"" << first_name <<"\"."<<std::endl;
          
          Nx=sim_param->Nx;
          Ny=sim_param->Ny;
          Nz=sim_param->Nz;
          break;
        }
      case parameter_status::unset:
        fprintf(stderr,"The parameter first_output_file_number was not found in your input file/n Please, review your input file/n");
        file_number=0;
        break;
    }
}

// Evaluate the n vector, along with the S and P parameters, and print in a ".csv" file
void container::write_state(double t , const double   * Qij,const int * pt)
{

  FILE * snapshot;
  char time_stamp[100];
  char fn[FILE_NUMBER_DIGITS];
  int i,j,k;
  double d0, d1, n[3], l[3], data[9];
  
  sprintf(fn, "%d", file_number);
  std::string output_name=output_fname;
  std::size_t found = output_name.rfind("$$");
  if (found!=std::string::npos)
    output_name.replace (found,2,fn);
    
  sprintf(time_stamp,"%s/%s", output_folder, output_name.c_str());



  
  snapshot=fopen(time_stamp,"w");

  

  fprintf(snapshot,"x,y,z,nx,ny,nz,lx,ly,lz,S,P,pt\n");


  

  for(k= 0; k< Nz; k++)
    {
        
      for(j= 0; j< Ny; j++)
        {

                  
          for(i= 0; i< Nx; i++)
            {

              data[0]= Qij[5*(Nx*(Ny*k+j)+i)+0];
              data[1]= Qij[5*(Nx*(Ny*k+j)+i)+1];
              data[2]= Qij[5*(Nx*(Ny*k+j)+i)+2];
              data[3]= Qij[5*(Nx*(Ny*k+j)+i)+1];
              data[4]= Qij[5*(Nx*(Ny*k+j)+i)+3];
              data[5]= Qij[5*(Nx*(Ny*k+j)+i)+4];
              data[6]= Qij[5*(Nx*(Ny*k+j)+i)+2];
              data[7]= Qij[5*(Nx*(Ny*k+j)+i)+4];
              data[8]= (-Qij[5*(Nx*(Ny*k+j)+i)+0]-Qij[5*(Nx*(Ny*k+j)+i)+3]);
                
              gsl_matrix_view m= gsl_matrix_view_array(data, 3, 3);
              gsl_vector *eval= gsl_vector_alloc(3);
              gsl_matrix *evec= gsl_matrix_alloc(3, 3);
              gsl_eigen_symmv_workspace * w= gsl_eigen_symmv_alloc(3);
              gsl_eigen_symmv(&m.matrix, eval, evec, w);
              gsl_eigen_symmv_free(w);
              gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
                
              d0= gsl_vector_get(eval, 0);
              d1= gsl_vector_get(eval, 1);
                
              n[0]= gsl_matrix_get(evec, 0, 0);
              n[1]= gsl_matrix_get(evec, 1, 0);
              n[2]= gsl_matrix_get(evec, 2, 0);
              
              l[0]= gsl_matrix_get(evec, 0, 1);
              l[1]= gsl_matrix_get(evec, 1, 1);
              l[2]= gsl_matrix_get(evec, 2, 1);
                
              gsl_vector_free(eval);
              gsl_matrix_free(evec);
              
              fprintf(snapshot,"%d,%d,%d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%d\n", i,  j,  k, n[0], n[1], n[2], l[0], l[1], l[2],d0, 2.0*d1+d0,pt[(Nx*(Ny*k+j)+i)]);

                
            }     
        }
    }

  fclose(snapshot);
  std::cout<<"Snapshot "<<file_number<<": "<< t << std::endl;
  file_number++;
}


      

     


