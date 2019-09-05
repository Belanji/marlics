#ifndef LC_DRIVER_
#define LC_DRIVER_

#include <iostream>
#define ERROr std::cout<<"cannot allocate memory\n"; exit(4);
#include <vector>
#include <string>
#include <map>

enum class parameter_status {unset=0, set=1};
enum class elastic_const_status {unset=0, Ls=1, Ks=2};
enum class chirality_status {unset=0, p0=1, q0=2};
enum class viscosity_status {unset=0, mu=1,gamma=2};

struct Simulation_Parameters
{

  //Domain Paramaters:
  
  int Nx;              /* grid size*/
  int Ny;              /* grid size*/
  int Nz;              /* grid size*/
  parameter_status grid_flag[3]={parameter_status::unset,
				    parameter_status::unset,
				    parameter_status::unset};
    
  double dx;           /* 10^-9 m */
  double dy;           /* 10^-9 m */
  double dz;           /* 10^-9 m */
  parameter_status grid_spacing_flag[3]={parameter_status::unset,
				    parameter_status::unset,
				    parameter_status::unset};

  
  char geometry[200];
  parameter_status geometry_flag=parameter_status::unset;
  
  /*Just for use in Sphere geometry */
  int R_in;            /* 10^-9 m */
  int R_out;           /* 10^-9 m */

  
  //Lc Parameters:

  double a=0;            /* 10^6 J/Km^3 */
  double B=0;            /* 10^6 J/m^3  */
  double C=0;            /* 10^6 J/m^3  */
  double L1=0;           /* 10^-12 N    */
  double L2=0;           /* 10^-12 N    */
  double L3=0;           /* 10^-12 N    */
  double Lq=0;           /* 10^-12 N    */
  double Ls=0;           /* 10^-12 N    */
  double k11=0, k22=0;      /* 10^-12 N    */
  double k33=0, k24=0;     /* 10^-12 N    */
  double p0=0;          /* 10^-9 m     */
  double q0=0;          /* 10^9 1/m    */
  double T=-1;            /* (T-T*)  */
  elastic_const_status elastic_flag=elastic_const_status::unset;
  chirality_status chirality_flag=chirality_status::unset;
  parameter_status thermal_flag[4]={parameter_status::unset,
				    parameter_status::unset,
				    parameter_status::unset,
				    parameter_status::unset};
  
  double mu_1=0;         /* Pa s    */
  double mu_1_s=0;       /* Pa s    */

  double gamma_1=0;
  double gamma_1_s=0;
  viscosity_status viscosity_flag[2]={viscosity_status::unset,
				      viscosity_status::unset};

  
  //Anchoring paramters:
  
  std::map<int,double> theta_0;      /* degrees */
  std::map<int,double> phi_0;        /* degrees */  
  std::map<int,double> Wo1;         /* 10^-3 m */
  std::map<int,double> Wo2;         /* 10^-3 m */

  std::map<int,std::string> anchoring_type;
  
  

  //Initial Conditions:

  
  char initial_conditions[200];
  char ic_file_name[200];
  double theta_i;       /* degrees */
  double phi_i;         /* degrees */ 
  parameter_status ic_flag[4]={parameter_status::unset,
			       parameter_status::unset,
			       parameter_status::unset,
			       parameter_status::unset};
  


  //Integrator Paramters:

  char integrator_type[200];  
  double Atol=0.001;
  double Rtol=0.001;
  double prefac=0.8;
  double facmin=0.4;
  double facmax=3.0;
  parameter_status integrator_flag=parameter_status::unset;
  int integrator_parameters_flag=0;
  
  //Output name parameters:
  
  char output_folder[200]=".";
  char output_fname[200]="director_field";
  parameter_status output_name_status[2]={parameter_status::unset, parameter_status::unset};
  
  //Execution Parameters:
  double timeprint;    /*10^-6 s */
  double timeprint_increase_factor; /*linear=10^-6 s, logarithmic=admensional */
  char print_time_type[200];
  int first_output_file_number;  
  parameter_status timeprint_status[4]={parameter_status::unset,
					parameter_status::unset,
					parameter_status::unset,
					parameter_status::unset};

  
  double ti;           /*10^-6 s */
  double tf;           /*10^-6 s */
  double dt;           /*10^-6 s */
  parameter_status time_status[3]={parameter_status::unset,
				   parameter_status::unset,
				   parameter_status::unset};
  
  double S_eq;
  double * Q_00;

  
} ;




  class driver
  {

  public:

  
    struct Simulation_Parameters sim_param;
    class GEOMETRY * LcS_Geometry;
    class Integrator * LcS_Integrator;
    class container * Data_Container;
    double * Qij;
 
    void update_timeprint (void);
    void (*next_time_print)  (double *, double);
    void setup_LC ( void);
    void setup_Simulation ( void);
    int parse_input_file(char[]);
    void error_check(int error_handler, char parser[]);

    driver(void);
    ~driver(void){};

  };

  void init_container(struct Simulation_Parameters * sim_param);
  void write_state(double , const double  *, const int *);
  void finish_container( void );


  void log_next_time_print(double *, double );
  void linear_next_time_print( double *, double );







#endif
