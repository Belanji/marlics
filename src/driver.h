#ifndef LC_DRIVER_
#define LC_DRIVER_

#define ERROr printf("cannot allocate memory\n"); exit(0);


  enum parameter_status {unset=0, set=1};

  struct Simulation_Parameters
  {

    //Domain Paramaters:
  
    int Nx;              /* grid size*/
    int Ny;              /* grid size*/
    int Nz;              /* grid size*/

    double dx;           /* 10^-9 m */
    double dy;           /* 10^-9 m */
    double dz;           /* 10^-9 m */
  
    char geometry[200];
    parameter_status geometry_flag=unset;
  
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
    int thermal_flag=0;
    int elastic_flag=0;
    int quirality_flag=0;
  
    double mu_1;         /* Pa s    */
    double mu_1_s;       /* Pa s    */
    double T=-1;            /* (T-T*)  */
    double gamma_1;
    double gamma_1_s;
    int viscosity_flag=0;
    int temperature_flag=0;

  
    //Anchoring paramters:
  
    double theta_u;      /* degrees */
    double phi_u;        /* degrees */
  
    double Wo1_u;         /* 10^-3 m */
    double Wo1_b;

    double theta_d;       /* degrees */
    double phi_d;         /* degrees */

    double theta_i;       /* degrees */
    double phi_i;         /* degrees */
 
    double Wo1_d;        /* 10^-3 m */     

    char upper_anc_type[200];
    char bottom_anc_type[200];  
    int anchoring_flag=0;
  

    //Initial Conditions:


    char ic_file_name[200];
    char initial_conditions[200];
    double disclination_number;
    int ic_file_flag=0;  


    //Integrator Paramters:

    char integrator_type[200];
    int integrator_flag=0;
    double Atol;
    double Rtol;
    double prefac;
    double facmin;
    double facmax;

  
    //Execution Parameters:
  

    char print_time_type[200];
    int firt_output_file_number;  
    int calls_for_director_output;

    double dt;           /*10^-6 s */
    double ti;           /*10^-6 s */
    double tf;           /*10^-6 s */
    double timeprint;    /*10^-6 s */
    double timeprint_increase_factor; /*linear=10^-6 s, logarithmic=admensional */

    double S_eq;
    double * Q_00;

  
  } ;




  class driver
  {

  public:

  
    struct Simulation_Parameters sim_param;
    class GEOMETRY * LcS_Geometry;
    class Integrator * LcS_Integrator;

    double * Qij;
 
    void update_timeprint (void);
    void (*next_time_print)  (double *, double);
    void setup_LC ( void);
    int parse_input_file(void);
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
