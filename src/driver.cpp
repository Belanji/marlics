#include "driver.h"
#include "container.h"
#include "geometry.h"
#include "geometry_slab.h"
#include <cstdio>           
#include <cstdlib>          
#include <cstring>
#include <cmath>            
#define MIN(a,b) ((a) < (b) ? a : b)



driver::driver(void)
{

     
}


void driver::setup_LC(void)
{

  //Setting the chirality parameter:
  switch(sim_param.chirality_flag)
    {
    case chirality_status::unset:

      printf("\nNo chirality set using q0=0 (nematic liquid crystal.\n\n");
      break;

    case chirality_status::p0:

      sim_param.q0=2.0*M_PI/sim_param.p0;
      break;

    }


  //Setting up the thermal LDG parameters:
  if(sim_param.thermal_flag[0]==parameter_status::set &&
     sim_param.thermal_flag[1]==parameter_status::set &&
     sim_param.thermal_flag[2]==parameter_status::set &&
     sim_param.thermal_flag[3]==parameter_status::set)
    {
      
      sim_param.S_eq=((-sim_param.B+sqrt(sim_param.B*sim_param.B-24.0*sim_param.a*sim_param.T*sim_param.C)))/(6.0*sim_param.C);
      
    }
  else
    {

      printf("You must define all the four thermal energy LDG parameters: temperaure, a,B and C (also C can not be 0).\n");
      printf("Please, review your input file");
      printf("Aborting the program.\n");
      exit(0);
           
    };
      
  
  //Setting up the elastic LDG parameters:  
  switch(sim_param.elastic_flag)
    {
    case elastic_const_status::Ls:

      break;
	
    case elastic_const_status::Ks:

      printf("Converting elastic constants Kii to Li.\n");

      double S_eq=sim_param.S_eq;
      double k11=sim_param.k11;
      double k22=sim_param.k22;
      double k33=sim_param.k33;
      double k24=sim_param.k24;

  
      sim_param.L1=2.0*(k33-k11+3.0*k22)/(27.0*S_eq*S_eq);
      sim_param.L2=4.0*(k11-k22-k24)/(9.0*S_eq*S_eq);
      sim_param.L3=4.0*(k33-k11)/(27.0*S_eq*S_eq*S_eq);
      sim_param.Lq=2.0*(k22)/(9.0*S_eq*S_eq);
      sim_param.Ls=4*(k24)/(9.0*S_eq*S_eq);
      break;

    case elastic_const_status::unset:

      printf("You must define at least one elastic constant K/L in your input file.\n");
      printf("Please, define one or check your input file to see if they are mispelled.\n");
      printf("Aborting the program");
      exit(0);
      break;
	
    }


  switch(sim_param.viscosity_flag[0])
    {
    case viscosity_status::unset:
      printf("You must define the liquid cristal viscosity (gamma_1 or mu_1) in your input file.\n");
      printf("Please, define one or check your input file to see if they are mispelled.\n");
      printf("Aborting the program");
      exit(0);
      break;
	
    case viscosity_status::gamma:

      sim_param.mu_1=sim_param.gamma_1/sim_param.S_eq;
      break;

    }

  switch(sim_param.viscosity_flag[1])
    {
    case viscosity_status::unset:
      printf("You must define the liquid cristal surface viscosity (gamma_1_s or mu_1_s) in your input file.\n");
      printf("Please, define one or check your input file to see if they are mispelled.\n");
      printf("Aborting the program");
      exit(0);
      break;
	
    case viscosity_status::gamma:

	
      sim_param.mu_1_s=sim_param.gamma_1_s/sim_param.S_eq;

      break;

    }

    
  //Checking and Setting up the grid:
  if(sim_param.grid_spacing_flag[0]==parameter_status::unset ||
     sim_param.grid_spacing_flag[1]==parameter_status::unset ||
     sim_param.grid_spacing_flag[2]==parameter_status::unset )

    {
	
      printf("You must define the grid spacing parameters {dx,dy,dz} in your input file.\n");
      printf("Please, define them or check your input file to see if they are mispelled.\n");
      printf("Aborting the program");
      exit(0);

    }
    
  if(sim_param.grid_flag[0]==parameter_status::set &&
     sim_param.grid_flag[1]==parameter_status::set &&
     sim_param.grid_flag[2]==parameter_status::set )
    {
      if((sim_param.Q_00= (double *)calloc(5*sim_param.Nx*sim_param.Ny*sim_param.Nz, sizeof(double)))==NULL){ERROr};  
      if((Qij= (double *)calloc(5*sim_param.Nx*sim_param.Ny*sim_param.Nz, sizeof(double)))==NULL){ERROr};
    }
  else
    {

      printf("You must define all three grid sizes {Nx,Ny,Nz} in your input file.\n");
      printf("Please, define them or check your input file to see if you mispelled one of them.\n");
      printf("Aborting the program");
	
    }

  

      
}
  

void driver::setup_Simulation(void)
{


  switch(sim_param.time_status[0])
    {
    case parameter_status::set:

      printf("ti=%lf\n", sim_param.ti);
      break;

    case parameter_status::unset:


      sim_param.ti=0;

      printf("Simulation parameter \"ti\" not set.");
      printf("Using standard value ti=0.\n");
      printf("ti=%lf\n",sim_param.ti);
      break;
    }

      
    switch(sim_param.time_status[1])
    {
    case parameter_status::set:

      printf("tf=%lf\n", sim_param.tf);
      break;

    case parameter_status::unset:

      printf("Simulation parameter \"tf\" not set.");
      printf("You mut define a value for \"tf\".\n");
      printf("Aborting the program.\n");
      exit(0);
      break;
    }



    switch(sim_param.time_status[2])
    {
    case parameter_status::set:

      printf("dt=%lf\n\n", sim_param.dt);
      break;

    case parameter_status::unset:


      sim_param.ti=0;

      printf("Simulation parameter \"dt\" not set.");
      printf("Using standard value dt=tf/1e6.\n");
      printf("dt=%lf\n\n",sim_param.tf/1e6);
      break;
    }



    
  switch (sim_param.timeprint_status[0])
    {
	
    case parameter_status::unset:

      next_time_print=linear_next_time_print;
      printf("No time printing time chosen.\n");
      printf("Using the standard type:linear.\n");
      break;

    case parameter_status::set:

      if ( strcasecmp(sim_param.print_time_type,"linear") == 0 )
	{
	  next_time_print=linear_next_time_print;
	  printf("Snapshots frequency type: %s\n",sim_param.print_time_type);
	  
	}
      else if ( strcasecmp(sim_param.print_time_type,"logarithmic") == 0 )
	{

	  printf("Snapshots frequency type: %s\n",sim_param.print_time_type);
	  next_time_print=log_next_time_print;
      
	}
      else
	{

	  printf("No time printing type named %s.\n",sim_param.print_time_type);
	  printf("Aborting the program.\n");
	  exit(0);
	};      
      break;
	  
    }

    

  
  switch (sim_param.timeprint_status[1])
    {
	
    case parameter_status::unset:

      
      printf("No time printing time chosen.\n");
      printf("Using the standard timeprint parameter:tf/20.\n");
      sim_param.timeprint=sim_param.tf/20;
      printf("First snapshot taken at:  %lf\n", sim_param.timeprint);
      break;
      
      case parameter_status::set:
	


	printf("First snapshot taken at:  %lf\n", sim_param.timeprint);
	break;
    }



  switch (sim_param.timeprint_status[2])
    {
	
    case parameter_status::unset:

      
      printf("No time printing frequency chosen.\n");
      printf("Using the standard time_print_frequency parameter:tf/20.\n");
      sim_param.timeprint_increase_factor=sim_param.tf/20;
      printf("Snapshot print frequency: %lf\n\n",sim_param.timeprint_increase_factor);

      break;
      
      case parameter_status::set:


	printf("Snapshot print frequency: %lf\n\n",sim_param.timeprint_increase_factor);

	break;
    }



  Data_Container= new container(& sim_param);

  
  switch(sim_param.geometry_flag)
    {
    case parameter_status::set:
      
      if ( strcasecmp(sim_param.geometry,"slab") == 0 )
	{

	  
	  LcS_Geometry= new slab( & sim_param);
  
	}
  
      else 
	{
	
	  printf("No geometry named %s is defined.\nAborting the program.\n\n",sim_param.geometry);
	  exit(0);
  
	};

      break;

    case parameter_status::unset:

      printf("Simulation geometry not set.");
      printf("You mut define a geometry for your simulation.\n");
      printf("Aborting the program.\n");
      exit(0);
      break;

    }

      
    LcS_Geometry->ic( & sim_param, Qij );
  //  LcS_Geometry->boundary_init( & sim_param );  


  //  if ( strcasecmp(sim_param.integrator_type,"DP5") == 0 || strcasecmp(sim_param.integrator_type,"Dormand-Prince") == 0
  //       || strcasecmp(sim_param.integrator_type,"Rk54") == 0 || strcasecmp(sim_param.integrator_type,"Rk5") == 0 )
  //    {
  //
  //      LcS_Integrator= new DP5(LcS_Geometry, & sim_param );
  //
  //    }
  //
  //  else
  //    {
  //      
  //      printf("No integrator named %s is defined.\nAborting the program.\n\n",sim_param.integrator_type);
  //      exit(0);
  //
  //      
  //    };
  



  
}
void driver::error_check(int error_handler,
		  char parser[])

{

  	  if (error_handler <= 0 )
	    {
	    printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	  printf("Please review your input file.\n Aborting the program\n");
	  exit(0);
	    }
	  
}

void driver::update_timeprint(void)
{


next_time_print(& sim_param.timeprint, sim_param.timeprint_increase_factor);
sim_param.timeprint=MIN(sim_param.timeprint,sim_param.tf);



};


void log_next_time_print(double *time_print  , double factor )
{

  *time_print*=factor;

};

void linear_next_time_print(double *time_print , double factor )
{

  *time_print+=factor;
  
};




int driver::parse_input_file(void)
{

  char parser[80];
  char garbage[400];
  int error_handler;

 
  while (   scanf("%200s",parser) != EOF )
    {

      if ( strcasecmp(parser,"Geometry") == 0 )
	{
	  sim_param.geometry_flag=parameter_status::set;
	  error_handler=scanf("%200s",&sim_param.geometry);
	  error_check(error_handler,parser);
		
	  printf("Geometry Used:  %s\n\n", sim_param.geometry);		  
	  fgets(garbage,400,stdin);

	}            
      else if ( strcasecmp(parser,"initial_conditions") == 0 )
	{
	  sim_param.ic_flag[0]=parameter_status::set;
	  error_handler=scanf("%200s",&sim_param.initial_conditions);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);
	  

	}
      else if ( strcasecmp(parser,"ic_file") == 0 || strcasecmp(parser,"initial_conditions_file") == 0 || strcasecmp(parser,"input_initial_conditions") == 0 || strcasecmp(parser,"input_initial_conditions_file") == 0)
	{
	  sim_param.ic_flag[1]=parameter_status::set;
	  error_handler=scanf("%200s",&sim_param.ic_file_name);

	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"integrator") == 0 || strcasecmp(parser,"integrator_type") == 0 )
	{

	  error_handler=scanf("%200s",&sim_param.integrator_type);
	  
	  error_check(error_handler,parser);
		

	  fgets(garbage,400,stdin);

	  
	}
      else if ( strcasecmp(parser,"time_print_type") == 0 || strcasecmp(parser,"snapshot_print_frequency_type") == 0 || strcasecmp(parser,"print_time_type") == 0 || strcasecmp(parser,"snap_print_freq_type") == 0)
	{
	  sim_param.timeprint_status[0]=parameter_status::set;
	  error_handler=scanf("%200s", &sim_param.print_time_type);
	  
	  error_check(error_handler,parser);
		

	  fgets(garbage,400,stdin);


	  
	}
      else if ( strcasecmp(parser,"timeprint") == 0 || strcasecmp(parser,"first_snapshot") ==0 )
	{

	  sim_param.timeprint_status[1]=parameter_status::set;
	  error_handler=scanf("%lf", &sim_param.timeprint);	  
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}
         else if ( strcasecmp(parser,"timeprint_increase_factor") == 0 )
	{

	  sim_param.timeprint_status[2]=parameter_status::set;
	  error_handler=scanf("%lf", &sim_param.timeprint_increase_factor);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}
	 else if ( strcasecmp(parser,"initial_output_file_number") == 0 ||
		strcasecmp(parser,"first_output_file_number") == 0 )
	{

	  sim_param.timeprint_status[3]=parameter_status::set;
	  error_handler=scanf("%d",&(sim_param.firt_output_file_number));
	  error_check(error_handler,parser);
		
	  
	  fgets(garbage,400,stdin);


	}
	 else if ( strcasecmp(parser,"L1") == 0 )
	{
	  sim_param.elastic_flag=elastic_const_status::Ls;
	  error_handler=scanf("%lf",&sim_param.L1);

	  error_check(error_handler,parser);
		

	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"L2")== 0 )
	{
	  
	  sim_param.elastic_flag=elastic_const_status::Ls;
	  error_handler=scanf("%lf",&sim_param.L2);
	  error_check(error_handler,parser);

	  fgets(garbage,400,stdin);
	  
	  
      
	}
      else if ( strcasecmp(parser,"L3") == 0)
	{

	  sim_param.elastic_flag=elastic_const_status::Ls;
	  error_handler=scanf("%lf",&sim_param.L3);
	  error_check(error_handler,parser);


	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"Lq")== 0 || strcasecmp(parser,"L_q")== 0 )
	{
	  sim_param.elastic_flag=elastic_const_status::Ls;
	  error_handler=scanf("%lf",&sim_param.Lq);

	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);
	   
      
	}
      else if ( strcasecmp(parser,"Ls")== 0 || strcasecmp(parser,"L_s")== 0 )
	{
	  sim_param.elastic_flag=elastic_const_status::Ls;
	  error_handler=scanf("%lf",&sim_param.Ls);
	  error_check(error_handler,parser);

	  fgets(garbage,400,stdin);
	}
      else if ( strcasecmp(parser,"k11") == 0  )
	{

	  sim_param.elastic_flag=elastic_const_status::Ks;
	  error_handler=scanf("%lf",&sim_param.k11);

	  error_check(error_handler,parser);
		

	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"k22")== 0 )
	{
	  sim_param.elastic_flag=elastic_const_status::Ks;	  
	  error_handler=scanf("%lf",&sim_param.k22);


	  error_check(error_handler,parser);


	  fgets(garbage,400,stdin);
	  
	  
      
	}
      else if ( strcasecmp(parser,"k33") == 0 )
	{

	  sim_param.elastic_flag=elastic_const_status::Ks;
	  error_handler=scanf("%lf",&sim_param.k33);
	  error_check(error_handler,parser);

	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"k24") == 0 )
	{

	  sim_param.elastic_flag=elastic_const_status::Ks;
	  error_handler=scanf("%lf",&sim_param.k24);
	  error_check(error_handler,parser);

	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"a") == 0 )
	{

	  sim_param.thermal_flag[0]=parameter_status::set;
	  error_handler=scanf("%lf",&sim_param.a);

	  error_check(error_handler,parser);
		
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"B")== 0 )
	{
	  sim_param.thermal_flag[1]=parameter_status::set;
	  error_handler=scanf("%lf",&sim_param.B);


	  error_check(error_handler,parser);		  	  
	  fgets(garbage,400,stdin);
	  
	  
      
	}
      else if ( strcasecmp(parser,"C") == 0 )
	{

	  sim_param.thermal_flag[2]=parameter_status::set;
	  error_handler=scanf("%lf",&sim_param.C);

	  error_check(error_handler,parser);
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"T")== 0 || strcasecmp(parser,"temperature")== 0 )
	{

	  sim_param.thermal_flag[3]=parameter_status::set;;
	  error_handler=scanf("%lf",&sim_param.T);
	  error_check(error_handler,parser);		  	  
	  fgets(garbage,400,stdin);
	  
	  
      
	}
      else if ( strcasecmp(parser,"dx") == 0  )
	{

	  sim_param.grid_spacing_flag[0]=parameter_status::set;
	  error_handler=scanf("%lf",&sim_param.dx);

	  error_check(error_handler,parser);
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"dy") == 0  )
	{
	  sim_param.grid_spacing_flag[1]=parameter_status::set;
	  error_handler=scanf("%lf",&sim_param.dy);

	  error_check(error_handler,parser);
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"dz") == 0  )
	{

	  sim_param.grid_spacing_flag[2]=parameter_status::set;
	  error_handler=scanf("%lf",&sim_param.dz);

	  error_check(error_handler,parser);
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"p0") == 0 )
	{
	  sim_param.chirality_flag=chirality_status::p0;
	  error_handler=scanf("%lf",&sim_param.p0);
	  error_check(error_handler,parser);

	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"q0") == 0 )
	{
	  sim_param.chirality_flag=chirality_status::q0;

	  error_handler=scanf("%lf",&sim_param.q0);
	  error_check(error_handler,parser);

	  
	  fgets(garbage,400,stdin);


	}
      else if (strcasecmp(parser,"tf") == 0 || strcasecmp(parser,"run_time") ==0 )
	{

	  sim_param.time_status[1]=parameter_status::set;
	  error_handler=scanf("%lf",&sim_param.tf);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}
      else if (strcasecmp(parser,"ti") == 0 || strcasecmp(parser,"start_time") ==0 )
	{

	  sim_param.time_status[0]=parameter_status::set;
	  error_handler=scanf("%lf",&sim_param.ti);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}
      else if (strcasecmp(parser,"dt") == 0 || strcasecmp(parser,"timestep") ==0 )
	{

	  sim_param.time_status[2]=parameter_status::set;
	  error_handler=scanf("%lf",&sim_param.dt);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}      
      else if ( strcasecmp(parser,"bulk_viscosity") == 0  || strcasecmp(parser,"bviscosity") == 0 || strcasecmp(parser,"bulk_visc") == 0 || strcasecmp(parser,"gamma1") == 0 || strcasecmp(parser,"gamma_1") == 0 )
	{

	  sim_param.viscosity_flag[0]=viscosity_status::gamma;
	  error_handler=scanf( "%lf",&sim_param.gamma_1 );
	  error_check(error_handler,parser);	  
	  fgets(garbage,400,stdin);

	}
      else if ( strcasecmp(parser,"surface_viscosity") == 0  || strcasecmp(parser,"sviscosity") == 0 || strcasecmp(parser,"surface_visc") == 0 || strcasecmp(parser,"gamma_1_s") == 0 || strcasecmp(parser,"gamma_s") == 0 )
	{

	  sim_param.viscosity_flag[1]=viscosity_status::gamma;
	  error_handler=scanf( "%lf",&sim_param.gamma_1_s );
	  error_check(error_handler,parser);	  
	  fgets(garbage,400,stdin);
	  
	}
      else if (strcasecmp(parser,"mu_1") == 0  || strcasecmp(parser,"mu1") == 0 )
	{

	  sim_param.viscosity_flag[0]=viscosity_status::mu;
	  error_handler=scanf( "%lf",&sim_param.mu_1 );
	  error_check(error_handler,parser);	  
	  fgets(garbage,400,stdin);	  
	}
      else if (strcasecmp(parser,"mu_1_s") == 0  || strcasecmp(parser,"mu1_s") == 0 )
	{

	  sim_param.viscosity_flag[1]=viscosity_status::mu;
	  error_handler=scanf( "%lf",&sim_param.mu_1_s );
	  error_check(error_handler,parser);	  
	  fgets(garbage,400,stdin);
	  
	}
      else if ( strcasecmp(parser,"Nx") == 0)
	{
	  sim_param.grid_flag[0]=parameter_status::set;
	  error_handler=scanf("%i",&sim_param.Nx);
	  error_check(error_handler,parser);
		
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"Ny")== 0 )
	{
	  sim_param.grid_flag[1]=parameter_status::set;
	  error_handler=scanf("%i",&sim_param.Ny);

	  error_check(error_handler,parser);		  	  
	  fgets(garbage,400,stdin);
	  
	  
      
	}
      else if ( strcasecmp(parser,"Nz") == 0 )
	{
	  sim_param.grid_flag[2]=parameter_status::set;
	  error_handler=scanf("%i",&sim_param.Nz);
	  error_check(error_handler,parser);
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"Wo1_u") == 0 || strcasecmp(parser,"Upper_wo1")== 0 )
	{
	  error_handler=scanf("%lf",&sim_param.Wo1_u);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);
	}
            else if ( strcasecmp(parser,"Wo1_b") == 0 ||strcasecmp(parser,"lower_wo1")== 0 )
	{

	  error_handler=scanf("%lf",&sim_param.Wo1_b);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}      
      else if ( strcasecmp(parser,"phi_u") == 0 )
	{

	  error_handler=scanf("%lf",&sim_param.phi_u);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"theta_i") == 0 )
	{

	  error_handler=scanf("%lf",&sim_param.theta_i);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"phi_i") == 0 )
	{

	  error_handler=scanf("%lf",&sim_param.phi_i);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"theta_u") == 0 )
	{

	  error_handler=scanf("%lf",&sim_param.theta_u);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"upper_anchoring_type") == 0 || strcasecmp(parser,"top_anchoring_type") == 0 || strcasecmp(parser,"upper_anchoring") == 0 || strcasecmp(parser,"top_anchoring") == 0)
	{

	  error_handler=scanf("%200s",&sim_param.upper_anc_type);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"bottom_anchoring_type") == 0 || strcasecmp(parser,"lower_anchoring_type") == 0 || strcasecmp(parser,"bottom_anchoring") == 0 || strcasecmp(parser,"lower_anchoring") == 0)
	{

	  error_handler=scanf("%200s",&sim_param.bottom_anc_type);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"atol") == 0  || strcasecmp(parser,"Absolute_Tolerance") == 0 )
	{

	  error_handler=scanf("%lf",&sim_param.Atol);

	  error_check(error_handler,parser);
		
	  
	  fgets(garbage,400,stdin);


	}
            else if ( strcasecmp(parser,"rtol") == 0  || strcasecmp(parser,"Relative_Tolerance") == 0 )
	{

	  error_handler=scanf("%lf",&sim_param.Rtol);

	  error_check(error_handler,parser);
		
	  
	  fgets(garbage,400,stdin);


	}
            else if ( strcasecmp(parser,"facmax") == 0 || strcasecmp(parser,"maximum_timestep_increase") == 0  )
	{

	  error_handler=scanf("%lf",&sim_param.facmax);

	  error_check(error_handler,parser);
		
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"facmin") == 0 || strcasecmp(parser,"maximum_timestep_decrease") == 0  )
	{

	  error_handler=scanf("%lf",&sim_param.facmin);

	  error_check(error_handler,parser);
		
	  
	  fgets(garbage,400,stdin);


	}
       else if ( strcasecmp(parser,"prefac") == 0 || strcasecmp(parser,"timestep_control_pre_factor") == 0  )
	{

	  error_handler=scanf("%lf",&sim_param.prefac);

	  error_check(error_handler,parser);
		
	  
	  fgets(garbage,400,stdin);


	}
      else if (strcasecmp(parser,"run") == 0)
	{

	  return 1;

	}	  
      else if (parser[0]=='#')
	{

	  fgets(garbage,400,stdin);

	}	  
      else
	{

	  printf("The parser did not recognize the option '%s'. Please review your input file\n", parser);
	  printf("Aborting the program\n");
	  exit(0);
	};

    }
  return 1;
}


//sim_param.R_out=(sim_param.Nx-sim_param.Nx/2)+1.0;
  //sim_param.R_in=(sim_param.Nx-sim_param.Nx/2)-0.1;
  
