#ifndef initial_condtions__
#define initial_condtions__


 //Initial conditions routines:
void apply_initial_conditions(struct Simulation_Parameters *, double *, const class GEOMETRY *);
void read_from_file_ic(  struct Simulation_Parameters * lc_droplet,double * Qij , const class GEOMETRY *);
void random_ic( struct Simulation_Parameters * lc_droplet,double * Qij, const class GEOMETRY *);
void homogeneous_ic( struct Simulation_Parameters * lc_droplet,double * Qij, const class GEOMETRY *);
void homogeneous_boundary( struct Simulation_Parameters * sim_param,double * Qij, const class GEOMETRY * );
void homogeneous_easy_axis_ic( struct Simulation_Parameters * sim_param,double * Qij, const class GEOMETRY * );
void random_bulk_homogeneous_easy_axis_ic( struct Simulation_Parameters * sim_param,double * Qij, const class GEOMETRY * );
void random_hemis_ic( struct Simulation_Parameters * sim_param,double * Qij, const GEOMETRY * geometry );
void cholesteric_ic( struct Simulation_Parameters * sim_param,double * Qij, const class GEOMETRY * );
void DSS_asatz( struct Simulation_Parameters * sim_param,double * Qij, const class GEOMETRY * );
void RSS_asatz( struct Simulation_Parameters * sim_param,double * Qij, const class GEOMETRY * );
void homeotropic_asatz( struct Simulation_Parameters * sim_param,double * Qij, const GEOMETRY * geometry );
void read_check(int read_status, int line);

#endif 

