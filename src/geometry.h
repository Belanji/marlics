#ifndef geometry_
#define geometry_
#include <vector>

class GEOMETRY
{

 public:
  
  friend class BOUNDARY;

  std::vector<class BOUNDARY *> bc_conditions;

  const int Nx;
  const int Ny;
  const int Nz;

  //Initial conditions routines:
  virtual void ic(struct Simulation_Parameters *, double *) = 0;
  void read_from_file_ic(  struct Simulation_Parameters * lc_droplet,double * Qij );
  void random_ic( struct Simulation_Parameters * lc_droplet,double * Qij);
  void homogeneous_ic( struct Simulation_Parameters * lc_droplet,double * Qij);
  void homogeneous_easy_axis_ic( struct Simulation_Parameters * sim_param,double * Qij );
  void read_check(int , int );
  
  virtual void boundary_init( struct Simulation_Parameters * ) = 0;
  virtual ~GEOMETRY() {};

  virtual void fill_ki(double * k_i, const  double * Qij, const int i,const int j,const int k) const = 0;

  const int *point_type;
  
 protected:


  virtual int * fill_point_type( void ) const = 0 ;  
  
  const double sigma;
  const double a;
  const double bb;
  const double cc;
  const double dx_1;
  const double dy_1;
  const double dz_1;
  const double L1;
  const double L2;
  const double L3;
  const double Lq;
  const double Ls;
  const double Lq_tilde;
  const double Lambda;
  const double Lambda_s;
  const double Nx_1;
  const double p0;
  const double q0;
  const double  S_eq;    
  double theta_i, phi_i;

  
  double dt;

  int number_of_boundaries;  
  GEOMETRY * geometry_pointer;
  GEOMETRY(const struct Simulation_Parameters * );

  virtual  double bulk_00(const double  QN[5],const double  dQ[15],const double  ddQ[30]) const;
  virtual  double bulk_01(const double  QN[5],const double  dQ[15],const double  ddQ[30]) const;
  virtual  double bulk_02(const double  QN[5],const double  dQ[15],const double  ddQ[30]) const;
  virtual  double bulk_11(const double  QN[5],const double  dQ[15],const double  ddQ[30]) const;
  virtual  double bulk_12(const double  QN[5],const double  dQ[15],const double  ddQ[30]) const;

  
};

#endif