#include "boundary_strong.h"
#include "geometry.h"
#include <iostream>

Strong_Boundary::Strong_Boundary(const Simulation_Parameters * sim_param, int boundary_id) : BOUNDARY ( sim_param, boundary_id )
{

  condition_name="strong";
  std::cout << std::endl;
  
  
} ;

double Strong_Boundary::surface_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const { return 0.0; }  ;
double Strong_Boundary::surface_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const { return 0.0; }  ;
double Strong_Boundary::surface_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const { return 0.0; }  ;
double Strong_Boundary::surface_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const { return 0.0; }  ;
double Strong_Boundary::surface_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const { return 0.0; }  ;

Strong_Boundary::~Strong_Boundary() {};

