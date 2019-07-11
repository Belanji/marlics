#include "boundary.h"
#include "driver.h"


BOUNDARY::BOUNDARY(const  Simulation_Parameters *  sim_par, int boundary_id):
  sigma(sim_par->T*sim_par->a), 
  a(sim_par->a),
  bb(sim_par->B),
  cc(sim_par->C),
  dx_1(1/sim_par->dx),
  dy_1(1/sim_par->dy),
  dz_1(1/sim_par->dz),
  L1(sim_par->L1),
  L2(sim_par->L2),
  L3(sim_par->L3),
  Lq(sim_par->Lq),
  Ls(sim_par->Ls),
  Lq_tilde(2.0*sim_par->Lq*sim_par->q0),
  Lambda(1/sim_par->mu_1)    ,
  Lambda_s(1/sim_par->mu_1_s),
  S_eq(sim_par->S_eq),
  q0(sim_par->q0)  
  { };
