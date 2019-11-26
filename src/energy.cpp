#include "driver.h"
#include "energy.h"


ENERGY::ENERGY(const struct Simulation_Parameters * lc) :
    Lambda(1/lc->mu_1),
    Lambda_s(1/(lc->mu_1_s)),
    a(lc->a),
    sigma(lc->a*lc->T),
    bb(lc->B),
    cc(lc->C),
    L1(lc->L1),
    L2(lc->L2),
    L3(lc->L3),
    Lq(lc->Lq),
    Ls(lc->Ls),
    Lq_tilde(lc->Lq*2.0*lc->q0),
    S_eq(lc->S_eq),
    p0(lc->p0),
    q0(lc->q0),
    deltaepslon(8.8541878176e-3*lc->deltaepslon),
    elecfieldx(lc->elecfieldx),
    elecfieldy(lc->elecfieldy),
    elecfieldz(lc->elecfieldz)
{ }


//  void  assert_parameter_is_set(bool parameter, std::string parameter_name, bool has_standard_value =false);
