#include "geometry.h"	
#include "boundary.h"



BOUNDARY::BOUNDARY(const GEOMETRY * calling_object) :
  sigma(calling_object->sigma), 
    a(calling_object->a),
    bb(calling_object->bb),
    cc(calling_object->cc),
    dx_1(calling_object->dx_1),
    L1(calling_object->L1),
    L2(calling_object->L2),
    L3(calling_object->L3),
    Lq(calling_object->Lq),
    Ls(calling_object->Ls),
    Lq_tilde(calling_object->Lq_tilde),
    Lambda(calling_object->Lambda)    ,
    Lambda_s(calling_object->Lambda_s),
    S_eq(calling_object->S_eq),
    q0(calling_object->q0)  
    { };
