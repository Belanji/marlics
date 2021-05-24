#ifndef ENERGY_LDG__
#define ENERGY_LDG__

class landau_de_gennes : public ENERGY
{
 public:

  landau_de_gennes(const struct Simulation_Parameters * lc);

    virtual double functional_derivative_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
    virtual double functional_derivative_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
    virtual double functional_derivative_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
    virtual double functional_derivative_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
    virtual double functional_derivative_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const ;
    virtual double energy_calculation      (const double  QN[5],const double  dQ[15],const double v[3]) const ;

};



#endif
