#ifndef Force_LDG__
#define Force_LDG__

class force_ldg : public FORCE
{

 public:

  force_ldg(const struct Simulation_Parameters * lc);
    double dx1, dy1, dz1;

    virtual double force_00(const double  QN[27*5], const double dQ[]) const ;
    virtual double force_01(const double  QN[27*5], const double dQ[]) const ;
    virtual double force_02(const double  QN[27*5], const double dQ[]) const ;
    virtual double force_11(const double  QN[27*5], const double dQ[]) const ;
    virtual double force_12(const double  QN[27*5], const double dQ[]) const ;

};



#endif
