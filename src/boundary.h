#ifndef BOUNDARY_
#define BOUNDARY_



class BOUNDARY
{
  
 protected:

    
   const double sigma;
   const double a;
   const double bb;
   const double cc;
   const double dx_1;
   const double L1;
   const double L2;
   const double L3;
   const double Lq;
   const double Ls;
   const double Lq_tilde;
   const double Lambda, Lambda_s;
   const double  S_eq;
   const double q0;
   const double * Wo;
    double theta_0, phi_0;
   
   public:

   BOUNDARY(const GEOMETRY * calling_object);
   virtual ~BOUNDARY() {};


  virtual double surface_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  = 0;
  virtual double surface_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  = 0;
  virtual double surface_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  = 0;
  virtual double surface_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  = 0;
  virtual double surface_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const  = 0;
  
};


#endif
