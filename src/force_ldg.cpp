#include "driver.h"
#include "force.h"
#include "force_ldg.h"
#include <stdio.h>
#include <math.h>

force_ldg::force_ldg(const struct Simulation_Parameters * lc) : FORCE(lc)
{
        
  std::cout << "L1=" << L1 << " \n";
  std::cout << "L2=" << L2 << " \n";
  std::cout << "L3=" << L3 << " \n";
  std::cout << "Lq=" << Lq << " \n";
  std::cout << "Ls=" << Ls << " \n\n";
  std::cout << "p0=" << p0 << " \n";
  std::cout << "dx=" << lc->dx << " \n";
  std::cout << "dy=" << lc->dy << " \n";
  std::cout << "dz=" << lc->dz << " \n";
  std::cout << "S_eq=" << S_eq << " \n\n";
        
  std::cout << "a=" << a << " \n";
  std::cout << "b=" << bb << " \n";
  std::cout << "c=" << cc << " \n";
  std::cout << "T=" << lc->T << " \n\n";

  std::cout << "Lambda=" << 1/Lambda << " \n";
  std::cout << "Lambda_s=" << 1/Lambda_s << " \n\n";
  
  dx1=1.0/lc->dx ;
  dy1=1.0/lc->dy ;
  dz1=1.0/lc->dz ;

};
#define QN00(i,j,k) QN[0+5*(-3+i+2*j+3*k)] 
#define QN01(i,j,k) QN[1+5*(-3+i+2*j+3*k)]
#define QN02(i,j,k) QN[2+5*(-3+i+2*j+3*k)]
#define QN11(i,j,k) QN[3+5*(-3+i+2*j+3*k)]
#define QN12(i,j,k) QN[4+5*(-3+i+2*j+3*k)]

#define dQ00_0(i,j,k) dQ[0+3*(0+5*(-3+i+2*j+3*k))] 
#define dQ01_0(i,j,k) dQ[0+3*(1+5*(-3+i+2*j+3*k))]
#define dQ02_0(i,j,k) dQ[0+3*(2+5*(-3+i+2*j+3*k))]
#define dQ11_0(i,j,k) dQ[0+3*(3+5*(-3+i+2*j+3*k))]
#define dQ12_0(i,j,k) dQ[0+3*(4+5*(-3+i+2*j+3*k))]

#define dQ00_1(i,j,k) dQ[1+3*(0+5*(-3+i+2*j+3*k))] 
#define dQ01_1(i,j,k) dQ[1+3*(1+5*(-3+i+2*j+3*k))]
#define dQ02_1(i,j,k) dQ[1+3*(2+5*(-3+i+2*j+3*k))]
#define dQ11_1(i,j,k) dQ[1+3*(3+5*(-3+i+2*j+3*k))]
#define dQ12_1(i,j,k) dQ[1+3*(4+5*(-3+i+2*j+3*k))]

#define dQ00_2(i,j,k) dQ[2+3*(0+5*(-3+i+2*j+3*k))] 
#define dQ01_2(i,j,k) dQ[2+3*(1+5*(-3+i+2*j+3*k))]
#define dQ02_2(i,j,k) dQ[2+3*(2+5*(-3+i+2*j+3*k))]
#define dQ11_2(i,j,k) dQ[2+3*(3+5*(-3+i+2*j+3*k))]
#define dQ12_2(i,j,k) dQ[2+3*(4+5*(-3+i+2*j+3*k))]


double inline force_ldg::force_00(const double  QN[],const double  dQ[]) const 
{ 
 return Lambda*((3*sigma*QN00(1,1,1) + 6*cc*QN00(1,1,1)*(pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) + pow(QN02(1,1,1),2) + QN00(1,1,1)*QN11(1,1,1) + pow(QN11(1,1,1),2) + pow(QN12(1,1,1),2)) + bb*(pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) + pow(QN02(1,1,1),2) - 2*QN00(1,1,1)*QN11(1,1,1) - 2*(pow(QN11(1,1,1),2) + pow(QN12(1,1,1),2))))/3.
+L1*(dx1*dQ00_0(0,1,1) - dx1*dQ00_0(2,1,1) + dy1*dQ00_1(1,0,1) - dy1*dQ00_1(1,2,1) + dz1*dQ00_2(1,1,0) - dz1*dQ00_2(1,1,2))
+0.5*Lq_tilde*(dz1*QN01(1,1,0) - dz1*QN01(1,1,2) - dy1*QN02(1,0,1) + dy1*QN02(1,2,1) - dQ01_2(1,1,1) + dQ02_1(1,1,1)));
}
double inline force_ldg::force_01(const double  QN[],const double  dQ[]) const 
{ 
 return Lambda*(QN01(1,1,1)*(sigma + bb*(QN00(1,1,1) + QN11(1,1,1)) + 2*cc*(pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) + pow(QN02(1,1,1),2) + QN00(1,1,1)*QN11(1,1,1) + pow(QN11(1,1,1),2))) + bb*QN02(1,1,1)*QN12(1,1,1) + 2*cc*QN01(1,1,1)*pow(QN12(1,1,1),2)
+L1*(dx1*dQ01_0(0,1,1) - dx1*dQ01_0(2,1,1) + dy1*dQ01_1(1,0,1) - dy1*dQ01_1(1,2,1) + dz1*dQ01_2(1,1,0) - dz1*dQ01_2(1,1,2))
+(0.5*Lq_tilde*(dx1*QN02(0,1,1) - dx1*QN02(2,1,1) + dz1*(-QN00(1,1,0) + QN00(1,1,2) + QN11(1,1,0) - QN11(1,1,2)) - dy1*QN12(1,0,1) + dy1*QN12(1,2,1) + dQ00_2(1,1,1) - dQ02_0(1,1,1) - dQ11_2(1,1,1) + dQ12_1(1,1,1)))/2.);
}
double inline force_ldg::force_02(const double  QN[],const double  dQ[]) const 
{ 
 return Lambda*(QN02(1,1,1)*(sigma - bb*QN11(1,1,1) + 2*cc*(pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) + pow(QN02(1,1,1),2) + QN00(1,1,1)*QN11(1,1,1) + pow(QN11(1,1,1),2))) + bb*QN01(1,1,1)*QN12(1,1,1) + 2*cc*QN02(1,1,1)*pow(QN12(1,1,1),2)
+L1*(dx1*dQ02_0(0,1,1) - dx1*dQ02_0(2,1,1) + dy1*dQ02_1(1,0,1) - dy1*dQ02_1(1,2,1) + dz1*dQ02_2(1,1,0) - dz1*dQ02_2(1,1,2))
+-(0.5*Lq_tilde*(dx1*QN01(0,1,1) - dx1*QN01(2,1,1) + dy1*(-2*QN00(1,0,1) + 2*QN00(1,2,1) - QN11(1,0,1) + QN11(1,2,1)) - dz1*QN12(1,1,0) + dz1*QN12(1,1,2) + 2*dQ00_1(1,1,1) - dQ01_0(1,1,1) + dQ11_1(1,1,1) + dQ12_2(1,1,1)))/2.);
}
double inline force_ldg::force_11(const double  QN[],const double  dQ[]) const 
{ 
 return Lambda*((bb*(-2*pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) - 2*pow(QN02(1,1,1),2) - 2*QN00(1,1,1)*QN11(1,1,1) + pow(QN11(1,1,1),2) + pow(QN12(1,1,1),2)) + 3*QN11(1,1,1)*(sigma + 2*cc*(pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) + pow(QN02(1,1,1),2) + QN00(1,1,1)*QN11(1,1,1) + pow(QN11(1,1,1),2) + pow(QN12(1,1,1),2))))/3.
+L1*(dx1*dQ11_0(0,1,1) - dx1*dQ11_0(2,1,1) + dy1*dQ11_1(1,0,1) - dy1*dQ11_1(1,2,1) + dz1*dQ11_2(1,1,0) - dz1*dQ11_2(1,1,2))
+0.5*Lq_tilde*(-(dz1*QN01(1,1,0)) + dz1*QN01(1,1,2) + dx1*QN12(0,1,1) - dx1*QN12(2,1,1) + dQ01_2(1,1,1) - dQ12_0(1,1,1)));
}
double inline force_ldg::force_12(const double  QN[],const double  dQ[]) const 
{ 
 return Lambda*(bb*(QN01(1,1,1)*QN02(1,1,1) - QN00(1,1,1)*QN12(1,1,1)) + QN12(1,1,1)*(sigma + 2*cc*(pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) + pow(QN02(1,1,1),2) + QN00(1,1,1)*QN11(1,1,1) + pow(QN11(1,1,1),2) + pow(QN12(1,1,1),2)))
+L1*(dx1*dQ12_0(0,1,1) - dx1*dQ12_0(2,1,1) + dy1*dQ12_1(1,0,1) - dy1*dQ12_1(1,2,1) + dz1*dQ12_2(1,1,0) - dz1*dQ12_2(1,1,2))
+(0.5*Lq_tilde*(dy1*QN01(1,0,1) - dy1*QN01(1,2,1) - dz1*QN02(1,1,0) + dz1*QN02(1,1,2) + dx1*(-QN00(0,1,1) + QN00(2,1,1) - 2*QN11(0,1,1) + 2*QN11(2,1,1)) + dQ00_0(1,1,1) - dQ01_1(1,1,1) + dQ02_2(1,1,1) + 2*dQ11_0(1,1,1)))/2.);
}
