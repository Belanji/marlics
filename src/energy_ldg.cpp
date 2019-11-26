#include "driver.h"
#include "energy.h"
#include "energy_ldg.h"
#include <stdio.h>

landau_de_gennes::landau_de_gennes(const struct Simulation_Parameters * lc) : ENERGY(lc)
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
  

};

#define QN00 QN[0] 
#define QN01 QN[1]
#define QN02 QN[2]
#define QN11 QN[3]
#define QN12 QN[4]

#define  Q_00_0 dQ[0]
#define  Q_01_0 dQ[1]
#define  Q_02_0 dQ[2]
#define  Q_11_0 dQ[3]
#define  Q_12_0 dQ[4]

#define  Q_00_1 dQ[5]
#define  Q_01_1 dQ[6]
#define  Q_02_1 dQ[7]
#define  Q_11_1 dQ[8]
#define  Q_12_1 dQ[9]

#define  Q_00_2 dQ[10]
#define  Q_01_2 dQ[11]
#define  Q_02_2 dQ[12]
#define  Q_11_2 dQ[13]
#define  Q_12_2 dQ[14]


#define  Q_00_00 ddQ[0]
#define  Q_01_00 ddQ[1]
#define  Q_02_00 ddQ[2]
#define  Q_11_00 ddQ[3]
#define  Q_12_00 ddQ[4]
                
#define  Q_00_01 ddQ[5]
#define  Q_01_01 ddQ[6]
#define  Q_02_01 ddQ[7]
#define  Q_11_01 ddQ[8]
#define  Q_12_01 ddQ[9]
                
#define  Q_00_02 ddQ[10]
#define  Q_01_02 ddQ[11]
#define  Q_02_02 ddQ[12]
#define  Q_11_02 ddQ[13]
#define  Q_12_02 ddQ[14]

#define  Q_00_11 ddQ[15]
#define  Q_01_11 ddQ[16]
#define  Q_02_11 ddQ[17]
#define  Q_11_11 ddQ[18]
#define  Q_12_11 ddQ[19]

#define  Q_00_12 ddQ[20]
#define  Q_01_12 ddQ[21]
#define  Q_02_12 ddQ[22]
#define  Q_11_12 ddQ[23]
#define  Q_12_12 ddQ[24]

#define  Q_00_22 ddQ[25]
#define  Q_01_22 ddQ[26]
#define  Q_02_22 ddQ[27]
#define  Q_11_22 ddQ[28]
#define  Q_12_22 ddQ[29]

double landau_de_gennes::functional_derivative_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const 
  { 
    return Lambda*((-deltaepslon*(-2*(elecfieldx*elecfieldx) + (elecfieldy*elecfieldy) + (elecfieldz*elecfieldz)) - 3*sigma*QN00 - 6*cc*QN00*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12)) - bb*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) - 2*QN00*QN11 - 2*((QN11*QN11) + (QN12*QN12))) + 6*Lq_tilde*(Q_01_2 - Q_02_1) + (3*L1 + 2*(L2 + Ls))*Q_00_00 + 3*L1*Q_00_11 + 3*L1*Q_00_22 + L2*Q_00_22 + Ls*Q_00_22 + L3*((Q_00_0*Q_00_0) + (Q_00_1*Q_00_1) - 2*(Q_00_2*Q_00_2) - 2*(Q_01_0*Q_01_0) + (Q_01_1*Q_01_1) + (Q_01_2*Q_01_2) - 2*(Q_02_0*Q_02_0) + (Q_02_1*Q_02_1) + (Q_02_2*Q_02_2) + Q_00_0*(3*(Q_01_1 + Q_02_2) - 2*Q_11_0) - 2*(Q_11_0*Q_11_0) + (Q_11_1*Q_11_1) + (Q_11_2*Q_11_2) - 2*(Q_12_0*Q_12_0) + (Q_12_1*Q_12_1) + Q_00_2*(3*Q_02_0 - 2*Q_11_2 + 3*Q_12_1) + (Q_12_2*Q_12_2) + Q_00_1*(3*Q_01_0 + 4*Q_11_1 + 3*Q_12_2) + 3*QN00*Q_00_00 + 6*QN01*Q_00_01 + 6*QN02*Q_00_02 + 3*QN11*Q_00_11 + 6*QN12*Q_00_12 - 3*(QN00 + QN11)*Q_00_22) + L2*Q_01_01 + Ls*Q_01_01 + L2*Q_02_02 + Ls*Q_02_02 - L2*Q_11_11 - Ls*Q_11_11 + L2*Q_11_22 + Ls*Q_11_22 - 2*(L2 + Ls)*Q_12_12)/3.);
  }

double landau_de_gennes::functional_derivative_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const 
  { 
    return Lambda*((-2*Q_00_2*(Lq_tilde + L3*Q_01_2) - L3*Q_00_0*(2*Q_00_1 - 2*Q_01_0 + Q_11_1) + 2*Lq_tilde*(Q_02_0 + Q_11_2 - Q_12_1) - 2*(-deltaepslon*elecfieldx*elecfieldy + sigma*QN01 + bb*QN01*(QN00 + QN11) + bb*QN02*QN12 + 2*cc*QN01*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12)) - L1*(Q_01_00 + Q_01_11 + Q_01_22)) + L3*(-2*Q_02_0*Q_02_1 + 2*Q_01_0*(Q_01_1 + Q_02_2) - Q_11_0*(Q_00_1 + 2*Q_11_1) + 2*Q_01_2*(Q_02_0 - Q_11_2 + Q_12_1) + 2*(-(Q_12_0*Q_12_1) + Q_01_1*(Q_11_1 + Q_12_2) + QN00*Q_01_00 + 2*QN01*Q_01_01 + 2*QN02*Q_01_02 + QN11*Q_01_11 + 2*QN12*Q_01_12 - (QN00 + QN11)*Q_01_22)) + (L2 + Ls)*(Q_00_01 + Q_01_00 + Q_01_11 + Q_02_12 + Q_11_01 + Q_12_02))/2.);
  }

double landau_de_gennes::functional_derivative_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const 
  { 
    return Lambda*((4*Lq_tilde*Q_00_1 + 2*Lq_tilde*(-Q_01_0 + Q_11_1 + Q_12_2) - 2*(-deltaepslon*elecfieldx*elecfieldz + sigma*QN02 - bb*QN02*QN11 + bb*QN01*QN12 + 2*cc*QN02*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12)) - L1*(Q_02_00 + Q_02_11 + Q_02_22)) + L3*(-(Q_00_2*(2*Q_02_2 + Q_11_0)) - Q_00_0*(2*Q_00_2 - 2*Q_02_0 + Q_11_2) + 2*(Q_01_1*Q_02_0 + Q_02_0*Q_02_2 - (Q_02_2 + Q_11_0)*Q_11_2 + Q_02_2*Q_12_1) + 2*(Q_01_0*(-Q_01_2 + Q_02_1) - Q_12_0*Q_12_2 + Q_02_1*(Q_11_1 + Q_12_2) + QN00*Q_02_00 + 2*QN01*Q_02_01 + 2*QN02*Q_02_02 + QN11*Q_02_11 + 2*QN12*Q_02_12 - (QN00 + QN11)*Q_02_22)) + (L2 + Ls)*(Q_01_12 + Q_02_00 + Q_02_22 - Q_11_02 + Q_12_01))/2.);
  }

double landau_de_gennes::functional_derivative_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const 
  { 
    return Lambda*((-deltaepslon*((elecfieldx*elecfieldx) - 2*(elecfieldy*elecfieldy) + (elecfieldz*elecfieldz)) + bb*(2*(QN00*QN00) - (QN01*QN01) + 2*(QN02*QN02) + 2*QN00*QN11 - (QN11*QN11) - (QN12*QN12)) - 3*QN11*(sigma + 2*cc*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12))) + 6*Lq_tilde*(-Q_01_2 + Q_12_0) - (L2 + Ls)*Q_00_00 + 3*L1*(Q_11_00 + Q_11_11 + Q_11_22) + L3*((Q_00_0*Q_00_0) - 2*(Q_00_1*Q_00_1) + (Q_00_2*Q_00_2) + (Q_01_0*Q_01_0) - 2*(Q_01_1*Q_01_1) + (Q_01_2*Q_01_2) + (Q_02_0*Q_02_0) - 2*(Q_02_1*Q_02_1) + (Q_02_2*Q_02_2) + (4*Q_00_0 + 3*(Q_01_1 + Q_02_2))*Q_11_0 + (Q_11_0*Q_11_0) - 2*Q_00_1*Q_11_1 + 3*Q_01_0*Q_11_1 + (Q_11_1*Q_11_1) - 2*Q_00_2*Q_11_2 + 3*Q_02_0*Q_11_2 - 2*(Q_11_2*Q_11_2) + (Q_12_0*Q_12_0) + 3*Q_11_2*Q_12_1 - 2*(Q_12_1*Q_12_1) + 3*Q_11_1*Q_12_2 + (Q_12_2*Q_12_2) + 3*QN00*Q_11_00 + 6*QN01*Q_11_01 + 6*QN02*Q_11_02 + 3*QN11*Q_11_11 + 6*QN12*Q_11_12 - 3*(QN00 + QN11)*Q_11_22) + (L2 + Ls)*(Q_00_22 + Q_01_01 - 2*Q_02_02 + 2*Q_11_11 + Q_11_22 + Q_12_12))/3.);
  }

double landau_de_gennes::functional_derivative_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const
  { 
    return Lambda*((2*Lq_tilde*(Q_01_1 - Q_02_2 - 2*Q_11_0) - 2*Q_00_0*(Lq_tilde - L3*Q_12_0) - (L2 + Ls)*Q_00_12 + (L2 + Ls)*(Q_01_02 + Q_02_01 + Q_12_11 + Q_12_22) - 2*(-deltaepslon*elecfieldy*elecfieldz + bb*QN01*QN02 + (sigma - bb*QN00 + 2*cc*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11)))*QN12 + 2*cc*(QN12*QN12*QN12) - L1*(Q_12_00 + Q_12_11 + Q_12_22)) + L3*(-(Q_00_1*(2*Q_00_2 + Q_11_2)) + 2*Q_02_0*Q_12_2 - 2*(Q_11_2 - Q_12_1)*(Q_11_1 + Q_12_2) - Q_00_2*(Q_11_1 + 2*Q_12_2) + 2*(Q_01_1*(-Q_01_2 + Q_12_0) + Q_02_2*(-Q_02_1 + Q_12_0) + Q_01_0*Q_12_1 + QN00*Q_12_00 + 2*QN01*Q_12_01 + 2*QN02*Q_12_02 + QN11*Q_12_11 + 2*QN12*Q_12_12 - (QN00 + QN11)*Q_12_22)))/2.);
  }
