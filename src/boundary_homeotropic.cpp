#include "driver.h"
#include "geometry.h"
#include "boundary.h"
#include "boundary_homeotropic.h"
#include <math.h>
#include <iostream>
#include <vector>



Boundary_Homeotropic::Boundary_Homeotropic(const Simulation_Parameters * sim_param, int boundary_number) : BOUNDARY ( sim_param, boundary_number)
{                      
  
  condition_name="Homeotropic";

  //Asserting anchoring energy is set and get its value:
  try
    {
      Wo1= sim_param->Wo1.at(boundary_id);
      std::cout << "Wo1= " << Wo1 <<".\n";
      
    }
  catch(std::out_of_range dummy_var )
    {

      assert_parameter_is_set(false, "Wo1");

    }
  dx=sim_param->dx;
  dy=sim_param->dy;
  dz=sim_param->dz;

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

//Evaluate dQ00 at boundary
double Boundary_Homeotropic::functional_derivative_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const
{
  double  Q0_00= (0.5*S_eq*(3.0*v[0]*v[0]-1.0));
  
  return Lambda_s*((-3*Wo1*QN00 + 3*Wo1*Q0_00 - 3*L3*QN01*Q_00_1*v[0] - 3*L3*QN02*Q_00_2*v[0] - 2*L2*Q_01_1*v[0] + Ls*Q_01_1*v[0] - 2*L2*Q_02_2*v[0] + Ls*Q_02_2*v[0] + 3*Lq_tilde*QN02*v[1] - 3*L1*Q_00_1*v[1] - 3*L3*QN11*Q_00_1*v[1] - 3*L3*QN12*Q_00_2*v[1] + L2*Q_01_0*v[1] - 2*Ls*Q_01_0*v[1] + L2*Q_11_1*v[1] + Ls*Q_11_1*v[1] + L2*Q_12_2*v[1] + Ls*Q_12_2*v[1] - (3*Lq_tilde*QN01 + 3*L3*QN12*Q_00_1 + (3*L1 + L2 + Ls)*Q_00_2 - 3*L3*(QN00 + QN11)*Q_00_2 - (L2 - 2*Ls)*Q_02_0 + (L2 + Ls)*(Q_11_2 - Q_12_1))*v[2] - Q_00_0*(3*L1*v[0] + 2*L2*v[0] + 2*Ls*v[0] + 3*L3*QN00*v[0] + 3*L3*QN01*v[1] + 3*L3*QN02*v[2]))/3.);
}
//Evaluate dQ01 at boundary
double Boundary_Homeotropic::functional_derivative_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double  v[3]) const
{
    double  Q0_01= (0.5*S_eq*(3.0*v[0]*v[1]));

  return Lambda_s*((-2*Wo1*QN01 + 2*Wo1*Q0_01 - Lq_tilde*QN02*v[0] - Ls*Q_00_1*v[0] - 2*L3*QN01*Q_01_1*v[0] - 2*L3*QN02*Q_01_2*v[0] - L2*Q_11_1*v[0] - L2*Q_12_2*v[0] + Lq_tilde*QN12*v[1] - L2*Q_00_0*v[1] - 2*L1*Q_01_1*v[1] - L2*Q_01_1*v[1] - Ls*Q_01_1*v[1] - 2*L3*QN11*Q_01_1*v[1] - 2*L3*QN12*Q_01_2*v[1] - L2*Q_02_2*v[1] - Ls*Q_11_0*v[1] - (Lq_tilde*(-QN00 + QN11) + 2*L3*QN12*Q_01_1 + 2*(L1 - L3*(QN00 + QN11))*Q_01_2 + Ls*(Q_02_1 + Q_12_0))*v[2] - Q_01_0*(2*L1*v[0] + L2*v[0] + Ls*v[0] + 2*L3*QN00*v[0] + 2*L3*QN01*v[1] + 2*L3*QN02*v[2]))/2.);
}
//Evaluate dQ02 at boundary
double Boundary_Homeotropic::functional_derivative_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double  v[3]) const
{
    double Q0_02= (0.5*S_eq*(3.0*v[0]*v[2]));
  
 return Lambda_s*((-2*Wo1*QN02 + 2*Wo1*Q0_02 + Lq_tilde*QN01*v[0] + (L2 - Ls)*Q_00_2*v[0] - 2*L3*QN01*Q_02_1*v[0] - 2*L3*QN02*Q_02_2*v[0] + L2*Q_11_2*v[0] - L2*Q_12_1*v[0] - (Lq_tilde*(2*QN00 + QN11) + Ls*Q_01_2)*v[1] - 2*L1*Q_02_1*v[1] - 2*L3*QN11*Q_02_1*v[1] - 2*L3*QN12*Q_02_2*v[1] - Ls*Q_12_0*v[1] - ((L2 - Ls)*Q_00_0 + L2*Q_01_1 + QN12*(Lq_tilde + 2*L3*Q_02_1) + 2*L1*Q_02_2 + (L2 + Ls - 2*L3*(QN00 + QN11))*Q_02_2 - Ls*Q_11_0)*v[2] - Q_02_0*(2*L1*v[0] + L2*v[0] + Ls*v[0] + 2*L3*QN00*v[0] + 2*L3*QN01*v[1] + 2*L3*QN02*v[2]))/2.);
}

//Evaluate dQ11 at boundary
double Boundary_Homeotropic::functional_derivative_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double  v[3]) const
{  
  double  Q0_11= (0.5*S_eq*(3.0*v[1]*v[1]-1.0));

  return Lambda_s*((-3*Wo1*QN11 + 3*Wo1*Q0_11 - 3*Lq_tilde*QN12*v[0] + L2*Q_00_0*v[0] + Ls*Q_00_0*v[0] + L2*Q_01_1*v[0] - 2*Ls*Q_01_1*v[0] + L2*Q_02_2*v[0] + Ls*Q_02_2*v[0] - 3*L1*Q_11_0*v[0] - 3*L3*QN00*Q_11_0*v[0] - 3*L3*QN01*Q_11_1*v[0] - 3*L3*QN02*Q_11_2*v[0] - 2*L2*Q_01_0*v[1] + Ls*Q_01_0*v[1] - 3*L3*QN01*Q_11_0*v[1] - 3*L1*Q_11_1*v[1] - 2*L2*Q_11_1*v[1] - 2*Ls*Q_11_1*v[1] - 3*L3*QN11*Q_11_1*v[1] - 3*L3*QN12*Q_11_2*v[1] - 2*L2*Q_12_2*v[1] + Ls*Q_12_2*v[1] + (3*Lq_tilde*QN01 - (L2 + Ls)*Q_00_2 + (L2 + Ls)*Q_02_0 - 3*L3*(QN02*Q_11_0 + QN12*Q_11_1) - (3*L1 + L2 + Ls)*Q_11_2 + 3*L3*(QN00 + QN11)*Q_11_2 + (L2 - 2*Ls)*Q_12_1)*v[2])/3.);
}
//Evaluate dQ12 at boundary
double Boundary_Homeotropic::functional_derivative_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double  v[3]) const
{
  double Q0_12= (0.5*S_eq*(3.0*v[1]*v[2]));
  
 return Lambda_s*((2*Wo1*(-QN12 + Q0_12) + Lq_tilde*(QN00 + 2*QN11)*v[0] - Ls*Q_01_2*v[0] - Ls*Q_02_1*v[0] - 2*L1*Q_12_0*v[0] - 2*L3*QN00*Q_12_0*v[0] - 2*L3*QN01*Q_12_1*v[0] - 2*L3*QN02*Q_12_2*v[0] - Lq_tilde*QN01*v[1] + L2*Q_00_2*v[1] - L2*Q_02_0*v[1] + L2*Q_11_2*v[1] - Ls*Q_11_2*v[1] - 2*L3*QN01*Q_12_0*v[1] - 2*L1*Q_12_1*v[1] - L2*Q_12_1*v[1] - Ls*Q_12_1*v[1] - 2*L3*QN11*Q_12_1*v[1] - 2*L3*QN12*Q_12_2*v[1] - (-(Lq_tilde*QN02) - Ls*Q_00_1 + L2*Q_01_0 + L2*Q_11_1 - Ls*Q_11_1 + 2*L3*QN02*Q_12_0 + 2*L3*QN12*Q_12_1 + (2*L1 + L2 + Ls - 2*L3*(QN00 + QN11))*Q_12_2)*v[2])/2.);
}

double Boundary_Homeotropic::energy_calculation(const double  QN[5],const double  dQ[15],const double  v[3]) const
{  
  double Q0_00= (0.5*S_eq*(3.0*v[0]*v[0]));
  double Q0_01= (0.5*S_eq*(3.0*v[0]*v[1]));
  double Q0_02= (0.5*S_eq*(3.0*v[0]*v[2]));
  double Q0_11= (0.5*S_eq*(3.0*v[1]*v[1]));
  double Q0_12= (0.5*S_eq*(3.0*v[1]*v[2]));
  return 0.5*Wo1*(pow(QN00 - Q0_00,2) + 2*pow(QN01 - Q0_01,2) + 2*pow(QN02 - Q0_02,2) + pow(QN11 - Q0_11,2) + pow(QN00 + QN11 - Q0_00 - Q0_11,2) + 2*pow(QN12 - Q0_12,2));
}

#undef QN00
#undef QN01
#undef QN02
#undef QN11
#undef QN12
#undef dQ00_0
#undef dQ01_0
#undef dQ02_0
#undef dQ11_0
#undef dQ12_0
#undef dQ00_1
#undef dQ01_1
#undef dQ02_1
#undef dQ11_1
#undef dQ12_1
#undef dQ00_2
#undef dQ01_2
#undef dQ02_2
#undef dQ11_2
#undef dQ12_2

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

#define dx1 (1.0/dx)
#define dy1 (1.0/dy)
#define dz1 (1.0/dz)

double inline Boundary_Homeotropic::force_00(const double  QN[27*5], const double dQ[], const double  v[3]) const 
{ 
  double  Q0_00= (0.5*S_eq*(3.0*v[0]*v[0]-1.0));
  double  Q0_11= (0.5*S_eq*(3.0*v[1]*v[1]-1.0));
  return Lambda*(Wo1*(2*QN00(1,1,1) + QN11(1,1,1) - 2*Q0_00 - Q0_11)
+(3*sigma*QN00(1,1,1) + 6*cc*QN00(1,1,1)*(pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) + pow(QN02(1,1,1),2) + QN00(1,1,1)*QN11(1,1,1) + pow(QN11(1,1,1),2) + pow(QN12(1,1,1),2)) + bb*(pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) + pow(QN02(1,1,1),2) - 2*QN00(1,1,1)*QN11(1,1,1) - 2*(pow(QN11(1,1,1),2) + pow(QN12(1,1,1),2))))/3.
+L1*(v[0]*dx1*dQ00_0(0,1,1) + v[0]*dx1*dQ00_0(2,1,1) + v[1]*dy1*dQ00_1(1,0,1) + v[1]*dy1*dQ00_1(1,2,1) + v[2]*dz1*dQ00_2(1,1,0) + v[2]*dz1*dQ00_2(1,1,2))
+0.5*Lq_tilde*(v[2]*dz1*QN01(1,1,0) + v[2]*dz1*QN01(1,1,2) - v[1]*dy1*QN02(1,0,1) - v[1]*dy1*QN02(1,2,1) - dQ01_2(1,1,1) + dQ02_1(1,1,1)));
}
double inline Boundary_Homeotropic::force_01(const double  QN[27*5], const double dQ[], const double  v[3]) const 
{ 
  double  Q0_01= (0.5*S_eq*(3.0*v[0]*v[1]));
  return Lambda*(2*Wo1*(QN01(1,1,1) - Q0_01)
+QN01(1,1,1)*(sigma + bb*(QN00(1,1,1) + QN11(1,1,1)) + 2*cc*(pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) + pow(QN02(1,1,1),2) + QN00(1,1,1)*QN11(1,1,1) + pow(QN11(1,1,1),2))) + bb*QN02(1,1,1)*QN12(1,1,1) + 2*cc*QN01(1,1,1)*pow(QN12(1,1,1),2)
+L1*(v[0]*dx1*dQ01_0(0,1,1) + v[0]*dx1*dQ01_0(2,1,1) + v[1]*dy1*dQ01_1(1,0,1) + v[1]*dy1*dQ01_1(1,2,1) + v[2]*dz1*dQ01_2(1,1,0) + v[2]*dz1*dQ01_2(1,1,2))
+(0.5*Lq_tilde*(v[0]*dx1*QN02(0,1,1) + v[0]*dx1*QN02(2,1,1) + v[2]*dz1*(-QN00(1,1,0) - QN00(1,1,2) + QN11(1,1,0) + QN11(1,1,2)) - v[1]*dy1*QN12(1,0,1) - v[1]*dy1*QN12(1,2,1) + dQ00_2(1,1,1) - dQ02_0(1,1,1) - dQ11_2(1,1,1) + dQ12_1(1,1,1)))/2.);
}
double inline Boundary_Homeotropic::force_02(const double  QN[27*5], const double dQ[], const double  v[3]) const 
{ 
  double Q0_02= (0.5*S_eq*(3.0*v[0]*v[2]));
  return Lambda*(2*Wo1*(QN02(1,1,1) - Q0_02)
+QN02(1,1,1)*(sigma - bb*QN11(1,1,1) + 2*cc*(pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) + pow(QN02(1,1,1),2) + QN00(1,1,1)*QN11(1,1,1) + pow(QN11(1,1,1),2))) + bb*QN01(1,1,1)*QN12(1,1,1) + 2*cc*QN02(1,1,1)*pow(QN12(1,1,1),2)
+L1*(v[0]*dx1*dQ02_0(0,1,1) + v[0]*dx1*dQ02_0(2,1,1) + v[1]*dy1*dQ02_1(1,0,1) + v[1]*dy1*dQ02_1(1,2,1) + v[2]*dz1*dQ02_2(1,1,0) + v[2]*dz1*dQ02_2(1,1,2))
+-(0.5*Lq_tilde*(v[0]*dx1*QN01(0,1,1) + v[0]*dx1*QN01(2,1,1) + v[1]*dy1*(-2*QN00(1,0,1) - 2*QN00(1,2,1) - QN11(1,0,1) - QN11(1,2,1)) - v[2]*dz1*QN12(1,1,0) - v[2]*dz1*QN12(1,1,2) + 2*dQ00_1(1,1,1) - dQ01_0(1,1,1) + dQ11_1(1,1,1) + dQ12_2(1,1,1)))/2.);
}
double inline Boundary_Homeotropic::force_11(const double  QN[27*5], const double dQ[], const double  v[3]) const 
{ 
  double  Q0_00= (0.5*S_eq*(3.0*v[0]*v[0]-1.0));
  double  Q0_11= (0.5*S_eq*(3.0*v[1]*v[1]-1.0));
  return Lambda*(Wo1*(QN00(1,1,1) + 2*QN11(1,1,1) - Q0_00 - 2*Q0_11)
+(bb*(-2*pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) - 2*pow(QN02(1,1,1),2) - 2*QN00(1,1,1)*QN11(1,1,1) + pow(QN11(1,1,1),2) + pow(QN12(1,1,1),2)) + 3*QN11(1,1,1)*(sigma + 2*cc*(pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) + pow(QN02(1,1,1),2) + QN00(1,1,1)*QN11(1,1,1) + pow(QN11(1,1,1),2) + pow(QN12(1,1,1),2))))/3.
+L1*(v[0]*dx1*dQ11_0(0,1,1) + v[0]*dx1*dQ11_0(2,1,1) + v[1]*dy1*dQ11_1(1,0,1) + v[1]*dy1*dQ11_1(1,2,1) + v[2]*dz1*dQ11_2(1,1,0) + v[2]*dz1*dQ11_2(1,1,2))
+0.5*Lq_tilde*(-(v[2]*dz1*QN01(1,1,0)) - v[2]*dz1*QN01(1,1,2) + v[0]*dx1*QN12(0,1,1) + v[0]*dx1*QN12(2,1,1) + dQ01_2(1,1,1) - dQ12_0(1,1,1)));
}
double inline Boundary_Homeotropic::force_12(const double  QN[27*5], const double dQ[], const double  v[3]) const 
{ 
  double Q0_12= (0.5*S_eq*(3.0*v[1]*v[2]));
  return Lambda*(2*Wo1*(QN12(1,1,1) - Q0_12)
+bb*(QN01(1,1,1)*QN02(1,1,1) - QN00(1,1,1)*QN12(1,1,1)) + QN12(1,1,1)*(sigma + 2*cc*(pow(QN00(1,1,1),2) + pow(QN01(1,1,1),2) + pow(QN02(1,1,1),2) + QN00(1,1,1)*QN11(1,1,1) + pow(QN11(1,1,1),2) + pow(QN12(1,1,1),2)))
+L1*(v[0]*dx1*dQ12_0(0,1,1) + v[0]*dx1*dQ12_0(2,1,1) + v[1]*dy1*dQ12_1(1,0,1) + v[1]*dy1*dQ12_1(1,2,1) + v[2]*dz1*dQ12_2(1,1,0) + v[2]*dz1*dQ12_2(1,1,2))
+(0.5*Lq_tilde*(v[1]*dy1*QN01(1,0,1) + v[1]*dy1*QN01(1,2,1) - v[2]*dz1*QN02(1,1,0) - v[2]*dz1*QN02(1,1,2) + v[0]*dx1*(-QN00(0,1,1) - QN00(2,1,1) - 2*QN11(0,1,1) - 2*QN11(2,1,1)) + dQ00_0(1,1,1) - dQ01_1(1,1,1) + dQ02_2(1,1,1) + 2*dQ11_0(1,1,1)))/2.);
}
