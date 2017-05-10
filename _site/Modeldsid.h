// modeldsid.h
//   A Continuum Damage Model used in FLAC3D.
// History:
// 2016/09/25  Hao Xu   First release 

#ifndef MODELDSID_H_
#define MODELDSID_H_

//#include "../mathLib/arithmetic.h"
#include "../mathLib/r1Tensor.h"
#include "../mathLib/r2Tensor.h"
#include "../mathLib/r3Tensor.h"
//#include "../mathLib/errInfo.h"
//#include "../mathLib/gaussj.h"
//#include "../mathLib/eigen.h"

class Modeldsid {
    public:
      // Creators
      Modeldsid(double E0, double Poisson0, double a1,    double a2, double a3,
                double a4, double C0, double C1, double alpha, double Debug);
      void          run(r1Tensor<double> & stnE_, r1Tensor<double> & stnS_, r1Tensor<double> & stnO_, double & fd0); 
      void          printVariables(); 

    private:
      double Bulk_,Shear_,E0_,Poisson0_,a1_,a2_,a3_,a4_,C0_,C1_;
      double alpha_,Debug_,Omega_00_,Omega_11_,Omega_22_;
      double Omega_01_,Omega_12_,Omega_20_,Epsid_00_,Epsid_11_;
      double Epsid_22_,Epsid_01_,Epsid_12_,Epsid_20_;

      r2Tensor<double> Matdom;
      r1Tensor<double> Omega,Epsid,dstran,Stress,dSig;
};

#endif 
