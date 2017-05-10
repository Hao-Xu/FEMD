// modeldsid.cpp
// subroutine:
//   A Continuum Damage Model used in FLAC3D.
// History:
// 2016/09/25  Hao Xu   First release 

#include "Modeldsid.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <assert.h>
#include "../mathLib/arithmetic.h"
#include "../mathLib/r1Tensor.h"
#include "../mathLib/r2Tensor.h"
#include "../mathLib/r3Tensor.h"
#include "../mathLib/errInfo.h"
#include "../mathLib/gaussj.h"
#include "../mathLib/eigen.h"



void effectiveStiffness(r2Tensor<double> &Matdom, const r1Tensor<double> &Omega, const double & E0_,
                        const double & Poisson0_, const double & a1_, const double & a2_,
                        const double & a3_, const double & a4_, const double & C0_,
                        const double & C1_, const double & alpha_);

void damageFunction(double & fd, const r1Tensor<double> &Sigma, const r1Tensor<double> &Omega, const double & E0_,
                   const double & Poisson0_, const double & a1_, const double & a2_,
                   const double & a3_, const double & a4_, const double & C0_,
                   const double & C1_, const double & alpha_, const int & ioptfd);

void matP1(r2Tensor<double> &P1, const r1Tensor<double> &Sigma);

void matP2(r2Tensor<double> &P2, const r1Tensor<double> &Sigma);

void dY_dSigFunction(const r1Tensor<double> &Sigma, r2Tensor<double> &dY_dSig,const double & E0_,
             const double & Poisson0_, const double & a1_, const double & a2_,
             const double & a3_, const double & a4_, const double & C0_,
             const double & C1_, const double & alpha_);

// SUB FD_LAM
void cuttingPlaneMethod(const r1Tensor<double> &Omega, const r1Tensor<double> &Sigma, r2Tensor<double> &Matdom,const double & E0_,
             const double & Poisson0_, const double & a1_, const double & a2_,
             const double & a3_, const double & a4_, const double & C0_,
             const double & C1_, const double & alpha_, double & H0,
             double & HP, r1Tensor<double> &dG_dY, r1Tensor<double> &df_dSig,
             r1Tensor<double> &temp, const int & iopt);

void Mat_dS_dOmega(r3Tensor<double> &dS_dO, const double & E0_,
                   const double & Poisson0_, const double & a1_, const double & a2_,
                   const double & a3_, const double & a4_, const double & C0_,
                   const double & C1_, const double & alpha_);



    Modeldsid::Modeldsid(double E0, double Poisson0, double a1,    double a2, double a3, double a4,
                         double C0, double C1,       double alpha, double Debug): 
               E0_(E0),Poisson0_(Poisson0),a1_(a1),
               a2_(a2),a3_(a3),a4_(a4),C0_(C0),C1_(C1),
               alpha_(alpha),Debug_(Debug),Omega_00_(0.0),Omega_11_(0.0),Omega_22_(0.0),
               Omega_01_(0.0),Omega_12_(0.0),Omega_20_(0.0),Epsid_00_(0.0),Epsid_11_(0.0),
               Epsid_22_(0.0),Epsid_01_(0.0),Epsid_12_(0.0),Epsid_20_(0.0),Matdom(6,6,0.),
               Omega(6,0.),Epsid(6,0.),dstran(6,0.),Stress(6,0.),dSig(6,0.) {   }


    //--------------------------------------------------------------------------
    void Modeldsid::run(r1Tensor<double> & stnE_, r1Tensor<double> & stnS_ , r1Tensor<double> & stnO_, double & fd0) {

        /* --- trial elastic stresses --- */
        int ntens = Omega.size();
        double zero = 0.;
        double deps = 0.;
        double fd1, fd2, XL;
        int iopt, ioptfd;
        double H0, HP;
        r1Tensor<double> dG_dY(ntens), df_dSig(ntens), temp(ntens);
        const double tol = 1e-6, ITmax = 250, tiny = 1e-3;
        r1Tensor<double> Omega0(ntens), Epsid0(ntens);
        Omega[0] = Omega_00_;
        Omega[1] = Omega_11_;
        Omega[2] = Omega_22_;
        Omega[3] = Omega_01_;
        Omega[4] = Omega_12_;
        Omega[5] = Omega_20_;
        Omega0 = Omega;
        effectiveStiffness(Matdom, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                           C0_, C1_, alpha_);
        Epsid[0] = Epsid_00_;
        Epsid[1] = Epsid_11_;
        Epsid[2] = Epsid_22_;
        Epsid[3] = Epsid_01_;
        Epsid[4] = Epsid_12_;
        Epsid[5] = Epsid_20_;
        Epsid0 = Epsid;

        dstran[0] = stnE_[0];
        dstran[1] = stnE_[1];
        dstran[2] = stnE_[2];
        dstran[3] = stnE_[3];
        dstran[4] = stnE_[4];
        dstran[5] = stnE_[5];

        Stress[0] = stnS_[0];
        Stress[1] = stnS_[1];
        Stress[2] = stnS_[2];
        Stress[3] = stnS_[3];
        Stress[4] = stnS_[4];
        Stress[5] = stnS_[5];
        for (int i=0; i<ntens; i++) {deps += dstran[i]*dstran[i];}
        deps = sqrt(deps);
        for (int i=0; i<ntens; i++) {
            dSig[i] = 0.;
            for (int j=0; j<ntens; j++) {
                if (j<3) {
                    dSig[i] += Matdom[i][j]*dstran[j];
                } else {
                    dSig[i] += 2.* Matdom[i][j]*dstran[j];
                }
            }
        }
        for (int i=0; i<ntens; i++) {
            Stress[i] += dSig[i];
        }
        stnS_[0] = Stress[0];
        stnS_[1] = Stress[1];
        stnS_[2] = Stress[2];
        stnS_[3] = Stress[3];
        stnS_[4] = Stress[4];
        stnS_[5] = Stress[5];

        ioptfd = 1;
        damageFunction(fd1, Stress, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                       C0_, C1_, alpha_, ioptfd);
        if (fd1 > tol) {
            int Inc = 0;
            double fdt = fd1;
            while(fdt>0. && (fdt/fd1)>tol && Inc < ITmax) {
                iopt = 0;
                if (Inc == 0) iopt =0;
                cuttingPlaneMethod(Omega, Stress, Matdom, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                                   C0_, C1_, alpha_, H0, HP, dG_dY, df_dSig, temp, iopt);
                XL = fdt/(H0-HP);
                if (deps>tiny) XL/=1.5;
                for (int i=0; i<ntens; i++) {
                    Omega[i] += XL*dG_dY[i];
                    Epsid[i] += XL*df_dSig[i];
                    Stress[i] -= XL*temp[i];
                }
                ioptfd = 2;
                damageFunction(fdt, Stress, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                       C0_, C1_, alpha_, ioptfd);
                ++Inc;
            }
            if (Inc>=ITmax) throw std::runtime_error("DSID: no convergence");
        }
            ioptfd = 3;
            damageFunction(fd2, Stress, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                           C0_, C1_, alpha_, ioptfd);
            fd0 = fd2; 
            stnS_[0] = Stress[0];
            stnS_[1] = Stress[1];
            stnS_[2] = Stress[2];
            stnS_[3] = Stress[3];
            stnS_[4] = Stress[4];
            stnS_[5] = Stress[5];

            Omega_00_ = Omega[0];
            Omega_11_ = Omega[1];
            Omega_22_ = Omega[2];
            Omega_01_ = Omega[3];
            Omega_12_ = Omega[4];
            Omega_20_ = Omega[5];
            Epsid_00_ = Epsid[0];
            Epsid_11_ = Epsid[1];
            Epsid_22_ = Epsid[2];
            Epsid_01_ = Epsid[3];
            Epsid_12_ = Epsid[4];
            Epsid_20_ = Epsid[5];

            stnO_[0] = Omega[0];
            stnO_[1] = Omega[1];
            stnO_[2] = Omega[2];
            stnO_[3] = Omega[3];
            stnO_[4] = Omega[4];
            stnO_[5] = Omega[5];
            /*
            Omega[0] = Omega_00_;
            Omega[1] = Omega_11_;
            Omega[2] = Omega_22_;
            Omega[3] = Omega_01_;
            Omega[4] = Omega_12_;
            Omega[5] = Omega_20_;
            */
            //effectiveStiffness(Matdom, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
            //                   C0_, C1_, alpha_);
            

    }

    void Modeldsid::printVariables(){
        cout << " E0       = " << E0_       << endl;
        cout << " Poisson0 = " << Poisson0_ << endl;
        cout << " a1       = " << a1_       << endl;
        cout << " a2       = " << a2_       << endl;
        cout << " a3       = " << a3_       << endl;
        cout << " a4       = " << a4_       << endl;
        cout << " C0       = " << C0_       << endl;
        cout << " C1       = " << C1_       << endl;
        cout << " alpha    = " << alpha_    << endl;
        cout << " Debug    = " << Debug_    << endl;

    }

void effectiveStiffness(r2Tensor<double> &Matdom, const r1Tensor<double> &Omega, const double & E0_,
                        const double & Poisson0_, const double & a1_, const double & a2_,
                        const double & a3_, const double & a4_, const double & C0_,
                        const double & C1_, const double & alpha_){
     int ntens = Matdom.dim1();
     double zero = 0;
     double trOmega, b1, b2, coe1, coe2;
     r2Tensor<double> MatS(ntens,ntens,zero);
     Matdom = MatS;

       trOmega = Omega[0] + Omega[1] + Omega[2];
       b1 = (1. + Poisson0_)/E0_/2.;
       b2 = Poisson0_/E0_;

       coe1 = 2.;
       coe2 = 4.;

       MatS[0][0]=2.*b1-b2+2.*a1_*trOmega+2.*(a2_+a3_)*Omega[0]+2.*a4_*trOmega;
       MatS[0][1]=-b2+2.*a1_*trOmega+a3_*(Omega[0]+Omega[1]);
       MatS[0][2]=-b2+2.*a1_*trOmega+a3_*(Omega[2]+Omega[0]);
       MatS[0][3]=coe1*(a2_*Omega[3]+a3_*Omega[3]);

       MatS[1][0]=MatS[0][1];
       MatS[1][1]=2.*b1-b2+2.*a1_*trOmega+2.*(a2_+a3_)*Omega[1]+2.*a4_*trOmega;
       MatS[1][2]=-b2+2.*a1_*trOmega+a3_*(Omega[2]+Omega[1]);
       MatS[1][3]=coe1*(a2_*Omega[3]+a3_*Omega[3]);

       MatS[2][0]=MatS[0][2];
       MatS[2][1]=MatS[1][2];
       MatS[2][2]=2.*b1-b2+2.*a1_*trOmega+2.*(a2_+a3_)*Omega[2]+2.*a4_*trOmega;
       MatS[2][3]=coe1*(a3_*Omega[3]);

       MatS[3][0]=MatS[0][3];
       MatS[3][1]=MatS[1][3];      
       MatS[3][2]=MatS[2][3];       
       MatS[3][3]=coe2*(b1+0.5*a2_*(Omega[0]+Omega[1])+a4_*trOmega);
 
       
//       IF (NTENS.EQ.6] THEN
         MatS[0][4]=coe1*(a3_*Omega[4]);
         MatS[0][5]=coe1*(a2_*Omega[5]+a3_*Omega[5]);

         MatS[1][4]=coe1*(a2_*Omega[4]+a3_*Omega[4]);
         MatS[1][5]=coe1*(a3_*Omega[5]);

         MatS[2][4]=coe1*(a2_*Omega[4]+a3_*Omega[4]);
         MatS[2][5]=coe1*(a2_*Omega[5]+a3_*Omega[5]);

         MatS[3][4]=coe2*0.5*a2_*Omega[5];
         MatS[3][5]=coe2*0.5*a2_*Omega[4];

         MatS[4][0]=MatS[0][4];
         MatS[4][1]=MatS[1][4];      
         MatS[4][2]=MatS[2][4];       
         MatS[4][3]=MatS[3][4];
         MatS[4][4]=coe2*(b1+0.5*a2_*(Omega[2]+Omega[1])+a4_*trOmega);
         MatS[4][5]=coe2*0.5*a2_*Omega[3];

         MatS[5][0]=MatS[0][5];
         MatS[5][1]=MatS[1][5];       
         MatS[5][2]=MatS[2][5];       
         MatS[5][3]=MatS[3][5];
         MatS[5][4]=MatS[4][5];
         MatS[5][5]=coe2*(b1+0.5*a2_*(Omega[2]+Omega[0])+a4_*trOmega);


//       ENDIF

        gaussj(MatS,Matdom);
     
}

void damageFunction(double & fd, const r1Tensor<double> &Sigma, const r1Tensor<double> &Omega, const double & E0_,
                   const double & Poisson0_, const double & a1_, const double & a2_,
                   const double & a3_, const double & a4_, const double & C0_,
                   const double & C1_, const double & alpha_, const int & ioptfd) {

     int ntens = Sigma.size();
     double zero = 0;
     double trSigma, trOmega, trSigSig, trY, SS;
     r1Tensor<double> SigSig(ntens,zero),P1Y(ntens,zero),sij(ntens,zero),e1(ntens,zero),yd1(ntens,zero);
     r2Tensor<double> P1(ntens,ntens,zero);
     for (int i=0; i<3; i++) e1[i]=1;
     
     trSigma = Sigma[0] + Sigma[1] + Sigma[2];
     trOmega = Omega[0] + Omega[1] + Omega[2];
     //cout << " trSigma = " << trSigma << endl;
     //cout << " " << Omega[0] << " " << Omega[1] << " " << Omega[2] <<  endl;
     Aik_Bkj(Sigma,Sigma,SigSig);

     trSigSig = SigSig[0] + SigSig[1] + SigSig[2];

     for (int i=0; i<ntens; i++)
          yd1[i]= a1_*trSigma*trSigma*e1[i]+a2_*SigSig[i]+a3_*trSigma*Sigma[i]+a4_*trSigSig*e1[i];
    /* 
    cout << " yd1= ";
    for (int i=0; i<6; i++) {
        cout << "  " << yd1[i];
    }
    cout << endl;
    */

    matP1(P1,Sigma);
    /*
    cout << " P1= ";
    for (int i=0; i<6; i++) {
        for (int j=0; j<6; j++) {
            cout << "  " << P1[i][j];
        }
    }
    cout << endl;
    */
     for (int i=0; i<ntens; ++i) {
            P1Y[i] = 0.;
         for (int j=0; j<ntens; j++) {
             if (j<3) {
                P1Y[i] += P1[i][j]*yd1[j];
             } else {
                P1Y[i] += 2.*P1[i][j]*yd1[j];
             }  
         }
     }

     trY = P1Y[0] + P1Y[1] + P1Y[2];

     for (int i=0; i<ntens; ++i)
          sij[i] = P1Y[i]-1./3.*trY*e1[i];


     SS=0.;
     for (int i=0; i<ntens; ++i) {
         if (i<3) {
            SS += sij[i]*sij[i];
         } else {
            SS += 2.*sij[i]*sij[i];
         }  
     }

    /*
    fd = alpha_*trY-C0_-C1_*trOmega;
    if (fd>0.) {
        throwout('DSID:  Stresses are in tension');
        cout << "IOPT = " << IOPTfd << endl;
        cout << "Sigma = " << Sigma[1] << " " << Sigma[2] << " " << Sigma[3] << " " << Sigma[4] << endl;
    } 
    */
    //THE SIGN BEFORE alpha_ IS "+" DUE TO MECHANICAL CONVENTION
    fd = sqrt(0.5*SS)+alpha_*trY-C0_-C1_*trOmega;   
}

void matP1(r2Tensor<double> &P1, const r1Tensor<double> &Sigma) {
    int ntens = Sigma.size();
    int m=3;
    r1Tensor<double> s(m), anan1(ntens), anan2(ntens), anan3(ntens);
    r2Tensor<double> sigma1(m,m), sigma0(m,m);
    const double one=1.,two=2.,tol=1e-6;
    vectorToTensor(Sigma, sigma1, one);
    sigma0 = sigma1;
    Unsymmeig h(sigma0);
    for (int i=0; i<m; i++) {
        if (abs(h.wri[i].real())<tol) h.wri[i].real(0.);
    }
    for (int i=0; i<m; i++) {
        if (h.wri[i].real()>=0){
            s[i] = 1.;
        } else {
            s[i] = -1.;
        }
    }
    for (int i=0; i<m; i++) {
        anan1[i] = h.zz[i][0]*h.zz[i][0];
        anan2[i] = h.zz[i][1]*h.zz[i][1];
        anan3[i] = h.zz[i][2]*h.zz[i][2];
    }

    anan1[3] = h.zz[0][0]*h.zz[1][0];
    anan2[3] = h.zz[0][1]*h.zz[1][1];
    anan3[3] = h.zz[0][2]*h.zz[1][2];

    anan1[4] = h.zz[1][0]*h.zz[2][0];
    anan2[4] = h.zz[1][1]*h.zz[2][1];
    anan3[4] = h.zz[1][2]*h.zz[2][2];

    anan1[5] = h.zz[2][0]*h.zz[0][0];
    anan2[5] = h.zz[2][1]*h.zz[0][1];
    anan3[5] = h.zz[2][2]*h.zz[0][2];


    for (int i=0; i<ntens; i++) {
        for (int j=0; j<ntens; j++) {
            P1[i][j] = s[0] * anan1[i]*anan1[j] +
                       s[1] * anan2[i]*anan2[j] +
                       s[2] * anan3[i]*anan3[j];
        }
    }
}

void matP2(r2Tensor<double> &P2, const r1Tensor<double> &Sigma) {
    int ntens = Sigma.size();
    int m=3;
    double res;
    r1Tensor<double> s(m), anan1(ntens), anan2(ntens), anan3(ntens);
    r2Tensor<double> sigma1(m,m), sigma0(m,m);
    const double one=1.,two=2.,tol=1e-6;
    vectorToTensor(Sigma, sigma1, one);
    sigma0 = sigma1;
    Unsymmeig h(sigma0);
    //for (int i=0; i<m; i++) {
    //    if (abs(h.wri[i].real())<tol) h.wri[i].real()=0.;
    //}
    for (int i=0; i<m; i++) {
        res = h.wri[i].real()- MIN(h.wri[0].real(),MIN(h.wri[1].real(),h.wri[2].real()));
        if (res>0.){
            s[i] = 1.;
        } else {
            s[i] = 0.;
        }
        if (abs(res)<tol) s[i] = 0.;
    }
    for (int i=0; i<m; i++) {
        anan1[i] = h.zz[i][0]*h.zz[i][0];
        anan2[i] = h.zz[i][1]*h.zz[i][1];
        anan3[i] = h.zz[i][2]*h.zz[i][2];
    }

    anan1[3] = h.zz[0][0]*h.zz[1][0];
    anan2[3] = h.zz[0][1]*h.zz[1][1];
    anan3[3] = h.zz[0][2]*h.zz[1][2];

    anan1[4] = h.zz[1][0]*h.zz[2][0];
    anan2[4] = h.zz[1][1]*h.zz[2][1];
    anan3[4] = h.zz[1][2]*h.zz[2][2];

    anan1[5] = h.zz[2][0]*h.zz[0][0];
    anan2[5] = h.zz[2][1]*h.zz[0][1];
    anan3[5] = h.zz[2][2]*h.zz[0][2];

    for (int i=0; i<ntens; i++) {
        for (int j=0; j<ntens; j++) {
            P2[i][j] = s[0] * anan1[i]*anan1[j] +
                       s[1] * anan2[i]*anan2[j] +
                       s[2] * anan3[i]*anan3[j];
        }
    }
}

void dY_dSigFunction(const r1Tensor<double> &Sigma, r2Tensor<double> &dY_dSig,const double & E0_,
             const double & Poisson0_, const double & a1_, const double & a2_,
             const double & a3_, const double & a4_, const double & C0_,
             const double & C1_, const double & alpha_) {
    int ntens = Sigma.size();
    double trSig;
    trSig = Sigma[0]+Sigma[1]+Sigma[2];
    dY_dSig[0][0]=2.*a1_*trSig+2.*a2_*Sigma[0]+a3_*(Sigma[0]+trSig)+2.*a4_*Sigma[0];
    dY_dSig[0][1]=2.*a1_*trSig+a3_*Sigma[0]+2.*a4_*Sigma[1];
    dY_dSig[0][2]=2.*a1_*trSig+a3_*Sigma[0]+2.*a4_*Sigma[2];
    dY_dSig[0][3]=a2_*Sigma[3]+2.*a4_*Sigma[3];
                                                                                   
    dY_dSig[1][0]=2.*a1_*trSig+a3_*Sigma[1]+2.*a4_*Sigma[0];
    dY_dSig[1][1]=2.*a1_*trSig+2.*a2_*Sigma[1]+a3_*(Sigma[1]+trSig)+2.*a4_*Sigma[1];
    dY_dSig[1][2]=2.*a1_*trSig+a3_*Sigma[1]+2.*a4_*Sigma[2];
    dY_dSig[1][3]=a2_*Sigma[3]+2.*a4_*Sigma[3];
                                                                                   
    dY_dSig[2][0]=2.*a1_*trSig+a3_*Sigma[2]+2.*a4_*Sigma[0];
    dY_dSig[2][1]=2.*a1_*trSig+a3_*Sigma[2]+2.*a4_*Sigma[1];
    dY_dSig[2][2]=2.*a1_*trSig+2.*a2_*Sigma[2]+a3_*(Sigma[2]+trSig)+2.*a4_*Sigma[2];
    dY_dSig[2][3]=2.*a4_*Sigma[3];
                                                                                   
    dY_dSig[3][0]=a2_*Sigma[3]+a3_*Sigma[3];
    dY_dSig[3][1]=a2_*Sigma[3]+a3_*Sigma[3];
    dY_dSig[3][2]=a3_*Sigma[3];
    dY_dSig[3][3]=0.5*a2_*(Sigma[0]+Sigma[1])+0.5*a3_*trSig;
                                                                                   
    dY_dSig[0][4]=2.*a4_*Sigma[4];
    dY_dSig[0][5]=a2_*Sigma[5]+2.*a4_*Sigma[5];
                                                                                   
    dY_dSig[1][4]=a2_*Sigma[4]+2.*a4_*Sigma[4];
    dY_dSig[1][5]=2.*a4_*Sigma[5];
                                                                                   
    dY_dSig[2][4]=a2_*Sigma[4]+2.*a4_*Sigma[4];
    dY_dSig[2][5]=a2_*Sigma[5]+2.*a4_*Sigma[5];
                                                                                   
    dY_dSig[3][4]=0.5*a2_*Sigma[5];
    dY_dSig[3][5]=0.5*a2_*Sigma[4];
                                                                                   
    dY_dSig[4][0]=a3_*Sigma[4];
    dY_dSig[4][1]=a2_*Sigma[4]+a3_*Sigma[4];
    dY_dSig[4][2]=a2_*Sigma[4]+a3_*Sigma[4];
    dY_dSig[4][3]=dY_dSig[3][4];
    dY_dSig[4][4]=0.5*a2_*(Sigma[1]+Sigma[2])+0.5*a3_*trSig;
    dY_dSig[4][5]=0.5*a2_*Sigma[3];
                                                                                   
    dY_dSig[5][0]=a2_*Sigma[5]+a3_*Sigma[5];
    dY_dSig[5][1]=a3_*Sigma[5];
    dY_dSig[5][2]=a2_*Sigma[5]+a3_*Sigma[5];
    dY_dSig[5][3]=dY_dSig[3][5];
    dY_dSig[5][4]=dY_dSig[4][5];
    dY_dSig[5][5]=0.5*a2_*(Sigma[2]+Sigma[0])+0.5*a3_*trSig;
}

// SUB FD_LAM
void cuttingPlaneMethod(const r1Tensor<double> &Omega, const r1Tensor<double> &Sigma, r2Tensor<double> &Matdom,const double & E0_,
             const double & Poisson0_, const double & a1_, const double & a2_,
             const double & a3_, const double & a4_, const double & C0_,
             const double & C1_, const double & alpha_, double & H0,
             double & HP, r1Tensor<double> &dG_dY, r1Tensor<double> &df_dSig,
             r1Tensor<double> &temp, const int & iopt){
    int ntens = Omega.size();
    double zero=0.;
    double trSigma, trOmega, trSigSig, P1yd1e, f1f1, f2f2;
    r1Tensor<double> e(ntens), yd1(ntens), SigSig(ntens), f1ij(ntens),f2ij(ntens);
    r1Tensor<double> P1yd1(ntens), f2p2(ntens), ep1(ntens),df_dY(ntens), df_dOmega(ntens);
    r1Tensor<double> f1p3(ntens),zeros(ntens,zero),temp1(ntens); 
    r2Tensor<double> P1(ntens,ntens), dY_dSig(ntens,ntens),eep1(ntens,ntens);
    r2Tensor<double> P3(ntens,ntens), P2(ntens,ntens), temp2(ntens,ntens);
    r3Tensor<double> dS_dOmega(ntens,ntens,ntens,0.); 

    //DATA ONE,TWO,HALF / 1.0D0,2.0D0,0.5D0 /

    for (int i=0; i<ntens; i++) {
        if (i<3) {
            e[i]=1.;
        } else {
            e[i]=0.;
        }
    }

    trSigma = Sigma[0]+Sigma[1]+Sigma[2];
    trOmega = Omega[0]+Omega[1]+Omega[2];

    Aik_Bkj(Sigma,Sigma,SigSig);

    trSigSig =SigSig[0]+SigSig[1]+SigSig[2];

    for (int i=0; i<ntens; i++) {
        yd1[i]= a1_*trSigma*trSigma*e[i]+a2_*SigSig[i]+
              a3_*trSigma*Sigma[i]+a4_*trSigSig*e[i];
    }

    dY_dSigFunction(Sigma,dY_dSig,
            E0_,Poisson0_,a1_,a2_,a3_,a4_,C0_,C1_,alpha_);
    matP1(P1,Sigma);

    matP2(P2,Sigma);
    for (int i=0;i<ntens;i++){
        P1yd1[i]=0.;
        for (int j=0;j<ntens;j++){
            if (j<3) {
                P1yd1[i]+=P1[i][j]*yd1[j];
            } else {
                P1yd1[i]+=2.*P1[i][j]*yd1[j];
            }
         }   
    }
    P1yd1e=P1yd1[0]+P1yd1[1]+P1yd1[2];
    //cout << " P1yd1e = " << P1yd1e << endl;
    //cout << " yd1= ";
    //cout << endl;
    for (int i=0; i<ntens; i++) {
        f2ij[i]=0.;
        for (int j=0; j<ntens; j++) {
            if (j<3) {
                f2ij[i]+=P2[i][j]*yd1[j];
            } else {
                f2ij[i]+=2.*P2[i][j]*yd1[j];
            }
         }
     }
    /*
    for (int i=0; i<ntens; i++){
        for(int j=0; j<ntens; j++){
            cout << " " << P2[i][j];
        }
    }
    cout << endl;
    */
    for (int i=0; i<ntens; i++){
        f2p2[i]=0.;
        for (int j=0;j<ntens;j++){
            if (j<3) {   
                f2p2[i]+=P2[j][i]*f2ij[j];
            } else {
                f2p2[i]+=2.*P2[j][i]*f2ij[j];
            }  
         }
     }
      
     f2f2=0.;
     for (int i=0; i<ntens; i++) {
         if (i<3) {
             f2f2+=f2ij[i]*f2ij[i];
         } else {
             f2f2+=2.*f2ij[i]*f2ij[i];
         }  
     }

    for (int i=0; i<ntens; i++) {
        ep1[i]=0.;
        for (int j=0; j<ntens; j++) {
            if (j<3) {
                ep1[i]+=P1[j][i]*e[j];
            } else {
                ep1[i]+=2.*P1[j][i]*e[j];
            }  
        }
    }
     
    for (int i=0; i<ntens; i++) {
        for (int j=0; j<ntens; j++) {
            eep1[i][j]=e[i]*ep1[j];
        }
    }

    for (int i=0; i<ntens; i++) {
        for (int j=0; j<ntens; j++) {
            P3[i][j]=P1[i][j]-1./3.*eep1[i][j];
        }
    }

    for (int i=0; i<ntens; i++) {
        f1ij[i]=P1yd1[i]-1./3.*P1yd1e*e[i];
        df_dOmega[i]=-C1_*e[i];
        if (f2f2 == 0.) {
            dG_dY[i]=0.;
        } else {  
            dG_dY[i]=f2p2[i]/sqrt(2.0*f2f2);
        }
    }
    //cout << " dG_dY = " << dG_dY[0] << " " << dG_dY[1] << " " << dG_dY[2] << endl;
    //cout << " f2ij = " << f2ij[0] << " " << f2ij[1] << " " << f2ij[2] << endl;
    //cout << " f2f2 = " << f2f2 << endl;

      
    for (int i=0; i<ntens; i++) {
        f1p3[i]=0.;
        for (int j=0; j<ntens; j++) {
            if (j<3) {
                f1p3[i]+=f1ij[j]*P3[j][i];
            } else {
                f1p3[i]+=2.*f1ij[j]*P3[j][i];
            }  
        }
    }
      
    f1f1=0.;
    for (int i=0; i<ntens; i++) {
        if (i<3) {
            f1f1+=f1ij[i]*f1ij[i];
        } else {
            f1f1+=2.*f1ij[i]*f1ij[i];
        }  
    }

    for (int i=0; i<ntens; i++) {
        df_dY[i]=f1p3[i]/sqrt(2.*f1f1)+alpha_*ep1[i];
    }
     
    for (int i=0; i<ntens; i++) {
        df_dSig[i]=0.;
        for (int j=0; j<ntens; j++) {
            if (j<3) {
                df_dSig[i]+=df_dY[j]*dY_dSig[j][i];
            } else {
                df_dSig[i]+=2.*df_dY[j]*dY_dSig[j][i];
            }  
        }
    }

      
    HP=0.;
    for (int i=0; i<ntens; i++) {
        if(i<3){
            HP+=df_dOmega[i]*dG_dY[i];
        } else {
            HP+=2.*df_dOmega[i]*dG_dY[i];
        } 
    }

    H0=0.;

    Mat_dS_dOmega(dS_dOmega,E0_,Poisson0_,a1_,a2_,a3_,a4_,C0_,C1_,alpha_);

    for (int i=0; i<ntens; i++) {
        for (int j=0; j<ntens; j++) {
            temp2[i][j]=0.;
            for (int k=0; k<ntens; k++) {
                if (k<3) {
                    temp2[i][j]=temp2[i][j]+Sigma[k]*dS_dOmega[k][i][j];
                } else {
                    temp2[i][j]=temp2[i][j]+2.*Sigma[k]*dS_dOmega[k][i][j];
                }
            }
        }
    }

    for (int i=0; i<ntens; i++) {
        temp1[i]=0.;
        for (int j=0; j<ntens; j++) {
            if (j<3){
                temp1[i]+=temp2[i][j]*dG_dY[j];
            } else {
                temp1[i]+=2.0*temp2[i][j]*dG_dY[j];
            }
        }
    }

    for (int i=0; i<ntens; i++) {
         temp1[i]+=df_dSig[i];
    }
     
    if(iopt == 1) {
      effectiveStiffness(Matdom,Omega,
               E0_,Poisson0_,a1_,a2_,a3_,a4_,C0_,C1_,alpha_);
    }
    for (int i=0; i<ntens; i++) {
        temp[i]=0.;
        for (int j=0; j<ntens; j++) {
            if(j<3){
                temp[i]+=Matdom[i][j]*temp1[j];
            } else {
                temp[i]+=2.*Matdom[i][j]*temp1[j];
            } 
        }
    }
    for (int i=0; i<ntens; i++) {
        if(i<3){
            H0+=df_dSig[i]*temp[i];
        } else {
            H0+=2.*df_dSig[i]*temp[i];
        } 
    }
}

void Mat_dS_dOmega(r3Tensor<double> &dS_dO, const double & E0_,
                   const double & Poisson0_, const double & a1_, const double & a2_,
                   const double & a3_, const double & a4_, const double & C0_,
                   const double & C1_, const double & alpha_){

// ------- page 1 -------
      dS_dO[0][0][0]=2.*a1_+2.*a2_+2.*a3_+2.*a4_;  
      dS_dO[0][1][0]=2.*a1_+a3_                 ;
      dS_dO[0][2][0]=2.*a1_+a3_                 ;
      dS_dO[0][3][0]=0.                         ;
                                                  
      dS_dO[1][0][0]=dS_dO[0][1][0]             ;
      dS_dO[1][1][0]=2.*a1_+2.*a4_              ;
      dS_dO[1][2][0]=2.*a1_                     ;
      dS_dO[1][3][0]=0.                         ;
                                                  
      dS_dO[2][0][0]=dS_dO[0][2][0]             ;
      dS_dO[2][1][0]=dS_dO[1][2][0]             ;
      dS_dO[2][2][0]=2.*a1_+2.*a4_              ;
      dS_dO[2][3][0]=0.                         ;
                                                  
      dS_dO[3][0][0]=dS_dO[0][3][0]             ;
      dS_dO[3][1][0]=dS_dO[1][3][0]             ;
      dS_dO[3][2][0]=dS_dO[2][3][0]             ;
      dS_dO[3][3][0]=0.5*a2_+a4_                ;
                                                  
//      IF[NTENS.EQ.6]THEN                        
      dS_dO[0][4][0]=0.                         ;
      dS_dO[0][5][0]=0.                         ;
      dS_dO[1][4][0]=0.                         ;
      dS_dO[1][5][0]=0.                         ;
      dS_dO[2][4][0]=0.                         ;
      dS_dO[2][5][0]=0.                         ;
      dS_dO[3][4][0]=0.                         ;
      dS_dO[3][5][0]=0.                         ;
      dS_dO[4][0][0]=dS_dO[0][4][0]             ;
      dS_dO[4][1][0]=dS_dO[1][4][0]             ;
      dS_dO[4][2][0]=dS_dO[2][4][0]             ;
      dS_dO[4][3][0]=dS_dO[3][4][0]             ;
      dS_dO[4][4][0]=a4_                        ;
      dS_dO[4][5][0]=0.                         ;
                                                  
      dS_dO[5][0][0]=dS_dO[0][5][0]             ;
      dS_dO[5][1][0]=dS_dO[1][5][0]             ;
      dS_dO[5][2][0]=dS_dO[2][5][0]             ;
      dS_dO[5][3][0]=dS_dO[3][5][0]             ;
      dS_dO[5][4][0]=dS_dO[4][5][0]             ;
      dS_dO[5][5][0]=0.5*a2_+a4_                ;
//      ENDIF                                     
                                                  
// ------- page 2 -------                       --
      dS_dO[0][0][1]=2.*a1_+2.*a4_              ;
      dS_dO[0][1][1]=2.*a1_+a3_                 ;
      dS_dO[0][2][1]=2.*a1_                     ;
      dS_dO[0][3][1]=0.                         ;
                                                  
      dS_dO[1][0][1]=dS_dO[0][1][1]             ;
      dS_dO[1][1][1]=2.*a1_+2.*a2_+2.*a3_+2.*a4_;
      dS_dO[1][2][1]=2.*a1_+a3_                 ;
      dS_dO[1][3][1]=0.                         ;
                                                  
      dS_dO[2][0][1]=dS_dO[0][2][1]             ;
      dS_dO[2][1][1]=dS_dO[1][2][1]             ;
      dS_dO[2][2][1]=2.*a1_+2.*a4_              ;
      dS_dO[2][3][1]=0.                         ;
                                                  
      dS_dO[3][0][1]=dS_dO[0][3][1]             ;
      dS_dO[3][1][1]=dS_dO[1][3][1]             ;
      dS_dO[3][2][1]=dS_dO[2][3][1]             ;
      dS_dO[3][3][1]=0.5*a2_+a4_                ;
                                                  
//      IF[NTENS.EQ.6]THEN                        
      dS_dO[0][4][1]=0.                         ;
      dS_dO[0][5][1]=0.                         ;
      dS_dO[1][4][1]=0.                         ;
      dS_dO[1][5][1]=0.                         ;
      dS_dO[2][4][1]=0.                         ;
      dS_dO[2][5][1]=0.                         ;
      dS_dO[3][4][1]=0.                         ;
      dS_dO[3][5][1]=0.                         ;
      dS_dO[4][0][1]=dS_dO[0][4][1]             ;
      dS_dO[4][1][1]=dS_dO[1][4][1]             ;
      dS_dO[4][2][1]=dS_dO[2][4][1]             ;
      dS_dO[4][3][1]=dS_dO[3][4][1]             ;
      dS_dO[4][4][1]=0.5*a2_+a4_                ;
      dS_dO[4][5][1]=0.                         ;
                                                  
      dS_dO[5][0][1]=dS_dO[0][5][1]             ;
      dS_dO[5][1][1]=dS_dO[1][5][1]             ;
      dS_dO[5][2][1]=dS_dO[2][5][1]             ;
      dS_dO[5][3][1]=dS_dO[3][5][1]             ;
      dS_dO[5][4][1]=dS_dO[4][5][1]             ;
      dS_dO[5][5][1]=0.5*a2_+a4_                ;
//      ENDIF                                     
                                                  
// ------- page 3 -------                       --
      dS_dO[0][0][2]=2.*a1_+2.*a4_              ;
      dS_dO[0][1][2]=2.*a1_                     ;
      dS_dO[0][2][2]=2.*a1_+a3_                 ;
      dS_dO[0][3][2]=0.                         ;
                                                  
      dS_dO[1][0][2]=dS_dO[0][1][2]             ;
      dS_dO[1][1][2]=2.*a1_+2.*a4_              ;
      dS_dO[1][2][2]=2.*a1_+a3_                 ;
      dS_dO[1][3][2]=0.                         ;
                                                  
      dS_dO[2][0][2]=dS_dO[0][2][2]             ;
      dS_dO[2][1][2]=dS_dO[1][2][2]             ;
      dS_dO[2][2][2]=2.*a1_+2.*a2_+2.*a3_+2.*a4_;
      dS_dO[2][3][2]=0.                         ;
                                                  
      dS_dO[3][0][2]=dS_dO[0][3][2]             ;
      dS_dO[3][1][2]=dS_dO[1][3][2]             ;
      dS_dO[3][2][2]=dS_dO[2][3][2]             ;
      dS_dO[3][3][2]=a4_                        ;
                                                  
//      IF[NTENS.EQ.6]THEN                        
      dS_dO[0][4][2]=0.                         ;
      dS_dO[0][5][2]=0.                         ;
      dS_dO[1][4][2]=0.                         ;
      dS_dO[1][5][2]=0.                         ;
      dS_dO[2][4][2]=0.                         ;
      dS_dO[2][5][2]=0.                         ;
      dS_dO[3][4][2]=0.                         ;
      dS_dO[3][5][2]=0.                         ;
      dS_dO[4][0][2]=dS_dO[0][4][2]             ;
      dS_dO[4][1][2]=dS_dO[1][4][2]             ;
      dS_dO[4][2][2]=dS_dO[2][4][2]             ;
      dS_dO[4][3][2]=dS_dO[3][4][2]             ;
      dS_dO[4][4][2]=0.5*a2_+a4_                ;
      dS_dO[4][5][2]=0.                         ;
                                                  
      dS_dO[5][0][2]=dS_dO[0][5][2]             ;
      dS_dO[5][1][2]=dS_dO[1][5][2]             ;
      dS_dO[5][2][2]=dS_dO[2][5][2]             ;
      dS_dO[5][3][2]=dS_dO[3][5][2]             ;
      dS_dO[5][4][2]=dS_dO[4][5][2]             ;
      dS_dO[5][5][2]=0.5*a2_+a4_                ;
//      ENDIF                                     
                                                  
// ------- page 4 -------                       --
      dS_dO[0][0][3]=0.                         ;
      dS_dO[0][1][3]=0.                         ;
      dS_dO[0][2][3]=0.                         ;
      dS_dO[0][3][3]=0.5*a2_+0.5*a3_            ;
                                                  
      dS_dO[1][0][3]=dS_dO[0][1][3]             ;
      dS_dO[1][1][3]=0.                         ;
      dS_dO[1][2][3]=0.                         ;
      dS_dO[1][3][3]=0.5*a2_+0.5*a3_            ;
                                                  
      dS_dO[2][0][3]=dS_dO[0][2][3]             ;
      dS_dO[2][1][3]=dS_dO[1][2][3]             ;
      dS_dO[2][2][3]=0.                         ;
      dS_dO[2][3][3]=0.5*a3_                    ;
                                                  
      dS_dO[3][0][3]=dS_dO[0][3][3]             ;
      dS_dO[3][1][3]=dS_dO[1][3][3]             ;
      dS_dO[3][2][3]=dS_dO[2][3][3]             ;
      dS_dO[3][3][3]=0.                         ;
                                                  
//      IF[NTENS.EQ.6]THEN                        
      dS_dO[0][4][3]=0.                         ;
      dS_dO[0][5][3]=0.                         ;
      dS_dO[1][4][3]=0.                         ;
      dS_dO[1][5][3]=0.                         ;
      dS_dO[2][4][3]=0.                         ;
      dS_dO[2][5][3]=0.                         ;
      dS_dO[3][4][3]=0.                         ;
      dS_dO[3][5][3]=0.                         ;
      dS_dO[4][0][3]=dS_dO[0][4][3]             ;
      dS_dO[4][1][3]=dS_dO[1][4][3]             ;
      dS_dO[4][2][3]=dS_dO[2][4][3]             ;
      dS_dO[4][3][3]=dS_dO[3][4][3]             ;
      dS_dO[4][4][3]=0.                         ;
      dS_dO[4][5][3]=0.25*a2_                   ;
                                                  
      dS_dO[5][0][3]=dS_dO[0][5][3]             ;
      dS_dO[5][1][3]=dS_dO[1][5][3]             ;
      dS_dO[5][2][3]=dS_dO[2][5][3]             ;
      dS_dO[5][3][3]=dS_dO[3][5][3]             ;
      dS_dO[5][4][3]=dS_dO[4][5][3]             ;
      dS_dO[5][5][3]=0.                         ;
//      ENDIF                                     
                                                  
//      IF[NTENS.EQ.6]THEN                        
// ------- page 5 -------                       --
      dS_dO[0][0][4]=0.                         ;
      dS_dO[0][1][4]=0.                         ;
      dS_dO[0][2][4]=0.                         ;
      dS_dO[0][3][4]=0.                         ;
      dS_dO[0][4][4]=0.5*a3_                    ;
      dS_dO[0][5][4]=0.                         ;
                                                  
      dS_dO[1][0][4]=dS_dO[0][1][4]             ;
      dS_dO[1][1][4]=0.                         ;
      dS_dO[1][2][4]=0.                         ;
      dS_dO[1][3][4]=0.                         ;
      dS_dO[1][4][4]=0.5*a2_+0.5*a3_            ;
      dS_dO[1][5][4]=0.                         ;
                                                  
      dS_dO[2][0][4]=dS_dO[0][2][4]             ;
      dS_dO[2][1][4]=dS_dO[1][2][4]             ;
      dS_dO[2][2][4]=0.                         ;
      dS_dO[2][3][4]=0.                         ;
      dS_dO[2][4][4]=0.5*a2_+0.5*a3_            ;
      dS_dO[2][5][4]=0.                         ;
                                                  
      dS_dO[3][0][4]=dS_dO[0][3][4]             ;
      dS_dO[3][1][4]=dS_dO[1][3][4]             ;
      dS_dO[3][2][4]=dS_dO[2][3][4]             ;
      dS_dO[3][3][4]=0.                         ;
      dS_dO[3][4][4]=0.                         ;
      dS_dO[3][5][4]=0.25*a2_                   ;
                                                  
      dS_dO[4][0][4]=dS_dO[0][4][4]             ;
      dS_dO[4][1][4]=dS_dO[1][4][4]             ;
      dS_dO[4][2][4]=dS_dO[2][4][4]             ;
      dS_dO[4][3][4]=dS_dO[3][4][4]             ;
      dS_dO[4][4][4]=0.                         ;
      dS_dO[4][5][4]=0.                         ;
                                                  
      dS_dO[5][0][4]=dS_dO[0][5][4]             ;
      dS_dO[5][1][4]=dS_dO[1][5][4]             ;
      dS_dO[5][2][4]=dS_dO[2][5][4]             ;
      dS_dO[5][3][4]=dS_dO[3][5][4]             ;
      dS_dO[5][4][4]=dS_dO[4][5][4]             ;
      dS_dO[5][5][4]=0.                         ;
                                                  
// ------- page 6 -------                       --
      dS_dO[0][0][5]=0.                         ;
      dS_dO[0][1][5]=0.                         ;
      dS_dO[0][2][5]=0.                         ;
      dS_dO[0][3][5]=0.                         ;
      dS_dO[0][4][5]=0.                         ;
      dS_dO[0][5][5]=0.5*a2_+0.5*a3_            ;
                                                  
      dS_dO[1][0][5]=dS_dO[0][1][5]             ;
      dS_dO[1][1][5]=0.                         ;
      dS_dO[1][2][5]=0.                         ;
      dS_dO[1][3][5]=0.                         ;
      dS_dO[1][4][5]=0.                         ;
      dS_dO[1][5][5]=0.5*a3_                    ;
                                                  
      dS_dO[2][0][5]=dS_dO[0][2][5]             ;
      dS_dO[2][1][5]=dS_dO[1][2][5]             ;
      dS_dO[2][2][5]=0.                         ;
      dS_dO[2][3][5]=0.                         ;
      dS_dO[2][4][5]=0.                         ;
      dS_dO[2][5][5]=0.5*a2_+0.5*a3_            ;
                                                  
      dS_dO[3][0][5]=dS_dO[0][3][5]             ;
      dS_dO[3][1][5]=dS_dO[1][3][5]             ;
      dS_dO[3][2][5]=dS_dO[2][3][5]             ;
      dS_dO[3][3][5]=0.                         ;
      dS_dO[3][4][5]=0.25*a2_                   ;
      dS_dO[3][5][5]=0.                         ;
                                                  
      dS_dO[4][0][5]=dS_dO[0][4][5]             ;
      dS_dO[4][1][5]=dS_dO[1][4][5]             ;
      dS_dO[4][2][5]=dS_dO[2][4][5]             ;
      dS_dO[4][3][5]=dS_dO[3][4][5]             ;
      dS_dO[4][4][5]=0.                         ;
      dS_dO[4][5][5]=0.                         ;
                                                  
      dS_dO[5][0][5]=dS_dO[0][5][5]             ;
      dS_dO[5][1][5]=dS_dO[1][5][5]             ;
      dS_dO[5][2][5]=dS_dO[2][5][5]             ;
      dS_dO[5][3][5]=dS_dO[3][5][5]             ;
      dS_dO[5][4][5]=dS_dO[4][5][5]             ;
      dS_dO[5][5][5]=0.                         ;
//      ENDIF

}


// EOF
