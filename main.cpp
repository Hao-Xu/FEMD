#include<iostream>
#include<fstream>
#include<iomanip>
#include"Modeldsid.h"
#include "../mathLib/r1Tensor.h"

using namespace std;

int main() {
    double E0, Poisson0, a1, a2, a3, a4, C0, C1, alpha, Debug,fd0;
    r1Tensor<double> stnE(6,0.), stnS(6,0.), stnO(6,0.),stran(6,0.);
    ofstream fstrm;
    E0 = 6.8e10;
    Poisson0 = 0.21;
    a1 = 1.26E-13;
    a2 = 3.94E-11;
    a3 = -1.26E-12;
    a4 = 2.51E-12;
    C0 = 0.E6;
    C1 = 2.2E6;
    alpha = 0.231;
    Debug = 0.;
    fd0   = 0.;
    Modeldsid dsid1(E0, Poisson0, a1, a2, a3, a4, C0, C1, alpha, Debug);
    //dsid1.printVariables();
    fstrm.open("result.out");
    fstrm << "# STEP  STRAN(1)   STRAN(2)   STRAN(3)  STRESS(1)  STRESS(2)  STRESS(1)"
          <<"    OMEGA(1)   OMEGA(2)   OMEGA(3)      fd0" << endl;
    fstrm << setprecision(4) << scientific;
    stnE[2] = 0.00002;
    for (int j=1; j<=1000; j++){
        dsid1.run(stnE,stnS,stnO,fd0);
        //cout << " load " << j << endl;
        fstrm << setw(5) << j;
        for (int i=0; i<stran.size();i++) {
            stran[i] += stnE[i];
        }
        for (int i=0; i<3; i++) {
            fstrm << " " << setw(10) << stran[i] ;
        }
        for (int i=0; i<3; i++) {
            fstrm << " " << setw(10) << stnS[i] ;
        }
        for (int i=0; i<3; i++) {
            fstrm << " " << setw(10) << stnO[i] ;
        }
        fstrm << " " << setw(10) << fd0;
        fstrm << endl;
    }
    return 0;
}
