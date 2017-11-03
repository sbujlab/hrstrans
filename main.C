#include "THRSTrans.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv) 
{
    THRSTrans *trans = new THRSTrans(0.099950, -0.132890, -0.171751, 0.050178, 0.037056, THRSTrans::kPREX);

    TMatrixD tovdc= *(trans->GetTransport(14));

    TMatrixD optics = *(trans->GetOptics(14));

    TMatrixD tofp = *(trans->GetTransport(14));
    TMatrixD tosen = *(trans->GetTransport(2));
    TMatrixD tosex = *(trans->GetTransport(3));
    TMatrixD q2ex = *(trans->GetTransport(9));
    TMatrixD den = *(trans->GetTransport(10));
    TMatrixD q3en = *(trans->GetTransport(12));
    TMatrixD q3ex = *(trans->GetTransport(13));

    cout << "Result of THRSTrans "<<endl;
    cout << "Show Output () -----------------------------"<<endl;
    trans->ShowOutput();
    printf("\n\n q3ex  Output \n");
    trans->ShowOutput(13);
    printf("\n\nq3ex  again,  Output \n");
    q3ex.Print();
    cout << "\nShow Acc -----------------------------"<<endl;
    trans->ShowAcc();
    cout << "\nShow Focal Lengths -----------------------------"<<endl;
    trans->ShowFocalLengths();
    cout << "Optics Print -----------------------------"<<endl;
    trans->GetOptics()->Print();
    cout << "----------------- end ----------------- "<<endl;

    TVectorD iv(20),v(20);

    double xtg[100000], ytg[100000], thtg[100000], phtg[100000], dp[100000];

    Int_t ntrk = 4;

    xtg[0] = -1.3e-5; ytg[0] = 1.2e-4; thtg[0] = -0.1176; phtg[0] = 0.0731; dp[0] = -8.68e-5;
    xtg[1] = 5.94e-4; ytg[1] = -7.9e-4; thtg[1] = -0.0914; phtg[1] = -0.04; dp[1] = -0.00607;
    xtg[2] = 2.26e-4; ytg[2] = -4.84e-4; thtg[2] = 0.1438; phtg[2] = 0.0864; dp[2] = -0.00247;
    xtg[3] = -3.69e-4; ytg[3] = -5.54e-4; thtg[3] = 0.0308; phtg[3] = -0.0872; dp[3] = -0.00169;

    for (Int_t itrk = 0; itrk < ntrk; itrk++) {    

      iv[THRSTrans::kX]  =  xtg[itrk];
      iv[THRSTrans::kTh] = thtg[itrk];
      iv[THRSTrans::kY]  =  ytg[itrk];
      iv[THRSTrans::kPh] = phtg[itrk];
      iv[THRSTrans::kd]  =   dp[itrk];

      trans->fillvector(iv);

      v = (tofp)*iv;

      cout << "\n\nTrack number "<<itrk<<endl;
      cout << "input vector << "<<iv[THRSTrans::kX]<<"  "<<iv[THRSTrans::kTh] << "  "<<iv[THRSTrans::kY]<<"  "<<iv[THRSTrans::kPh]<<"  "<<iv[THRSTrans::kd]<<endl;
      cout << "output vector << "<<v[THRSTrans::kX]<<"  "<<v[THRSTrans::kTh] << "  "<<v[THRSTrans::kY]<<"  "<<v[THRSTrans::kPh]<<"  "<<v[THRSTrans::kd]<<endl;

    }


    return 1;

}



