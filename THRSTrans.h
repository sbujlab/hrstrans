#ifndef __THRSTRANS_H
#define __THRSTRANS_H

#include "TMatrixD.h"
class TH1F;
class TGraph;

class THRSTrans {
    public:
    THRSTrans( double, double, double, double );
    ~THRSTrans();

    TMatrixD *swapxy(double sign = 1.0);
    TMatrixD *makedrift( double l );
    TMatrixD *makequad( double B0, double a, double l );
    TMatrixD *makedip( double theta, double rho, double n, double beta, double beta2, double psi  );

    void DoTransport();
    void PrintSimple(TMatrixD *);

    enum trans_t {kX, kTh, kd, kX2, kXTh, kXd, kTh2, kThd, kd2, kY2, kYPh, kPh2, kY, kPh, kXY, kXPh, kThY, kThPh, kYd, kPhd };

    void ShowOutput();

    TMatrixD *GetTransport(){ return trans[nelm]; }

    TGraph *GetXAcc(){ return gacc[0][0]; }
    TGraph *GetThAcc(){ return gacc[1][0]; }
    TGraph *GetYAcc(){ return gacc[2][0]; }
    TGraph *GetPhAcc(){ return gacc[3][0]; }
    TGraph *GetDpAcc(){ return gacc[4][0]; }

    private:
        double Bq1, Bq2, Bq3,psi;
        TMatrixD *trans[20], *chain[20];

        int nelm;

        int nbin;
        TH1F   *hacc[5][2][20];
        TGraph *gacc[5][20];

        double larr[20];
        double x_x0[20];
        double x_th[20];
        double x_d[20];
        double y_ph[20];
        double y_y0[20];
        double ph_ph[20];
        double ph_y0[20];

        int ntrk;
        double ytg[100000], xtg[100000], thtg[100000], phtg[100000], dp[100000];
        bool acc[100000];

        void setcrossterms(TMatrixD *);
        void fillvector(TVectorD &);

         TGraph *gx_th, *gx_x0, *gx_d, *gy_ph, *gy_y0 , *gph_y0, *gph_ph;
};
#endif//__THRSTRANS_H
