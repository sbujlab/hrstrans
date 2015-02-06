#ifndef __THRSTRANS_H
#define __THRSTRANS_H

#include "TMatrixD.h"
class TH1F;
class TGraph;

class THRSTrans {
    public:
    enum tune_t {kStd, kPREX, kCREX, kAPEX };

    THRSTrans( double, double, double, double, double, tune_t tune = kStd, int n = 100000 );
    ~THRSTrans();

    TMatrixD *swapxy(double sign = 1.0);
    TMatrixD *makedrift( double l );
    TMatrixD *makequad( double B0, double a, double l );
    TMatrixD *makedip( double theta, double rho, double n, double beta, double beta2, double ps1, double ps2);

    void DoTransport();
    void PrintSimple(TMatrixD *, double l = 0.0);

    enum trans_t {kX, kTh, kd, kX2, kXTh, kXd, kTh2, kThd, kd2, kY2, kYPh, kPh2, kY, kPh, kXY, kXPh, kThY, kThPh, kYd, kPhd };

    void ShowOutput(int i = -1, double l = 0.0);
    void ShowFocalLengths();
    void ShowAcc();

    TMatrixD *GetTransport(){ return trans[nelm]; }
    TMatrixD *GetTransport(int i){ return trans[i]; }
    TMatrixD *GetElement(int i){ return chain[i]; }
    TMatrixD *GetOptics(int idx = -1);

    TGraph *GetXAcc(){ return gacc[0][0]; }
    TGraph *GetThAcc(){ return gacc[1][0]; }
    TGraph *GetYAcc(){ return gacc[2][0]; }
    TGraph *GetPhAcc(){ return gacc[3][0]; }
    TGraph *GetDpAcc(){ return gacc[4][0]; }


    void SetSeptumPsi( double k1, double k2 ){ sept_K1 = k1; sept_K2 = k2; DoTransport(); }

    private:
        tune_t fTune;

        double Bq1, Bq2, Bq3,K1, K2;
        double sept_K1, sept_K2;
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

        double fxtg_rng;
        double fytg_rng;
        double fthtg_rng;
        double fphtg_rng;
        double fdp_rng;
        double fxtg_lim;
        double fytg_lim;
        double fthtg_lim;
        double fphtg_lim;
        double fdp_lim;

        int ntrk;
        double ytg[100000], xtg[100000], thtg[100000], phtg[100000], dp[100000];
        bool acc[100000];

        void setcrossterms(TMatrixD *);
        void fillvector(TVectorD &);

         TGraph *gx_th, *gx_x0, *gx_d, *gy_ph, *gy_y0 , *gph_y0, *gph_ph;
};
#endif//__THRSTRANS_H
