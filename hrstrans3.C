#include "THRSTrans.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TCanvas.h"



void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    THRSTrans *trans = new THRSTrans(par[0], par[1], par[2], par[3]);

    double chi2 = 0.0;

    TMatrixD res = *(trans->GetTransport());

    chi2 += pow(res[THRSTrans::kX][THRSTrans::kX] + 2.48, 2.0);  
    chi2 += pow(res[THRSTrans::kTh][THRSTrans::kX] + 0.15, 2.0); 
    chi2 += 1000*pow(res[THRSTrans::kX][THRSTrans::kTh], 2.0);  // This should be zero
    chi2 += pow(res[THRSTrans::kTh][THRSTrans::kTh] + 0.40, 2.0); 
//    chi2 += 1.0/pow( res[THRSTrans::kX][THRSTrans::kd]/res[THRSTrans::kX][THRSTrans::kX], 2.0);  // This should be maximized
//
    // Force a couple matrix elements

    chi2 += 10*pow(res[THRSTrans::kY][THRSTrans::kY] + 0.4 , 2.0);  // This should be zero
    chi2 += 10*pow(res[THRSTrans::kPh][THRSTrans::kPh] + 0.78 , 2.0);  // This should be zero

    /*
    chi2 += 1000*pow(res[THRSTrans::kX][THRSTrans::kX], 2.0);  // This should be zero
    chi2 += 1000*pow(res[THRSTrans::kX][THRSTrans::kTh], 2.0);  // This should be zero
    chi2 += pow(res[THRSTrans::kY][THRSTrans::kY], 2.0);  // This should be zero
    chi2 += pow(res[THRSTrans::kY][THRSTrans::kPh], 2.0);  // This should be zero
    */

    // Maximize specified acceptance
    int i;

    double thsum = 0.0;
    // Th
    for( i = 0; i < trans->GetThAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetThAcc()->GetX())[i] ) < 0.03 ){
            thsum += pow( (trans->GetThAcc()->GetY())[i], 2.0 );
        }
    }
    if( thsum > 0 ){
    //    chi2 += 1.0/thsum;
    }

    double ysum = 0.0;
    // Y
    for( i = 0; i < trans->GetYAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetYAcc()->GetX())[i] ) < 0.05 ){
            ysum += pow( (trans->GetYAcc()->GetY())[i], 2.0 );
        }
    }
    if( ysum > 0 ){
    //chi2 += 1.0/ysum;
    }
    // Ph 
    double phsum = 0.0;
    for( i = 0; i < trans->GetPhAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetPhAcc()->GetX())[i] ) < 0.014 ){
            phsum += pow( (trans->GetPhAcc()->GetY())[i], 2.0 );
        }
    }
    if( phsum > 0 ){
    //chi2 += 1.0/phsum;
    }
    // d
    double dsum = 0.0;
    for( i = 0; i < trans->GetDpAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetDpAcc()->GetX())[i] ) < 0.045 ){
            dsum += pow( (trans->GetDpAcc()->GetY())[i], 2.0 );
        }
    }
    if( dsum > 0 ){
    //chi2 += 1.0/dsum;
    }


    f = chi2;

    delete trans;

    printf("%f %f %f -> %f\n",par[0], par[1], par[2], f);
}




void hrstrans3(){

    THRSTrans *trans = new THRSTrans( 1.86122e-01, -1.41506e-01, -1.37479e-01, 2.07500e-02, THRSTrans::kStd);
   trans->ShowOutput();
   return;

    TMinuit *gMinuit = new TMinuit(4);
    gMinuit->SetFCN(fcn);
    Double_t arglist[10];
    Int_t ierflg = 0;

    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    Double_t vstart[4] = {
    // Standard tune
    1.86122e-01, -1.41506e-01, -1.37479e-01, 2.07500e-02 };
    //PREX Tune

    Double_t step[4] = {0.001 , 0.001 , 0.001, 0.001};

    gMinuit->mnparm(0, "q1", vstart[0], step[0], 0,0,ierflg);
    gMinuit->mnparm(1, "q2", vstart[1], step[1], 0,0,ierflg);
    gMinuit->mnparm(2, "q3", vstart[2], step[2], 0, 0,ierflg);
    gMinuit->mnparm(3, "psi", vstart[3], step[3], 0,0,ierflg);


    arglist[0] = 5000;
    arglist[1] = 1.;
    gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    double q1, q2, q3, psi, e1;

    gMinuit->GetParameter(0, q1, e1);
    gMinuit->GetParameter(1, q2, e1);
    gMinuit->GetParameter(2, q3, e1);
    gMinuit->GetParameter(3, psi, e1);

    THRSTrans *result = new THRSTrans( q1, q2, q3, psi);

    result->ShowOutput();

    return;
}


