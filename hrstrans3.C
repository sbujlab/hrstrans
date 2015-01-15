#include "THRSTrans.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TCanvas.h"



void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    THRSTrans *trans = new THRSTrans(par[0], par[1], par[2], par[3], par[4], THRSTrans::kStd, 10);

    double chi2 = 0.0;

    TMatrixD res = *(trans->GetTransport());
    TMatrixD optics = *(trans->GetOptics());

    // PVDIS Optics DB from Kiads directory
    chi2 += pow(optics[0][0], 2.0);  
    chi2 += pow(optics[0][1] + 2.307403, 2.0);  
    chi2 += pow(optics[0][2], 2.0);  
    chi2 += pow(optics[0][3], 2.0);  

    chi2 += pow(optics[1][0] - 7.981594e-04, 2.0);  
    chi2 += pow(optics[1][1] + 1.853464e-02, 2.0);  
    chi2 += pow(optics[1][2] + 1.149558e+00, 2.0);  
    chi2 += pow(optics[1][3] - 8.035791e-01, 2.0);  

    chi2 += pow(optics[2][0] - 5.098623e-04, 2.0);  
    chi2 += pow(optics[2][1] , 2.0);  
    chi2 += pow(optics[2][2] + 3.104328e-01, 2.0);  
    chi2 += pow(optics[2][3] + 6.248311e-01, 2.0);  

    chi2 += pow(optics[3][0] - 8.370861e-02, 2.0);  
    chi2 += pow(optics[3][1] + 1.737617e-02, 2.0);  
    chi2 += pow(optics[3][2], 2.0);  
    chi2 += pow(optics[3][3], 2.0);  


    //chi2 += 100*pow(res[THRSTrans::kX][THRSTrans::kTh] + 0.2, 2.0);  // This should be absolutely zero


    /*
    chi2 += pow(res[THRSTrans::kX][THRSTrans::kX] + 2.48, 2.0);  
    chi2 += pow(res[THRSTrans::kTh][THRSTrans::kX] + 0.15, 2.0); 
    chi2 += pow(res[THRSTrans::kTh][THRSTrans::kTh] + 0.40, 2.0); 

    chi2 += 1000*pow(res[THRSTrans::kX][THRSTrans::kTh], 2.0);  // This should be zero
//    chi2 += 1.0/pow( res[THRSTrans::kX][THRSTrans::kd]/res[THRSTrans::kX][THRSTrans::kX], 2.0);  // This should be maximized
//
    // Force a couple matrix elements

    chi2 += 10*pow(res[THRSTrans::kY][THRSTrans::kY] + 0.4 , 2.0);  // This should be zero
    chi2 += 10*pow(res[THRSTrans::kPh][THRSTrans::kPh] + 0.78 , 2.0);  // This should be zero

    // PREX Criteria
//    chi2 += 1000*pow(res[THRSTrans::kX][THRSTrans::kTh], 2.0);  // This should be zero
//    chi2 += 1000*pow(res[THRSTrans::kY][THRSTrans::kPh], 2.0);  // This should be zero

    // Maximize specified acceptance
    int i;

    double thsum = 0.0;
    // Th
    for( i = 0; i < trans->GetThAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetThAcc()->GetX())[i] ) < 0.03 ){
            thsum += pow( (trans->GetThAcc()->GetY())[i], 2.0 );
        }
    }

    double ysum = 0.0;
    // Y
    for( i = 0; i < trans->GetYAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetYAcc()->GetX())[i] ) < 0.05 ){
            ysum += pow( (trans->GetYAcc()->GetY())[i], 2.0 );
        }
    }
    // Ph 
    double phsum = 0.0;
    for( i = 0; i < trans->GetPhAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetPhAcc()->GetX())[i] ) < 0.014 ){
            phsum += pow( (trans->GetPhAcc()->GetY())[i], 2.0 );
        }
    }
    // d
    double dsum = 0.0;
    for( i = 0; i < trans->GetDpAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetDpAcc()->GetX())[i] ) < 0.045 ){
            dsum += pow( (trans->GetDpAcc()->GetY())[i], 2.0 );
        }
    }
    */

    /*
    if( dsum > 0 ){
    chi2 += 1.0/dsum;
    }
    if( phsum > 0 ){
    chi2 += 1.0/phsum;
    }
    if( ysum > 0 ){
    chi2 += 1.0/ysum;
    }
    if( thsum > 0 ){
        chi2 += 1.0/thsum;
    }
    */


    f = chi2;

    delete trans;

    printf("%f %f %f -> %f\n",par[0], par[1], par[2], f);
}




void hrstrans3(){

//    LeRose SNAKE scaled by 1.4
    THRSTrans *trans = new THRSTrans( 0.1746, -0.1385, -0.1281, 0.050178, 0.037056, THRSTrans::kStd);
//
//
//    // Fit to data optics
//    THRSTrans *trans = new THRSTrans(  0.178259, -0.137781, -0.128674, 0.050178, 0.037056, THRSTrans::kStd);

    /*

//    trans->GetElement(5)->Print();
   trans->ShowOutput();
   trans->ShowAcc();
   trans->ShowFocalLengths();
   return;

   */

    // LeRose Translated Numbers
    //THRSTrans *trans = new THRSTrans( 0.1861, -0.1415, -0.1375, 0.5, 0.5, THRSTrans::kStd);
//   trans->ShowOutput();
  // return;
    // LeRose Numbers Fit to Transport matrix elements, PREX
//    THRSTrans *trans = new THRSTrans( 1.34761e-01, -1.38523e-01, -2.13742e-01, 0.02075 };
//
//

//     THRSTrans *trans = new THRSTrans(  1.344e-01, -1.431e-01, -1.855e-01, 0.005763, 0.000014, THRSTrans::kPREX);
//     trans->SetSeptumPsi(  1.22603, 0.929424);


    // HRS Database fit forcing 
//    THRSTrans *trans = new THRSTrans(     1.38796e-01, -1.39940e-01, -1.73034e-01, 0.02075, THRSTrans::kPREX);

//     trans->ShowOutput(4);

//     trans->ShowOutput();
//     trans->ShowFocalLengths();
//     trans->ShowAcc();
//   return;


    TMinuit *gMinuit = new TMinuit(5);
    gMinuit->SetFCN(fcn);
    Double_t arglist[10];
    Int_t ierflg = 0;

    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    Double_t vstart[5] = {
    // Standard tune from data
    0.178615, -0.137474, -0.128384, 0.043667, 0.036814};

    //PREX Tune
       

    Double_t step[5] = {0.001 , 0.001 , 0.001, 0.001, 0.001};

    gMinuit->mnparm(0, "q1", vstart[0], step[0], 0,0,ierflg);
    gMinuit->mnparm(1, "q2", vstart[1], step[1], 0,0,ierflg);
    gMinuit->mnparm(2, "q3", vstart[2], step[2], 0, 0,ierflg);
    gMinuit->mnparm(3, "K1", vstart[3], step[3], 0.0,1.5,ierflg);
    gMinuit->mnparm(4, "K2", vstart[4], step[4], 0.0,1.5,ierflg);


    arglist[0] = 5000;
    arglist[1] = 1.;
    gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    gMinuit->mnimpr();
    gMinuit->mnimpr();
    double q1, q2, q3, psi1, psi2, e1;

    gMinuit->GetParameter(0, q1, e1);
    gMinuit->GetParameter(1, q2, e1);
    gMinuit->GetParameter(2, q3, e1);
    gMinuit->GetParameter(3, psi1, e1);
    gMinuit->GetParameter(4, psi2, e1);


    THRSTrans *result = new THRSTrans( q1, q2, q3, psi1, psi2);

    result->ShowOutput();
    result->GetOptics()->Print();
    result->ShowAcc();
    result->ShowFocalLengths();

    return;
}


