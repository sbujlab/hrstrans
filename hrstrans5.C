#include "THRSTrans.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TCanvas.h"


//  Try to fit to all intermediate matrices

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    THRSTrans *trans = new THRSTrans(par[0], par[1], par[2], 0.050178, 0.037056, THRSTrans::kAPEX, 50000);


    TMatrixD tovdc= *(trans->GetTransport(14));

    TMatrixD optics = *(trans->GetOptics(14));
;

    TMatrixD tofp = *(trans->GetTransport(14));
    double chi2 = 0.0;

    TMatrixD tosen = *(trans->GetTransport(2));
    TMatrixD tosex = *(trans->GetTransport(3));

    TMatrixD q2ex = *(trans->GetTransport(9));
    TMatrixD den = *(trans->GetTransport(10));
    TMatrixD q3en = *(trans->GetTransport(12));
    TMatrixD q3ex = *(trans->GetTransport(13));


    int i;
    double thsum = 0.0;
    // Th
    for( i = 0; i < trans->GetThAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetThAcc()->GetX())[i] ) < 0.60 ){
            thsum += pow( (trans->GetThAcc()->GetY())[i], 2.0 );
        }
    }

    double ysum = 0.0;
    // Y
    // APEX target is 50cm long at six deg
    for( i = 0; i < trans->GetYAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetYAcc()->GetX())[i] ) <  0.25*sin(6.0*3.14159/180) ){
            ysum += pow( (trans->GetYAcc()->GetY())[i], 2.0 );
        }
    }
    // Ph 
    double phsum = 0.0;
    for( i = 0; i < trans->GetPhAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetPhAcc()->GetX())[i] ) < 0.030 ){
            phsum += pow( (trans->GetPhAcc()->GetY())[i], 2.0 );
        }
    }
    // d
    double dsum = 0.0;
    for( i = 0; i < trans->GetDpAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetDpAcc()->GetX())[i] ) < 0.045 ){
            dsum += pow( (trans->GetDpAcc()->GetY())[i], 2.0 );
        //    printf("Adding (%f) %f\n", (trans->GetDpAcc()->GetX())[i],    (trans->GetDpAcc()->GetY())[i] );
        }
    }

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

    chi2 += 100.*pow(optics[2][0], 2.0);  
    chi2 += 500.*pow(optics[2][1], 2.0);  
    chi2 += 100.*pow(optics[2][2], 2.0);  
    chi2 += 500.*pow(optics[2][3], 2.0);  

    f = chi2;

    delete trans;

    printf("%f %f %f -> %f\n",par[0], par[1], par[2], f);
}




void hrstrans5(){

    // LeRose SNAKE scaled by 1.4 (in agreement with standard tune)
//    THRSTrans *trans = new THRSTrans( 0.1746, -0.1385, -0.1281, 0.050178, 0.037056, THRSTrans::kAPEX);


    // Just set Q1 focus further down
 //   THRSTrans *trans = new THRSTrans( 0.09, -0.1385, -0.1281, 0.050178, 0.037056, THRSTrans::kAPEX);
    // First shot at optmization
    THRSTrans *trans = new THRSTrans( 0.148284, -0.126489, -0.162333, 0.050178, 0.037056, THRSTrans::kAPEX);
  //  Tune from Logbook
//  THRSTrans *trans = new THRSTrans( 0.1183, -0.1310, -0.1554, 0.050178, 0.037056, THRSTrans::kAPEX);


    trans->ShowOutput();
    trans->ShowAcc();
    trans->ShowFocalLengths();
    trans->GetOptics()->Print();

    return;

    TMinuit *gMinuit = new TMinuit(3);
    gMinuit->SetFCN(fcn);
    Double_t arglist[10];
    Int_t ierflg = 0;

    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    Double_t vstart[5] = {
    // Standard tune
    //PREX Tune
    1.38796e-01, -1.39940e-01, -1.73034e-01};

    Double_t step[5] = {0.001 , 0.001 , 0.001};

    gMinuit->mnparm(0, "q1", vstart[0], step[0], 0,0,ierflg);
    gMinuit->mnparm(1, "q2", vstart[1], step[1], 0,0,ierflg);
    gMinuit->mnparm(2, "q3", vstart[2], step[2], 0, 0,ierflg);


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

    THRSTrans *result = new THRSTrans( q1, q2, q3, 0.050178, 0.037056, THRSTrans::kAPEX);

    /*
    result->ShowOutput(9);
    result->ShowOutput(12);
    printf("q3ex\n");
    result->ShowOutput(13);

    printf("fp + 0.75\n");
    result->ShowOutput(-1, 0.75);
    printf("fp + 1.43\n");
    result->ShowOutput(-1, 1.43);
    */

    //result->ShowOutput();
    result->ShowOutput();
    result->GetOptics()->Print();
    result->ShowAcc();
    result->ShowFocalLengths();

    return;
}


