#include "THRSTrans.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TCanvas.h"


//  Try to fit to all intermediate matrices

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    THRSTrans *trans = new THRSTrans(par[0], par[1], par[2], 0.050178, 0.037056, THRSTrans::kPREX, 50000);


    TMatrixD tovdc= *(trans->GetTransport(14));

    TMatrixD optics = *(trans->GetOptics(14));
;

    TMatrixD tofp = *(trans->GetTransport());
    double chi2 = 0.0;

    TMatrixD tosen = *(trans->GetTransport(2));
    TMatrixD tosex = *(trans->GetTransport(3));

    TMatrixD q2ex = *(trans->GetTransport(9));
    TMatrixD den = *(trans->GetTransport(10));
    TMatrixD q3en = *(trans->GetTransport(12));
    TMatrixD q3ex = *(trans->GetTransport(13));


    /*
      // LeRose Predicted PREX elements
    chi2 += pow(q2ex[THRSTrans::kX][THRSTrans::kTh] - 3.715, 2.0);  
    chi2 += pow(q2ex[THRSTrans::kTh][THRSTrans::kX] - 0.476, 2.0);  
    chi2 += pow(q2ex[THRSTrans::kTh][THRSTrans::kTh] - 0.936, 2.0);  

    chi2 += pow(q2ex[THRSTrans::kY][THRSTrans::kY] - 1.842, 2.0);  
    chi2 += pow(q2ex[THRSTrans::kY][THRSTrans::kPh] - 7.898, 2.0);  
    chi2 += pow(q2ex[THRSTrans::kPh][THRSTrans::kY] + 0.387, 2.0);  
    chi2 += pow(q2ex[THRSTrans::kPh][THRSTrans::kPh] + 1.110, 2.0);  

    chi2 += pow(q2ex[THRSTrans::kY][THRSTrans::kd] - 0.627, 2.0);  

    ////////
    chi2 += pow(q3en[THRSTrans::kX][THRSTrans::kX] + 2.072, 2.0);  
    chi2 += pow(q3en[THRSTrans::kX][THRSTrans::kTh] - 3.707, 2.0);  
    chi2 += pow(q3en[THRSTrans::kTh][THRSTrans::kX] - 0.494, 2.0);  
    chi2 += pow(q3en[THRSTrans::kTh][THRSTrans::kTh] + 1.380, 2.0);  

    chi2 += pow(q3en[THRSTrans::kY][THRSTrans::kY] + 2.622, 2.0);  
    chi2 += pow(q3en[THRSTrans::kY][THRSTrans::kPh] + 1.583, 2.0);  
    chi2 += pow(q3en[THRSTrans::kPh][THRSTrans::kY] - 0.602, 2.0);  
    chi2 += pow(q3en[THRSTrans::kPh][THRSTrans::kPh] - 0.748, 2.0);  

    chi2 += pow(q3en[THRSTrans::kX][THRSTrans::kd] - 2.770, 2.0);  
    chi2 += pow(q3en[THRSTrans::kTh][THRSTrans::kd] - 0.455, 2.0);  

    */

    /*
    chi2 += pow(q3ex[THRSTrans::kX][THRSTrans::kX] + 1.854, 2.0);  
    chi2 += pow(q3ex[THRSTrans::kX][THRSTrans::kTh] - 1.393, 2.0);  
    chi2 += pow(q3ex[THRSTrans::kTh][THRSTrans::kX] + 0.341, 2.0);  
    chi2 += pow(q3ex[THRSTrans::kTh][THRSTrans::kTh] - 0.279, 2.0);  

    chi2 += pow(q3ex[THRSTrans::kY][THRSTrans::kY] + 2., 2.0);  
    chi2 += pow(q3ex[THRSTrans::kY][THRSTrans::kPh] + 1.892, 2.0);  
    chi2 += pow(q3ex[THRSTrans::kPh][THRSTrans::kY] - 0.955, 2.0);  
    chi2 += pow(q3ex[THRSTrans::kPh][THRSTrans::kPh] - 0.399, 2.0);  

    chi2 += pow(q3ex[THRSTrans::kX][THRSTrans::kd] - 6.992, 2.0);  
    chi2 += pow(q3ex[THRSTrans::kTh][THRSTrans::kd] - 2.543, 2.0);  
    chi2 += pow(q3ex[THRSTrans::kY][THRSTrans::kd] - 0.169, 2.0);  
    chi2 += pow(q3ex[THRSTrans::kPh][THRSTrans::kd] + 0.137, 2.0);  
    */

    /*

//    chi2 += pow(den[THRSTrans::kX][THRSTrans::kX] - 0.8, 2.0);  
    chi2 += pow(den[THRSTrans::kX][THRSTrans::kTh] + 2.198, 2.0);  
    chi2 += pow(den[THRSTrans::kTh][THRSTrans::kX] + 0.482, 2.0);  
    chi2 += pow(den[THRSTrans::kTh][THRSTrans::kTh] - 0.962, 2.0);  

    chi2 += pow(den[THRSTrans::kY][THRSTrans::kY] - 0.346, 2.0);  
    chi2 += pow(den[THRSTrans::kY][THRSTrans::kPh] - 3.628, 2.0);  
    chi2 += pow(den[THRSTrans::kPh][THRSTrans::kY] + 1.366, 2.0);  
    chi2 += pow(den[THRSTrans::kPh][THRSTrans::kPh] + 3.534, 2.0);  

    chi2 += pow(den[THRSTrans::kY][THRSTrans::kd] - 0.376, 2.0);  
    chi2 += pow(den[THRSTrans::kPh][THRSTrans::kd] - 0.176, 2.0);  


    */
    // Elements from database
    
    /*
    chi2 += pow(optics[0][0] - 0.0392, 2.0);  
    chi2 += pow(optics[0][1] + 3.135, 2.0);  
    chi2 += pow(optics[0][2] + 0.03369, 2.0);  
    chi2 += pow(optics[0][3] - 0.3608, 2.0);  

    chi2 += pow(optics[1][0] + 0.0192, 2.0);  
    chi2 += pow(optics[1][1] , 2.0);  
    chi2 += pow(optics[1][2] - 1, 2.0);  
    chi2 += pow(optics[1][3] -1.361, 2.0);  

    chi2 += pow(optics[2][0], 2.0);  
    chi2 += pow(optics[2][1] + 0.1618, 2.0);  
    chi2 += pow(optics[2][2] + 0.9916 , 2.0);  
    chi2 += pow(optics[2][3] - 1.121, 2.0);  

    chi2 += pow(optics[3][0] -0.06876, 2.0);  
    chi2 += pow(optics[3][1] -0.07035, 2.0);  
    chi2 += pow(optics[3][2] + 0.0042, 2.0);  
    chi2 += pow(optics[3][3] - 0.0311, 2.0);  
    */

    int i;
    double thsum = 0.0;
    // Th
    for( i = 0; i < trans->GetThAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetThAcc()->GetX())[i] ) < 0.98 ){
            thsum += pow( (trans->GetThAcc()->GetY())[i], 2.0 );
        }
    }

    double ysum = 0.0;
    // Y
    for( i = 0; i < trans->GetYAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetYAcc()->GetX())[i] ) < 0.01 ){
            ysum += pow( (trans->GetYAcc()->GetY())[i], 2.0 );
        }
    }
    // Ph 
    double phsum = 0.0;
    for( i = 0; i < trans->GetPhAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetPhAcc()->GetX())[i] ) < 0.020 ){
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


    chi2 += 100*pow(tofp[THRSTrans::kX][THRSTrans::kTh], 2.0);  
    chi2 += 100*pow(tofp[THRSTrans::kY][THRSTrans::kPh], 2.0);  

    f = chi2;

    delete trans;

    printf("%f %f %f -> %f\n",par[0], par[1], par[2], f);
}




void hrstrans4(){

    // LeRose SNAKE scaled by 1.4 (in agreement with standard tune)
    THRSTrans *trans = new THRSTrans( 0.1288, -0.1393, -0.1736, 0.050178, 0.037056, THRSTrans::kPREX);

    // Fit to LeRose
//    THRSTrans *trans = new THRSTrans( 0.098786, -0.132014, -0.172415, 0.050178, 0.037056, THRSTrans::kPREX);

    printf("LeRose focal plane\n");
    trans->ShowOutput(14, 0.75);
    printf("True focal plane\n");
    trans->ShowOutput();
    trans->ShowAcc();

    trans->GetOptics(14)->Print();

    trans->ShowFocalLengths();

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

    THRSTrans *result = new THRSTrans( q1, q2, q3, 0.050178, 0.037056, THRSTrans::kPREX);

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
    printf("LeRose Focal Plane\n");
    result->ShowOutput(14, 0.75);
    printf("True Focal Plane\n");
    result->ShowOutput();
    printf("Optics Matrix\n");
    result->GetOptics()->Print();
    result->ShowAcc();
    result->ShowFocalLengths();

    return;
}


