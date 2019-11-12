#include "THRSTrans.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TCanvas.h"

//Try to fit to all intermediate matrices 

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
   THRSTrans *trans = new THRSTrans(par[0], par[1], par[2], 0.5,0.5, 0.050178, 0.037056, THRSTrans::kPREX, THRSTrans::kRHRS, 50000);


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

 /////////
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

 /* chi2 += pow(q3ex[THRSTrans::kX][THRSTrans::kX] + 1.854, 2.0);  
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
  //Elements from database
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
 //Th
   for( i = 0; i < trans->GetThAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetThAcc()->GetX())[i] ) < 0.98 ){
            thsum += pow( (trans->GetThAcc()->GetY())[i], 2.0 );
        }
    }
  double ysum = 0.0;
  //Y
   for( i = 0; i < trans->GetYAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetYAcc()->GetX())[i] ) < 0.01 ){
            ysum += pow( (trans->GetYAcc()->GetY())[i], 2.0 );
        }
    }
  //Ph
  double phsum = 0.0;
    for( i = 0; i < trans->GetPhAcc()->GetN(); i++ ){
        if(  fabs( (trans->GetPhAcc()->GetX())[i] ) < 0.020 ){
            phsum += pow( (trans->GetPhAcc()->GetY())[i], 2.0 );
        }
    } 
  //d 
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

 chi2 += 500*pow(tofp[THRSTrans::kX][THRSTrans::kTh], 2.0);  
    chi2 += 100*pow(tofp[THRSTrans::kY][THRSTrans::kPh], 2.0);
/*
   chi2 += 500*pow(tovdc[THRSTrans::kY][THRSTrans::kPh] + 1.4, 2.0);  
   chi2 += 500*pow(tovdc[THRSTrans::kPh][THRSTrans::kPh], 2.0);  
        */

  f = chi2;

    delete trans;

    printf("%f %f %f -> %f\n",par[0], par[1], par[2], f);
}

void hrstrans_prexii_tuneP(){


    // LeRose SNAKE scaled by 1.4 (in agreement with standard tune)
   //THRSTrans *trans = new THRSTrans( 0.1288, -0.1393, -0.1736, 0.050178, 0.037056, THRSTrans::kPREX);
 
//for(int i = 0; i < 1; i++){

//double z = 1.0 + i*0.01;
  
//printf("%.2f",z);

THRSTrans *trans = new THRSTrans(0.122396*0.94, -0.136543, -0.171633*0.93, 0.5, 0.5, 0.050178, 0.037056, THRSTrans::kPREX, THRSTrans::kLHRS);

/*
     printf("Q1 Entrance Transport\n");
     trans->ShowOutput(6,0);

     printf("Q1 Exit Transport\n");
     trans->ShowOutput(7,0);

     printf("Q2 Entrance Transport\n");
     trans->ShowOutput(8,0);

     printf("Q2 Exit Transport\n");
     trans->ShowOutput(9,0);
   
     printf("Dipole entrance Transport\n");
     trans->ShowOutput(10, 0);

      printf("Dipole exit Transport\n");
      trans->ShowOutput(11, 0);

      printf("Q3 Entrance Transport\n");
      trans->ShowOutput(12,0);

      printf("Q3 Exit Transport\n");
      trans->ShowOutput(13,0);
*/


  

//         printf("Septum Ent\n");
  //       trans->ShowOutput(2,0);
 
    //     printf("Septum Exit\n");
      //   trans->ShowOutput(3,0);
   
      //   printf("VDC \n");
       //  trans->ShowOutput(14,0);



      //   trans->Matri();


//         printf("Collimator\n");
  //       trans->ShowOutput(5,0);
//         trans->ShowFocalLengths();

     printf("Focal Plane\n");
     trans->ShowOutput();
//}

//    printf("Optics Matrix\n");
  //  trans->GetOptics()->Print();

    return;
}


