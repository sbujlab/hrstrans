#include "THRSTrans.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom.h"
#endif
#include "TGraph.h"
#include "TCanvas.h"
#include <iostream>
#include <stdio.h>

using namespace std;

int main(int argc, char **argv) 
{
 
   Int_t debug_out = 2;

 // Initialize root and output.  
   TROOT scalana("hrsroot","Hall A spectrometer optics analysis");
   TFile hfile("hrs.root","RECREATE","HRSTrans analysis");

// Define the ntuple here. 
//               0   1    2    3   4    5    6    7      8   9    10    11     
   char ctrans[]="evt:inx:inth:iny:inph:indp:outx:outth:outy:outph:outdp:hamcx:hamcth:hamcy:hamcph:hamcdp:wsx:wsth:wsy:wsph:wsdp";
//12  13    14      15    16  17  18  19    20

   Int_t nlen = strlen(ctrans);
   char *string_ntup = new char[nlen+1];
   strcpy(string_ntup, ctrans);

   TNtuple *hrt = new TNtuple("hrt","HRS Trans analysis",string_ntup);
   Float_t* farray_ntup = new Float_t[nlen+1];  
// end of ntuple definition

   THRSTrans *trans = new THRSTrans(0.096255, -0.131739, -0.170480, 0.050178, 0.037056, THRSTrans::kPREX);

     // new THRSTrans(0.099950, -0.132890, -0.171751, 0.050178, 0.037056, THRSTrans::kPREX);

    trans->ShowOutput();

    // For a check use the 1st order warm septum transport to focal plane
    TMatrixD *wS = new TMatrixD(20,20);

    (*wS)[THRSTrans::kX][THRSTrans::kX]  = -2.507;
    (*wS)[THRSTrans::kX][THRSTrans::kTh] = -0.0148;
    (*wS)[THRSTrans::kX][THRSTrans::kY]  = -0.089;
    (*wS)[THRSTrans::kX][THRSTrans::kPh] = 0;
    (*wS)[THRSTrans::kX][THRSTrans::kd]  = 14.2;

    (*wS)[THRSTrans::kTh][THRSTrans::kX]  = -0.309;
    (*wS)[THRSTrans::kTh][THRSTrans::kTh] = -0.401;
    (*wS)[THRSTrans::kTh][THRSTrans::kY]  = 0;
    (*wS)[THRSTrans::kTh][THRSTrans::kPh] = 0.002;
    (*wS)[THRSTrans::kTh][THRSTrans::kd]  = 2.501;

    (*wS)[THRSTrans::kY][THRSTrans::kX]  = 0.008;
    (*wS)[THRSTrans::kY][THRSTrans::kTh] = 0.013;
    (*wS)[THRSTrans::kY][THRSTrans::kY]  = 0.321;
    (*wS)[THRSTrans::kY][THRSTrans::kPh] = -2.181;
    (*wS)[THRSTrans::kY][THRSTrans::kd]  = -0.351;

    (*wS)[THRSTrans::kPh][THRSTrans::kX]  = 0.007;
    (*wS)[THRSTrans::kPh][THRSTrans::kTh] = 0.008;
    (*wS)[THRSTrans::kPh][THRSTrans::kY]  = 0.610;
    (*wS)[THRSTrans::kPh][THRSTrans::kPh] = -1.030;
    (*wS)[THRSTrans::kPh][THRSTrans::kd]  = -0.270;

    (*wS)[THRSTrans::kd][THRSTrans::kX]  = 0;
    (*wS)[THRSTrans::kd][THRSTrans::kTh] = 0;
    (*wS)[THRSTrans::kd][THRSTrans::kY]  = 0;
    (*wS)[THRSTrans::kd][THRSTrans::kPh] = 0;
    (*wS)[THRSTrans::kd][THRSTrans::kd]  = 1;

    TVectorD iv(20),v(20),v2(20),vT(20);

   // Read in the tracks from hamc

    char strin[150];
    Float_t d1,d2,d3,d4,d5,d6,d7,d8,d9,d10;
    Float_t transOut[5];
      
    FILE *fd = fopen("tracks_out.txt","r");
    if (fd == NULL) {
      cout << "No input data !  Missing file ! Quitting. "<<endl;
      exit(0);
    }

    cout << "\n\n--------------------------------------"<<endl<<endl;

    Int_t iev=0;

    while (fgets(strin,150,fd) != NULL)  {

      sscanf(strin,"HRSTR %f %f %f %f %f %f %f %f %f %f",&d1,&d2,&d3,&d4,&d5,&d6,&d7,&d8,&d9,&d10);

      transOut[0] = d6;
      transOut[1] = d7;
      transOut[2] = d8;
      transOut[3] = d9;
      transOut[4] = d10;

// I think sign convention for X, Y are opposite for THRSTrans but same for Warm-Septum matrix

      iv[THRSTrans::kX]  = d2;
      iv[THRSTrans::kTh] = d3;
      iv[THRSTrans::kY]  = d4;
      iv[THRSTrans::kPh] = d5;
      iv[THRSTrans::kd]  = d10;

      v2 = (*wS)*iv;

      iv[THRSTrans::kX]  = 1.0*d2;
      iv[THRSTrans::kY]  = 1.0*d4;

      if (debug_out==2) {
         cout << "\n\n--------------------------------------"<<endl<<endl;
         cout << "input string "<<strin<<endl;
         cout << "d:  "<<d1<<"  "<<d2<<"  "<<d3<<"  "<<d4<<"  "<<d5<<"  "<<d6<<"  "<<d7<<"  "<<d8<<"  "<<d9<<"  "<<d10<<endl;
      }

      trans->fillvector(iv);  // add 2nd order

      TMatrixD tofp= *(trans->GetTransport());  // At the focal plane

      v = (tofp)*iv;

      farray_ntup[0] = iev++;

      farray_ntup[1] = iv[THRSTrans::kX];
      farray_ntup[2] = iv[THRSTrans::kTh];
      farray_ntup[3] = iv[THRSTrans::kY];
      farray_ntup[4] = iv[THRSTrans::kPh];
      farray_ntup[5] = iv[THRSTrans::kd];

      farray_ntup[6] = v[THRSTrans::kX];
      farray_ntup[7] = v[THRSTrans::kTh];
      farray_ntup[8] = v[THRSTrans::kY];
      farray_ntup[9] = v[THRSTrans::kPh];
      farray_ntup[10] = v[THRSTrans::kd];

      farray_ntup[11] = transOut[0];
      farray_ntup[12] = transOut[1];
      farray_ntup[13] = transOut[2];
      farray_ntup[14] = transOut[3];
      farray_ntup[15] = transOut[4];

      farray_ntup[16] = v2[THRSTrans::kX];
      farray_ntup[17] = v2[THRSTrans::kTh];
      farray_ntup[18] = v2[THRSTrans::kY];
      farray_ntup[19] = v2[THRSTrans::kPh];
      farray_ntup[20] = v2[THRSTrans::kd];

      hrt->Fill(farray_ntup);

      if (debug_out) {
        cout << "\n\nEvent number "<<iev<<endl;
        cout << "input vector "<<iv[THRSTrans::kX]<<"  "<<iv[THRSTrans::kTh] << "  "<<iv[THRSTrans::kY]<<"  "<<iv[THRSTrans::kPh]<<"  "<<iv[THRSTrans::kd]<<endl;
        cout << "output vector  "<<v[THRSTrans::kX]<<"  "<<v[THRSTrans::kTh] << "  "<<v[THRSTrans::kY]<<"  "<<v[THRSTrans::kPh]<<"  "<<v[THRSTrans::kd]<<endl;
        cout << "from wS  "<<v2[THRSTrans::kX]<<"  "<<v2[THRSTrans::kTh] << "  "<<v2[THRSTrans::kY]<<"  "<<v2[THRSTrans::kPh]<<"  "<<v2[THRSTrans::kd]<<endl;
        cout << "from hamc  "<<transOut[0]<<"  "<<transOut[1]<<"  "<<transOut[2]<<"  "<<transOut[3]<<"  "<<transOut[4]<<endl;
      }

    }

    hfile.Write();
    hfile.Close();

    return 1;

}



