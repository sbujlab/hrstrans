void runmin(){
    gROOT->ProcessLine(".L THRSTrans.C+");

//   return;
//    gROOT->ProcessLine(".x hrstrans_scan.C+");

//    gROOT->ProcessLine(".x hrstrans3.C+");
//   gROOT->ProcessLine(".x hrstrans4.C+");
//   gROOT->ProcessLine(".x hrstrans4.C+");
//   gROOT->ProcessLine(".x hrstrans_prexii.C+");
   gROOT->ProcessLine(".x hrstrans_prexii_tuneP.C+");
}
