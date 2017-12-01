{
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatH(0.3);
   gStyle->SetStatW(0.3);
   gStyle->SetTitleH(0.07);
   gStyle->SetTitleW(0.9);
   gStyle->SetLabelSize(0.045,"x");
   gStyle->SetLabelSize(0.045,"y");
   gROOT->ForceStyle();

   // 0 = hamc vs HRSTRANS
   // 1 = hamc vs LeRose's warmseptum 1st order

   Int_t which=0; 

   if (which == 0) {
     TH2F *h2x = new TH2F("h2x","X (HRSTrans) vs X (LeRose fcn)",100,-1.0,0.08,100,-1.0,0.08);
     TH2F *h3x = new TH2F("h3x","X (HRSTrans) vs X (LeRose fcn)",100,-1.0,0.08,100,-1.0,0.08);
     TH2F *h2th = new TH2F("h2th","tan(th) HRSTrans vs LeRose fcn",100,-0.08,0.04,100,-0.08,0.04);
     TH2F *h3th = new TH2F("h3th","tan(th) HRSTrans vs LeRose fcn",100,-0.08,0.04,100,-0.08,0.04);
     TH2F *h2y = new TH2F("h2y","Y (HRSTrans) vs Y (LeRose fcn)",100,-0.08,0.08,100,-0.08,0.08);
     TH2F *h3y = new TH2F("h3y","Y (HRSTrans) vs Y (LeRose fcn)",100,-0.08,0.08,100,-0.08,0.08);
     TH2F *h2ph = new TH2F("h2ph","tan(phi) HRSTrans vs LeRose fcn",100,-0.04,0.04,100,-0.04,0.04);
     TH2F *h3ph = new TH2F("h3ph","tan(phi) HRSTrans vs LeRose fcn",100,-0.04,0.04,100,-0.04,0.04);
   } else {
     TH2F *h2x = new TH2F("h2x","X (LeRose fcn) vs X (LeRose 1st Order)",100,-1.0,0.08,100,-1.0,0.08);
     TH2F *h3x = new TH2F("h3x","X (LeRose fcn) vs X (LeRose 1st Order)",100,-1.0,0.08,100,-1.0,0.08);
     TH2F *h2th = new TH2F("h2th","tan(th) LeRose fcn vs LeRose 1st Order",100,-0.08,0.04,100,-0.08,0.04);
     TH2F *h3th = new TH2F("h3th","tan(th) LeRose fcn vs LeRose 1st Order",100,-0.08,0.04,100,-0.08,0.04);
     TH2F *h2y = new TH2F("h2y","Y (LeRose fcn) vs Y (LeRose 1st Order)",100,-0.08,0.08,100,-0.08,0.08);
     TH2F *h3y = new TH2F("h3y","Y (LeRose fcn) vs Y (LeRose 1st Order)",100,-0.08,0.08,100,-0.08,0.08);
     TH2F *h2ph = new TH2F("h2ph","tan(phi) LeRose fcn vs LeRose 1st Order",100,-0.04,0.04,100,-0.04,0.04);
     TH2F *h3ph = new TH2F("h3ph","tan(phi) LeRose fcn vs LeRose 1st Order",100,-0.04,0.04,100,-0.04,0.04);
   }

   Float_t zloc = -1.4;
   char strplot[100];
 
   if (which==0) {
     sprintf(strplot,"outx:hamcx-%f*hamcth>>h2x",zloc);
     hrt->Draw(strplot);
     sprintf(strplot,"outx:hamcx-%f*hamcth>>h3x",zloc);
     hrt->Draw(strplot,"wsdp>-0.002");
     hrt->Draw("outth:hamcth>>h2th");
     hrt->Draw("outth:hamcth>>h3th","wsdp>-0.002");
     sprintf(strplot,"outy:hamcy-%f*hamcph>>h2y",zloc);
     hrt->Draw(strplot);
     sprintf(strplot,"outy:hamcy-%f*hamcph>>h3y",zloc);
     hrt->Draw(strplot,"wsdp>-0.002");
     hrt->Draw("outph:hamcph>>h2ph");
     hrt->Draw("outph:hamcph>>h3ph","wsdp>-0.002");
   } else {
     sprintf(strplot,"hamcx:wsx-%f*wsth>>h2x",zloc);
     hrt->Draw(strplot);
     sprintf(strplot,"hamcx:wsx-%f*wsth>>h3x",zloc);
     hrt->Draw(strplot,"wsdp>-0.002");
     hrt->Draw("hamcth:wsth>>h2th");
     hrt->Draw("hamcth:wsth>>h3th","wsdp>-0.002");
     sprintf(strplot,"hamcy:wsy-%f*wsph>>h2y",zloc);
     hrt->Draw(strplot);
     sprintf(strplot,"hamcy:wsy-%f*wsph>>h3y",zloc);
     hrt->Draw(strplot,"wsdp>-0.002");
     hrt->Draw("hamcph:wsph>>h2ph");
     hrt->Draw("hamcph:wsph>>h3ph","wsdp>-0.002");
  }

   h3x->SetMarkerColor(2);
   h3th->SetMarkerColor(2);
   h3y->SetMarkerColor(2);
   h3ph->SetMarkerColor(2);
 
   TCanvas c1;
   c1->Divide(2,2);

   c1->cd(1);

   h2x->Draw();
   h3x->Draw("same");

   c1->cd(2);

   h2th->Draw();
   h3th->Draw("same");

   c1->cd(3);

   h2y->Draw();
   h3y->Draw("same");

   c1->cd(4);

   h2ph->Draw();
   h3ph->Draw("same");

   c1->Update();







}
