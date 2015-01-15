#include "THRSTrans.h"
#include "TMatrixD.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVectorD.h"


THRSTrans::THRSTrans(double bq1, double bq2, double bq3, double k1, double k2, tune_t tune, int nt ):
    Bq1(bq1),Bq2(bq2),Bq3(bq3), K1(k1), K2(k2),fTune(tune), ntrk(nt) {
        TH1::AddDirectory(kFALSE);
        int i,j;



        double q_l[3] = {0.9, 1.8, 1.8};
        double apps[3] = {0.3, 0.6, 0.6};

        double lchain[20];
        char pname[20][255];

        ///////////////////////////////////////////////////////////////////////////////////////

        if( fTune == kStd ){
            double dr_l[6] = {1.59, 1.25, 4.42, 1.50, 3.57};

            //  Standard Tune
            nelm = 9;
            chain[0] = makedrift( dr_l[0] );
            chain[1] = makequad( Bq1, apps[0], q_l[0] );
            chain[2] = makedrift( dr_l[1] );
            chain[3] = makequad( Bq2, apps[1], q_l[1] );
            chain[4] = makedrift( dr_l[2] );
            chain[5] = makedip( 45.0, 8.4, -1.25, -30.0, -30.0, K1, K2 );
            chain[6] = makedrift( dr_l[3] );
            chain[7] = makequad( Bq3, apps[2], q_l[2] );
            chain[8] = makedrift( dr_l[4] );

            double thischain[20] = {
                dr_l[0], q_l[0], dr_l[1], q_l[1], dr_l[2],
                8.4*TMath::Pi()/4, dr_l[3], q_l[2], dr_l[4]
            };



            char thisname[20][255] = { "Origin",
                "Q1 ent", "Q1 ext", "Q2 ent", "Q2 ext", "Dip ent", "Dip ext",
                "Q3 ent", "Q3 ext", "Focal Plane"
            };

            for( i = 0; i < 20; i++){
                lchain[i] = thischain[i];
                strcpy(pname[i], thisname[i]);
            }

        }

        ///////////////////////////////////////////////////////////////////////////////////////
        double th0 = -5.0*3.14159/180; // PREX central angle
        double th0_apex = -6.0*3.14159/180; // APEX central angle
        double th1 = -12.5*3.14159/180; // HRS angle

        double l_sept = 0.95; //  Matches SNAKE

        double sept_rho = (l_sept/cos(th0))/
            ( sin(th1-th0) - sin(th0) + cos(th1-th0)*sin(th0) );

        // Septum center is 1.7538 downstream of target

        double sept_face = 1.7538 - l_sept/2;
        double sept_exit = 1.7538 + l_sept/2;
        double sept_face_apex = 1.46 - l_sept/2;
        double sept_exit_apex = 1.46 + l_sept/2;
        double Q1_z = 2.3956;// From survey +
        double Q1coll_to_bore = 0.31 ; //  From HRS survey assuming spectrometers at 12.5 deg
                                       //  and survey

        double sept_exit_to_q1 = (Q1_z-sept_exit)/cos(th1);
        double sept_exit_to_q1_apex = (Q1_z-sept_exit_apex)/cos(th1);

        sept_K1 = 0.5;
        sept_K2 = 0.5;

        if( fTune == kPREX ){
            nelm = 15;
            double dr_l[8] = {sept_face, sept_exit_to_q1, Q1coll_to_bore, 1.25, 4.42, 1.50, 3.57, 1.43 };

            // Septum tune
            chain[0] = makedrift( dr_l[0] );
            chain[1] = swapxy();
            chain[2] = makedip( (th1-th0)*180/3.14159, sept_rho, 0, -th0*180/3.14159, th1*180/3.14159, sept_K1, sept_K2);
            chain[3] = swapxy(-1.0);
            chain[4] = makedrift( dr_l[1] );
            chain[5] = makedrift( dr_l[2] );

            chain[6] = makequad( Bq1, apps[0], q_l[0] );
            chain[7] = makedrift( dr_l[3] );
            chain[8] = makequad( Bq2, apps[1], q_l[1] );
            chain[9] = makedrift( dr_l[4] );
            chain[10] = makedip( 45.0, 8.4, -1.25, -30.0, -30.0, K1, K2 );
            chain[11] = makedrift( dr_l[5] );
            chain[12] = makequad( Bq3, apps[2], q_l[2] );
            chain[13] = makedrift( dr_l[6] );
            chain[14] = makedrift( dr_l[7] );


            double thischain[20] = {
                dr_l[0], 0.0, sept_rho*(th1-th0), 0.0,
                dr_l[1], dr_l[2],

                q_l[0], dr_l[3], q_l[1], dr_l[4],
                8.4*TMath::Pi()/4, dr_l[5], q_l[2], dr_l[6], dr_l[7]
            };

            char thisname[20][255] = { "Origin", "Septum Ent", "Rot 1", "Rot 2", "Septum Ext",
                "Q1 Coll",
                "Q1 ent", "Q1 ext", "Q2 ent", "Q2 ext", "Dip ent", "Dip ext",
                "Q3 ent", "Q3 ext", "Focal Plane"
            };

            for( i = 0; i < 20; i++){
                lchain[i] = thischain[i];
                strcpy(pname[i], thisname[i]);
            }
        }

        if( fTune == kAPEX ){
            nelm = 14;
            double dr_l[8] = {sept_face_apex, sept_exit_to_q1_apex, Q1coll_to_bore, 1.25, 4.42, 1.50, 3.57 };

            // Septum tune
            chain[0] = makedrift( dr_l[0] );
            chain[1] = swapxy();
            chain[2] = makedip( (th1-th0_apex)*180/3.14159, sept_rho, 0, -th0_apex*180/3.14159, th1*180/3.14159, sept_K1, sept_K2);
            chain[3] = swapxy(-1.0);
            chain[4] = makedrift( dr_l[1] );
            chain[5] = makedrift( dr_l[2] );

            chain[6] = makequad( Bq1, apps[0], q_l[0] );
            chain[7] = makedrift( dr_l[3] );
            chain[8] = makequad( Bq2, apps[1], q_l[1] );
            chain[9] = makedrift( dr_l[4] );
            chain[10] = makedip( 45.0, 8.4, -1.25, -30.0, -30.0, K1, K2 );
            chain[11] = makedrift( dr_l[5] );
            chain[12] = makequad( Bq3, apps[2], q_l[2] );
            chain[13] = makedrift( dr_l[6] );


            double thischain[20] = {
                dr_l[0], 0.0, sept_rho*(th1-th0_apex), 0.0,
                dr_l[1], dr_l[2],

                q_l[0], dr_l[3], q_l[1], dr_l[4],
                8.4*TMath::Pi()/4, dr_l[5], q_l[2], dr_l[6], dr_l[7]
            };

            char thisname[20][255] = { "Origin", "Septum Ent", "Rot 1", "Rot 2", "Septum Ext",
                "Q1 Coll",
                "Q1 ent", "Q1 ext", "Q2 ent", "Q2 ext", "Dip ent", "Dip ext",
                "Q3 ent", "Q3 ext", "Focal Plane"
            };

            for( i = 0; i < 20; i++){
                lchain[i] = thischain[i];
                strcpy(pname[i], thisname[i]);
            }
        }



        ///////////////////////////////////////////////////////////////////////////////////////




        larr[0] = 0.0;
        for( i = 0; i < nelm; i++ ){
            larr[i+1] = larr[i]+lchain[i];
        }

        fdp_rng= 0.1;
        fytg_rng = 0.1;
        fxtg_rng = 0.005;
        fthtg_rng = 0.15;
        fphtg_rng = 0.08;

        fdp_lim= 0.045;
        fytg_lim = 0.05;
        fxtg_lim = 0.005;
        fthtg_lim = 0.06;
        fphtg_lim = 0.028;

        gRandom->SetSeed(100);
        for( i = 0; i < ntrk; i++ ){
            ytg[i]  = gRandom->Uniform(-fytg_rng,fytg_rng);
            xtg[i]  = gRandom->Uniform(-fxtg_rng,fxtg_rng);
            thtg[i] = gRandom->Uniform(-fthtg_rng, fthtg_rng);
            phtg[i] = gRandom->Uniform(-fphtg_rng, fphtg_rng);
            dp[i]   = gRandom->Uniform(-fdp_rng,fdp_rng);
            acc[i] = true;
        }


        nbin = 40;

        for( i = 0; i < nelm+1; i++ ){
            hacc[0][0][i] = new TH1F(Form("hacc_%d_%d_%d",0, i+1,0), Form("X %s - All", pname[i]), nbin, -1.2, 1.2);
            hacc[0][1][i] = new TH1F(Form("hacc_%d_%d_%d",0, i+1,1), Form("X %s - Accepted", pname[i]), nbin, -1.2, 1.2);
            hacc[0][1][i]->GetXaxis()->SetTitle("x [m]");

            hacc[1][0][i] = new TH1F(Form("hacc_%d_%d_%d",1, i+1,0), Form("#theta %s - All", pname[i]), nbin, -0.2, 0.2);
            hacc[1][1][i] = new TH1F(Form("hacc_%d_%d_%d",1, i+1,1), Form("#theta %s - Accepted", pname[i]), nbin, -0.2, 0.2);
            hacc[1][1][i]->GetXaxis()->SetTitle("#theta");

            hacc[2][0][i] = new TH1F(Form("hacc_%d_%d_%d",2, i+1,0), Form("Y %s - All", pname[i]), nbin, -0.2, 0.2);
            hacc[2][1][i] = new TH1F(Form("hacc_%d_%d_%d",2, i+1,1), Form("Y %s - Accepted", pname[i]), nbin, -0.2, 0.2);
            hacc[2][1][i]->GetXaxis()->SetTitle("y [m]");

            hacc[3][0][i] = new TH1F(Form("hacc_%d_%d_%d",3, i+1,0), Form("#phi %s - All", pname[i]), nbin, -0.1, 0.1);
            hacc[3][1][i] = new TH1F(Form("hacc_%d_%d_%d",3, i+1,1), Form("#phi %s - Accepted", pname[i]), nbin, -0.1, 0.1);
            hacc[3][1][i]->GetXaxis()->SetTitle("#phi");

            hacc[4][0][i] = new TH1F(Form("hacc_%d_%d_%d",4, i+1,0), Form("#delta %s - All", pname[i]), nbin, -0.11, 0.11);
            hacc[4][1][i] = new TH1F(Form("hacc_%d_%d_%d",4, i+1,1), Form("#delta %s - Accepted", pname[i]), nbin, -0.11, 0.11);
            hacc[4][1][i]->GetXaxis()->SetTitle("#delta");

            for( j = 0; j < 5; j++ ){
                 hacc[j][1][i]->SetLineColor(kRed);
            }
        }
        
        DoTransport();
}

THRSTrans::~THRSTrans(){
    int i,j;
    for( i = 0; i < 5; i++ ){
        for( j = 0; j < nelm+1; j++ ){
            delete hacc[i][0][j];
            delete hacc[i][1][j];
            delete gacc[i][j];
        }
    }

    for( i = 0; i < nelm; i++ ){
        delete chain[i];
    }

    for( i = 0; i < nelm+1; i++ ){
        delete trans[i];
    }

    delete gx_th;
    delete gx_x0;
    delete gx_d;
    delete gy_ph;
    delete gy_y0;

    delete gph_y0;
    delete gph_ph;

}

void THRSTrans::DoTransport(){
    int i,j,k,bidx;
        double apps[3] = {0.3, 0.6, 0.6};

        double accdat[5][20][1000];
        double accpnt[5][20][1000];

   
        // Make all the transport matrices for each element

        trans[0] = new TMatrixD(20,20);
        for( j = 0; j < 20; j++ ){
            for( k = 0; k < 20; k++ ){
                if( j == k ){
                    (*trans[0])[j][k] = 1.0;
                } else {
                    (*trans[0])[j][k] = 0.0;
                }
            }
        }

        for( i = 0; i < nelm; i++ ){
            trans[i+1] =  new TMatrixD( (*chain[i])*(*trans[i]) ); 
        }


        for( i = 0; i < nelm+1; i++ ){
            // Pull out matrix elements for plotting
            x_th[i] = (*trans[i])[kX][kTh];
            x_x0[i] = (*trans[i])[kX][kX];
            x_d[i] = (*trans[i])[kX][kd];

            y_ph[i] = (*trans[i])[kY][kPh];
            y_y0[i] = (*trans[i])[kY][kY];

            ph_ph[i] = (*trans[i])[kPh][kPh];
            ph_y0[i] = (*trans[i])[kPh][kY];

            // Propagate base rays and test for acceptance

            for( j = 0; j < ntrk; j++ ){
                TVectorD iv(20), v(20);

                iv[kX]  =  xtg[j];
                iv[kTh] = thtg[j];
                iv[kY]  =  ytg[j];
                iv[kPh] = phtg[j];
                iv[kd]  =   dp[j];

                fillvector(iv);

                v = (*trans[i])*iv;

                hacc[0][0][i]->Fill(v[kX]);
                hacc[1][0][i]->Fill(v[kTh]);
                hacc[2][0][i]->Fill(v[kY]);
                hacc[3][0][i]->Fill(v[kPh]);
                hacc[4][0][i]->Fill(v[kd]);

                double xvdc, yvdc, zint;

                // Test acceptance
                if( fTune == kStd ){
                    //standard
                    switch(i){
                        case 1:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[0] ){
                                acc[j] = false;
                            }
                            break;
                        case 2:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[0] ){
                                acc[j] = false;
                            }
                            break;
                        case 3:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[1] ){
                                acc[j] = false;
                            }
                            break;
                        case 4:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[1] ){
                                acc[j] = false;
                            }
                            break;

                        case 5:
                            if( fabs(v[kY]) > 0.25/2  || fabs(v[kX])> 0.4  ){
                                acc[j] = false;
                            }
                            break;
                        case 6:
                            if( fabs(v[kY]) > 0.25/2  || fabs(v[kX])> 0.4  ){
                                acc[j] = false;
                            }
                            break;
                        case 7:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[2] ){
                                acc[j] = false;
                            }
                            break;
                        case 8:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[2] ){
                                acc[j] = false;
                            }
                            break;
                        case 9: // VDC
                            xvdc = v[kX]/(1.0 - v[kTh]);
                            zint = xvdc;
                            yvdc = v[kY] + v[kPh]*zint;

                            if( fabs(xvdc)>2.118/2 || fabs(yvdc)>0.288/2  ){
                                acc[j] = false;
                            }
                            break;
                        default:
                            break;
                    }
                }

                if( fTune == kPREX ){


                    // Septum

                    switch(i){
                        case 6:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[0] ){
                                acc[j] = false;
                            }
                            break;
                        case 7:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[0] ){
                                acc[j] = false;
                            }
                            break;
                        case 8:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[1] ){
                                acc[j] = false;
                            }
                            break;
                        case 9:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[1] ){
                                acc[j] = false;
                            }
                            break;

                        case 10:
                            if( fabs(v[kY]) > 0.25/2  || fabs(v[kX])> 0.4  ){
                                acc[j] = false;
                            }
                            break;
                        case 11:
                            if( fabs(v[kY]) > 0.25/2  || fabs(v[kX])> 0.4  ){
                                acc[j] = false;
                            }
                            break;
                        case 12:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[2] ){
                                acc[j] = false;
                            }
                            break;
                        case 13:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[2] ){
                                acc[j] = false;
                            }
                            break;
                        case 14: // VDC
                            xvdc = v[kX]/(1.0 - v[kTh]);
                            zint = xvdc;
                            yvdc = v[kY] + v[kPh]*zint;

                            if( fabs(xvdc)>2.118/2 || fabs(yvdc)>0.288/2  ){
                                acc[j] = false;
                            }
                            break;

                        default:
                            break;
                    }
                }


                if( fTune == kAPEX ){


                    // Septum

                    switch(i){
                        case 6:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[0] ){
                                acc[j] = false;
                            }
                            break;
                        case 7:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[0] ){
                                acc[j] = false;
                            }
                            break;
                        case 8:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[1] ){
                                acc[j] = false;
                            }
                            break;
                        case 9:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[1] ){
                                acc[j] = false;
                            }
                            break;

                        case 10:
                            if( fabs(v[kY]) > 0.25/2  || fabs(v[kX])> 0.4  ){
                                acc[j] = false;
                            }
                            break;
                        case 11:
                            if( fabs(v[kY]) > 0.25/2  || fabs(v[kX])> 0.4  ){
                                acc[j] = false;
                            }
                            break;
                        case 12:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[2] ){
                                acc[j] = false;
                            }
                            break;
                        case 13:
                            if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[2] ){
                                acc[j] = false;
                            }
                            break;
                        case 14: // VDC
                            xvdc = v[kX]/(1.0 - v[kTh]);
                            zint = xvdc;
                            yvdc = v[kY] + v[kPh]*zint;

                            if( fabs(xvdc)>2.118/2 || fabs(yvdc)>0.288/2  ){
                                acc[j] = false;
                            }
                            break;

                        default:
                            break;
                    }
                }


            }
        }

        for( i = 0; i < nelm+1; i++ ){
            for( j = 0; j < ntrk; j++ ){
                TVectorD iv(20), v(20);

                iv[kX]  =  xtg[j];
                iv[kTh] = thtg[j];
                iv[kY]  =  ytg[j];
                iv[kPh] = phtg[j];
                iv[kd]  =   dp[j];

                fillvector(iv);

                v = (*trans[i])*iv;

                if( acc[j]) {
                    hacc[0][1][i]->Fill(v[kX]);
                    hacc[1][1][i]->Fill(v[kTh]);
                    hacc[2][1][i]->Fill(v[kY]);
                    hacc[3][1][i]->Fill(v[kPh]);
                    hacc[4][1][i]->Fill(v[kd]);
                }

            }

            for( k = 0; k < 5; k++ ){
                for( bidx = 0; bidx < nbin; bidx++ ){
                    if( hacc[k][0][i]->GetBinContent(bidx+1) > 0 ){
                        accdat[k][i][bidx] =  hacc[k][1][i]->GetBinContent(bidx+1)/hacc[k][0][i]->GetBinContent(bidx+1);
                    } else {
                        accdat[k][i][bidx] = 0.0;
                    }
                    accpnt[k][i][bidx] = hacc[k][0][i]->GetBinCenter(bidx+1);
                }


                gacc[k][i] = new TGraph(nbin, accpnt[k][i], accdat[k][i]);
                gacc[k][i]->SetTitle(hacc[k][1][i]->GetTitle());
                gacc[k][i]->GetXaxis()->SetTitle(hacc[k][1][i]->GetXaxis()->GetTitle());
                gacc[k][i]->GetXaxis()->CenterTitle();
                if( 3 != k ){
                    gacc[k][i]->SetMaximum(0.30);
                } else {
                    gacc[k][i]->SetMaximum(0.65);
                }
            }


        }

        gx_th = new TGraph(nelm+1, larr, x_th);
        gx_th->SetLineColor(kBlue);
        gx_x0 = new TGraph(nelm+1, larr, x_x0);
        gx_d  = new TGraph(nelm+1, larr, x_d);
        gx_d->SetLineColor(kRed);
        gy_ph = new TGraph(nelm+1, larr, y_ph);
        gy_ph->SetLineColor(kBlue);
        gy_y0 = new TGraph(nelm+1, larr, y_y0);

        gph_y0 = new TGraph(nelm+1, larr, ph_y0);
        gph_ph = new TGraph(nelm+1, larr, ph_ph);
        gph_ph->SetLineColor(kMagenta);

        return;
}

void THRSTrans::ShowOutput(int idx, double l){
    int i;

    TMultiGraph *mgx = new TMultiGraph();
    mgx->Add(gx_th, "L");
    mgx->Add(gx_x0, "L");
    mgx->Add(gx_d, "L");


    TLegend *leg = new TLegend(0.14, 0.64, 0.43, 0.86);
    leg->AddEntry(gx_x0, "(x|x_{tg})", "l");
    leg->AddEntry(gx_th, "(x|#theta_{tg})", "l");
    leg->AddEntry(gx_d, "(x|#delta)", "l");
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0);

    TMultiGraph *mgy = new TMultiGraph();
    mgy->Add(gy_ph, "L");
    mgy->Add(gy_y0, "L");
//    mgy->Add(gph_ph, "L");


    TLegend *legy = new TLegend(0.64, 0.70, 0.87, 0.85);
    legy->AddEntry(gy_y0, "(y|y_{tg})", "l");
    legy->AddEntry(gy_ph, "(y|#phi_{tg})", "l");
    legy->SetFillColor(kWhite);
    legy->SetBorderSize(0);



    TMultiGraph *mgph = new TMultiGraph();
    mgph->Add(gph_ph, "L");
    mgph->Add(gph_y0, "L");

    TCanvas *cx =new TCanvas();
    cx->SetGridx();
    cx->SetGridy();
    mgx->Draw("A");
    mgx->SetTitle("x Couplings");
    mgx->GetXaxis()->SetTitle("t [m]");
    mgx->GetXaxis()->CenterTitle();
    mgx->Draw("A");

    leg->Draw();


    TCanvas *cy = new TCanvas();
    cy->SetGridx();
    cy->SetGridy();
    mgy->Draw("A");
    mgy->SetTitle("y Couplings");
    mgy->GetXaxis()->SetTitle("t [m]");
    mgy->GetXaxis()->CenterTitle();
    mgy->Draw("A");
    legy->Draw();

    TCanvas *c = new TCanvas();
    /*
    c->Divide(5,2);
    for( i = 0; i < 5; i++ ){
        c->cd(i+1);
        hacc[i][0][0]->Draw();
        hacc[i][1][0]->Draw("same");
        c->cd(i+6);
        gacc[i][0]->Draw("AC");
    }
    */
    c->Divide(2,2);
    for( i = 1; i < 5; i++ ){
        c->cd(i);
//        hacc[i][0][0]->Draw();
//        hacc[i][1][0]->Draw("same");
        gacc[i][0]->Draw("AC");
    }


    if( idx == -1 ){
        PrintSimple(trans[nelm], l);
    } else {
        PrintSimple(trans[idx], l);
    }

    printf("%f, %f, %f, %f, %f\n", Bq1, Bq2, Bq3, K1, K2);
}

void THRSTrans::ShowFocalLengths(){
    double kq1 = sqrt(fabs(Bq1)/0.3);
    double kq2 = sqrt(fabs(Bq2)/0.3);
    double kq3 = sqrt(fabs(Bq3)/0.3);
    printf("Q1 ");
    if( Bq1 > 0 ){
        printf("x Foc ");
    } else {
        printf("x Def ");
    }
    printf("%5.2f m\n", 1.0/(kq1*sin(kq1*0.8)) );

    printf("Q2 ");
    if( Bq2 > 0 ){
        printf("x Foc ");
    } else {
        printf("x Def ");
    }
    printf("%5.2f m\n", 1.0/(kq2*sin(kq2*0.8)) );

    printf("Q3 ");
    if( Bq3 > 0 ){
        printf("x Foc ");
    } else {
        printf("x Def ");
    }
    printf("%5.2f m\n", 1.0/(kq3*sin(kq3*0.8)) );

}

void THRSTrans::ShowAcc(){
    double dsigma = 4.0*fthtg_rng*fphtg_rng/ntrk;

    dsigma *= (fytg_rng/fytg_lim)*(fdp_rng/fdp_lim);

    double sigmasum = 0.0;
    int i;
    for( i = 0; i < ntrk; i++ ){
        if( acc[i] && fabs(ytg[i])<fytg_lim && fabs(dp[i])<fdp_lim ){
            sigmasum += dsigma;
        }
    }

    printf("Accepted %f mrad\n", sigmasum*1e3);
}

TMatrixD *THRSTrans::makedrift(double l ){
    TMatrixD *m = new TMatrixD(20,20);
   (*m)[kX][kX] = 1.0;
    (*m)[kX][kTh] = l;
    (*m)[kTh][kTh] = 1.0;
    (*m)[kd][kd] = 1.0;

    (*m)[kY][kY] = 1.0;
    (*m)[kY][kPh] = l;
    (*m)[kPh][kPh] = 1.0;

    setcrossterms(m);

    return m;
}

TMatrixD *THRSTrans::makequad( double B0, double a, double l   ){
    TMatrixD *m = new TMatrixD(20,20);

    int i;

    double kq2 = fabs(B0)/a;
    double kq = sqrt(kq2);

    if( B0 > 0 ){
        (*m)[kX][kX]   = cos(kq*l);
        (*m)[kX][kTh]  = sin(kq*l)/kq;
        (*m)[kTh][kX]  = -kq*sin(kq*l);
        (*m)[kTh][kTh] = cos(kq*l);

        (*m)[kY][kY]   = cosh(kq*l);
        (*m)[kY][kPh]  = sinh(kq*l)/kq;
        (*m)[kPh][kY]  = kq*sinh(kq*l);
        (*m)[kPh][kPh] = cosh(kq*l);



        (*m)[kX][kXd] = 0.5*kq*l*sin(kq*l);
        (*m)[kX][kThd] = 0.5*sin(kq*l)/kq - 0.5*l*cos(kq*l);

        (*m)[kTh][kXd] = 0.5*kq*( kq*l*cos(kq*l) + sin(kq*l) );
        (*m)[kTh][kThd] = 0.5*kq*l*sin(kq*l);

        (*m)[kY][kYd] = -0.5*kq*l*sinh(kq*l);
        (*m)[kY][kPhd] = 0.5*(sinh(kq*l)/kq - l*cosh(kq*l));

        (*m)[kPh][kYd] =  -kq*0.5*(cosh(kq*l)*kq *l + sinh(kq*l));
        (*m)[kPh][kPhd] = -0.5*kq*l*sinh(kq*l);
    } else {
        (*m)[kY][kY] = cos(kq*l);
        (*m)[kY][kPh] = sin(kq*l)/kq;

        (*m)[kPh][kY] = -kq*sin(kq*l);
        (*m)[kPh][kPh] = cos(kq*l);

        (*m)[kX][kX] = cosh(kq*l);
        (*m)[kX][kTh] = sinh(kq*l)/kq;

        (*m)[kTh][kX] =  kq*sinh(kq*l);
        (*m)[kTh][kTh] = cosh(kq*l);

        (*m)[kY][kYd] = 0.5*kq*l*sin(kq*l);
        (*m)[kY][kPhd] = 0.5*sin(kq*l)/kq - 0.5*l*cos(kq*l);

        (*m)[kPh][kYd] = 0.5*kq*( kq*l*cos(kq*l) + sin(kq*l) );
        (*m)[kPh][kPhd] = 0.5*kq*l*sin(kq*l);

        (*m)[kX][kXd] = -0.5*kq*l*sinh(kq*l);
        (*m)[kX][kThd] = 0.5*(sinh(kq*l)/kq - l*cosh(kq*l));

        (*m)[kTh][kXd] =  -kq*0.5*(cosh(kq*l)*kq *l + sinh(kq*l));
        (*m)[kTh][kThd] = -0.5*kq*l*sinh(kq*l);
    }

    (*m)[kd][kd] = 1.0;

    setcrossterms(m);

   
    return m;
}
TMatrixD *THRSTrans::makedip( double theta, double rho, double n, double beta1, double beta2, double ps1, double ps2){
    double b1 = beta1*3.14159/180;
    double b2 = beta2*3.14159/180;


    int i,j;

    double h = 1./rho;
    double kx2 = (1.0-n)*h*h;
    double ky2 = n*h*h;

    double psi1 = ps1*h*0.25*(1.0+sin(b1)*sin(b1))/cos(b1);
    double psi2 = ps2*h*0.25*(1.0+sin(b2)*sin(b2))/cos(b2);

    double kx = sqrt(kx2);
    double ky = sqrt(-ky2);

    double l = theta*rho*3.14159/180;

    // En fringe
    TMatrixD *men = new TMatrixD(20,20);
    (*men)[kX][kX] = 1.0;
    (*men)[kX][kX2] = -h*tan(b1)*tan(b1)/2.0;
    (*men)[kX][kY2] = h/(cos(b1)*cos(b1))/2.0;
    (*men)[kTh][kX] = h*tan(b1);
    (*men)[kTh][kTh] = 1.0;
    (*men)[kTh][kX2] = -n*h*h*tan(b1);
    (*men)[kTh][kXTh] = h*tan(b1)*tan(b1);
    (*men)[kTh][kXd] = -h*tan(b1);
    (*men)[kTh][kY2] =  h*h*(n+0.5+tan(b1)*tan(b1))*tan(b1);
    (*men)[kTh][kYPh] = -h*tan(b1)*tan(b1);

    (*men)[kY][kY] = 1.0;
    (*men)[kY][kXY] = h*tan(b1)*tan(b1);

    (*men)[kPh][kY] = -h*tan(b1-ps1);
    (*men)[kPh][kPh] = 1.0;

    (*men)[kPh][kXY] = 2*h*h*n*tan(b1);
    (*men)[kPh][kXPh] = -h*tan(b1)*tan(b1);
    (*men)[kPh][kThY] = -h/(cos(b1)*cos(b1));

    (*men)[kPh][kYd] = h*tan(b1)-h*ps1/(cos(b1-ps1)*cos(b1-ps1));

    (*men)[kd][kd]= 1.0;


    setcrossterms(men);

    // Ex fringe
    TMatrixD *mex = new TMatrixD(20,20);
    (*mex)[kX][kX] = 1.0;
    (*mex)[kX][kX2] = h*tan(b2)*tan(b2)/2.0;
    (*mex)[kX][kY2] = -h/(cos(b2)*cos(b2))/2.0;
    (*mex)[kTh][kX] = h*tan(b2);
    (*mex)[kTh][kTh] = 1.0;
    (*mex)[kTh][kX2] = -h*h*(n+0.5*tan(b2)*tan(b2))*tan(b2);

    (*mex)[kTh][kXTh] = -h*tan(b2)*tan(b2);
    (*mex)[kTh][kXd] = -h*tan(b2);
    (*mex)[kTh][kY2] =  h*h*(n-0.5*tan(b2)*tan(b2))*tan(b2);
    (*mex)[kTh][kYPh] = h*tan(b2)*tan(b2);

    (*mex)[kY][kY] = 1.0;
    (*mex)[kY][kXY] = -h*tan(b2)*tan(b2);

    (*mex)[kPh][kY] = -h*tan(b2-ps2);
    (*mex)[kPh][kPh] = 1.0;

    (*mex)[kPh][kXY] = h*h*(2*n+1.0/(cos(b2)*cos(b2)))*tan(b2);

    (*mex)[kPh][kXPh] = h*tan(b2)*tan(b2);
    (*mex)[kPh][kThY] = h/(cos(b2)*cos(b2));

    (*mex)[kPh][kYd] = h*tan(b2)-h*ps2/(cos(b2-ps2)*cos(b2-ps2));

    (*mex)[kd][kd] = 1.0;

    setcrossterms(mex);



    //   Main dipole field
    TMatrixD *mdip = new TMatrixD(20,20);

    double cxt, sxt, cpxt, spxt, dxt, dpxt, cyt, syt;
    double cpyt, spyt;

    if( kx2 > 0 ){
        cxt = cos(kx*l);
        sxt = sin(kx*l)/kx;
    } else {
        if( kx > 0 ){
            cxt = cosh(kx*l);
            sxt = sinh(kx*l)/kx;
        } else {
            cxt = 1.0;
            sxt = l;
        }
    }

    if( ky2 > 0 ){
        cyt = cos(ky*l);
        syt = sin(ky*l)/ky;
    } else {
        if( ky > 0 ){
            cyt = cosh(ky*l);
            syt = sinh(ky*l)/ky;
        } else {
            cyt = 1.0;
            syt = l;
        }
    }

    cpxt = -kx2*sxt;
    spxt = cxt;
   
    cpyt = -ky2*syt;
    spyt = cyt;

    dxt = h*(1.0-cxt)/kx2;
    dpxt = h*sxt;

    double I10 = dxt/h;
    double I11 = 0.5*l*sxt;
    double I12 = (sxt-l*cxt)/(2.0*kx2);
    double I16 = h*(I10 - I11)/kx2;

    double I111 = (sxt*sxt + I10)/3.0;
    double I112 = sxt*I10/3.0;
    double I116 = h*(I11 - I111)/kx2;
    double I122 = (2*I10 - sxt*sxt)/(3.0*kx2);
    double I126 = h*(I12 - I112)/kx2;
    double I166 = h*h*( 4.0*I10/3.0 + sxt*sxt/3.0 - l*sxt)/kx2/kx2;

    double I20 = sxt;
    double I21 = 0.5*(sxt+l*cxt);
    double I22 = I11;
    double I26 = 0.5*h*(sxt - l*cxt)/kx2;

    double I211 = sxt*(1.0+2.0*cxt)/3.0;
    double I212 = (2.0*sxt*sxt - dxt/h)/3.0;
    double I216 = h*(0.5*l*cxt + sxt/6.0 - 2.0*sxt*cxt/3.0)/kx2;
    double I222 = 2.0*sxt*dxt/(h*3.0);
    double I226 = h*(0.5*l*sxt - 2.0*sxt*sxt/3.0 + dxt/(h*3.0))/kx2;
    double I266 = h*h*(sxt/3.0 + 2.0*sxt*cxt/3.0 - l*cxt)/kx2/kx2;

    double I33 = 0.5*l*syt;
    double I34;
    if( fabs(ky2)<1e-5 ){
        I34 = 0.5*l*l*l;
    } else {
        I34 = 0.5*(syt - l*cyt)/ky2;
    }

    double I313 = (kx2*cyt*dxt/h - 2.0*ky2*sxt*syt)/(kx2 - 4.0*ky2);
    double I314 = (2.0*sxt*cyt - syt*(1.0+cxt))/(kx2 - 4.0*ky2);
    double I323 = (2.0*ky2*syt*(1.0+cxt)/kx2 - sxt*cyt)/(kx2 - 4.0*ky2) +syt/kx2;
    double I324 = ( 2.0*cyt*dxt/h - sxt*syt )/(kx2 - 4*ky2);
    double I336 = h*(0.5*l*syt - ( cyt*(1.0-cxt) - 2.0*ky2*sxt*syt)/(kx2 - 4.0*ky2))/kx2;
    double I346 = h*( I34 - (2.0*sxt*cyt - syt*(1.0+cxt))/(kx2 - 4.0*ky2) )/kx2;
    
    double I43 = 0.5*(syt + l*cyt);
    double I44 = I33;

    double I413 = ( (kx2 - 2.0*ky2)*sxt*cyt - ky2*syt*(1.0+cxt))/(kx2 - 4.0*ky2);
    double I414 = ( (kx2 - 2.0*ky2)*sxt*syt - cyt*(1.0-cxt))/(kx2 - 4.0*ky2);
    double I423 = ( 2.0*ky2*cyt*(1.0+cxt)/kx2 - cxt*cyt - ky2*sxt*syt )/(kx2 - 4.0*ky2);
    double I424 = ( cyt*sxt - cxt*syt - ky2*syt*dxt/h  )/(kx2 - 4.0*ky2);
    double I436 = h*(0.5*l*cyt + 0.5*syt + (ky2*syt*(1.0+cxt) - (kx2-2.0*ky2)*sxt*cyt)/(kx2-4.0*ky2) )/kx2;
    double I446 = h*(0.5*l*syt - ( (kx2-2.0*ky2)*sxt*syt - cyt*(1.0-cxt)  )/(kx2-4.0*ky2) )/kx2;


    (*mdip)[kX][kX] = cxt;
    (*mdip)[kX][kTh] = sxt;
    (*mdip)[kX][kd] = dxt;

    (*mdip)[kX][kX2] = (2.0*n-1)*h*h*h*I111 + 0.5*kx2*kx2*h*I122;
    (*mdip)[kX][kXTh] = h*sxt + 2.0*(2.0*n-1)*h*h*h*I112 - kx2*h*I112;
    (*mdip)[kX][kXd] = (2.0-n)*h*h*I11 + 2.0*(2.0*n-1)*h*h*h*I116 - kx2*h*h*I122;
    (*mdip)[kX][kTh2] =  (2.0*n-1)*h*h*h*I122 + 0.5*h*I111;
    (*mdip)[kX][kThd] =  (2.0-n)*h*h*I12 + 2.0*(2.0*n-1)*h*h*h*I126 + h*h*I112;

    (*mdip)[kX][kd2] =  -h*I10 + (2.0-n)*h*h*I16 + (2.0*n-1)*h*h*h*I166 + 0.5*h*h*h*I122;
    (*mdip)[kX][kY2] =  -0.5*ky2*h*I10;
    (*mdip)[kX][kPh2] =  -0.5*h*I10;

    (*mdip)[kTh][kX] = cpxt;
    (*mdip)[kTh][kTh] = cxt;
    (*mdip)[kTh][kd] = h*sxt;

    (*mdip)[kTh][kX2] = (2.0*n-1)*h*h*h*I211 + 0.5*kx2*kx2*h*I222 - h*cxt*cpxt;
    (*mdip)[kTh][kXTh] = h*spxt + 2.0*(2.0*n-1)*h*h*h*I212 - kx2*h*I212 - h*(cxt*spxt + cpxt*sxt);
    (*mdip)[kTh][kXd] = (2.0-n)*h*h*I21 + 2.0*(2.0*n-1)*h*h*h*I216 -kx2*h*h*I222 - h*(cxt*dpxt + cpxt*dxt);
    (*mdip)[kTh][kTh2] = (2.0*n-1)*h*h*h*I222 + 0.5*h*I211 - h*sxt*spxt;
    (*mdip)[kTh][kThd] = (2.0-n)*h*h*I22 + 2.0*(2.0*n-1)*h*h*h*I226 + h*h*I212 - h*(sxt*dpxt + spxt*dxt);

    (*mdip)[kTh][kd2] = -h*I20 + (2.0-n)*h*h*I26 + (2.0*n-1)*h*h*h*I266 + 0.5*h*h*h*I222 - h*dxt*dpxt;
    (*mdip)[kTh][kY2] = -0.5*ky2*h*I20; 
    (*mdip)[kTh][kPh2] = -0.5*h*I20;

    (*mdip)[kY][kY] = cyt;
    (*mdip)[kY][kPh] = syt;
    (*mdip)[kY][kXY] = -2.0*n*h*h*h*I313+kx2*ky2*h*I324;
    (*mdip)[kY][kXPh] = h*syt - 2.0*n*h*h*h*I314 - kx2*h*I323;
    (*mdip)[kY][kThY] = -2.0*n*h*h*h*I323 - ky2*h*I314;
    (*mdip)[kY][kThPh] = -2.0*n*h*h*h*I324 + h*I313;
    (*mdip)[kY][kYd] = ky2*I33 - 2.0*n*h*h*h*I336 - ky2*h*h*I324;
    (*mdip)[kY][kPhd] = ky2*I34 - 2.0*n*h*h*h*I346 + h*h*I323;


    (*mdip)[kPh][kY] = cpyt;
    (*mdip)[kPh][kPh] = spyt;
    (*mdip)[kPh][kXY] = -2.0*n*h*h*h*I413+kx2*ky2*h*I424 - h*cxt*cpyt;
    (*mdip)[kPh][kXPh] = h*spyt - 2.0*n*h*h*h*I414 - kx2*h*I423 - h*cxt*spyt;
    (*mdip)[kPh][kThY] = -2.0*n*h*h*h*I423 - ky2*h*I414 - h*sxt*cpyt;
    (*mdip)[kPh][kThPh] = -2.0*n*h*h*h*I424 + h*I413 - h*sxt*spyt;
    (*mdip)[kPh][kYd] = ky2*I43 - 2.0*n*h*h*h*I436 - ky2*h*h*I424 - h*dxt*cpyt;
    (*mdip)[kPh][kPhd] = ky2*I44 - 2.0*n*h*h*h*I446 + h*h*I423 - h*dxt*spyt;

    (*mdip)[kd][kd] = 1.0;

    setcrossterms(mdip);

    TMatrixD *m = new TMatrixD((*mex)*(*mdip)*(*men));

    delete mex;
    delete mdip;
    delete men;

    return m;
}



void THRSTrans::PrintSimple(TMatrixD *m, double l){
    int i, j;

    int map[5] = {kX, kTh, kY, kPh, kd };

    TMatrixD *newm  = new TMatrixD(5,5);
    TMatrixD *drift = new TMatrixD(5,5);

    for( i = 0; i < 5; i++ ){
        for( j = 0; j < 5; j++ ){
            (*newm)[i][j] = (*m)[map[i]][map[j]];
            if( i == j ){
                (*drift)[i][j] = 1.0;
            } else {
                (*drift)[i][j] = 0.0;
            }
        }
    }
    (*drift)[0][1] = l;
    (*drift)[3][4] = l;

    (*newm) = (*drift)*(*newm);

    newm->Print("f=%6.2f");

    delete newm;
    delete drift;
    return;
}

TMatrixD *THRSTrans::GetOptics(int idx){
    int i, j;

    TMatrixD *m;
    if( idx != -1 ){
        m = trans[idx];
    } else {
        m = trans[nelm];
    }

    int map[5] = {kX, kTh, kY, kPh, kd };

    TMatrixD *newm  = new TMatrixD(4,4);

    for( i = 0; i < 4; i++ ){
        for( j = 1; j < 5; j++ ){
            (*newm)[i][j-1] = (*m)[map[i]][map[j]];
        }
    }

    newm->Invert();

    return newm;
}




TMatrixD *THRSTrans::swapxy(double sign){
    // Sign + (RHRS convention)
    // x-> y
    // y->-x

    if( sign == 0.0 ) exit(1);

    sign /= fabs(sign); // normalize to 1

    TMatrixD *m = new TMatrixD(20,20);
    (*m)[kX][kY] =   -1.0*sign;
    (*m)[kTh][kPh] = -1.0*sign;
    (*m)[kY][kX]   = sign;
    (*m)[kPh][kTh] = sign;
    (*m)[kd][kd]   = 1.0;

    (*m)[kX2][kY2]   = 1.0;
    (*m)[kXTh][kYPh] = 1.0;
    (*m)[kXd][kYd]   = -1.0*sign;
    (*m)[kTh2][kPh2]   = 1.0;
    (*m)[kThd][kPhd]   = -1.0*sign;
    (*m)[kd2][kd2]   = 1.0;
    (*m)[kY2][kX2]   = 1.0;
    (*m)[kYPh][kXTh]   = 1.0;
    (*m)[kPh2][kTh2]   = 1.0;
    (*m)[kXY][kXY]     = -1.0;
    (*m)[kXPh][kThY]   = -1.0;
    (*m)[kThY][kXPh]   = -1.0;
    (*m)[kThPh][kThPh] = -1.0;
    (*m)[kYd][kXd]     = sign;
    (*m)[kPhd][kThd]   = sign;

    return m;
}



void THRSTrans::fillvector(TVectorD &v){
    v[kX2] = v[kX]*v[kX];
    v[kXTh] = v[kX]*v[kTh];
    v[kXd] = v[kX]*v[kd];
    v[kTh2] = v[kTh]*v[kTh];
    v[kThd] = v[kTh]*v[kd];
    v[kd2] = v[kd]*v[kd];
    v[kY2] = v[kY]*v[kY];
    v[kYPh] = v[kY]*v[kPh];
    v[kPh2] = v[kPh]*v[kPh];
    v[kXY] = v[kX]*v[kY];
    v[kXPh] = v[kX]*v[kPh];
    v[kThY] = v[kTh]*v[kY];
    v[kThPh] = v[kTh]*v[kPh];
    v[kYd] = v[kY]*v[kd];
    v[kPhd] = v[kPh]*v[kd];
    return;
}

void THRSTrans::setcrossterms(TMatrixD *m){
    (*m)[kX2][kX2] = (*m)[kX][kX]*(*m)[kX][kX];
    (*m)[kX2][kXTh] = 2.0*(*m)[kX][kX]*(*m)[kX][kTh];
    (*m)[kX2][kXd] = 2.0*(*m)[kX][kX]*(*m)[kX][kd];
    (*m)[kX2][kTh2] = (*m)[kX][kTh]*(*m)[kX][kTh];
    (*m)[kX2][kThd] = 2.0*(*m)[kX][kTh]*(*m)[kX][kd];
    (*m)[kX2][kd2] = (*m)[kX][kd]*(*m)[kX][kd];
    (*m)[kXTh][kX2] = (*m)[kX][kX]*(*m)[kTh][kX];
    (*m)[kXTh][kXTh] = (*m)[kX][kX]*(*m)[kTh][kTh] +
        (*m)[kX][kTh]*(*m)[kTh][kX];
    (*m)[kXTh][kXd] = (*m)[kX][kX]*(*m)[kTh][kd] +
        (*m)[kX][kd]*(*m)[kTh][kX];
    (*m)[kXTh][kTh2] = (*m)[kX][kTh]*(*m)[kTh][kTh];
    (*m)[kXTh][kThd] = (*m)[kX][kTh]*(*m)[kTh][kd] +
        (*m)[kX][kd]*(*m)[kTh][kTh];
    (*m)[kXTh][kd2] = (*m)[kX][kd]*(*m)[kTh][kd];
    (*m)[kXd][kX2] = (*m)[kX][kX]*(*m)[kd][kX];
    (*m)[kXd][kXTh] = (*m)[kX][kX]*(*m)[kd][kTh] +
        (*m)[kX][kTh]*(*m)[kd][kX];
    (*m)[kXd][kXd] = (*m)[kX][kX]*(*m)[kd][kd] +
        (*m)[kX][kd]*(*m)[kd][kX];
    (*m)[kXd][kTh2] = (*m)[kX][kTh]*(*m)[kd][kTh];

    (*m)[kXd][kThd] = (*m)[kX][kTh]*(*m)[kd][kd] +
        (*m)[kX][kd]*(*m)[kd][kTh];
    (*m)[kXd][kd2] = (*m)[kX][kd]*(*m)[kd][kd];
    (*m)[kTh2][kX2] = (*m)[kTh][kX]*(*m)[kTh][kX];
    (*m)[kTh2][kXTh] = 2.0*(*m)[kTh][kX]*(*m)[kTh][kTh];
    (*m)[kTh2][kXd] = 2.0*(*m)[kTh][kX]*(*m)[kTh][kd];
    (*m)[kTh2][kTh2] = (*m)[kTh][kTh]*(*m)[kTh][kTh];
    (*m)[kTh2][kThd] = 2.0*(*m)[kTh][kTh]*(*m)[kTh][kd];
    (*m)[kTh2][kd2] = (*m)[kTh][kd]*(*m)[kTh][kd];
    (*m)[kThd][kX2] = (*m)[kTh][kX]*(*m)[kd][kX];
    (*m)[kThd][kXTh] = (*m)[kTh][kX]*(*m)[kd][kTh] +
        (*m)[kTh][kTh]*(*m)[kd][kX];
    (*m)[kThd][kXd] = (*m)[kTh][kX]*(*m)[kd][kd] +
                      (*m)[kTh][kd]*(*m)[kd][kX];

    (*m)[kThd][kTh2] = (*m)[kTh][kTh]*(*m)[kd][kTh];

    (*m)[kThd][kThd] = (*m)[kTh][kTh]*(*m)[kd][kd] +
        (*m)[kTh][kd]*(*m)[kd][kTh];

    (*m)[kThd][kd2] = (*m)[kTh][kd]*(*m)[kd][kd];
    (*m)[kd2][kX2] = (*m)[kd][kX]*(*m)[kd][kX];
    (*m)[kd2][kXTh] = 2.0*(*m)[kd][kX]*(*m)[kd][kTh];

    (*m)[kd2][kXd] = 2.0*(*m)[kd][kX]*(*m)[kd][kd];

    (*m)[kd2][kTh2] = (*m)[kd][kTh]*(*m)[kd][kTh];
    (*m)[kd2][kThd] = 2.0*(*m)[kd][kTh]*(*m)[kd][kd];

    (*m)[kd2][kd2] = (*m)[kd][kd]*(*m)[kd][kd];
    (*m)[kY2][kY2] = (*m)[kY][kY]*(*m)[kY][kY];

    (*m)[kY2][kYPh] = 2.0*(*m)[kY][kY]*(*m)[kY][kPh];
    (*m)[kY2][kPh2] = (*m)[kY][kPh]*(*m)[kY][kPh];

    (*m)[kYPh][kY2] = (*m)[kY][kY]*(*m)[kPh][kY];

    (*m)[kYPh][kYPh] = (*m)[kY][kY]*(*m)[kPh][kPh] +
                       (*m)[kY][kPh]*(*m)[kPh][kY];
    (*m)[kYPh][kPh2] = (*m)[kY][kPh]*(*m)[kPh][kPh];

    (*m)[kPh2][kY2] = (*m)[kPh][kY]*(*m)[kPh][kY];
    (*m)[kPh2][kYPh] = 2.0*(*m)[kPh][kY]*(*m)[kPh][kPh];
    
    (*m)[kPh2][kPh2] = (*m)[kPh][kPh]*(*m)[kPh][kPh];

    (*m)[kXY][kXY] = (*m)[kX][kX]*(*m)[kY][kY];

    (*m)[kXY][kXPh] = (*m)[kX][kX]*(*m)[kY][kPh];
    (*m)[kXY][kThY] = (*m)[kX][kTh]*(*m)[kY][kY];
    (*m)[kXY][kThPh] = (*m)[kX][kTh]*(*m)[kY][kPh];
    (*m)[kXY][kYd] = (*m)[kX][kd]*(*m)[kY][kY];
    (*m)[kXY][kPhd] = (*m)[kX][kd]*(*m)[kY][kPh];
    (*m)[kXPh][kXY] = (*m)[kX][kX]*(*m)[kPh][kY];
    (*m)[kXPh][kXPh] = (*m)[kX][kX]*(*m)[kPh][kPh];
    (*m)[kXPh][kThY] = (*m)[kX][kTh]*(*m)[kPh][kY];
    (*m)[kXPh][kThPh] = (*m)[kX][kTh]*(*m)[kPh][kPh];
    (*m)[kXPh][kYd] = (*m)[kX][kd]*(*m)[kPh][kY];
    (*m)[kXPh][kPhd] = (*m)[kX][kd]*(*m)[kPh][kPh];
    (*m)[kThY][kXY] = (*m)[kTh][kX]*(*m)[kY][kY];
    (*m)[kThY][kXPh] = (*m)[kTh][kX]*(*m)[kY][kPh];
    (*m)[kThY][kThY] = (*m)[kTh][kTh]*(*m)[kY][kY];
    (*m)[kThY][kThPh] = (*m)[kTh][kTh]*(*m)[kY][kPh];
    (*m)[kThY][kYd] = (*m)[kTh][kd]*(*m)[kY][kY];
    (*m)[kThY][kPhd] = (*m)[kTh][kd]*(*m)[kY][kPh];
    (*m)[kThPh][kXY] = (*m)[kTh][kX]*(*m)[kPh][kY];
    (*m)[kThPh][kXPh] = (*m)[kTh][kX]*(*m)[kPh][kPh];
    (*m)[kThPh][kThY] = (*m)[kTh][kTh]*(*m)[kPh][kY];
    (*m)[kThPh][kThPh] = (*m)[kTh][kTh]*(*m)[kPh][kPh];
    (*m)[kThPh][kYd] = (*m)[kTh][kd]*(*m)[kPh][kY];
    (*m)[kThPh][kPhd] = (*m)[kTh][kd]*(*m)[kPh][kPh];
    (*m)[kYd][kXY] = (*m)[kd][kX]*(*m)[kY][kY];
    (*m)[kYd][kXPh] = (*m)[kd][kX]*(*m)[kY][kPh];
    (*m)[kYd][kThY] = (*m)[kd][kTh]*(*m)[kY][kY];
    (*m)[kYd][kThPh] = (*m)[kd][kTh]*(*m)[kY][kPh];
    (*m)[kYd][kYd] = (*m)[kd][kd]*(*m)[kY][kY];
    (*m)[kYd][kPhd] = (*m)[kd][kd]*(*m)[kY][kPh];
    (*m)[kPhd][kXY] = (*m)[kd][kX]*(*m)[kPh][kY];
    (*m)[kPhd][kXPh] = (*m)[kd][kX]*(*m)[kPh][kPh];
    (*m)[kPhd][kThY] = (*m)[kd][kTh]*(*m)[kPh][kY];
    (*m)[kPhd][kThPh] = (*m)[kd][kTh]*(*m)[kPh][kPh];
    (*m)[kPhd][kYd] = (*m)[kd][kd]*(*m)[kPh][kY];
    (*m)[kPhd][kPhd] = (*m)[kd][kd]*(*m)[kPh][kPh];
}




