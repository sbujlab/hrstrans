#include "THRSTrans.h"
#include "TMatrixD.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVectorD.h"


THRSTrans::THRSTrans(double bq1, double bq2, double bq3, double ps ):
    Bq1(bq1),Bq2(bq2),Bq3(bq3),psi(ps) {
        TH1::AddDirectory(kFALSE);
        int i,j;

        double dr_l[6] = {1.59, 1.25, 4.42, 1.50, 3.57, 1.16 };
        double q_l[3] = {0.9, 1.8, 1.8};
        double apps[3] = {0.3, 0.6, 0.6};

        nelm = 9;

        chain[0] = makedrift( dr_l[0] );
        chain[1] = makequad( Bq1, apps[0], q_l[0] );
        chain[2] = makedrift( dr_l[1] );
        chain[3] = makequad( Bq2, apps[1], q_l[1] );
        chain[4] = makedrift( dr_l[2] );
        chain[5] = makedip( 45.0, 8.4, -1.25, -30.0, -30.0, psi );
        chain[6] = makedrift( dr_l[3] );
        chain[7] = makequad( Bq3, apps[2], q_l[2] );
        chain[8] = makedrift( dr_l[4] );

        double lchain[20] = {
            dr_l[0], q_l[0], dr_l[1], q_l[1], dr_l[2],
            8.4*TMath::Pi()/4, dr_l[3], q_l[2], dr_l[4]
        };


        larr[0] = 0.0;
        for( i = 0; i < nelm; i++ ){
            larr[i+1] = larr[i]+lchain[i];
        }


        ntrk = 40000;
        gRandom->SetSeed(100);
        for( i = 0; i < ntrk; i++ ){
            ytg[i]  = gRandom->Uniform(-0.2/2.0, 0.2/2.0);
            xtg[i]  = gRandom->Uniform(-0.1,0.1);
            thtg[i] = gRandom->Uniform(-0.1, 0.1);
            phtg[i] = gRandom->Uniform(-0.08, 0.08);
            dp[i]   = gRandom->Uniform(-0.1, 0.1);
            acc[i] = true;
        }

        char pname[20][255] = { "Origin",
            "Q1 ent", "Q1 ext", "Q2 ent", "Q2 ext", "Dip ent", "Dip ext",
            "Q3 ent", "Q3 ext", "Focal Plane"
        };

        nbin = 40;

        for( i = 0; i < nelm+1; i++ ){
            hacc[0][0][i] = new TH1F(Form("hacc_%d_%d_%d",0, i+1,0), Form("X %s - All", pname[i]), nbin, -1.2, 1.2);
            hacc[0][1][i] = new TH1F(Form("hacc_%d_%d_%d",0, i+1,1), Form("X %s - Accepted", pname[i]), nbin, -1.2, 1.2);

            hacc[1][0][i] = new TH1F(Form("hacc_%d_%d_%d",1, i+1,0), Form("#theta %s - All", pname[i]), nbin, -0.2, 0.2);
            hacc[1][1][i] = new TH1F(Form("hacc_%d_%d_%d",1, i+1,1), Form("#theta X %s - Accepted", pname[i]), nbin, -0.2, 0.2);

            hacc[2][0][i] = new TH1F(Form("hacc_%d_%d_%d",2, i+1,0), Form("Y %s - All", pname[i]), nbin, -0.2, 0.2);
            hacc[2][1][i] = new TH1F(Form("hacc_%d_%d_%d",2, i+1,1), Form("Y %s - Accepted", pname[i]), nbin, -0.2, 0.2);

            hacc[3][0][i] = new TH1F(Form("hacc_%d_%d_%d",3, i+1,0), Form("#phi %s - All", pname[i]), nbin, -0.2, 0.2);
            hacc[3][1][i] = new TH1F(Form("hacc_%d_%d_%d",3, i+1,1), Form("#phi X %s - Accepted", pname[i]), nbin, -0.2, 0.2);

            hacc[4][0][i] = new TH1F(Form("hacc_%d_%d_%d",4, i+1,0), Form("#delta %s - All", pname[i]), nbin, -0.1, 0.1);
            hacc[4][1][i] = new TH1F(Form("hacc_%d_%d_%d",4, i+1,1), Form("#delta X %s - Accepted", pname[i]), nbin, -0.1, 0.1);

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

                // Test acceptance
                switch(i){
                    case 0:
                        if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[0] ){
                            acc[j] = false;
                        }
                    case 1:
                        if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[0] ){
                            acc[j] = false;
                        }
                    case 2:
                        if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[1] ){
                            acc[j] = false;
                        }
                    case 3:
                        if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[1] ){
                            acc[j] = false;
                        }

                    case 4:
                        if( fabs(v[kY]) > 0.25/2  || fabs(v[kX])> 0.4  ){
                            acc[j] = false;
                        }
                    case 5:
                        if( fabs(v[kY]) > 0.25/2  || fabs(v[kX])> 0.4  ){
                            acc[j] = false;
                        }
                    case 6:
                        if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[2] ){
                            acc[j] = false;
                        }
                    case 7:
                        if( sqrt(v[kX]*v[kX] + v[kY]*v[kY]) > apps[2] ){
                            acc[j] = false;
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

void THRSTrans::ShowOutput(){
    int i;

    TMultiGraph *mgx = new TMultiGraph();
    mgx->Add(gx_th, "L");
    mgx->Add(gx_x0, "L");
    mgx->Add(gx_d, "L");


    TMultiGraph *mgy = new TMultiGraph();
    mgy->Add(gy_ph, "L");
    mgy->Add(gy_y0, "L");
    mgy->Add(gph_ph, "L");

    TMultiGraph *mgph = new TMultiGraph();
    mgph->Add(gph_ph, "L");
    mgph->Add(gph_y0, "L");

    TCanvas *cx =new TCanvas();
    cx->SetGridx();
    cx->SetGridy();
    mgx->Draw("A");


    TCanvas *cy = new TCanvas();
    cy->SetGridx();
    cy->SetGridy();
    mgy->Draw("A");

    TCanvas *c = new TCanvas();
    c->Divide(5,2);
    for( i = 0; i < 5; i++ ){
        c->cd(i+1);
        hacc[i][0][0]->Draw();
        hacc[i][1][0]->Draw("same");
        c->cd(i+6);
        gacc[i][0]->Draw("AC");
    }


    PrintSimple(trans[nelm]);

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
        (*m)[kX][kTh]  = sin(kq*l);
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
TMatrixD *THRSTrans::makedip( double theta, double rho, double n, double beta1, double beta2, double psi  ){
    double b1 = beta1*3.14159/180;
    double b2 = beta2*3.14159/180;

    int i,j;

    double h = 1./rho;
    double kx2 = (1.0-n)*h*h;
    double ky2 = n*h*h;


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

    (*men)[kPh][kY] = -h*tan(b1-psi);
    (*men)[kPh][kPh] = 1.0;

    (*men)[kPh][kXY] = 2*h*h*n*tan(b1);
    (*men)[kPh][kXPh] = -h*tan(b1)*tan(b1);
    (*men)[kPh][kThY] = -h/(cos(b1)*cos(b1));

    (*men)[kPh][kYd] = h*tan(b1)-h*psi/(cos(b1-psi)*cos(b1-psi));

    (*men)[kd][kd]= 1.0;


    setcrossterms(men);

    // Ex fringe
    TMatrixD *mex = new TMatrixD(20,20);
    (*mex)[kX][kX] = 1.0;
    (*mex)[kX][kX2] = h*tan(b2)*tan(b2)/2.0;
    (*mex)[kX][kY2] = -h/(cos(b2)*cos(b2))/2.0;
    (*mex)[kTh][kX] = h*tan(b2);
    (*mex)[kTh][kTh] = 1.0;
    (*mex)[kTh][kX2] = -h*h*(n+0.5*tan(b2));

    (*mex)[kTh][kXTh] = -h*tan(b2)*tan(b2);
    (*mex)[kTh][kXd] = -h*tan(b2);
    (*mex)[kTh][kY2] =  h*h*(n-0.5*tan(b2)*tan(b2))*tan(b2);
    (*mex)[kTh][kYPh] = h*tan(b2)*tan(b2);

    (*mex)[kY][kY] = 1.0;
    (*mex)[kY][kXY] = -h*tan(b2)*tan(b2);

    (*mex)[kPh][kY] = -h*tan(b2-psi);
    (*mex)[kPh][kPh] = 1.0;

    (*mex)[kPh][kXY] = h*h*(2*n+1.0/(cos(b2)*cos(b2)))*tan(b2);

    (*mex)[kPh][kXPh] = h*tan(b2)*tan(b2);
    (*mex)[kPh][kThY] = h/(cos(b2)*cos(b2));

    (*mex)[kPh][kYd] = h*tan(b2)-h*psi/(cos(b2-psi)*cos(b2-psi));

    (*mex)[kd][kd] = 1.0;

    setcrossterms(mex);



    //   Main dipole field
    TMatrixD *mdip = new TMatrixD(20,20);

    double cxt = cos(kx*l);
    double sxt = sin(kx*l)/kx;
    double cpxt = -kx2*sxt;
    double spxt = cos(kx*l);
    double dxt = h*(1.0-cxt)/kx2;
    double dpxt = h*sxt;
    double cyt = cosh(ky*l);
    double syt = sinh(ky*l)/ky;
    double cpyt = -ky2*syt;
    double spyt = cyt;

    double I10 = dxt/h;
    double I11 = 0.5*l*dxt;
    double I12 = (dxt-l*cxt)/(2.0*kx2);
    double I16 = h*(I10 - I11)/kx2;

    double I111 = (sxt*sxt + I10)/3.0;
    double I112 = sxt*I10/3.0;
    double I116 = h*(I11 - I111)/kx2;
    double I122 = (2*I10 - sxt*sxt)/(3.0*kx2);
    double I166 = h*h*( 4.0*I10/3.0 + sxt*sxt/3.0 - l*sxt)/kx2/kx2;

    double I20 = sxt;
    double I21 = 0.5*(dxt+l*cxt);
    double I22 = I11;
    double I26 = 0.5*h*(sxt - l*cxt)/kx2;

    double I211 = sxt*(1.0-2.0*cxt)/3.0;
    double I212 = (2.0*sxt*sxt - dxt/h)/3.0;
    double I216 = h*(0.5*l*cxt + sxt/6.0 - 2.0*sxt*cxt/3.0)/kx2;
    double I222 = 2.0*sxt*dxt/(h*3.0);
    double I226 = h*(0.5*l*sxt - 2.0*sxt*sxt/3.0 + dxt/(h*3.0))/kx2;
    double I266 = h*h*(sxt/3.0 + 2.0*sxt*cxt/3.0 - l*cxt)/kx2/kx2;

    double I33 = 0.5*l*syt;
    double I34 = 0.5*(syt - l*cyt)/ky2;

    double I313 = (kx2*cyt*dxt/h - 2.0*ky2*sxt*syt)/(kx2 - 4.0*ky2);
    double I314 = (2.0*sxt*cyt - syt*(1.0+cxt))/(kx2 - 4.0*ky2);
    double I323 = (2.0*ky2*syt*(1.0+cxt)/kx2 - sxt*cyt)/(kx2 - 4.0*ky2) +syt/kx2;
    double I324 = ( 2.0*cyt*dxt/h - sxt*syt )/(kx2 - 4*ky2);
    double I336 = h*(0.5*l*cyt - ( cyt*(1.0-cxt) - 2.0*ky2*sxt*syt)/(kx2 - 4.0*ky2))/kx2;
    double I346 = h*( 0.5*(syt-l*cyt)/ky2 - (2.0*sxt*cyt - syt*(1.0+cxt))/(kx2 - 4.0*ky2) );
    
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
    (*mdip)[kX][kXTh] = 2.0*(2.0*n-1)*h*h*h*I112 - kx2*h*I112;
    (*mdip)[kX][kXd] = (2.0-n)*h*h*I11 + 2.0*(2.0*n-1)*h*h*h*I116 - kx2*h*h*I122;
    (*mdip)[kX][kTh2] =  (2.0*n-1)*h*h*h*I122 + 0.5*h*I111;
    (*mdip)[kX][kThd] =  (2.0-n)*h*h*I12 + 2.0*(2.0*n-1)*h*h*h*I122 + h*h*I112;

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
    (*mdip)[kPh][kThY] = -2.0*n*h*h*h*I423 - ky2*h*I414 - h*dxt*cpyt;
    (*mdip)[kPh][kThPh] = -2.0*n*h*h*h*I424 + h*I413 - h*dxt*spyt;
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



void THRSTrans::PrintSimple(TMatrixD *m){
    int i, j;

    int map[5] = {kX, kTh, kY, kPh, kd };

    TMatrixD *newm = new TMatrixD(5,5);

    for( i = 0; i < 5; i++ ){
        for( j = 0; j < 5; j++ ){
            (*newm)[i][j] = (*m)[map[i]][map[j]];
        }
    }

    newm->Print("f=%6.2f");

    delete newm;
    return;
}

void THRSTrans::setcrossterms(TMatrixD *m){
    (*m)[kX2][kX2] =  (*m)[kX][kX]*(*m)[kX][kX]; 
    (*m)[kX2][kXTh] = 2.0*(*m)[kX][kX]*(*m)[kX][kTh]; 
    (*m)[kX2][kXd] =  2.0*(*m)[kX][kX]*(*m)[kX][kd]; 
    (*m)[kX2][kTh2] = (*m)[kX][kTh]*(*m)[kX][kTh]; 
    (*m)[kX2][kThd] = 2.0*(*m)[kX][kTh]*(*m)[kX][kd]; 
    (*m)[kX2][kd2] =  (*m)[kX][kd]*(*m)[kX][kd]; 


    (*m)[kXTh][kX2] =   (*m)[kX][kX]*(*m)[kTh][kX]; 
    (*m)[kXTh][kXTh] =  (*m)[kX][kX]*(*m)[kTh][kTh] +
                            (*m)[kX][kTh]*(*m)[kTh][kX]; 
    (*m)[kXTh][kXd] =   (*m)[kX][kX]*(*m)[kTh][kd] +
                            (*m)[kX][kd]*(*m)[kTh][kX]; 
    (*m)[kXTh][kTh2] =  (*m)[kX][kTh]*(*m)[kTh][kTh] +
                            (*m)[kX][kTh]*(*m)[kTh][kTh]; 
    (*m)[kXTh][kThd] =  (*m)[kX][kTh]*(*m)[kTh][kd] +
                            (*m)[kX][kd]*(*m)[kTh][kTh]; 
    (*m)[kXTh][kd2] =   (*m)[kX][kd]*(*m)[kTh][kd] +
                            (*m)[kX][kd]*(*m)[kTh][kd]; 

    (*m)[kXd][kX2] =    (*m)[kX][kX]*(*m)[kd][kX]; 
    (*m)[kXd][kXTh] =   (*m)[kX][kX]*(*m)[kd][kTh] +
                            (*m)[kX][kTh]*(*m)[kd][kX]; 
    (*m)[kXd][kXd] =    (*m)[kX][kX]*(*m)[kd][kd] +
                            (*m)[kX][kd]*(*m)[kd][kX]; 
    (*m)[kXd][kTh2] =   (*m)[kX][kTh]*(*m)[kd][kTh] +
                            (*m)[kX][kTh]*(*m)[kd][kTh]; 
    (*m)[kXd][kThd] =   (*m)[kX][kTh]*(*m)[kd][kd] +
                            (*m)[kX][kd]*(*m)[kd][kTh]; 
    (*m)[kXd][kd2] =    (*m)[kX][kd]*(*m)[kd][kd] +
                            (*m)[kX][kd]*(*m)[kd][kd]; 

    (*m)[kTh2][kX2] =    (*m)[kTh][kX]*(*m)[kTh][kX]; 
    (*m)[kTh2][kXTh] =   (*m)[kTh][kX]*(*m)[kTh][kTh] +
                            (*m)[kTh][kTh]*(*m)[kTh][kX]; 
    (*m)[kTh2][kXd] =    (*m)[kTh][kX]*(*m)[kTh][kd] +
                            (*m)[kTh][kd]*(*m)[kTh][kX]; 
    (*m)[kTh2][kTh2] =   (*m)[kTh][kTh]*(*m)[kTh][kTh] +
                            (*m)[kTh][kTh]*(*m)[kTh][kTh]; 
    (*m)[kTh2][kThd] =   (*m)[kTh][kTh]*(*m)[kTh][kd] +
                            (*m)[kTh][kd]*(*m)[kTh][kTh]; 
    (*m)[kTh2][kd2] =    (*m)[kTh][kd]*(*m)[kTh][kd] +
                            (*m)[kTh][kd]*(*m)[kTh][kd]; 

    (*m)[kThd][kX2] =    (*m)[kTh][kX]*(*m)[kd][kX]; 
    (*m)[kThd][kXTh] =   (*m)[kTh][kX]*(*m)[kd][kTh] +
                            (*m)[kTh][kTh]*(*m)[kd][kX]; 
    (*m)[kThd][kXd] =    (*m)[kTh][kX]*(*m)[kd][kd] +
                            (*m)[kTh][kd]*(*m)[kd][kX]; 
    (*m)[kThd][kTh2] =   (*m)[kTh][kTh]*(*m)[kd][kTh] +
                            (*m)[kTh][kTh]*(*m)[kd][kTh]; 
    (*m)[kThd][kThd] =   (*m)[kTh][kTh]*(*m)[kd][kd] +
                            (*m)[kTh][kd]*(*m)[kd][kTh]; 
    (*m)[kThd][kd2] =    (*m)[kTh][kd]*(*m)[kd][kd] +
                            (*m)[kTh][kd]*(*m)[kd][kd]; 


    (*m)[kd2][kX2] =    (*m)[kd][kX]*(*m)[kd][kX]; 
    (*m)[kd2][kXTh] =   (*m)[kd][kX]*(*m)[kd][kTh] +
                            (*m)[kd][kTh]*(*m)[kd][kX]; 
    (*m)[kd2][kXd] =    (*m)[kd][kX]*(*m)[kd][kd] +
                            (*m)[kd][kd]*(*m)[kd][kX]; 
    (*m)[kd2][kTh2] =   (*m)[kd][kTh]*(*m)[kd][kTh] +
                            (*m)[kd][kTh]*(*m)[kd][kTh]; 
    (*m)[kd2][kThd] =   (*m)[kd][kTh]*(*m)[kd][kd] +
                            (*m)[kd][kd]*(*m)[kd][kTh]; 
    (*m)[kd2][kd2] =    (*m)[kd][kd]*(*m)[kd][kd] +
                            (*m)[kd][kd]*(*m)[kd][kd]; 

    (*m)[kY2][kY2] =    (*m)[kY][kY]*(*m)[kY][kY]; 
    (*m)[kY2][kYPh] =   (*m)[kY][kY]*(*m)[kY][kPh] +
                            (*m)[kY][kPh]*(*m)[kY][kY]; 
    (*m)[kY2][kPh2] =   (*m)[kY][kPh]*(*m)[kY][kPh] +
                            (*m)[kY][kPh]*(*m)[kY][kPh]; 

    (*m)[kYPh][kY2] =   (*m)[kY][kY]*(*m)[kPh][kY]; 
    (*m)[kYPh][kYPh] =  (*m)[kY][kY]*(*m)[kPh][kPh] +
                            (*m)[kY][kPh]*(*m)[kPh][kY]; 
    (*m)[kYPh][kPh2] =  (*m)[kY][kPh]*(*m)[kPh][kPh] +
                            (*m)[kY][kPh]*(*m)[kPh][kPh]; 


    (*m)[kPh2][kY2] =   (*m)[kPh][kY]*(*m)[kPh][kY]; 
    (*m)[kPh2][kYPh] =  (*m)[kPh][kY]*(*m)[kPh][kPh] +
                            (*m)[kPh][kPh]*(*m)[kPh][kY]; 
    (*m)[kPh2][kPh2] =  (*m)[kPh][kPh]*(*m)[kPh][kPh] +
                            (*m)[kPh][kPh]*(*m)[kPh][kPh]; 

    (*m)[kXY][kXY] =    (*m)[kX][kX]*(*m)[kY][kY]; 
    (*m)[kXY][kXPh] =   (*m)[kX][kX]*(*m)[kY][kPh]; 
    (*m)[kXY][kThY] =   (*m)[kX][kTh]*(*m)[kY][kY]; 
    (*m)[kXY][kThPh] =  (*m)[kX][kTh]*(*m)[kY][kPh]; 
    (*m)[kXY][kYd] =    (*m)[kX][kd]*(*m)[kY][kY]; 
    (*m)[kXY][kPhd] =   (*m)[kX][kd]*(*m)[kY][kPh]; 

    (*m)[kXPh][kXY] =    (*m)[kX][kX]*(*m)[kPh][kY]; 
    (*m)[kXPh][kXPh] =   (*m)[kX][kX]*(*m)[kPh][kPh]; 
    (*m)[kXPh][kThY] =   (*m)[kX][kTh]*(*m)[kPh][kY]; 
    (*m)[kXPh][kThPh] =  (*m)[kX][kTh]*(*m)[kPh][kPh]; 
    (*m)[kXPh][kYd] =    (*m)[kX][kd]*(*m)[kPh][kY]; 
    (*m)[kXPh][kPhd] =   (*m)[kX][kd]*(*m)[kPh][kPh]; 

    (*m)[kThY][kXY] =    (*m)[kTh][kX]*(*m)[kY][kY]; 
    (*m)[kThY][kXPh] =   (*m)[kTh][kX]*(*m)[kY][kPh]; 
    (*m)[kThY][kThY] =   (*m)[kTh][kTh]*(*m)[kY][kY]; 
    (*m)[kThY][kThPh] =  (*m)[kTh][kTh]*(*m)[kY][kPh]; 
    (*m)[kThY][kYd] =    (*m)[kTh][kd]*(*m)[kY][kY]; 
    (*m)[kThY][kPhd] =   (*m)[kTh][kd]*(*m)[kY][kPh]; 

    (*m)[kThPh][kXY] =    (*m)[kTh][kX]*(*m)[kPh][kY]; 
    (*m)[kThPh][kXPh] =   (*m)[kTh][kX]*(*m)[kPh][kPh]; 
    (*m)[kThPh][kThY] =   (*m)[kTh][kTh]*(*m)[kPh][kY]; 
    (*m)[kThPh][kThPh] =  (*m)[kTh][kTh]*(*m)[kPh][kPh]; 
    (*m)[kThPh][kYd] =    (*m)[kTh][kd]*(*m)[kPh][kY]; 
    (*m)[kThPh][kPhd] =   (*m)[kTh][kd]*(*m)[kPh][kPh]; 

    (*m)[kYd][kXY] =    (*m)[kd][kX]*(*m)[kY][kY]; 
    (*m)[kYd][kXPh] =   (*m)[kd][kX]*(*m)[kY][kPh]; 
    (*m)[kYd][kThY] =   (*m)[kd][kTh]*(*m)[kY][kY]; 
    (*m)[kYd][kThPh] =  (*m)[kd][kTh]*(*m)[kY][kPh]; 
    (*m)[kYd][kYd] =    (*m)[kd][kd]*(*m)[kY][kY]; 
    (*m)[kYd][kPhd] =   (*m)[kd][kd]*(*m)[kY][kPh]; 

    (*m)[kPhd][kXY] =    (*m)[kd][kX]*(*m)[kPh][kY]; 
    (*m)[kPhd][kXPh] =   (*m)[kd][kX]*(*m)[kPh][kPh]; 
    (*m)[kPhd][kThY] =   (*m)[kd][kTh]*(*m)[kPh][kY]; 
    (*m)[kPhd][kThPh] =  (*m)[kd][kTh]*(*m)[kPh][kPh]; 
    (*m)[kPhd][kYd] =    (*m)[kd][kd]*(*m)[kPh][kY]; 
    (*m)[kPhd][kPhd] =   (*m)[kd][kd]*(*m)[kPh][kPh]; 
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

