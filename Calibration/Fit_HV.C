void Fit_HV(int opt = -1 ) {
    const int nset = 9;
    int run[nset] = {60102, 60103, 60104, 60105, 60106, 60107, 60108, 60109, 60110}; // Run number
    int hv[nset]  = {500, 550, 600, 650, 700, 750, 800, 850, 900};                   // HV of Run
    float intADC_val[8] = {7362.5, 32557.8, 48000, 41976, 27920, 16015, 8510, 0};    // Desired IntADC value

    TGraph* gHV_L[4][8] = {};
    TGraph* gHV_R[4][8] = {};

    // TGraph Initialize
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 8; ++c) {
            gHV_L[r][c] = new TGraph();
            gHV_R[r][c] = new TGraph();
        }
    }

    // Read TFile, Get Mean
    for (int i = 0; i < nset-1; ++i) {
        TFile* infile = new TFile(Form("./output/Run_%d_intADC.root", run[i]), "READ");
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 8; ++c) {

//	    if ( r != opt) continue;
    
	        TH1F* hL = (TH1F*)infile->Get(Form("h_adc_0_%d_%d", r, c));
                TH1F* hR = (TH1F*)infile->Get(Form("h_adc_1_%d_%d", r, c));

                if (hL) gHV_L[r][c]->SetPoint(i, hv[i], hL->GetMean());
                if (hR) gHV_R[r][c]->SetPoint(i, hv[i], hR->GetMean());
            }
        }
        infile->Close();
    }

if (opt == -1)    {
    // 3x8 Canvas
    TCanvas* cL = new TCanvas("cL", "Left", 1600, 800);
    TCanvas* cR = new TCanvas("cR", "Right", 1600, 800);
    cL->Divide(8, 4); cR->Divide(8, 4);
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 8; ++c) {
            TF1* fitL = new TF1("fitL", "exp([0]*x + [1])", 450, 850);
            TF1* fitR = new TF1("fitR", "exp([0]*x + [1])", 450, 850);

            cL->cd(r * 8 + c + 1);
            gHV_L[3-r][c]->Fit(fitL, "R0Q");
            gHV_L[3-r][c]->SetMarkerStyle(20);
            gHV_L[3-r][c]->SetTitle(Form("L (%d,%d);HV;intADC", 3-r, c));
            gHV_L[3-r][c]->Draw("AP");
            fitL->Draw("same");
/*
            double HV_L = fitL->GetX(intADC_val[c], 450, 850);
            TLatex* latexL = new TLatex(500, intADC_val[c]*0.8, Form("HV=%.1f", HV_L));
            latexL->SetTextSize(0.04);
            latexL->Draw();
*/
            cR->cd(r * 8 + c + 1);
            gHV_R[3-r][c]->Fit(fitR, "R0Q");
            gHV_R[3-r][c]->SetMarkerStyle(20);
            gHV_R[3-r][c]->SetTitle(Form("R (%d,%d);HV;intADC", 3-r, c));
            gHV_R[3-r][c]->Draw("AP");
            fitR->Draw("same");
/*
            double HV_R = fitR->GetX(intADC_val[c], 450, 850);
            TLatex* latexR = new TLatex(500, intADC_val[c]*0.8, Form("HV=%.1f", HV_R));
            latexR->SetTextSize(0.04);
            latexR->Draw();
*/
        }
    }
}    
if (opt != -1) {
    TCanvas* cL = new TCanvas("cL", "Left", 1200, 600);
    TCanvas* cR = new TCanvas("cR", "Right", 1200, 600);
    cL->Divide(4, 2); cR->Divide(4, 2);
        for (int c = 0; c < 7; ++c) {
            TF1* fitL = new TF1("fitL", "exp([0]*x + [1])", 450, 850);
            TF1* fitR = new TF1("fitR", "exp([0]*x + [1])", 450, 850);

            cL->cd(c + 1);
            gHV_L[opt][c]->Fit(fitL, "R0Q");
            gHV_L[opt][c]->SetMarkerStyle(20);
            gHV_L[opt][c]->SetTitle(Form("L (%d,%d);HV;intADC", opt, c));
            gHV_L[opt][c]->Draw("AP");
            fitL->Draw("same");

            double HV_L = fitL->GetX(intADC_val[c], 450, 850);
	    TArrow* arrow_L = new TArrow(HV_L, fitL->Eval(HV_L)-5000, HV_L, fitL->Eval(HV_L), 0.007, "|>");
	    arrow_L->SetLineColor(kBlue);
	    arrow_L->SetLineWidth(2);
	    arrow_L->SetFillColor(kBlue+2);         
	    arrow_L->Draw();
            TLatex* latexL = new TLatex(500, 20000+intADC_val[c], Form("HV=%.1f", HV_L));
            latexL->SetTextSize(0.06);
            latexL->Draw();

            cR->cd(c + 1);
            gHV_R[opt][c]->Fit(fitR, "R0Q");
            gHV_R[opt][c]->SetMarkerStyle(20);
            gHV_R[opt][c]->SetTitle(Form("R (%d,%d);HV;intADC", opt, c));
            gHV_R[opt][c]->Draw("AP");
            fitR->Draw("same");

            double HV_R = fitR->GetX(intADC_val[c], 450, 850);
            TArrow* arrow_R = new TArrow(HV_R, fitR->Eval(HV_R)-5000, HV_R, fitR->Eval(HV_R), 0.007, "|>");
            arrow_R->SetLineColor(kBlue);
            arrow_R->SetLineWidth(2);
            arrow_R->SetFillColor(kBlue+2);
            arrow_R->Draw();
	    TLatex* latexR = new TLatex(500,20000+intADC_val[c], Form("HV=%.1f", HV_R));
            latexR->SetTextSize(0.06);
            latexR->Draw();
        }
}	    
}
