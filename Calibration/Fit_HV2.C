void Fit_HV2(int opt = -1, int extra_run = -1, int extra_hv = -1)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  const int base_nset = 4;
  int run[base_nset] = {60103, 60104, 60105, 60106};  // HV = 550, 600, 650, 700
  int hv[base_nset]  = {550,   600,   650,   700};
  float intADC_val[8] = {7362.5, 32557.8, 48000, 41976, 27920, 16015, 8510, 0};

  if (extra_run == -1 || extra_hv == -1) {
    cout << "[ERROR] You must provide both extra_run and extra_hv!" << endl;
    return;
  }

  const int nset = base_nset + 1;

  vector<double> hv_all(hv, hv + base_nset);
  vector<TFile*> file_all;

  for (int i = 0; i < base_nset; i++) {
    file_all.push_back(new TFile(Form("../output/Run_%d_intADC.root", run[i]), "READ"));
  }

  file_all.push_back(new TFile(Form("../output/Run_%d_intADC.root", extra_run), "READ"));
  hv_all.push_back(extra_hv);

  // Sort by HV
  vector<pair<double, TFile*>> hv_file_pairs;
  for (int i = 0; i < nset; ++i)
    hv_file_pairs.push_back({hv_all[i], file_all[i]});
  sort(hv_file_pairs.begin(), hv_file_pairs.end());

  const int nCh = 8;
  TGraph* gHV_L[nCh];
  TGraph* gHV_R[nCh];

  TCanvas* cL = new TCanvas("cL", "Left Channels", 1200, 800);
  TCanvas* cR = new TCanvas("cR", "Right Channels", 1200, 800);
  cL->Divide(4, 2);
  cR->Divide(4, 2);

  for (int c = 0; c < nCh; ++c) {
    gHV_L[c] = new TGraph();
    gHV_R[c] = new TGraph();

    for (int i = 0; i < nset; ++i) {
      TFile* f = hv_file_pairs[i].second;
      double hv_val = hv_file_pairs[i].first;

      TH1F* hL = (TH1F*)f->Get(Form("h_adc_0_%d_%d", opt, c));
      TH1F* hR = (TH1F*)f->Get(Form("h_adc_1_%d_%d", opt, c));
      if (!hL || !hR) continue;

      gHV_L[c]->SetPoint(i, hv_val, hL->GetMean());
      gHV_R[c]->SetPoint(i, hv_val, hR->GetMean());
    }

    TF1* fitL = new TF1(Form("fitL_%d", c), "pol2", 480, 850);
    TF1* fitR = new TF1(Form("fitR_%d", c), "pol2", 480, 850);

    gHV_L[c]->Fit(fitL, "R0Q");
    gHV_R[c]->Fit(fitR, "R0Q");

    // ---------- Draw Left ----------
    cL->cd(c + 1);
    gHV_L[c]->SetMarkerStyle(20);
    gHV_L[c]->SetTitle(Form("Column %d Left;HV;IntADC Mean", c));
    gHV_L[c]->Draw("AP");
    fitL->SetLineColor(kRed);
    fitL->Draw("same");

    // Desired point marker
    double HV_L = fitL->GetX(intADC_val[c], 490, 850);
    TArrow* arrow_L = new TArrow(HV_L, fitL->Eval(HV_L) - 5000, HV_L, fitL->Eval(HV_L), 0.007, "|>");
    arrow_L->SetLineColor(kRed);
    arrow_L->SetLineWidth(2);
    arrow_L->SetFillColor(kRed + 2);
    arrow_L->Draw();

    TLatex* latexL = new TLatex(HV_L + 10, fitL->Eval(HV_L) - 3000, Form("HV=%.1f", HV_L));
    latexL->SetTextSize(0.05);
    latexL->SetTextColor(kRed + 1);
    latexL->Draw();

    // ---------- Draw Right ----------
    cR->cd(c + 1);
    gHV_R[c]->SetMarkerStyle(21);
    gHV_R[c]->SetMarkerColor(kBlue);
    gHV_R[c]->SetTitle(Form("Column %d Right;HV;IntADC Mean", c));
    gHV_R[c]->Draw("AP");
    fitR->SetLineColor(kBlue);
    fitR->Draw("same");

    double HV_R = fitR->GetX(intADC_val[c], 490, 850);
    TArrow* arrow_R = new TArrow(HV_R, fitR->Eval(HV_R) - 5000, HV_R, fitR->Eval(HV_R), 0.007, "|>");
    arrow_R->SetLineColor(kBlue);
    arrow_R->SetLineWidth(2);
    arrow_R->SetFillColor(kBlue + 2);
    arrow_R->Draw();

    TLatex* latexR = new TLatex(HV_R + 10, fitR->Eval(HV_R) - 3000, Form("HV=%.1f", HV_R));
    latexR->SetTextSize(0.05);
    latexR->SetTextColor(kBlue + 1);
    latexR->Draw();
  }

  for (auto f : file_all) f->Close();
}
