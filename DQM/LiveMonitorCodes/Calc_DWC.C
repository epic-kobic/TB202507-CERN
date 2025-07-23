#include <TGraph.h>
#include <TF1.h>

struct DWCPoint {
    double x;   // mm (ex. -30, 0, 30)
    double y;   // mm
    double deltaLR; // ns
    double deltaUD; // ns
};

void fitDWC(const std::vector<DWCPoint>& points, 
            double& horizontalSlope, double& horizontalOffset, 
            double& verticalSlope, double& verticalOffset,
            TString tag="DWC1")
{
    std::vector<double> posX, posY, dLR, dUD;
    for(const auto& p : points) {
        posX.push_back(p.x);
        posY.push_back(p.y);
        dLR.push_back(p.deltaLR);
        dUD.push_back(p.deltaUD);
    }

    // posX vs deltaLR
    TGraph* gLR = new TGraph(posX.size(), posX.data(), dLR.data());
    gLR->SetTitle(Form("%s: deltaLR vs posX", tag.Data()));
    TF1* fLR = new TF1("fLR", "pol1");
    gLR->Fit(fLR, "Q"); // quiet
    horizontalSlope  = fLR->GetParameter(1);
    horizontalOffset = fLR->GetParameter(0);

    // posY vs deltaUD
    TGraph* gUD = new TGraph(posY.size(), posY.data(), dUD.data());
    gUD->SetTitle(Form("%s: deltaUD vs posY", tag.Data()));
    TF1* fUD = new TF1("fUD", "pol1");
    gUD->Fit(fUD, "Q"); // quiet
    verticalSlope  = fUD->GetParameter(1);
    verticalOffset = fUD->GetParameter(0);

    // gLR->Draw("AP");
    // gUD->Draw("AP");

    delete fLR;
    delete gLR;
    delete fUD;
    delete gUD;
}

void Calc_DWC()
{
    std::vector<int> runList   = {, , }; //runnumber
    std::vector<double> posX   = {-30, 0, 30}; 
    std::vector<double> posY   = { 30, 0,-30}; 

    std::vector<DWCPoint> dwc1points, dwc2points;
    for (int i = 0; i < runList.size(); ++i) {
        TString fname = Form("path/DWC_calib_Run_%d.root", runList[i]);
        TFile* f = TFile::Open(fname, "READ");
        if (!f || f->IsZombie()) { std::cerr << fname << " open fail\n"; continue; }

        // 평균값 추출 (mean, peak 등 원하는대로)
        double mean_DWC1_LR = ((TH1F*)f->Get("h_deltaDWC1_LR"))->GetMean();
        double mean_DWC1_UD = ((TH1F*)f->Get("h_deltaDWC1_UD"))->GetMean();
        double mean_DWC2_LR = ((TH1F*)f->Get("h_deltaDWC2_LR"))->GetMean();
        double mean_DWC2_UD = ((TH1F*)f->Get("h_deltaDWC2_UD"))->GetMean();

        dwc1points.push_back({posX[i], posY[i], mean_DWC1_LR, mean_DWC1_UD});
        dwc2points.push_back({posX[i], posY[i], mean_DWC2_LR, mean_DWC2_UD});
        f->Close();
    }

    double dwc1horizontalSlope, dwc1horizontalOffset;
    double dwc1VerticalSlope,   dwc1VerticalOffset;
    double dwc2horizontalSlope, dwc2horizontalOffset;
    double dwc2VerticalSlope,   dwc2VerticalOffset;

    fitDWC(dwc1points, dwc1horizontalSlope, dwc1horizontalOffset,
                       dwc1VerticalSlope,   dwc1VerticalOffset, "DWC1");
    fitDWC(dwc2points, dwc2horizontalSlope, dwc2horizontalOffset,
                       dwc2VerticalSlope,   dwc2VerticalOffset, "DWC2");

    printf("DWC1 Horizontal: slope = %.5f, offset = %.5f\n", dwc1horizontalSlope, dwc1horizontalOffset);
    printf("DWC1 Vertical:   slope = %.5f, offset = %.5f\n", dwc1VerticalSlope, dwc1VerticalOffset);
    printf("DWC2 Horizontal: slope = %.5f, offset = %.5f\n", dwc2horizontalSlope, dwc2horizontalOffset);
    printf("DWC2 Vertical:   slope = %.5f, offset = %.5f\n", dwc2VerticalSlope, dwc2VerticalOffset);

}
