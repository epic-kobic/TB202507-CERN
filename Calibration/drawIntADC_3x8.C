#include "caloMap.h"
#include <map>
#include <tuple>
#include <utility>
#include <iostream>

void drawIntADC_3x8(int RunNo = 60102, int N = 1000, int rowN = -1, int colN = -1, TString mode = "draw" ) {

    int BinNum = 1000; int BinMax = 300000;


    // Mapping {(MID, ch)} → (left/right, col, row)
    std::map<std::pair<int,int>, std::tuple<int,int,int>> channelMap;


    auto fullMap = GetCaloChMap();
    
    for (const auto& [key, value] : fullMap) {
        int mid = key.first;
        int ch  = key.second;
        int lr  = value[0]; // left/right
        int col = value[2];
        int row = value[3];  

        channelMap[{mid, ch}] = std::make_tuple(lr, col, row);
    }


    // Histogram
    TH1F* hist_L[3][8] = {};
    TH1F* hist_R[3][8] = {};

    // Open ROOT file
    TFile* f = TFile::Open(Form("../input/Run_%d_waveform.root", RunNo));
    TTree* t = (TTree*)f->Get("event_build");

    std::vector<int>* MID = 0;
    std::vector<int>* ch = 0;
    std::vector<short>* waveform_total = 0;
    std::vector<int>* waveform_idx = 0;
    std::vector<int>* data_length = 0;

    t->SetBranchAddress("MID", &MID);
    t->SetBranchAddress("ch", &ch);
    t->SetBranchAddress("waveform_total", &waveform_total);
    t->SetBranchAddress("waveform_idx", &waveform_idx);
    t->SetBranchAddress("data_length", &data_length);

    int validCount = 0;
    int totalEntries = t->GetEntries();

std::map<std::pair<int,int>, double> eventSumMap;
    for (int evt = 0; evt < totalEntries; ++evt) {
        if (validCount >= N) break;

        t->GetEntry(evt);  eventSumMap.clear();
        bool eventValid = false;
        int nCh = MID->size();
        for (size_t i = 0; i < nCh; ++i) {
            if (data_length->size() != 92) continue;
    	    int mid_val = MID->at(i);
            int ch_val = ch->at(i);
            auto key = std::make_pair(mid_val, ch_val);
            if (mid_val != 41 && mid_val != 42)        continue;
	    if (channelMap.find(key) == channelMap.end()) continue;

            auto [lr, col, row] = channelMap[key];
            int nCh = MID->size();
            int idx = waveform_idx->at(i);
            int len = (i + 1 < nCh ? waveform_idx->at(i + 1) : waveform_total->size());
            int sum = 0;
            for (int j = idx; j < len; ++j) {
//  	      if(j%2==0)
	        sum += waveform_total->at(j);
            }
            auto key2 = std::make_pair(mid_val, ch_val);
    	    eventSumMap[key2] += sum;
	}

	for (auto &p : eventSumMap) {
            auto [mid, ch] = p.first;
            double val = p.second;
            auto [lr, col, row] = channelMap[{mid, ch}];
	    TH1F*& hist = (lr == 0) ? hist_L[row][col] : hist_R[row][col];
        if (!hist) {
            TString title = Form("(%d,%d)->(%d,%d,%d)", mid, ch, lr, row, col);
            hist = new TH1F(Form("h_adc_%d_%d_%d", lr, row, col), title, BinNum, 0, BinMax);
            }
            hist->Fill(val);
        }  
        validCount++ ;
    }

    // Canvas
    if (mode == "draw" && rowN == -1){  
	 TCanvas* cL = new TCanvas("cL", "(MID,ch)->(lr,row,col)", 2000, 800);
	 TCanvas* cR = new TCanvas("cR", "(MID,ch)->(lr,row,col)", 2000, 800);
         cL->Divide(8, 3);    cR->Divide(8, 3);    
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 8; ++c) {
                int pad = r * 8 + c + 1;
                if (hist_L[2-r][c]) {
                   cL->cd(pad); hist_L[2-r][c]->Draw(); cL->Update();
		   auto stats_L = (TPaveStats*)hist_L[2-r][c]->FindObject("stats");
	           stats_L->SetX1NDC(0.4); stats_L->SetX2NDC(0.94); stats_L->SetY1NDC(0.6);stats_L->SetY2NDC(0.93);
		   stats_L->SetTextSize(0.06);
		   stats_L->Draw("same");
                }
                if (hist_R[2-r][c]) {
                   cR->cd(pad); hist_R[2-r][c]->Draw(); cR->Update();
                   auto stats_R = (TPaveStats*)hist_R[2-r][c]->FindObject("stats");
                   stats_R->SetX1NDC(0.4); stats_R->SetX2NDC(0.94); stats_R->SetY1NDC(0.6);stats_R->SetY2NDC(0.93);
                   stats_R->SetTextSize(0.06);
                   stats_R->Draw("same");
		}
            }
        }
    }
    
    if (mode == "draw" && rowN != -1 && colN == -1){	
	TCanvas* cL = new TCanvas("cL", "(MID,ch)->(lr,row,col)", 1200, 600);
	TCanvas* cR = new TCanvas("cR", "(MID,ch)->(lr,row,col)", 1200, 600);
    	cL->Divide(4, 2);    cR->Divide(4, 2);
        for (int c = 0; c < 8; ++c) {
            int pad = c + 1;
            if (hist_L[rowN][c]) {
               cL->cd(pad); hist_L[rowN][c]->Draw(); cL->Update();
	       auto stats_L = (TPaveStats*)hist_L[rowN][c]->FindObject("stats");
	       stats_L->SetX1NDC(0.5);stats_L->SetX2NDC(0.94); stats_L->SetY1NDC(0.6);stats_L->SetY2NDC(0.93);
	       stats_L->SetTextSize(0.045);
	       stats_L->Draw("same");
	    }
	    if (hist_R[rowN][c]) {
               cR->cd(pad); hist_R[rowN][c]->Draw(); cR->Update();
	       auto stats_R = (TPaveStats*)hist_R[rowN][c]->FindObject("stats");
	       stats_R->SetX1NDC(0.5);stats_R->SetX2NDC(0.94); stats_R->SetY1NDC(0.6);stats_R->SetY2NDC(0.93);
 	       stats_R->SetTextSize(0.045);
   	       stats_R->Draw("same");
	    }
        }    
    }

    if (mode == "draw" && rowN != -1 && colN != -1 ){ 
        TCanvas* cR = new TCanvas("cR", "(MID,ch)->(lr,row,col)", 1000, 400);
        cR->Divide(2,1);
        if (hist_L[rowN][colN]) {cR->cd(1);
           hist_L[rowN][colN]->Draw(); cR->Update();
	   auto stats_L = (TPaveStats*)hist_L[rowN][colN]->FindObject("stats");
	   stats_L->SetX1NDC(0.6);stats_L->SetX2NDC(0.93); stats_L->SetY1NDC(0.6);stats_L->SetY2NDC(0.93);
	   stats_L->SetTextSize(0.04);
	   stats_L->Draw("same");

        }
        if (hist_R[rowN][colN]) {cR->cd(2);
           hist_R[rowN][colN]->Draw(); cR->Update();
  	   auto stats_h1 = (TPaveStats*)hist_R[rowN][colN]->FindObject("stats");
           stats_h1->SetX1NDC(0.6);stats_h1->SetX2NDC(0.93); stats_h1->SetY1NDC(0.6);stats_h1->SetY2NDC(0.93);
	   stats_h1->SetTextSize(0.04);
	   stats_h1->Draw("same");
        }
    }

    if (mode == "save"){
    TFile *outputFile = new TFile(Form("../output/Run_%d_intADC.root",RunNo), "RECREATE");
        for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 8; ++c) {
            if (hist_L[r][c]) {
               hist_L[r][c]->Write();
                }
            if (hist_R[r][c]) {  
	       hist_R[r][c]->Write();
             }
        }
        }
	outputFile->Close(); 
    }	
}
