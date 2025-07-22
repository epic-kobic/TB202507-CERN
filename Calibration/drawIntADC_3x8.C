void drawIntADC_3x8(int RunNo = 60102, int N = 1000, TString mode = "draw", int rowN = -1, int colN = -1 ) {
    const int cutPed = 0;

    // Mapping {(MID, ch)} → (left/right, col, row)
    std::map<std::pair<int,int>, std::tuple<int,int,int>> channelMap;

    std::vector<std::tuple<int,int,int>> mid41_map = {
        {0,7,0}, {1,7,0}, {0,7,1}, {1,7,1}, {1,4,0}, {1,4,1}, {1,4,2}, {1,4,3},
        {1,5,0}, {1,5,1}, {1,5,2}, {1,5,3}, {1,6,0}, {1,6,1}, {1,6,2}, {1,6,3},
        {0,7,2}, {1,7,2}, {0,7,3}, {1,7,3}, {0,4,0}, {0,4,1}, {0,4,2}, {0,4,3},
        {0,5,0}, {0,5,1}, {0,5,2}, {0,5,3}, {0,6,0}, {0,6,1}, {0,6,2}, {0,6,3}
    };
    for (int i = 0; i < 32; ++i) channelMap[{41, i+1}] = mid41_map[i];

    std::vector<std::tuple<int,int,int>> mid42_map = {
        {1,0,0}, {1,0,1}, {1,0,2}, {1,0,3}, {1,1,0}, {1,1,1}, {1,1,2}, {1,1,3},
        {1,2,0}, {1,2,1}, {1,2,2}, {1,2,3}, {1,3,0}, {1,3,1}, {1,3,2}, {1,3,3},
        {0,0,0}, {0,0,1}, {0,0,2}, {0,0,3}, {0,1,0}, {0,1,1}, {0,1,2}, {0,1,3},
        {0,2,0}, {0,2,1}, {0,2,2}, {0,2,3}, {0,3,0}, {0,3,1}, {0,3,2}, {0,3,3}
    };
    for (int i = 0; i < 32; ++i) channelMap[{42, i+1}] = mid42_map[i];

    // Histogram 선언
    TH1F* hist_L[4][8] = {};
    TH1F* hist_R[4][8] = {};

    // ROOT 파일 열기
    TFile* f = TFile::Open(Form("./input/Run_%d_waveform.root", RunNo));
//TFile* f = TFile::Open(Form("test.root", RunNo));
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

std::map<std::tuple<int,int,int>, double> eventSumMap;

    for (int evt = 0; evt < totalEntries; ++evt) {
        if (validCount >= N) break;

        t->GetEntry(evt);  eventSumMap.clear();
        bool eventValid = false;
        int nCh = MID->size();
        for (size_t i = 0; i < nCh; ++i) {
            int mid_val = MID->at(i);
            int ch_val = ch->at(i);
            auto key = std::make_pair(mid_val, ch_val);
            if (mid_val != 41 && mid_val != 42)        continue;
	    if (channelMap.find(key) == channelMap.end()) continue;

            auto [lr, col, row] = channelMap[key];

            // 옵션 처리
//            if (mode == "draw" && row != rowN) continue;

            int nCh = MID->size();
            int idx = waveform_idx->at(i);
            int len = (i + 1 < nCh ? waveform_idx->at(i + 1) : waveform_total->size());
	    int bin_start = 0; int bin_end = 10000;
            int sum = 0;
            for (int j = idx; j < len; ++j) {
            if (j >= bin_start && j < bin_end)
	        sum += waveform_total->at(j);
            }
            auto key2 = std::make_tuple(lr, row, col);
    	    eventSumMap[key2] += sum;
	}
        for (auto &p : eventSumMap) {
 	auto [lr, row, col] = p.first;
	double val = p.second;
	TH1F*& hist = (lr == 0) ? hist_L[row][col] : hist_R[row][col];
            if (!hist) {
                TString title = Form("L/R %d (%d,%d)", lr, row, col);
                hist = new TH1F(Form("h_adc_%d_%d_%d", lr, row, col), title, 1000, 0, 200000);
            }
            hist->Fill(val);
        }  
        validCount++ ;
    }
    // Canvas
    if (mode == "draw" && rowN == -1){  
	 TCanvas* cL = new TCanvas("cL", "Left (0)", 1600, 800);
	 TCanvas* cR = new TCanvas("cR", "Right (1)", 1600, 800);
         cL->Divide(8, 4);    cR->Divide(8, 4);    
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 8; ++c) {
                int pad = r * 8 + c + 1;
                if (hist_L[3-r][c]) {
                   cL->cd(pad); hist_L[3-r][c]->Draw(); cL->Update();
		   auto stats_L = (TPaveStats*)hist_L[3-r][c]->FindObject("stats");
	           stats_L->SetX1NDC(0.4); stats_L->SetX2NDC(0.94); stats_L->SetY1NDC(0.6);stats_L->SetY2NDC(0.93);
		   stats_L->SetTextSize(0.06);
		   stats_L->Draw("same");
                }
                if (hist_R[3-r][c]) {
                   cR->cd(pad); hist_R[3-r][c]->Draw(); cR->Update();
                   auto stats_R = (TPaveStats*)hist_R[3-r][c]->FindObject("stats");
                   stats_R->SetX1NDC(0.4); stats_R->SetX2NDC(0.94); stats_R->SetY1NDC(0.6);stats_R->SetY2NDC(0.93);
                   stats_R->SetTextSize(0.06);
                   stats_R->Draw("same");
		}
            }
        }
    }
    
    if (mode == "draw" && rowN != -1 && colN == -1){	
	TCanvas* cL = new TCanvas("cL", "Left (0)", 1200, 600);
	TCanvas* cR = new TCanvas("cR", "Right (1)", 1200, 600);
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
 //       TCanvas* cL = new TCanvas("cL", "Left (0)", 1000, 400);
        TCanvas* cR = new TCanvas("cR", "Left/Right,(0)/(1)", 1000, 400);
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
    TFile *outputFile = new TFile(Form("./output/Run_%d_intADC.root",RunNo), "RECREATE");
        for (int r = 0; r < 4; ++r) {
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
