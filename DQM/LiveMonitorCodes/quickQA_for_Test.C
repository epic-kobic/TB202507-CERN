#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>

#include "CaloMap_Mar.h"
using namespace std;


map<int,int> cutPed = {{41,0},{42,0}};
const float xCVS = 2.0f;
static map<pair<int,int>,vector<int>> chMapCalo = GetCaloChMap();


static vector<unsigned char> buf;
static vector<int>           wf;
void ensureBuffers(int maxLen)
{
    if ((int)buf.capacity() < maxLen*2) buf.reserve(maxLen*2);
    if ((int)wf.capacity()  < maxLen)   wf.reserve(maxLen);
}


TCanvas* CreateCanvas_LR(int RunNo, const char* side, int layer = -1) {
    int nlayers = (layer<0 ? 4 : 1);
    TCanvas* c = new TCanvas(
        Form("c_Run%d_%s", RunNo, side),
        Form("Run %d %s QA, %d floor", RunNo, side, layer),
        100, 100,
        int(120*8*xCVS), int(120*nlayers*xCVS)
    );
    c->Divide(8, nlayers);
    return c;
}


void resetAll(
  TH2F* hO[2][32],
  TProfile* hA[2][32],
  TH1F* hP[2][32],
  TH1F* hI[2][32]
){
  for(int lr=0; lr<2; ++lr){
    for(int pad=0; pad<32; ++pad){
      if (hO[lr][pad]) hO[lr][pad]->Reset();
      if (hA[lr][pad]) hA[lr][pad]->Reset();
      if (hP[lr][pad]) hP[lr][pad]->Reset();
      if (hI[lr][pad]) hI[lr][pad]->Reset();
    }
  }
}


void bic_daq_quickQA(int RunNo, int nEvt, const char* inPath, int mid, const char* plotType, TH2F* hO[2][32], TProfile* hA[2][32], TH1F* hP[2][32], TH1F* hI[2][32]){
    TString fn = Form("%s/Run_%d/Run_%d_MID_%d/bic_daq_%d_%d.dat",
                      inPath, RunNo, RunNo, mid, mid, RunNo);
    FILE* fp = fopen(fn.Data(),"rb");
    if(!fp){ cerr<<"[ERROR] "<<fn<<" open failed\n"; return; }

    cout<<"File: "<<fn<<" open" << endl;

    fseek(fp,0,SEEK_END);
    int fsize = ftell(fp);
    rewind(fp);
    // TH1F* H1_pulse = new TH1F(Form("H1_pulse_%d",mid), "", 240, 0, 240);

    ensureBuffers(240);
    buf.resize(0);
    wf.resize(0);
    const int nCh = 32;
    int readBytes=0, cnt=0, ped=cutPed[mid];
    char hdr[32];

    while(readBytes<fsize){
        // if (fread(hdr,1,32,fp)<32) break;
        fread(hdr, 1, 32, fp);
        int len = 0;
        for(int i=0;i<4;++i) len |= (hdr[i]&0xFF) << (8*i);
        cout << "DLen  " << len << endl;
        // int tcb_trigger_number = 0;
		// for (int a=0; a<4; a++) tcb_trigger_number += ((int)(header[a+7] & 0xFF) << 8*a);
        // if ( fabs(tcb_trigger_number-trigN)>10000 ) continue;

        // if ( trigN!=tcb_trigger_number ){
		// 	if ( trigT!=tcb_trigger_time ){
		// 		trigN = tcb_trigger_number;
		// 		trigT = tcb_trigger_time;
		// 	}else{
		// 		cout << "WARNNING! different trigger number but same trigger time!" << endl;
		// 	}
		// }
        if (len != 512) {cout << "WARN: suspicious data length! DLen=" << len << " at event " << cnt << endl;
            continue;
        }


        int ch = hdr[16]&0xFF;
        auto it = chMapCalo.find({mid,ch});
        if(ch<1||ch>32 || it==chMapCalo.end()){
            fseek(fp,len-32,SEEK_CUR);
            continue;
        }

        int nbytes = len-32;
        buf.resize(nbytes);
        fread(buf.data(),1,nbytes,fp);
        
        int wl = nbytes/2;
        if(wl<=0){ ++cnt; continue; }
        // if (cnt ==1) cout << "WaveLength  " << wl << endl;
        wf.resize(wl);
        bool valid=false;
        // H1_pulse->Reset();
        for(int i=0;i<wl;++i){
            int t_adc1 = (buf[2*i]&0xFF);
			int t_adc2 = (buf[2*i+1]&0xFF)<<8;
            int v = (short)(t_adc1 + t_adc2);
            wf[i]=v;
            // H1_pulse->SetBinContent(i+1, v);
            if (abs(v)>ped) valid=true;
        }
        if(!valid){ ++cnt; continue; }

        auto info = it->second;
        int LR=info[0], col=info[2], row=info[3];
        int pad = row*8 + col;
        int from = 100, to = min(wl,200);

        if (strcmp(plotType, "o") == 0) {
        for(int i = 0; i < wl; ++i)
                hO[LR][pad]->Fill(i, wf[i]);
        }
        else if (strcmp(plotType, "a") == 0) {
            for(int i = 0; i < wl; ++i)
                hA[LR][pad]->Fill(i, wf[i]);
        }
        else if (strcmp(plotType, "p") == 0) {
            if (to > from) {
                auto peak = *max_element(wf.begin() + from, wf.begin() + to);
                hP[LR][pad]->Fill(peak);
            }
        }
        else if (strcmp(plotType, "i") == 0) {
            if (to > from) {
                double sum = accumulate(wf.begin() + from, wf.begin() + to, 0.0);
                hI[LR][pad]->Fill(sum);
            }
        }

        ++cnt;
         readBytes = readBytes + len;
        cout << "Read byte  " << readBytes << endl;
        if (nEvt>0 && cnt >= nEvt) {
            cout << "Run end. " << cnt <<" events processed" << endl;
            break;
        }
    }
    fclose(fp);
}


void quickQA_for_Test(int RunNo=60292, int nEvt=50000, const char* inPath="/Users/jay.ryu/KoBIC/25KEKDATA", const char* plotType="i", int layer = 1)
{
    gStyle->SetOptStat(1111);
    gStyle->SetTitleSize(0.04f);


    TCanvas* cL = CreateCanvas_LR(RunNo,"LEFT", layer);
    TCanvas* cR = CreateCanvas_LR(RunNo,"RIGHT", layer);

    static TH2F*    hO[2][32];
    static TProfile*hA[2][32];
    static TH1F*    hP[2][32];
    static TH1F*    hI[2][32];
    static bool init=false;
    // if(!init){
    //   init=true;
    //   for(int lr=0;lr<2;++lr){
    //     const char* side = lr==0?"L":"R";
    //     for(int pad=0;pad<32;++pad){
    //       hO[lr][pad] = new TH2F(Form("O_%s_%02d",side,pad),Form("Overlay %s pad %02d",side,pad),240,0,240, 2000,-200,12000); 
    //       hO[lr][pad]->SetDirectory(0);

    //       hA[lr][pad] = new TProfile(Form("A_%s_%02d",side,pad),Form("AvgTime %s pad %02d",side,pad),240,0,240, -200,12000); 
    //       hA[lr][pad]->SetDirectory(0);

    //       hP[lr][pad] = new TH1F(Form("P_%s_%02d",side,pad),Form("PeakADC %s pad %02d",side,pad),10000,0,20000); 
    //       hP[lr][pad]->SetDirectory(0);

    //       hI[lr][pad] = new TH1F(Form("I_%s_%02d",side,pad),Form("IntADC %s pad %02d",side,pad),100,-100,150000); 
    //       hI[lr][pad]->SetDirectory(0);
    //     }
    //   }
    // }

     if(!init){
      init = true;
      for(int lr=0; lr<2; ++lr){
        const char* side = lr==0 ? "L" : "R";
        char sideChar = side[0];
        for(int pad=0; pad<32; ++pad){
          int row = pad/8, col = pad%8;

          TString mod41(""), mod42("");
          for(auto &kv: chMapCalo){
            int mid = kv.first.first;
            auto info = kv.second; // {LR, caloCh, row, col}
            if(info[0]==lr && info[2]==col && info[3]==row){
              TString name = Form("%c%d", sideChar, info[1]);
              if(mid==41) mod41 = name;
              if(mid==42) mod42 = name;
            }
          }
          TString mods = "";
          if(mod41.Length()) mods += mod41;
          if(mod42.Length()){
            if(mods.Length()) mods += " ";
            mods += mod42;
          }
          if(mods.IsNull()) mods = "NA";

          TString tO = Form("Overlay %s pad%02d  [%s]", side, pad, mods.Data());
          TString tA = Form("AvgTime  %s pad%02d  [%s]", side, pad, mods.Data());
          TString tP = Form("PeakADC  %s pad%02d  [%s]", side, pad, mods.Data());
          TString tI = Form("IntADC   %s pad%02d  [%s]", side, pad, mods.Data());

          hO[lr][pad] = new TH2F(Form("O_%s_%02d",side,pad), tO, 240,0,240, 2000,-200,12000);
          hO[lr][pad]->SetDirectory(0);

          hA[lr][pad] = new TProfile(Form("A_%s_%02d",side,pad), tA, 240,0,240, -200,12000);
          hA[lr][pad]->SetDirectory(0);

          hP[lr][pad] = new TH1F(Form("P_%s_%02d",side,pad), tP, 10000,0,20000);
          hP[lr][pad]->SetDirectory(0);

          hI[lr][pad] = new TH1F(Form("I_%s_%02d",side,pad), tI, 1000,-100,150000);
          hI[lr][pad]->SetDirectory(0);
        }
      }
    }

    resetAll(hO,hA,hP,hI);

    bic_daq_quickQA(RunNo,nEvt,inPath,41,plotType, hO,hA,hP,hI);
    bic_daq_quickQA(RunNo,nEvt,inPath,42,plotType, hO,hA,hP,hI);

    // LEFT
    if(layer<0){
        for(int i=0;i<32;++i){
          cL->cd(i+1);
          if      (!strcmp(plotType,"o")) hO[0][i]->Draw("COLZ");
          else if (!strcmp(plotType,"a")) hA[0][i]->Draw();
          else if (!strcmp(plotType,"p")) hP[0][i]->Draw("hist");
          else if (!strcmp(plotType,"i")) hI[0][i]->Draw("hist");
        }
    }else{
        for (int col=0; col<8; col++){
            cL->cd(col+1);
            int pad = layer*8 + col;
            if      (!strcmp(plotType,"o")) hO[0][pad]->Draw("COLZ");
            else if (!strcmp(plotType,"a")) hA[0][pad]->Draw();
            else if (!strcmp(plotType,"p")) hP[0][pad]->Draw("hist");
            else if (!strcmp(plotType,"i")) hI[0][pad]->Draw("hist");
        }
    }
    cL->Update();

    // RIGHT
    if(layer<0){
        for(int i=0;i<32;++i){
          cR->cd(i+1);
          if      (!strcmp(plotType,"o")) hO[1][i]->Draw("COLZ");
          else if (!strcmp(plotType,"a")) hA[1][i]->Draw();
          else if (!strcmp(plotType,"p")) hP[1][i]->Draw("hist");
          else if (!strcmp(plotType,"i")) hI[1][i]->Draw("hist");
        }
    }else{
        for (int col=0; col<8; col++){
            cR->cd(col+1);
            int pad = layer*8 + col;
            if      (!strcmp(plotType,"o")) hO[1][pad]->Draw("COLZ");
            else if (!strcmp(plotType,"a")) hA[1][pad]->Draw();
            else if (!strcmp(plotType,"p")) hP[1][pad]->Draw("hist");
            else if (!strcmp(plotType,"i")) hI[1][pad]->Draw("hist");
        }

    }
    cR->Update();

}
