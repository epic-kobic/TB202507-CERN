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

void ensureBuffers(int maxLen) {
    if ((int)buf.capacity() < maxLen*2) buf.reserve(maxLen*2);
    if ((int)wf.capacity()  < maxLen)   wf.reserve(maxLen);
}

TCanvas* CreateCanvas_LR(int RunNo, const char* side, int layer = -1) {
    int nlayers = (layer<0 ? 3 : 1);
    TCanvas* c = new TCanvas(
        Form("c_Run%d_%s", RunNo, side),
        Form("Run %d %s QA, %d floor", RunNo, side, layer), 
        100, 100,
        int(120*8*xCVS), int(120*nlayers*xCVS)); 
    c->Divide(8, nlayers);
    return c;
}

void resetAll(
  TH2F* hO[2][24],
  TProfile* hA[2][24],
  TH1F* hP[2][24],
  TH1F* hI[2][24]
){
  for(int lr=0; lr<2; ++lr){
    for(int pad=0; pad<24; ++pad){
      if (hO[lr][pad]) hO[lr][pad]->Reset();
      if (hA[lr][pad]) hA[lr][pad]->Reset();
      if (hP[lr][pad]) hP[lr][pad]->Reset();
      if (hI[lr][pad]) hI[lr][pad]->Reset();
    }
  }
}

void bic_daq_quickQA(int RunNo, int nEvt, const char* inPath, int mid, const char* plotType,
    TH2F* hO[2][24], TProfile* hA[2][24], TH1F* hP[2][24], TH1F* hI[2][24])
{
    TString fn = Form("%s/Run_%d/Run_%d_MID_%d/bic_daq_%d_%d.dat",
                      inPath, RunNo, RunNo, mid, mid, RunNo);
    FILE* fp = fopen(fn.Data(),"rb");
    if(!fp){ cerr<<"[ERROR] "<<fn<<" open failed\n"; return; }

    fseek(fp,0,SEEK_END);
    int fsize = ftell(fp);
    rewind(fp);

    ensureBuffers(240);
    buf.clear(); wf.clear();

    const int nCh = 32;
    int readBytes=0, cnt=0, ped=cutPed[mid];
    char hdr[32];

    while(readBytes < fsize){
        if(fread(hdr,1,32,fp) < 32) break;

        int len = 0;
        for(int i=0;i<4;++i) len |= (hdr[i]&0xFF) << (8*i);
        if (len != 512) {
            cout<<"WARN: suspicious data length! DLen="<<len<<" at evt="<<cnt<<"\n";
            continue;
        }

        int ch = hdr[16]&0xFF;
        auto it = chMapCalo.find({mid,ch});
        if(ch<1||ch>nCh || it==chMapCalo.end()){
            fseek(fp,len-32,SEEK_CUR);
            readBytes += len;
            ++cnt;
            continue;
        }

        int nbytes = len-32;
        buf.resize(nbytes);
        fread(buf.data(),1,nbytes,fp);

        int wl = nbytes/2;
        wf.resize(wl);
        bool valid=false;
        for(int i=0;i<wl;++i){
            int lo = buf[2*i]&0xFF;
            int hi = (buf[2*i+1]&0xFF)<<8;
            wf[i] = (short)(lo|hi);
            if(abs(wf[i]) > ped) valid=true;
        }
        if(!valid){ readBytes+=len; ++cnt; continue; }

        auto info = it->second;
        int LR=info[0], col=info[2], row=info[3];
        int pad = row*8 + col;
        int from = 100, to = min(wl,200);
        if (pad < 0 || pad >= 24) continue;
        if (strcmp(plotType,"o")==0) {
            for(int i=0;i<wl;++i) hO[LR][pad]->Fill(i,wf[i]);
        }
        else if (strcmp(plotType,"a")==0) {
            for(int i=0;i<wl;++i) hA[LR][pad]->Fill(i,wf[i]);
        }
        else if (strcmp(plotType,"p")==0) {
            if(to>from){
                int peak = *max_element(wf.begin()+from, wf.begin()+to);
                hP[LR][pad]->Fill(peak);
            }
        }
        else if (strcmp(plotType,"i")==0) {
            if(to>from){
                double sum = accumulate(wf.begin()+from, wf.begin()+to, 0.0);
                hI[LR][pad]->Fill(sum);
            }
        }

        readBytes += len;
        ++cnt;
        if(nEvt>0 && cnt>=nEvt) break;
    }
    fclose(fp);
}

void quickQA_for_July(int RunNo=60265, int nEvt=50000,
                      const char* inPath="/Users/jay.ryu/KoBIC/25KEKDATA",
                      const char* plotType="o", int layer = 1)
{
    gStyle->SetOptStat(1111);
    gStyle->SetTitleSize(0.04f);

    if (layer > 3){
        cout << "Error: Wrong layer number !!!" << endl;
        return;
    }

    TCanvas* cL = CreateCanvas_LR(RunNo,"LEFT", layer);
    TCanvas* cR = CreateCanvas_LR(RunNo,"RIGHT", layer);

    static TH2F*     hO[2][24];
    static TProfile* hA[2][24];
    static TH1F*     hP[2][24];
    static TH1F*     hI[2][24];
    static bool init=false;

    if(!init){
      init = true;
      for(int lr=0; lr<2; ++lr){
        const char* side = lr==0 ? "L" : "R";
        char sideChar = side[0];
        for(int pad=0; pad<24; ++pad){
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

    bic_daq_quickQA(RunNo,nEvt,inPath,41,plotType,hO,hA,hP,hI);
    bic_daq_quickQA(RunNo,nEvt,inPath,42,plotType,hO,hA,hP,hI);

    // LEFT
    if(layer<0){
        for(int i=0;i<24;++i){
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
        for(int i=0;i<24;++i){
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
