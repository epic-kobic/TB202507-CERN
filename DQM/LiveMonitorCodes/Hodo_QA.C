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
using namespace std;

map<int,int> GetHodoChMap(void) {
	std::map<int,int> chMap;
	// 1–16
	chMap.insert({ 1,  1}); chMap.insert({ 2,  2}); chMap.insert({ 3,  3}); chMap.insert({ 4,  4});
	chMap.insert({ 5,  8}); chMap.insert({ 6,  7}); chMap.insert({ 7,  6}); chMap.insert({ 8,  5});
	chMap.insert({ 9, 12}); chMap.insert({10, 11}); chMap.insert({11, 10}); chMap.insert({12,  9});
	chMap.insert({13, 13}); chMap.insert({14, 14}); chMap.insert({15, 15}); chMap.insert({16, 16});
	// 17–32
	for(int i=1;i<=16;++i) chMap.insert({i+16, chMap[i]+16});
	return chMap;
}

map<int,int> cutPed = {{0,0}}; 
static map<int,int> hodoMap = GetHodoChMap();
static vector<unsigned char> buf;
static vector<int> wf;

void ensureBuffers(int maxLen){
	if((int)buf.capacity()<maxLen*2) buf.reserve(maxLen*2);
	if((int)wf.capacity()<maxLen)   wf.reserve(maxLen);
}

TCanvas* CreateCanvas_Hodo(int RunNo, const char* side) {
	TCanvas* c = new TCanvas(
			Form("c_Hodo_Run%d_%s",RunNo,side),
			Form("Run %d Hodo %s QA",RunNo,side),
			100, 120, 1140, 780
			);
	c->Divide(4,4);
	return c;
}

void resetAll(TH2F* hO[2][16], TProfile* hA[2][16], TH1F* hP[2][16], TH1F* hI[2][16]){
	for(int lr=0;lr<2;++lr){
		for(int pad=0;pad<16;++pad){
			if(hO[lr][pad]) hO[lr][pad]->Reset();
			if(hA[lr][pad]) hA[lr][pad]->Reset();
			if(hP[lr][pad]) hP[lr][pad]->Reset();
			if(hI[lr][pad]) hI[lr][pad]->Reset();
		}
	}
}

void hodo_daq_quickQA(int RunNo, int nEvt, const char* inPath, char plotType,
		TH2F* hO[2][16], TProfile* hA[2][16], TH1F* hP[2][16], TH1F* hI[2][16])
{
	const int imid = 31;
	TString fn = Form("%s/Run_%d/Run_%d_MID_%d/jbnu_daq_%d_%d.dat", inPath, RunNo, RunNo, imid, imid, RunNo);
	FILE* fp = fopen(fn.Data(),"rb");
	if(!fp){ cerr<<"[ERROR] "<<fn<<" open failed\n"; return; }

	fseek(fp,0,SEEK_END);
	int fsize = ftell(fp);
	rewind(fp);

	ensureBuffers(10000);
	buf.clear(); wf.clear();

	static char dummy[65536];
	int readBytes=0, cnt=0;
	char hdr[32];

	while(readBytes < fsize){
		if(fread(hdr,1,32,fp) < 32) break;

		int data_length = 0;
		for(int a=0;a<4;++a) data_length |= (unsigned int)(hdr[a]&0xFF) << (8*a);

		if(data_length != 256){
			fread(dummy,1,256-32,fp);
			readBytes += 256;
			++cnt;
			continue;
		}

		int nbytes = data_length - 32;
		buf.resize(nbytes);
		if(fread(buf.data(),1,nbytes,fp) != (size_t)nbytes){
			cerr<<"[FATAL] fread failed at evt "<<cnt<<"\n"; break;
		}

		int wl = nbytes / 2;
		wf.resize(wl);
		for(int i=0;i<wl;++i){
			int lo = buf[2*i]&0xFF;
			int hi = (buf[2*i+1]&0xFF)<<8;
			wf[i] = (short)(lo+hi);
		}

		int channel = hdr[16]&0xFF;
		auto it = hodoMap.find(channel);
		if(it == hodoMap.end()){
			++cnt; readBytes += data_length;
			continue;
		}
		int mapped = it->second;
		int LR  = (mapped<=16?0:1);
		int pad = (mapped<=16? mapped-1 : mapped-17);
		if(pad<0 || pad>=16){ ++cnt; readBytes += data_length; continue; }

		int from = 0, to = min(wl,240);

		if(plotType=='o'){
			for(int i=0;i<wl;++i) hO[LR][pad]->Fill(i,wf[i]);
		} else if(plotType=='a'){
			for(int i=0;i<wl;++i) hA[LR][pad]->Fill(i,wf[i]);
		} else if(plotType=='p'){
			if(to>from){
				auto peak = *max_element(wf.begin()+from, wf.begin()+to);
				hP[LR][pad]->Fill(peak);
			}
		} else if(plotType=='i'){
			if(to>from){
				double sum = accumulate(wf.begin()+from, wf.begin()+to, 0.0);
				hI[LR][pad]->Fill(sum);
			}
		}

		++cnt;
		readBytes += data_length;
		if(nEvt>0 && cnt>=nEvt){
			cout<<"Run end. "<<cnt<<" events processed\n";
			break;
		}
	}

	fclose(fp);
}

void Hodo_QA(int RunNo=60264, int nEvt=100000, const char* inPath="/Users/jay.ryu/KoBIC/25KEKDATA", char plotType='o'){
	gStyle->SetOptStat(1111);
	gStyle->SetTitleSize(0.04f);

	TCanvas* cL = CreateCanvas_Hodo(RunNo,"1~16");
	TCanvas* cR = CreateCanvas_Hodo(RunNo,"17~32");

	static TH2F*     hO[2][16];
	static TProfile* hA[2][16];
	static TH1F*     hP[2][16];
	static TH1F*     hI[2][16];
	static bool init=false;
	if(!init){
		init=true;
		for(int lr=0;lr<2;++lr){
			const char* side = lr==0?"L":"R";
			for(int pad=0;pad<16;++pad){
				hO[lr][pad] = new TH2F(Form("O_Hodo_%02d",pad),Form("Overlay Hodo pad %02d",pad),120,0,120,1200,-200,1000); hO[lr][pad]->SetDirectory(0);
				hA[lr][pad] = new TProfile(Form("A_Hodo_%02d",pad),Form("Avg Hodo pad %02d",pad),120,0,120,-200,1200);      hA[lr][pad]->SetDirectory(0);
				hP[lr][pad] = new TH1F(Form("P_Hodo_%02d",pad),Form("Peak Hodo pad %02d",pad),10000,0,20000);              hP[lr][pad]->SetDirectory(0);
				hI[lr][pad] = new TH1F(Form("I_Hodo_%02d",pad),Form("Int Hodo pad %02d",pad),1000,-100,150000);            hI[lr][pad]->SetDirectory(0);
			}
		}
	}

	resetAll(hO,hA,hP,hI);
	hodo_daq_quickQA(RunNo,nEvt,inPath,plotType,hO,hA,hP,hI);

	cL->cd();
	for(int i=0;i<16;++i){
		cL->cd(i+1);
		if(plotType=='o')      hO[0][i]->Draw("COLZ");
		else if(plotType=='a') hA[0][i]->Draw();
		else if(plotType=='p') hP[0][i]->Draw("hist");
		else if(plotType=='i') hI[0][i]->Draw("hist");
	}
	cL->Update();

	for(int i=0;i<16;++i){
		cR->cd(i+1);
		if(plotType=='o')      hO[1][i]->Draw("COLZ");
		else if(plotType=='a') hA[1][i]->Draw();
		else if(plotType=='p') hP[1][i]->Draw("hist");
		else if(plotType=='i') hI[1][i]->Draw("hist");
	}
	cR->Update();
}
