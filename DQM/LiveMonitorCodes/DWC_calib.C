// DWC_calib.C - ROOT macro for DWC timing calibration

#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

int GetDataLength(const char* inFile) {
    std::ifstream in(inFile, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "GetDataLength: Cannot open " << inFile << std::endl;
        return 0;
    }
    unsigned char header[4] = {0};
    in.read(reinterpret_cast<char*>(header), 4);
    if (!in.good()) {
        std::cerr << "GetDataLength: Failed reading header from " << inFile << std::endl;
        in.close();
        return 0;
    }
    unsigned int data_length = 0;
    for (int a = 0; a < 4; ++a) {
        data_length += (static_cast<unsigned int>(header[a]) << (8 * a));
    }
    in.close();
    return static_cast<int>(data_length);
}

float getPed(const std::vector<short>& waveform) {
    if (waveform.size() < 101) return 0.0f;
    double sum = std::accumulate(
        waveform.begin() + 1,
        waveform.begin() + 101,
        0.0);
    return static_cast<float>(sum / 100.0);
}

std::vector<float> makePedCorrected(
    const std::vector<short>& raw,
    float ped)
{
    std::vector<float> corr;
    corr.reserve(raw.size());
    for (auto s : raw) corr.push_back(ped - static_cast<float>(s));
    return corr;
}

int findThresholdBin(
    const std::vector<float>& wf,
    float thrFactor)
{
    if (wf.size() < 2) return -1;
    float maxVal = *std::max_element(wf.begin(), wf.end());
    float thr    = maxVal * thrFactor;
    for (int i = 1; i < static_cast<int>(wf.size()); ++i) {
        if (wf[i] > thr) return i;
    }
    return -1;
}

double leadingEdgeTime(
    const std::vector<float>& wf,
    int bin,
    float thrFactor,
    double dt)
{
    if (bin <= 0 || bin >= static_cast<int>(wf.size())) return NAN;
    float maxVal = *std::max_element(wf.begin(), wf.end());
    float thr    = maxVal * thrFactor;
    double x0 = (bin - 1) * dt;
    double x1 = bin * dt;
    double y0 = wf[bin - 1];
    double y1 = wf[bin];
    return x0 + (thr - y0) * (x1 - x0) / (y1 - y0);
}

void processMID(int RunNo,
                int nEvtToRead,
                const char* inPath,
                int mid,
                TH1F* h_time[4],
                TH1F* h_delta[2],
                float thrFactor,
                double dt)
{
    const int nCh = 4;
    TString inFile = Form("%s/Run_%d/Run_%d_MID_%d/FADCData_%d_%i.dat",
                          inPath, RunNo, RunNo, mid, mid, RunNo);
    std::ifstream in(inFile.Data(), std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "Cannot open file: " << inFile << std::endl;
        return;
    }
    int DLen = GetDataLength(inFile.Data());
    if (DLen > 16384) DLen = 2048;
    int nADC = (DLen - 32) / 2;

    std::vector<char> buffer(nCh * DLen);
    int packets = 0;
    while (in.peek() != EOF && packets < nEvtToRead) {
        in.read(buffer.data(), buffer.size());
        if (!in.good()) break;

        // raw waveform 추출
        std::vector<std::vector<short>> rawWF(nCh, std::vector<short>(nADC));
        for (int ch = 0; ch < nCh; ++ch) {
            for (int samp = 0; samp < nADC; ++samp) {
                int idx = 32 + 2 * samp;
                int lo = buffer[ch + nCh * idx]   & 0xFF;
                int hi = buffer[ch + nCh * (idx+1)] & 0x0F;
                rawWF[ch][samp] = static_cast<short>(lo + (hi << 8));
            }
        }

        // 채널별 leading-edge 시간 계산
        double times[4];
        for (int ch = 0; ch < nCh; ++ch) {
            float ped = getPed(rawWF[ch]);
            auto wf = makePedCorrected(rawWF[ch], ped);
            int bin = findThresholdBin(wf, thrFactor);
            times[ch] = (bin < 0)
                ? NAN
                : leadingEdgeTime(wf, bin, thrFactor, dt);
            if (!std::isnan(times[ch])) {
                h_time[ch]->Fill(times[ch] * 1e9); // ns
            }
        }
        // delta UD (2-U,3-D), LR(0-R,1-L)
        if (!std::isnan(times[2]) && !std::isnan(times[3])) {
            h_delta[0]->Fill((times[2] - times[3]) * 1e9);
        }
        if (!std::isnan(times[0]) && !std::isnan(times[1])) {
            h_delta[1]->Fill((times[1] - times[0]) * 1e9);
        }
        ++packets;
    }
    in.close();
    std::cout << Form("MID %d processed %d events\n", mid, packets);
}

void DWC_calib(int RunNo, int nEvtToRead, const char* inPath) {
    float   thrFactor = 0.3f;
    double  dt        = 2e-9;
    TString outRoot = Form("%s/DWC_calib_Run_%d.root", inPath, RunNo);
    TFile* ofile = TFile::Open(outRoot, "RECREATE");

    TH1F* h_time_mid3[4]; TH1F* h_time_mid4[4];
    TH1F* h_delta_mid3[2]; TH1F* h_delta_mid4[2];
    const char* names[4] = {"R","L","U","D"};
    for (int ch = 0; ch < 4; ++ch) {
        h_time_mid3[ch] = new TH1F(Form("h_time_DWC1_%s", names[ch]),
                                   Form("DWC1 %s time [ns]", names[ch]),
                                   1000, 0, 100);
        h_time_mid4[ch] = new TH1F(Form("h_time_DWC2_%s", names[ch]),
                                   Form("DWC2 %s time [ns]", names[ch]),
                                   1000, 0, 100);
    }
    h_delta_mid3[0] = new TH1F("h_deltaDWC1_UD","DWC1 U-D [ns]", 500,-50,50);
    h_delta_mid3[1] = new TH1F("h_deltaDWC1_LR","DWC1 L-R [ns]", 500,-50,50);
    h_delta_mid4[0] = new TH1F("h_deltaDWC2_UD","DWC2 U-D [ns]", 500,-50,50);
    h_delta_mid4[1] = new TH1F("h_deltaDWC2_LR","DWC2 L-R [ns]", 500,-50,50);

    // DWC1 (mid=3)
    processMID(RunNo, nEvtToRead, inPath, 3,
               h_time_mid3, h_delta_mid3,
               thrFactor, dt);
    // DWC2 (mid=4)
    processMID(RunNo, nEvtToRead, inPath, 4,
               h_time_mid4, h_delta_mid4,
               thrFactor, dt);

    ofile->Write();
    ofile->Close();
    std::cout << "Output written to " << outRoot << std::endl;
}
