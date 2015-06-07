#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include <TCanvas.h>
#include <TLatex.h>
#include "TGraphErrors.h"
#include "TLegend.h"
#include <TPad.h>
#include <sstream>
#include "TVectorD.h"
#include "TGraph.h"

using std::string;
using std::cout;
using std::endl;
using std::istringstream;

void Test(){

gROOT->SetStyle("Plain");
gStyle->SetErrorX(0);

double ZGammaLL_SF = (2063*0.1*7.4*1000)/(6583032.0);

TCanvas c1("c1","Stacked Histogram",10,10,700,800);
TPad *p_2=new TPad("p_2", "p_2", 0, 0, 1, 0.35);
TPad *p_1=new TPad("p_1", "p_1", 0, 0.35, 1, 1);
p_1->SetBottomMargin(0.005);
p_1->SetFillStyle(4000);
p_1->SetFrameFillColor(0);
p_2->SetFillStyle(4000);
p_2->SetFrameFillColor(0);
p_1->SetLogy();
p_1->Draw();
p_2->Draw();
p_1->cd();

THStack hs("hs","PF MET");

Double_t rebin_array[8] = {0, 20, 40, 60, 80, 100, 160, 300};

TFile* file0 = TFile::Open("MC_OS/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
TH1F *h0 = (TH1F*)file0->Get("h_MET_ElMu");
h0->Rebin(7, "temp1", rebin_array);
temp1->SetBinContent(6, h0->GetBinContent(6)/3.0);
temp1->SetBinError(6, TMath::Sqrt(h0->GetBinContent(6)/3.0));
temp1->SetBinContent(7, h0->GetBinContent(7)/7.0);
temp1->SetBinError(7, TMath::Sqrt(h0->GetBinContent(7)/7.0));
temp1->GetXaxis()->SetRangeUser(0, 300);
temp1->Scale(ZGammaLL_SF*(116.688/2391.76)); //Loose OP
temp1->SetFillColor(kYellow);
hs.Add(temp1);

TFile* fileSS = TFile::Open("Data_SS/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
TH1F *hSS = (TH1F*)fileSS->Get("h_MET_ElMu");
hSS->Rebin(7, "temp2", rebin_array);
temp2->SetBinContent(6, hSS->GetBinContent(6)/3.0);
temp2->SetBinError(6, TMath::Sqrt(hSS->GetBinContent(6)/3.0));
temp2->SetBinContent(7, hSS->GetBinContent(7)/7.0);
temp2->SetBinError(7, TMath::Sqrt(hSS->GetBinContent(7)/7.0));
temp2->GetXaxis()->SetRangeUser(0, 300);
temp2->SetFillColor(kRed+3);
hs.Add(temp2);

TFile* file1a = TFile::Open("MC_OS/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_TauTau.root");
TH1F *h1a = (TH1F*)file1a->Get("h_MET_ElMu");
h1a->Rebin(7, "temp3", rebin_array);
temp3->SetBinContent(6, h1a->GetBinContent(6)/3.0);
temp3->SetBinError(6, TMath::Sqrt(h1a->GetBinContent(6)/3.0));
temp3->SetBinContent(7, h1a->GetBinContent(7)/7.0);
temp3->SetBinError(7, TMath::Sqrt(h1a->GetBinContent(7)/7.0));
temp3->GetXaxis()->SetRangeUser(0, 300);
temp3->Scale(ZGammaLL_SF);
temp3->SetFillColor(kBlue);
hs.Add(temp3);

hs.Draw("HIST");
hs.SetMaximum(100.0);
hs.SetMinimum(0.001);
hs.GetXaxis()->SetRangeUser(0, 300);
hs.GetXaxis()->SetTitle("MET (GeV)");
hs.GetXaxis()->SetTitle("");
hs.GetXaxis()->SetLabelSize(0.02);
hs.GetXaxis()->SetTitleSize(0.02);

c1.SaveAs("h_MET_ElMu_Test.png");
}
