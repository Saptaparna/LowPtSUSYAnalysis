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


void makeStackedHisto_MET(){

gROOT->SetStyle("Plain");
gStyle->SetErrorX(0);

TCanvas c1("c1","Stacked Histogram",10,10,600,400);
c1.SetLogy();

double ZGammaLL_SF = (1.2*132.6*7.3*1000)/6583032;
double ZGammaNU_SF = (1.2*123.9*7.3*1000)/3169096;
double WGamma_SF = (1.2*461.6*7.3*1000)/4802339;
double TTG_SF = (1.2*2.166*7.3*1000)/71328;
double DY_SF = (1.2*3532.8*7.3*1000)/27457743;
double TT_SF = (13.43*2.0*7.3*1000)/11843493;
double WWG_SF = (0.528*1.2*7.3*1000)/214538;
double WGSE_SF = (5.873*1.2*7.3*1000)/314539;
double WGSMu_SF = (1.914*1.2*7.3*1000)/302312;
double WGSTau_SF = (0.336*1.2*7.3*1000)/49981;
double ZZ_SF = (0.1769*1.2*7.3*1000)/18393019;
double Tbar_tW_SF = (11.1*1.2*7.3*1000)/492542; 
double T_tW_SF = (11.1*1.2*7.3*1000)/496681;
double DY_1050SF =  ((11050.0*1.3*0.069)*(7.3*1000))/7115972; 

THStack hs("hs","Stacked Plot of MET");

TFile* file5 = new TFile("MC_MuMu/Output_LowPtSUSY_Tree_DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_All.root");
TH1F *h5 = (TH1F*)file5->Get("h_MET");
h5->Rebin(50);
h5->GetXaxis()->SetRangeUser(0, 300);
h5->Scale(DY_SF);
h5->SetFillColor(kMagenta);
hs.Add(h5);

TFile* file1 = TFile::Open("MC_MuMu/Output_LowPtSUSY_Tree_DYJetsToLL_M_10To50filter_8TeV_madgraph_All.root");
TH1F *h1 = (TH1F*)file1->Get("h_MET");
h1->Rebin(50);
h1->GetXaxis()->SetRangeUser(0, 300);
h1->Scale(DY_1050SF);
h1->SetFillColor(kMagenta);
hs.Add(h1);

TFile* file2 = TFile::Open("MC_MuMu/Output_LowPtSUSY_Tree_ZGToLLG_v2_All.root");
TH1F *h2 = (TH1F*)file2->Get("h_MET");
h2->Rebin(50);
h2->GetXaxis()->SetRangeUser(0, 300);
h2->Scale(ZGammaLL_SF);
h2->SetFillColor(kRed);
hs.Add(h2);

TFile* file4 = TFile::Open("MC_MuMu/Output_LowPtSUSY_Tree_TTG.root");
TH1F *h4 = (TH1F*)file4->Get("h_MET");
h4->Rebin(50);
h4->GetXaxis()->SetRangeUser(0, 300);
h4->Scale(TTG_SF);
h4->SetFillColor(kCyan);
hs.Add(h4);

TFile* file6 = new TFile("MC_MuMu/Output_LowPtSUSY_Tree_TTbar_8TeV_madgraph_All.root");
TH1F *h6 = (TH1F*)file6->Get("h_MET");
h6->Rebin(50);
h6->GetXaxis()->SetRangeUser(0, 300);
h6->Scale(TT_SF);
h6->SetFillColor(kBlue+1);
hs.Add(h6);

TFile* file7 = new TFile("MC_MuMu/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star.root");
TH1F *h7 = (TH1F*)file7->Get("h_MET");
h7->Rebin(50);
h7->GetXaxis()->SetRangeUser(0, 300);
h7->Scale(WGSMu_SF);
h7->SetFillColor(kOrange);
hs.Add(h7);

TFile* file8 = new TFile("MC_MuMu/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_All.root");
TH1F *h8 = (TH1F*)file8->Get("h_MET");
h8->Rebin(50);
h8->GetXaxis()->SetRangeUser(0, 300);
h8->Scale(ZZ_SF);
h8->SetFillColor(kGreen);
hs.Add(h8);

TFile* file9 = new TFile("MC_MuMu/Output_LowPtSUSY_Tree_T_tW_All.root");
TH1F *h9 = (TH1F*)file9->Get("h_MET");
h9->Rebin(50);
h9->GetXaxis()->SetRangeUser(0, 300);
h9->Scale(T_tW_SF);
h9->SetFillColor(kGray);
hs.Add(h9);

TFile* file10 = new TFile("MC_MuMu/Output_LowPtSUSY_Tree_Tbar_tW_All.root");
TH1F *h10 = (TH1F*)file10->Get("h_MET");
h10->Rebin(50);
h10->GetXaxis()->SetRangeUser(0, 300);
h10->Scale(Tbar_tW_SF);
h10->SetFillColor(kGray);
hs.Add(h10);

hs.Draw("HIST");
//hs.SetMinimum(1.0);
hs.GetXaxis()->SetRangeUser(0, 300);

TFile* file11 = new TFile("Data_MuMu/Output_Data_ParkedDataFixed_All.root");
TH1F *h11 = (TH1F*)file11->Get("h_MET");
h11->Rebin(50);
h11->SetLineColor(kBlack);
h11->SetMarkerStyle(20);
h11->SetMarkerSize(1.0);
h11->GetXaxis()->SetRangeUser(0, 300);
h11->Draw("SAME E");

TLegend *leg1 = new TLegend(0.55,0.55,0.90,0.90,NULL,"brNDC");
leg1->SetBorderSize(0);
leg1->SetTextSize(0.03);
leg1->SetLineColor(1);
leg1->SetLineStyle(0);
leg1->SetLineWidth(1);
leg1->SetFillColor(10);
leg1->SetFillStyle(0);
leg1->AddEntry(h2,"Z#gamma #rightarrow ll#gamma","f");
leg1->AddEntry(h4,"t#bar{t}#gamma","f");
leg1->AddEntry(h5,"Drell Yan","f");
leg1->AddEntry(h6,"t#bar{t}","f");
leg1->AddEntry(h7,"WGStar","f");
leg1->AddEntry(h9,"Single T","f");
leg1->AddEntry(h11,"Data","lp");
leg1->Draw();


/*TFile* file5 = TFile::Open("Output_PULooseJetID.root");
TH1F *h5 = (TH1F*)file5->Get("h_MET");
h5->GetXaxis()->SetRangeUser(1, 450);
h5->SetMarkerSize(1.3);
h5->Draw("SAME e");
*/
c1.SaveAs("h_MET_Stacked.pdf");
c1.SaveAs("h_MET_Stacked.png");
}
