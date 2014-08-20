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


void makeStackedHisto_InvMassPh(){

gROOT->SetStyle("Plain");
gStyle->SetErrorX(0);

TCanvas c1("c1","Stacked Histogram",10,10,700,700);
TPad *p_2=new TPad("p_2", "p_2", 0, 0, 1, 0.35);
TPad *p_1=new TPad("p_1", "p_1", 0, 0.35, 1, 1);
p_1->SetBottomMargin(0.05);
p_1->SetFillStyle(4000);
p_1->SetFrameFillColor(0);
p_2->SetFillStyle(4000);
p_2->SetFrameFillColor(0);
p_1->Draw();
p_2->Draw();
p_1->cd();

double ZGammaLL_SF = (1.2*132.6*7.4*1000)/5750111;
double ZGammaNU_SF = (1.2*123.9*7.4*1000)/3169096;
double WGamma_SF = (1.2*461.6*7.4*1000)/4802339;
double TTG_SF = (1.2*2.166*7.4*1000)/71328;
double DY_SF = (1.2*3532.8*7.4*1000)/25283595;
double TT_SF = (13.43*2.0*7.4*1000)/11843493;
double WWG_SF = (0.528*1.2*7.4*1000)/218079;
double WGSE_SF = (5.873*1.2*7.4*1000)/314539;
double WGSMu_SF = (1.914*1.2*7.4*1000)/302312;
double WGSTau_SF = (0.336*1.2*7.4*1000)/49981;
double ZZ_SF = (0.1769*1.2*7.4*1000)/4781210;
double Tbar_tW_SF = (11.1*1.2*7.4*1000)/491809; 
double T_tW_SF = (11.1*1.2*7.4*1000)/496558;
double DY_1050SF =  ((11050.0*1.3*0.069)*(7.4*1000))/7115972; 
double WZ_SF = (1.0575*1.2*(7.4*1000))/2016285;

THStack hs("hs","Invariant Mass: M_{#mu#mu#gamma} (GeV)");

TFile* file1 = TFile::Open("MC_MuMu/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
TH1F *h1 = (TH1F*)file1->Get("h_InvariantMass_MuPh");
h1->Rebin(50);
h1->GetXaxis()->SetRangeUser(110, 185);
h1->Scale(ZGammaLL_SF);
h1->SetFillColor(kRed);
hs.Add(h1);

TFile* file2 = TFile::Open("MC_MuMu/Output_LowPtSUSY_Tree_TTGJets_8TeV_madgraph_CaloMET_All.root");
TH1F *h2 = (TH1F*)file2->Get("h_InvariantMass_MuPh");
h2->Rebin(50);
h2->GetXaxis()->SetRangeUser(110, 185);
h2->Scale(TTG_SF);
h2->SetFillColor(kCyan);
hs.Add(h2);

TFile* file3 = new TFile("MC_MuMu/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star_CaloMET_All.root");
TH1F *h3 = (TH1F*)file3->Get("h_InvariantMass_MuPh");
h3->Rebin(50);
h3->GetXaxis()->SetRangeUser(110, 185);
h3->Scale(WGSMu_SF);
h3->SetFillColor(kOrange);
hs.Add(h3);

TFile* file4 = new TFile("MC_MuMu/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_CaloMET_All.root");
TH1F *h4 = (TH1F*)file4->Get("h_InvariantMass_MuPh");
h4->Rebin(50);
h4->GetXaxis()->SetRangeUser(110, 185);
h4->Scale(ZZ_SF);
h4->SetFillColor(kGreen);
hs.Add(h4);

TFile* file5 = new TFile("MC_MuMu/Output_LowPtSUSY_Tree_T_tW_CaloMET_All.root");
TH1F *h5 = (TH1F*)file5->Get("h_InvariantMass_MuPh");
h5->Rebin(50);
h5->GetXaxis()->SetRangeUser(110, 185);
h5->Scale(T_tW_SF);
h5->SetFillColor(kGray);
hs.Add(h5);

TFile* file6 = new TFile("MC_MuMu/Output_LowPtSUSY_Tree_Tbar_tW_CaloMET_All.root");
TH1F *h6 = (TH1F*)file6->Get("h_InvariantMass_MuPh");
h6->Rebin(50);
h6->GetXaxis()->SetRangeUser(110, 185);
h6->Scale(Tbar_tW_SF);
h6->SetFillColor(kGray);
hs.Add(h6);

TFile* file7 = new TFile("MC_MuMu/Output_LowPtSUSY_Tree_DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraphCaloMET_All.root");
TH1F *h7 = (TH1F*)file7->Get("h_InvariantMass_MuPh");
h7->Rebin(50);
h7->GetXaxis()->SetRangeUser(110, 185);
h7->Scale(DY_SF);
h7->SetFillColor(kBlue);
hs.Add(h7);

TFile* file8 = new TFile("MC_MuMu/Output_LowPtSUSY_Tree_WZJetsTo3LNu_CaloMET_All.root");
TH1F *h8 = (TH1F*)file8->Get("h_InvariantMass_MuPh");
h8->Rebin(50);
h8->GetXaxis()->SetRangeUser(110, 185);
h8->Scale(WZ_SF);
h8->SetFillColor(kViolet+1);
hs.Add(h8);

TFile* file9 = new TFile("MC_MuMu/Output_LowPtSUSY_Tree_WWGJets_8TeV_madgraph_CaloMET_All.root");
TH1F *h9 = (TH1F*)file9->Get("h_InvariantMass_MuPh");
h9->Rebin(50);
h9->GetXaxis()->SetRangeUser(110, 185);
h9->Scale(WWG_SF);
h9->SetFillColor(kYellow);
hs.Add(h9);

hs.Draw("HIST");
hs.SetMaximum(15.0);
hs.SetMinimum(1.0);
hs.GetXaxis()->SetRangeUser(110, 185);
hs.GetXaxis()->SetTitle("Invariant Mass: M_{#mu#mu#gamma} (GeV)");
hs.GetXaxis()->SetLabelSize(0.02);
hs.GetXaxis()->SetTitleSize(0.02);

TFile* file_Data = new TFile("Data_MuMu/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
TH1F *h_Data = (TH1F*)file_Data->Get("h_InvariantMass_MuPh");
h_Data->Rebin(50);
h_Data->SetLineColor(kBlack);
h_Data->SetMarkerStyle(20);
h_Data->SetMarkerSize(1.0);
h_Data->GetXaxis()->SetRangeUser(110, 185);
h_Data->Draw("SAME E");

TLegend *leg1 = new TLegend(0.75,0.60,0.90,0.90,NULL,"brNDC"); 
leg1->SetBorderSize(0);
leg1->SetTextSize(0.03);
leg1->SetLineColor(1);
leg1->SetLineStyle(0);
leg1->SetLineWidth(1);
leg1->SetFillColor(10);
leg1->SetFillStyle(0);
leg1->AddEntry(h1,"Z#gamma #rightarrow ll#gamma","f");
leg1->AddEntry(h2,"t#bar{t}#gamma","f");
leg1->AddEntry(h3,"WGStar","f");
leg1->AddEntry(h4,"ZZ","f");
leg1->AddEntry(h5,"Single T","f");
leg1->AddEntry(h7,"Drell-Yan","f");
leg1->AddEntry(h8,"WZ","f");
leg1->AddEntry(h9,"WWG","f");
leg1->AddEntry(h_Data,"Data","lp");
leg1->Draw();

p_2->cd();
p_2->SetGridy();
TH1F *h_ratio=(TH1F*)h_Data->Clone("h_ratio");
h_ratio->SetTitle("; Invariant Mass: M_{#mu#mu#gamma} (GeV); Data/MC Ratio");
TH1F *h_AllMC=(TH1F*)h_Data->Clone("h_AllMC");
h_AllMC->Reset();
h_AllMC->Add(h1);
h_AllMC->Add(h2);
h_AllMC->Add(h3);
h_AllMC->Add(h4);
h_AllMC->Add(h5);
h_AllMC->Add(h6);
h_AllMC->Add(h7);
h_AllMC->Add(h8);
h_AllMC->Add(h9);
h_ratio->SetStats(kFALSE);
h_ratio->Divide(h_AllMC);
h_ratio->SetLineColor(1);
h_ratio->SetMarkerStyle(20);
h_ratio->SetMinimum(-0.5); 
h_ratio->SetMaximum(2.0); 
h_ratio->Draw();

c1.SaveAs("h_InvariantMass_MuPh_Stacked.pdf");
c1.SaveAs("h_InvariantMass_MuPh_Stacked.png");
}
