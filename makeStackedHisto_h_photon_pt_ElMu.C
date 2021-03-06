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


void makeStackedHisto_h_photon_pt_ElMu(){

gROOT->SetStyle("Plain");
gStyle->SetErrorX(0);

TCanvas c1("c1","Stacked Histogram",10,10,700,800);
TPad *p_2=new TPad("p_2", "p_2", 0, 0, 1, 0.25);
TPad *p_1=new TPad("p_1", "p_1", 0, 0.25, 1, 1);
p_1->SetBottomMargin(0.005);
p_1->SetFillStyle(4000);
p_1->SetFrameFillColor(0);
p_2->SetFillStyle(4000);
p_2->SetFrameFillColor(0);
p_1->SetLogy();
p_1->Draw();
p_2->Draw();
p_1->cd();


double ZGammaLL_SF = (2063*0.1*7.4*1000)/(6583032.0);
double ZGammaNU_SF = (1.2*123.9*7.4*1000)/3169096;
double WGamma_SF = (1.2*461.6*7.4*1000)/4802339;
double TTG_SF = (1.2*2.166*7.4*1000)/71328;
double DY_SF = (1.2*3532.8*7.4*1000)/25283595;
double TT_SF = (13.43*2.0*7.4*1000)/11843493;
double WWG_SF = (0.528*1.2*7.4*1000)/215994;
double WGSE_SF = (5.873*1.2*7.4*1000)/314539;
double WGSMu_SF = (1.914*1.2*7.4*1000)/280925;
double WGSTau_SF = (0.336*1.2*7.4*1000)/49981;
double ZZ_SF = (0.1769*1.2*7.4*1000)/4792929;
double Tbar_tW_SF = (11.1*1.2*7.4*1000)/473163;
double T_tW_SF = (11.1*1.2*7.4*1000)/496553;
double DY_1050SF =  ((11050.0*1.3*0.069)*(7.4*1000))/7115972;
double WZ_SF = (1.0575*1.2*(7.4*1000))/2015278;
double signal_stop_SF_100 = (1.66*0.3983*1000*7.4*10*0.33*0.33*2)/99986;
double signal_SF = (0.0004318*(7.4)*1000*1000*1.66)/99990;

THStack hs("hs","Leading Photon p_{T} [GeV]");

TFile* file0 = TFile::Open("MC_OS/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
TH1F *h0 = (TH1F*)file0->Get("h_photon_pt_ElMu");
h0->Rebin(20);
h0->GetXaxis()->SetRangeUser(0, 300);
h0->Scale(ZGammaLL_SF*(116.688/2391.76)); //Loose OP
h0->SetFillColor(kYellow);
hs.Add(h0);

TFile* fileSS = TFile::Open("Data_SS/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
TH1F *hSS = (TH1F*)fileSS->Get("h_photon_pt_ElMu");
hSS->Rebin(20);
hSS->GetXaxis()->SetRangeUser(0, 300);
hSS->SetFillColor(kRed+3);
hs.Add(hSS);
/*
TFile* file1 = TFile::Open("MC_OS/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_ElMu.root");
TH1F *h1 = (TH1F*)file1->Get("h_photon_pt_ElMu");
h1->Rebin(20);
h1->GetXaxis()->SetRangeUser(0, 300);
h1->Scale(ZGammaLL_SF);
h1->SetFillColor(kRed);
hs.Add(h1);
*/

TFile* file1a = TFile::Open("MC_OS/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_TauTau.root");
TH1F *h1a = (TH1F*)file1a->Get("h_photon_pt_ElMu");
h1a->Rebin(20);
h1a->GetXaxis()->SetRangeUser(0, 300);
h1a->Scale(ZGammaLL_SF);
h1a->SetFillColor(kBlue);
hs.Add(h1a);

TFile* file2 = TFile::Open("MC_OS/Output_LowPtSUSY_Tree_TTGJets_8TeV_madgraph_CaloMET_All.root");
TH1F *h2 = (TH1F*)file2->Get("h_photon_pt_ElMu");
h2->Rebin(20);
h2->GetXaxis()->SetRangeUser(0, 300);
h2->Scale(TTG_SF);
h2->SetFillColor(kCyan);
hs.Add(h2);

TFile* file3 = new TFile("MC_OS/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star_CaloMET_All.root");
TH1F *h3 = (TH1F*)file3->Get("h_photon_pt_ElMu");
h3->Rebin(20);
h3->GetXaxis()->SetRangeUser(0, 300);
h3->Scale(WGSMu_SF);
h3->SetFillColor(kOrange);
hs.Add(h3);

TFile* file4 = new TFile("MC_OS/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_CaloMET_All.root");
TH1F *h4 = (TH1F*)file4->Get("h_photon_pt_ElMu");
h4->Rebin(20);
h4->GetXaxis()->SetRangeUser(0, 300);
h4->Scale(ZZ_SF);
h4->SetFillColor(kGreen);
hs.Add(h4);

TFile* file5 = new TFile("MC_OS/Output_LowPtSUSY_Tree_T_tW_CaloMET_All.root");
TH1F *h5 = (TH1F*)file5->Get("h_photon_pt_ElMu");
h5->Rebin(20);
h5->GetXaxis()->SetRangeUser(0, 300);
h5->Scale(T_tW_SF);
h5->SetFillColor(kGray);
hs.Add(h5);

TFile* file6 = new TFile("MC_OS/Output_LowPtSUSY_Tree_Tbar_tW_CaloMET_All.root");
TH1F *h6 = (TH1F*)file6->Get("h_photon_pt_ElMu");
h6->Rebin(20);
h6->GetXaxis()->SetRangeUser(0, 300);
h6->Scale(Tbar_tW_SF);
h6->SetFillColor(kGray);
hs.Add(h6);

TFile* file8 = new TFile("MC_OS/Output_LowPtSUSY_Tree_WZJetsTo3LNu_CaloMET_All.root");
TH1F *h8 = (TH1F*)file8->Get("h_photon_pt_ElMu");
h8->Rebin(20);
h8->GetXaxis()->SetRangeUser(0, 300);
h8->Scale(WZ_SF);
h8->SetFillColor(kViolet+1);
hs.Add(h8);

TFile* file9 = new TFile("MC_OS/Output_LowPtSUSY_Tree_WWGJets_8TeV_madgraph_CaloMET_All.root");
TH1F *h9 = (TH1F*)file9->Get("h_photon_pt_ElMu");
h9->Rebin(20);
h9->GetXaxis()->SetRangeUser(0, 300);
h9->Scale(WWG_SF);
h9->SetFillColor(kMagenta);
hs.Add(h9);

hs.Draw("HIST");
hs.SetMaximum(1000.0);
hs.SetMinimum(0.1);
hs.GetXaxis()->SetRangeUser(0, 300);
hs.GetXaxis()->SetLabelSize(0.02);
hs.GetXaxis()->SetTitleSize(0.02);

TFile* file_Data = new TFile("Data_OS/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
TH1F *h_Data = (TH1F*)file_Data->Get("h_photon_pt_ElMu");
h_Data->Rebin(20);
h_Data->SetLineColor(kBlack);
h_Data->SetMarkerStyle(20);
h_Data->SetMarkerSize(1.0);
h_Data->GetXaxis()->SetRangeUser(0, 300);
h_Data->Draw("SAME E");

TFile* file_Signal = TFile::Open("Signal/Output_LowPtSUSY_Tree_Stop100_Chargino75_Neutralino30_Ntuple_Pta20.root");
TH1F *h_Signal = (TH1F*)file_Signal->Get("h_photon_pt_ElMu");
h_Signal->Rebin(20);
h_Signal->GetXaxis()->SetRangeUser(0, 300);
h_Signal->Scale(signal_stop_SF_100);
h_Signal->SetLineColor(kBlack);
h_Signal->SetLineWidth(3);
h_Signal->SetLineStyle(2);
h_Signal->Draw("SAME HIST");

TFile* file_Signal = TFile::Open("Signal/Output_LowPtSUSY_Tree_Chargino1.root");
TH1F *h_Signal_Stop = (TH1F*)file_Signal->Get("h_photon_pt_ElMu");
h_Signal_Stop->Rebin(20);
h_Signal_Stop->GetXaxis()->SetRangeUser(0, 300);
h_Signal_Stop->Scale(signal_SF);
h_Signal_Stop->SetLineColor(kRed);
h_Signal_Stop->SetLineWidth(3);
h_Signal_Stop->SetLineStyle(2);
h_Signal_Stop->Draw("SAME HIST");

TH1F *h_AllMC=(TH1F*)h_Data->Clone("h_AllMC");
h_AllMC->Reset();
h_AllMC->Add(hSS);
h_AllMC->Add(h0);
//h_AllMC->Add(h1);
h_AllMC->Add(h1a);
h_AllMC->Add(h2);
h_AllMC->Add(h3);
h_AllMC->Add(h4);
h_AllMC->Add(h5);
h_AllMC->Add(h6);
h_AllMC->Add(h8);
h_AllMC->Add(h9);

TH1F *h_AllMC_Unc=(TH1F*)h_AllMC->Clone("h_AllMC_Unc");
for (int ibin = 1; ibin < h_AllMC_Unc->GetNbinsX()+1; ibin++){
      
  double uncStat = 0;
  double uncSyst = 0;
  double uncTot  = 0;
  uncStat += hSS->GetBinError(ibin)*hSS->GetBinError(ibin) + h1a->GetBinError(ibin)*h1a->GetBinError(ibin) + h0->GetBinError(ibin)*h0->GetBinError(ibin) + h2->GetBinError(ibin)*h2->GetBinError(ibin) + h3->GetBinError(ibin)*h3->GetBinError(ibin) + h4->GetBinError(ibin)*h4->GetBinError(ibin) + h5->GetBinError(ibin)*h5->GetBinError(ibin) + h6->GetBinError(ibin)*h6->GetBinError(ibin) + h8->GetBinError(ibin)*h8->GetBinError(ibin) + h9->GetBinError(ibin)*h9->GetBinError(ibin);
  uncSyst += h1a->GetBinContent(ibin)*h1a->GetBinContent(ibin)*0.05*0.05 + h0->GetBinContent(ibin)*h0->GetBinContent(ibin)*0.50*0.50 + h2->GetBinContent(ibin)*h2->GetBinContent(ibin)*0.50*.50 + h3->GetBinContent(ibin)*h3->GetBinContent(ibin)*0.50*.50 + h4->GetBinContent(ibin)*h4->GetBinContent(ibin)*0.50*.50 + h5->GetBinContent(ibin)*h5->GetBinContent(ibin)*0.50*.50 + h6->GetBinContent(ibin)*h6->GetBinContent(ibin)*0.50*.50 + h8->GetBinContent(ibin)*h8->GetBinContent(ibin)*0.50*.50 + h9->GetBinContent(ibin)*h9->GetBinContent(ibin)*0.50*.50;

  uncTot = sqrt(uncStat + uncSyst);
  h_AllMC_Unc->SetBinError(ibin, uncTot);

}

gStyle->SetHatchesLineWidth(4);
gStyle->SetErrorX(0.5);
h_AllMC_Unc->SetFillStyle(3005);
h_AllMC_Unc->SetFillColor(1);
h_AllMC_Unc->SetMarkerStyle(1);
h_AllMC_Unc->Draw("SAME E2");

TLegend *leg1 = new TLegend(0.62,0.45,0.90,0.90,NULL,"brNDC"); 
leg1->SetBorderSize(0);
leg1->SetTextSize(0.03);
leg1->SetLineColor(1);
leg1->SetLineStyle(0);
leg1->SetLineWidth(1);
leg1->SetFillColor(10);
leg1->SetFillStyle(0);
leg1->AddEntry(h0,"Data driven jet #rightarrow #gamma","f");
//leg1->AddEntry(h1,"Z#gamma #rightarrow #mu#mu#gamma","f");
leg1->AddEntry(h1a,"Z#gamma #rightarrow #tau #tau#gamma","f");
leg1->AddEntry(h2,"t#bar{t}#gamma","f");
leg1->AddEntry(h3,"WGStar","f");
leg1->AddEntry(h4,"ZZ","f");
leg1->AddEntry(h5,"Single T","f");
leg1->AddEntry(h8,"WZ","f");
leg1->AddEntry(h9,"WWG","f");
leg1->AddEntry(hSS,"Data: Same signed","f");
leg1->AddEntry(h_Data,"Data","lp");
leg1->AddEntry(h_Signal,"Stop 100 X 10","lp");
leg1->AddEntry(h_Signal_Stop,"Chargino 100 X 1000","lp");
leg1->Draw();

TLatex* text1 = new TLatex(2.570061,23.08044,"CMS Preliminary, 7.4 fb^{-1} at #sqrt{s} = 8 TeV");
text1->SetNDC();
text1->SetTextAlign(13);
text1->SetX(0.15);
text1->SetY(0.85);
text1->SetTextFont(42);
text1->SetTextSizePixels(20);
text1->Draw();

p_2->cd();
p_2->SetGridy();

TH1F *h_ratio=(TH1F*)h_Data->Clone("h_ratio");
h_ratio->SetLabelSize(0.05);
h_ratio->SetTitleSize(0.05);
h_ratio->SetTitle("; #bf{Leading photon pT [GeV]}; Data/MC Ratio");
h_ratio->SetStats(kFALSE);
h_ratio->Divide(h_AllMC);
h_ratio->SetLineColor(kBlack);
h_ratio->SetMarkerStyle(20);
h_ratio->SetMinimum(-0.5); 
h_ratio->SetMaximum(3.0);
h_ratio->GetXaxis()->SetRangeUser(0, 300.0); 
h_ratio->Draw();

TF1 *fit_ratio = new TF1("fit_ratio","[0]*x + [1]",0.0, 300.0);
fit_ratio->SetParLimits(0,0.0,0.00001);
fit_ratio->SetParLimits(1,0.0, 1.2);
fit_ratio->SetLineColor(kRed);
fit_ratio->SetLineWidth(3);
h_ratio->Fit("fit_ratio", "", "", 0.0,300.0);

TH1F *h_ratio_Unc=(TH1F*)h_ratio->Clone("h_ratio_Unc");
for (int ibin = 1; ibin < h_ratio_Unc->GetNbinsX()+1; ibin++){
  h_ratio_Unc->SetBinError(ibin, (h_ratio->GetBinContent(ibin)*h_AllMC_Unc->GetBinError(ibin))/h_AllMC->GetBinContent(ibin));
}

gStyle->SetHatchesLineWidth(4);
gStyle->SetErrorX(0.5);
h_ratio_Unc->GetXaxis()->SetRangeUser(32, 300);
h_ratio_Unc->SetFillStyle(3005);
h_ratio_Unc->SetFillColor(1);
h_ratio_Unc->SetMarkerStyle(1);
h_ratio_Unc->Draw("SAME E2");


c1.SaveAs("h_photon_pt_ElMu_Stacked.pdf");
c1.SaveAs("h_photon_pt_ElMu_Stacked.png");
}
