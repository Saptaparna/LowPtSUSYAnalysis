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
#include "TROOT.h"

using std::string;
using std::cout;
using std::endl;
using std::istringstream;


void makeStackedHisto_MET_ElMu_Pretty(){

gROOT->SetStyle("Plain");
gStyle->SetErrorX(0);

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

THStack hs("hs","PF MET");

Double_t rebin_array[8] = {0, 20, 40, 60, 80, 100, 160, 300};

TFile* file0 = TFile::Open("MC_OS/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
TH1F *h0 = (TH1F*)file0->Get("h_MET_ElMu");
TH1F *temp1 = (TH1F*)h0->Rebin(7, "temp1", rebin_array);	
temp1->SetBinContent(6, h0->GetBinContent(6)/3.0);
temp1->SetBinError(6, TMath::Sqrt(h0->GetBinContent(6)/3.0));
temp1->SetBinContent(7, h0->GetBinContent(7)/7.0);
temp1->SetBinError(7, TMath::Sqrt(h0->GetBinContent(7)/7.0));
temp1->GetXaxis()->SetRangeUser(0, 300);
temp1->Scale(ZGammaLL_SF*(116.688/2391.76)); //Loose OP
temp1->SetFillColor(kYellow);
hs.Add(temp1);
cout << "temp1->GetNbinsX() = " << temp1->GetNbinsX() << endl;

TFile* fileSS = TFile::Open("Data_SS/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
TH1F *hSS = (TH1F*)fileSS->Get("h_MET_ElMu");
TH1F *temp2 = (TH1F*)hSS->Rebin(7, "temp2", rebin_array);     
temp2->SetBinContent(6, hSS->GetBinContent(6)/3.0);
temp2->SetBinError(6, TMath::Sqrt(hSS->GetBinContent(6)/3.0));
temp2->SetBinContent(7, hSS->GetBinContent(7)/7.0);
temp2->SetBinError(7, TMath::Sqrt(hSS->GetBinContent(7)/7.0));
temp2->GetXaxis()->SetRangeUser(0, 300);
temp2->SetFillColor(kRed+3);
hs.Add(temp2);
cout << "temp2->GetNbinsX() = " << temp2->GetNbinsX() << endl;

TFile* file1a = TFile::Open("MC_OS/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_TauTau.root");
TH1F *h1a = (TH1F*)file1a->Get("h_MET_ElMu");
TH1F *temp3 = (TH1F*)h1a->Rebin(7, "temp3", rebin_array);
temp3->SetBinContent(6, h1a->GetBinContent(6)/3.0);
temp3->SetBinError(6, TMath::Sqrt(h1a->GetBinContent(6)/3.0));
temp3->SetBinContent(7, h1a->GetBinContent(7)/7.0);
temp3->SetBinError(7, TMath::Sqrt(h1a->GetBinContent(7)/7.0));
temp3->GetXaxis()->SetRangeUser(0, 300);
temp3->Scale(ZGammaLL_SF);
temp3->SetFillColor(kBlue);
hs.Add(temp3);
cout << "temp3->GetNbinsX() = " << temp3->GetNbinsX() << endl;

TFile* file2 = TFile::Open("MC_OS/Output_LowPtSUSY_Tree_TTGJets_8TeV_madgraph_CaloMET_All.root");
TH1F *h2 = (TH1F*)file2->Get("h_MET_ElMu");
TH1F *temp4 = (TH1F*)h2->Rebin(7, "temp4", rebin_array);     
temp4->SetBinContent(6, h2->GetBinContent(6)/3.0);
temp4->SetBinError(6, TMath::Sqrt(h2->GetBinContent(6)/3.0));
temp4->SetBinContent(7, h2->GetBinContent(7)/7.0);
temp4->SetBinError(7, TMath::Sqrt(h2->GetBinContent(7)/7.0));
temp4->GetXaxis()->SetRangeUser(0, 300);
temp4->Scale(TTG_SF);
temp4->SetFillColor(kCyan);
hs.Add(temp4);
cout << "temp4->GetNbinsX() = " << temp4->GetNbinsX() << endl;

TFile* file3 = new TFile("MC_OS/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star_CaloMET_All.root");
TH1F *h3 = (TH1F*)file3->Get("h_MET_ElMu");
TH1F *temp5 = (TH1F*)h3->Rebin(7, "temp5", rebin_array);     
temp5->SetBinContent(6, h3->GetBinContent(6)/3.0);
temp5->SetBinError(6, TMath::Sqrt(h3->GetBinContent(6)/3.0));
temp5->SetBinContent(7, h3->GetBinContent(7)/7.0);
temp5->SetBinError(7, TMath::Sqrt(h3->GetBinContent(7)/7.0));
temp5->GetXaxis()->SetRangeUser(0, 300);
temp5->Scale(WGSMu_SF);
temp5->SetFillColor(kOrange);
hs.Add(temp5);
cout << "temp5->GetNbinsX() = " << temp5->GetNbinsX() << endl;

TFile* file4 = new TFile("MC_OS/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_CaloMET_All.root");
TH1F *h4 = (TH1F*)file4->Get("h_MET_ElMu");
TH1F *temp6 = (TH1F*)h4->Rebin(7, "temp6", rebin_array);     
temp6->SetBinContent(6, h4->GetBinContent(6)/3.0);
temp6->SetBinError(6, TMath::Sqrt(h4->GetBinContent(6)/3.0));
temp6->SetBinContent(7, h4->GetBinContent(7)/7.0);
temp6->SetBinError(7, TMath::Sqrt(h4->GetBinContent(7)/7.0));
temp6->GetXaxis()->SetRangeUser(0, 300);
temp6->Scale(ZZ_SF);
temp6->SetFillColor(kGreen);
hs.Add(temp6);
cout << "temp6->GetNbinsX() = " << temp6->GetNbinsX() << endl;

TFile* file5 = new TFile("MC_OS/Output_LowPtSUSY_Tree_T_tW_CaloMET_All.root");
TH1F *h5 = (TH1F*)file5->Get("h_MET_ElMu");
TH1F *temp7 = (TH1F*)h5->Rebin(7, "temp7", rebin_array);     
temp7->SetBinContent(6, h5->GetBinContent(6)/3.0);
temp7->SetBinError(6, TMath::Sqrt(h5->GetBinContent(6)/3.0));
temp7->SetBinContent(7, h5->GetBinContent(7)/7.0);
temp7->SetBinError(7, TMath::Sqrt(h5->GetBinContent(7)/7.0));
temp7->GetXaxis()->SetRangeUser(0, 300);
temp7->Scale(T_tW_SF);
temp7->SetFillColor(kGray);
hs.Add(temp7);
cout << "temp7->GetNbinsX() = " << temp7->GetNbinsX() << endl;

TFile* file6 = new TFile("MC_OS/Output_LowPtSUSY_Tree_Tbar_tW_CaloMET_All.root");
TH1F *h6 = (TH1F*)file6->Get("h_MET_ElMu");
TH1F *temp8 = (TH1F*)h6->Rebin(7, "temp8", rebin_array);     
temp8->SetBinContent(6, h6->GetBinContent(6)/3.0);
temp8->SetBinError(6, TMath::Sqrt(h6->GetBinContent(6)/3.0));
temp8->SetBinContent(7, h6->GetBinContent(7)/7.0);
temp8->SetBinError(7, TMath::Sqrt(h6->GetBinContent(7)/7.0));
temp8->GetXaxis()->SetRangeUser(0, 300);
temp8->Scale(Tbar_tW_SF);
temp8->SetFillColor(kGray);
hs.Add(temp8);
cout << "temp8->GetNbinsX() = " << temp8->GetNbinsX() << endl;

TFile* file8 = new TFile("MC_OS/Output_LowPtSUSY_Tree_WZJetsTo3LNu_CaloMET_All.root");
TH1F *h8 = (TH1F*)file8->Get("h_MET_ElMu");
TH1F *temp10 = (TH1F*)h8->Rebin(7, "temp10", rebin_array);     
temp10->SetBinContent(6, h8->GetBinContent(6)/3.0);
temp10->SetBinError(6, TMath::Sqrt(h8->GetBinContent(6)/3.0));
temp10->SetBinContent(7, h8->GetBinContent(7)/7.0);
temp10->SetBinError(7, TMath::Sqrt(h8->GetBinContent(7)/7.0));
temp10->GetXaxis()->SetRangeUser(0, 300);
temp10->Scale(WZ_SF);
temp10->SetFillColor(kViolet+1);
hs.Add(temp10);
cout << "temp10->GetNbinsX() = " << temp10->GetNbinsX() << endl;

TFile* file9 = new TFile("MC_OS/Output_LowPtSUSY_Tree_WWGJets_8TeV_madgraph_CaloMET_All.root");
TH1F *h9 = (TH1F*)file9->Get("h_MET_ElMu");
TH1F *temp11 = (TH1F*)h9->Rebin(7, "temp11", rebin_array);     
temp11->SetBinContent(6, h9->GetBinContent(6)/3.0);
temp11->SetBinError(6, TMath::Sqrt(h9->GetBinContent(6)/3.0));
temp11->SetBinContent(7, h9->GetBinContent(7)/7.0);
temp11->SetBinError(7, TMath::Sqrt(h9->GetBinContent(7)/7.0));
temp11->GetXaxis()->SetRangeUser(0, 300);
temp11->Scale(WWG_SF);
temp11->SetFillColor(kMagenta);
hs.Add(temp11);
cout << "temp11->GetNbinsX() = " << temp11->GetNbinsX() << endl;

hs.Draw("HIST");
hs.SetMaximum(100.0);
hs.SetMinimum(0.01);
hs.GetXaxis()->SetRangeUser(0, 300);
//hs.GetXaxis()->SetTitle("MET [GeV]");
hs.GetXaxis()->SetTitle("");
hs.GetXaxis()->SetLabelSize(0.02);
hs.GetXaxis()->SetTitleSize(0.02);


TFile* file_Data = new TFile("Data_OS/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
TH1F *h_Data = (TH1F*)file_Data->Get("h_MET_ElMu");
TH1F *temp12 = (TH1F*)h_Data->Rebin(7, "temp12", rebin_array);     
temp12->SetBinContent(6, h_Data->GetBinContent(6)/3.0);
temp12->SetBinError(6, TMath::Sqrt(h_Data->GetBinContent(6)/3.0));
temp12->SetBinContent(7, h_Data->GetBinContent(7)/7.0);
temp12->SetBinError(7, TMath::Sqrt(h_Data->GetBinContent(7)/7.0));
temp12->SetLineColor(kBlack);
temp12->SetMarkerStyle(20);
temp12->SetMarkerSize(1.0);
temp12->GetXaxis()->SetRangeUser(0, 300);
temp12->Draw("SAME E");
cout << "temp12->GetNbinsX() = " << temp12->GetNbinsX() << endl;

TFile* file_Signal = TFile::Open("Signal/Output_LowPtSUSY_Tree_Stop100_Chargino75_Neutralino30_Ntuple_Pta20.root");
TH1F *h_Signal = (TH1F*)file_Signal->Get("h_MET_ElMu");
TH1F *temp13   = (TH1F*)h_Signal->Rebin(7, "temp13", rebin_array);     
temp13->SetBinContent(6, h_Signal->GetBinContent(6)/3.0);
temp13->SetBinError(6, TMath::Sqrt(h_Signal->GetBinContent(6)/3.0));
temp13->SetBinContent(7, h_Signal->GetBinContent(7)/7.0);
temp13->SetBinError(7, TMath::Sqrt(h_Signal->GetBinContent(7)/7.0));
temp13->GetXaxis()->SetRangeUser(0, 300);
temp13->Scale(signal_stop_SF_100);
temp13->SetLineColor(kBlack);
temp13->SetLineWidth(3);
temp13->SetLineStyle(2);
temp13->Draw("SAME HIST");
cout << "temp13->GetNbinsX() = " << temp13->GetNbinsX() << endl;

TFile* file_Signal = TFile::Open("Signal/Output_LowPtSUSY_Tree_Chargino1.root");
TH1F *h_Signal_Stop = (TH1F*)file_Signal->Get("h_MET_ElMu");
TH1F *temp14 = (TH1F*)h_Signal_Stop->Rebin(7, "temp14", rebin_array);     
temp14->SetBinContent(6, h_Signal_Stop->GetBinContent(6)/3.0);
temp14->SetBinError(6, TMath::Sqrt(h_Signal_Stop->GetBinContent(6)/3.0));
temp14->SetBinContent(7, h_Signal_Stop->GetBinContent(7)/7.0);
temp14->SetBinError(7, TMath::Sqrt(h_Signal_Stop->GetBinContent(7)/7.0));
temp14->GetXaxis()->SetRangeUser(0, 300);
temp14->Scale(signal_SF);
temp14->SetLineColor(kRed);
temp14->SetLineWidth(3);
temp14->SetLineStyle(2);
temp14->Draw("SAME HIST");
cout << "temp14->GetNbinsX() = " << temp14->GetNbinsX() << endl;

TH1F *h_AllMC=(TH1F*)temp12->Clone("h_AllMC");
h_AllMC->Reset();
h_AllMC->Add(temp1);
h_AllMC->Add(temp2);
h_AllMC->Add(temp3);
h_AllMC->Add(temp4);
h_AllMC->Add(temp5);
h_AllMC->Add(temp6);
h_AllMC->Add(temp7);
h_AllMC->Add(temp8);
h_AllMC->Add(temp10);
h_AllMC->Add(temp11);

TH1F *h_AllMC_Unc=(TH1F*)h_AllMC->Clone("h_AllMC_Unc");
for (int ibin = 1; ibin < h_AllMC_Unc->GetNbinsX()+1; ibin++){
      
  double uncStat = 0;
  double uncSyst = 0;
  double uncTot  = 0;

  uncStat += temp1->GetBinError(ibin)*temp1->GetBinError(ibin) + temp2->GetBinError(ibin)*temp2->GetBinError(ibin) + temp3->GetBinError(ibin)*temp3->GetBinError(ibin) + temp4->GetBinError(ibin)*temp4->GetBinError(ibin) + temp5->GetBinError(ibin)*temp5->GetBinError(ibin) + temp6->GetBinError(ibin)*temp6->GetBinError(ibin) + temp7->GetBinError(ibin)*temp7->GetBinError(ibin) + temp8->GetBinError(ibin)*temp8->GetBinError(ibin) + temp10->GetBinError(ibin)*temp10->GetBinError(ibin) + temp11->GetBinError(ibin)*temp11->GetBinError(ibin);
  
  uncSyst += temp1->GetBinContent(ibin)*temp1->GetBinContent(ibin)*1.00*1.00 + temp2->GetBinContent(ibin)*temp2->GetBinContent(ibin)*1.00*1.00 + temp3->GetBinContent(ibin)*temp3->GetBinContent(ibin)*0.05*0.05 + temp4->GetBinContent(ibin)*temp4->GetBinContent(ibin)*0.50*0.50 + temp5->GetBinContent(ibin)*temp5->GetBinContent(ibin)*0.50*0.50 +
 temp6->GetBinContent(ibin)*temp6->GetBinContent(ibin)*0.50*0.50 + temp7->GetBinContent(ibin)*temp7->GetBinContent(ibin)*0.15*0.15 + temp8->GetBinContent(ibin)*temp8->GetBinContent(ibin)*0.15*0.15 + temp10->GetBinContent(ibin)*temp10->GetBinContent(ibin)*0.50*0.50 + temp11->GetBinContent(ibin)*temp11->GetBinContent(ibin)*0.50*0.50;
   
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
leg1->AddEntry(temp1,"Data driven jet #rightarrow #gamma","f");
//leg1->AddEntry(h1,"Z#gamma #rightarrow #mu#mu#gamma","f");
leg1->AddEntry(temp2,"Z#gamma #rightarrow #tau #tau#gamma","f");
leg1->AddEntry(temp3,"t#bar{t}#gamma","f");
leg1->AddEntry(temp4,"WGStar","f");
leg1->AddEntry(temp6,"ZZ","f");
leg1->AddEntry(temp7,"Single T","f");
leg1->AddEntry(temp10,"WZ","f");
leg1->AddEntry(temp11,"WWG","f");
leg1->AddEntry(temp2,"Data: Same signed","f");
leg1->AddEntry(temp12,"Data","lp");
leg1->AddEntry(h_Signal,"Stop 100 X 10","lp");
leg1->AddEntry(h_Signal_Stop,"Chargino 100 X 1000","lp");
leg1->Draw();

TLatex* text1 = new TLatex(2.570061,23.08044,"CMS Preliminary, 7.4 fb^{-1} at #sqrt{s} = 8 TeV");
text1->SetNDC();
text1->SetTextAlign(13);
text1->SetX(0.18);
text1->SetY(0.982);
text1->SetTextFont(42);
text1->SetTextSizePixels(20);
text1->Draw();

p_2->cd();
p_2->SetGridy();

TH1F *h_ratio=(TH1F*)temp12->Clone("h_ratio");
h_ratio->SetLabelSize(0.05);
h_ratio->SetTitleSize(0.05);
h_ratio->SetTitle("; #bf{MET (GeV)}; Data/MC Ratio");
h_ratio->SetStats(kFALSE);
h_ratio->Divide(h_AllMC);
h_ratio->SetLineColor(kBlack);
h_ratio->SetMarkerStyle(20);
h_ratio->SetMinimum(0.0); 
h_ratio->SetMaximum(4.0);
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
h_ratio_Unc->SetFillStyle(3005);
h_ratio_Unc->SetFillColor(1);
h_ratio_Unc->SetMarkerStyle(1);
h_ratio_Unc->Draw("SAME E2");


c1.SaveAs("h_MET_ElMu_Stacked_pretty.pdf");
c1.SaveAs("h_MET_ElMu_Stacked_pretty.png");
}
