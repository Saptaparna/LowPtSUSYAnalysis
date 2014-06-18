#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFrame.h>

void Fitting_Zpeak(std::string infile){
 

  int nbins = 9000/300;
 
  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);
  c1->SetGridx();
  c1->SetGridy();
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderMode(-1);
  c1->GetFrame()->SetBorderSize(5); 
 
  gStyle->SetErrorX(0);

  
  TFile* file = TFile::Open((infile+".root").c_str());
  TH1F *his = (TH1F*) file->Get("h_InvariantMass_Mu");
  his->Rebin(30);
  his->GetXaxis()->SetRangeUser(70, 110); 
  his->GetYaxis()->SetTitle("Events/1.0 GeV");
 
  TF1 *fSignal = new TF1("fSignal","gaus",70.0,110.0);
  TF1 *fBackground = new TF1("fBackground","pol1", 70.0,110.0);
  fBackground->SetLineColor(kBlue);
  fBackground->SetLineWidth(2);

  TF1 *fSpectrum = new TF1("fSpectrum","gaus+pol1(3)", 70.0,110.0);
  fSpectrum->SetLineColor(kRed);   
  fSpectrum->SetLineWidth(2);
  fSpectrum->SetParameters(40, his->GetMean(), his->GetRMS(), 1.0, 1.0);  
  his->SetMarkerStyle(20);
  his->SetMarkerSize(0.9);
  his->Fit("fSpectrum", "", "", 70.0, 110.0);
  Double_t param[5];
  fSpectrum->GetParameters(param);
  fSignal->SetParameters(&param[0]);
  fBackground->SetParameters(&param[3]);
  cout << "fSpectrum->GetParameter(0) = " << fSpectrum->GetParameter(0) << endl;
  cout << "fSpectrum->GetParameter(1) = " << fSpectrum->GetParameter(1) << endl;
  cout << "fSpectrum->GetParameter(2) = " << fSpectrum->GetParameter(2) << endl;
  cout << "fSpectrum->GetParameter(3) = " << fSpectrum->GetParameter(3) << endl;  
  cout << "fSpectrum->GetParameter(4) = " << fSpectrum->GetParameter(4) << endl;
 
  TH1F *hisSignal = new TH1F(*his);
  hisSignal->Sumw2();
  his->Draw("e");
  cout << "his->Integral() = " << his->Integral() << endl;
  fBackground->Draw("SAME");
  cout << "Sanity Check = " << fBackground->Integral(70, 110) + fSignal->Integral(70, 110) << " = Total = " << fSpectrum->Integral(70, 110) << endl; 
  cout << "Fitted result = " << fSpectrum->Integral(70, 110) << " Data = " << his->Integral() << endl;
  cout << "Denominator of THE Ratio = " << fSignal->Integral(70, 110) << endl; //corrected for binning 

  TF1 *fitResult = his->GetFunction("fSpectrum");
  cout << "fitResult->GetChisquare() = " << fitResult->GetChisquare() << endl;
  cout << "fitResult->GetNDF() = " << fitResult->GetNDF() << endl; 
  
  c1->SaveAs((infile+"_ZPeak.pdf").c_str());
  c1->SaveAs((infile+"_ZPeak.png").c_str());

}


void Fitting_JpsiPeak(std::string infile){

  int nbins = 9000/300;

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);
  c1->SetGridx();
  c1->SetGridy();
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderMode(-1);
  c1->GetFrame()->SetBorderSize(5);

  gStyle->SetErrorX(0);

  TFile* file = TFile::Open((infile+".root").c_str());
  TH1F *his = (TH1F*) file->Get("h_InvariantMass_Mu");
  his->GetXaxis()->SetRangeUser(2.5,3.5);
  his->GetYaxis()->SetTitle("Events/0.03 GeV");

  TF1 *fSignal = new TF1("fSignal","gaus",2.5,3.5);
  TF1 *fBackground = new TF1("fBackground","pol1", 2.5,3.5);
  TF1 *fSpectrum = new TF1("fSpectrum","gaus+pol1(3)", 2.5,3.5);
  fBackground->SetLineColor(kBlue);
  fBackground->SetLineWidth(2);

  fSpectrum->SetParameters(10, his->GetMean(), his->GetRMS(), 1.0, 1.0);
  his->SetMarkerStyle(20);
  his->SetMarkerSize(0.9);

  fSpectrum->SetLineColor(kRed);
  fSpectrum->SetLineWidth(2);
  
  his->Fit("fSpectrum", "", "", 2.5,3.5);
  Double_t param[5];
  fSpectrum->GetParameters(param);
  fSignal->SetParameters(&param[0]);
  fBackground->SetParameters(&param[3]);
  cout << "fSpectrum->GetParameter(0) = " << fSpectrum->GetParameter(0) << endl;
  cout << "fSpectrum->GetParameter(1) = " << fSpectrum->GetParameter(1) << endl;
  cout << "fSpectrum->GetParameter(2) = " << fSpectrum->GetParameter(2) << endl;
  cout << "fSpectrum->GetParameter(3) = " << fSpectrum->GetParameter(3) << endl;
  cout << "fSpectrum->GetParameter(4) = " << fSpectrum->GetParameter(4) << endl;


  TH1F *hisSignal = new TH1F(*his);
  hisSignal->Sumw2();
  his->Draw("e");
  cout << "his->Integral() = " << his->Integral() << endl;
  fBackground->Draw("SAME");

  TF1 *fitResult = his->GetFunction("fSpectrum");
  cout << "fitResult->GetChisquare() = " << fitResult->GetChisquare() << endl;
  cout << "fitResult->GetNDF() = " << fitResult->GetNDF() << endl;
  cout << "Sanity Check = " << fBackground->Integral(2.5, 3.5) + fSignal->Integral(2.5, 3.5) << " = Total = " << fSpectrum->Integral(2.5, 3.5) << endl;
  cout << "Fitted result = " << fSpectrum->Integral(2.5, 3.5)*nbins << " Data = " << his->Integral() << endl;//corrected for binning
  cout << "Numerator of THE Ratio = " << fSignal->Integral(2.5, 3.5)*nbins << endl; //corrected for binning

  c1->SaveAs((infile+"_JpsiPeak.pdf").c_str());
  c1->SaveAs((infile+"_JpsiPeak.png").c_str());
}

void Fitting_Upsilon(std::string infile){

  int nbins = 9000/300;

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);
  c1->SetGridx();
  c1->SetGridy();
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderMode(-1);
  c1->GetFrame()->SetBorderSize(5);

  TFile* file = TFile::Open((infile+".root").c_str());
  TH1F *his = (TH1F*) file->Get("h_InvariantMass_Mu");
  his->GetXaxis()->SetRangeUser(8.5,10.0);
  his->GetYaxis()->SetTitle("Events/0.03 GeV");

  TF1 *fSignal = new TF1("fSignal","gaus",8.5,10.0);
  TF1 *fBackground = new TF1("fBackground","pol1", 8.5,10.0);
  fBackground->SetLineColor(kBlue);
  fBackground->SetLineWidth(2);


  TF1 *fSpectrum = new TF1("fSpectrum","gaus+pol1(3)",8.5,10.0);

  fSpectrum->SetParameters(2.0, his->GetMean(), his->GetRMS(), 0.0, 0.0);
  his->Fit("fSpectrum", "LL", "", 8.5,10.0); //log likelihood fit
  his->SetMarkerStyle(20);
  his->SetMarkerSize(0.9);

  fSpectrum->SetLineColor(kRed);
  fSpectrum->SetLineWidth(2);


  Double_t param[5];
  fSpectrum->GetParameters(param);
  fSignal->SetParameters(&param[0]);
  fBackground->SetParameters(&param[3]);

  cout << "fSpectrum->GetParameter(0) = " << fSpectrum->GetParameter(0) << endl;
  cout << "fSpectrum->GetParameter(1) = " << fSpectrum->GetParameter(1) << endl;
  cout << "fSpectrum->GetParameter(2) = " << fSpectrum->GetParameter(2) << endl;
  cout << "fSpectrum->GetParameter(3) = " << fSpectrum->GetParameter(3) << endl;
  cout << "fSpectrum->GetParameter(4) = " << fSpectrum->GetParameter(4) << endl;

  TH1F *hisSignal = new TH1F(*his);
  hisSignal->Sumw2();
  his->Draw("e");
  cout << "his->Integral() = " << his->Integral() << endl;
  fBackground->Draw("SAME");
  cout << "fBackground->Integral(8.5,10.0) = " << fBackground->Integral(8.0,10.0) << endl;
  cout << "fSignal->Integral(8.5, 10.0) = " << fSignal->Integral(8.5, 10.0) << endl;
  cout << "fSpectrum->Integral(8.5, 10.0) = " << fSpectrum->Integral(8.5, 10.0) << endl;


  TF1 *fitResult = his->GetFunction("fSpectrum");
  cout << "fitResult->GetChisquare() = " << fitResult->GetChisquare() << endl;
  cout << "fitResult->GetNDF() = " << fitResult->GetNDF() << endl;
  cout << "Sanity Check = " << fBackground->Integral(8.5, 10.0) + fSignal->Integral(8.5, 10.0) << " = Total = " << fSpectrum->Integral(8.5, 10.0) << endl;
  cout << "Fitted result = " << fSpectrum->Integral(8.5, 10.0)*nbins << " Data = " << his->Integral() << endl;


  c1->SaveAs((infile+"_Upsilon_LL.pdf").c_str());
  c1->SaveAs((infile+"_Upsilon_LL.png").c_str());

}

void Fitting_MET_TurnOn(std::string infile){

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);

  gStyle->SetErrorX(0);

  TFile* file = TFile::Open((infile+".root").c_str());
  TH1F *h_MET_HLTPhoIdMet = (TH1F*) file->Get("h_MET_HLTPhoIdMet");
  h_MET_HLTPhoIdMet->Rebin(3);
  TH1F *h_MET_HLTPhoId = (TH1F*) file->Get("h_MET_HLTPhoId");
  h_MET_HLTPhoId->Rebin(3);
  TGraphAsymmErrors *Eff_MET = new TGraphAsymmErrors;
  Eff_MET->BayesDivide(h_MET_HLTPhoIdMet, h_MET_HLTPhoId, "");
  Eff_MET->SetLineColor(kBlue);
  Eff_MET->GetXaxis()->SetRangeUser(25, 800);
  Eff_MET->GetYaxis()->SetTitle("Efficiency");
  Eff_MET->GetXaxis()->SetTitle("Offline MET");
  Eff_MET->SetMinimum(0.0);
  Eff_MET->SetTitle("Efficiency: HLT_Photon30_R9Id90_Met25 w.r.t HLT_Photon30_R9Id90");

  Eff_MET->Draw("AP");
  
  TF1 *fit = new TF1("fit","0.5*[2]*(1+TMath::Erf((x-[0])/[1]))",25.0, 800.0); 
  fit->SetParameter(0, 20);
  fit->SetParameter(1, 40);
  fit->SetParameter(2, 0.8);
  fit->SetParLimits(2,0.0,1.0);
  Eff_MET->Fit("fit", "", "", 25.0,800.0); 

  fit->SetLineColor(kRed);
  fit->SetLineWidth(3);

  c1->SaveAs((infile+"_Efficiency_MET_TurnOn.pdf").c_str());
  c1->SaveAs((infile+"_Efficiency_MET_TurnOn.png").c_str());

}

void Fitting_Trigger_MET(std::string infile){

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);

  gStyle->SetErrorX(0);

  TFile* file = TFile::Open((infile+".root").c_str());
  TH1F *h_MET_HLTPhoIdMet = (TH1F*) file->Get("h_MET_HLTPhoIdMet");
  h_MET_HLTPhoIdMet->Rebin(3);
  TH1F *h_MET_HLTPhoId = (TH1F*) file->Get("h_MET_HLTPhoId");
  h_MET_HLTPhoId->Rebin(3);
  TGraphAsymmErrors *Eff_MET = new TGraphAsymmErrors;
  Eff_MET->BayesDivide(h_MET_HLTPhoIdMet, h_MET_HLTPhoId, "");
  Eff_MET->SetLineColor(kBlue);
  Eff_MET->GetXaxis()->SetRangeUser(0, 800); 
  Eff_MET->GetYaxis()->SetTitle("Efficiency");
  Eff_MET->GetXaxis()->SetTitle("Offline MET");
  Eff_MET->SetTitle("Efficiency: HLT_Photon30_R9Id90_Met25 w.r.t HLT_Photon30_R9Id90");
  
  Eff_MET->Draw("AP");

  TH1F *h_MET_HLTPhoIdMet_MuonVeto = (TH1F*) file->Get("h_MET_HLTPhoIdMet_MuonVeto");
  h_MET_HLTPhoIdMet_MuonVeto->Rebin(3);
  TH1F *h_MET_HLTPhoId_MuonVeto = (TH1F*) file->Get("h_MET_HLTPhoId_MuonVeto");
  h_MET_HLTPhoId_MuonVeto->Rebin(3);
  TGraphAsymmErrors *Eff_MET_MuonVeto = new TGraphAsymmErrors;
  Eff_MET_MuonVeto->BayesDivide(h_MET_HLTPhoIdMet_MuonVeto, h_MET_HLTPhoId_MuonVeto, ""); 
  Eff_MET_MuonVeto->SetLineColor(kRed);
  Eff_MET_MuonVeto->Draw("SAME");

  TLegend *leg1 = new TLegend(0.614094,0.1538462,0.8993289,0.3321678,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->AddEntry(Eff_MET,"Without muon veto","l");
  leg1->AddEntry(Eff_MET_MuonVeto,"With muon veto","l");
  leg1->Draw();

  c1->SaveAs((infile+"_Efficiency_Trigger_MET.pdf").c_str());
  c1->SaveAs((infile+"_Efficiency_Trigger_MET.png").c_str());

}

void Fitting_Trigger_MET_HT(std::string infile){

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);

  gStyle->SetErrorX(0);

  TFile* file = TFile::Open((infile+".root").c_str());
  TH1F *h_MET_HLTPhoIdMet_HighHT = (TH1F*) file->Get("h_MET_HLTPhoIdMet_HighHT");
  h_MET_HLTPhoIdMet_HighHT->Rebin(3);
  TH1F *h_MET_HLTPhoId_HighHT = (TH1F*) file->Get("h_MET_HLTPhoId_HighHT");
  h_MET_HLTPhoId_HighHT->Rebin(3);
  TGraphAsymmErrors *Eff_MET_HighHT = new TGraphAsymmErrors;
  Eff_MET_HighHT->BayesDivide(h_MET_HLTPhoIdMet_HighHT, h_MET_HLTPhoId_HighHT, "");
  Eff_MET_HighHT->GetXaxis()->SetRangeUser(0, 800);
  Eff_MET_HighHT->GetYaxis()->SetTitle("Efficiency");
  Eff_MET_HighHT->GetXaxis()->SetTitle("Offline MET");
  Eff_MET_HighHT->SetLineColor(kBlue);
  Eff_MET_HighHT->SetTitle("Efficiency: HLT_Photon30_R9Id90_Met25 w.r.t HLT_Photon30_R9Id90");
  Eff_MET_HighHT->Draw("AP");

  TH1F *h_MET_HLTPhoIdMet_LowHT = (TH1F*) file->Get("h_MET_HLTPhoIdMet_LowHT");
  h_MET_HLTPhoIdMet_LowHT->Rebin(3);
  TH1F *h_MET_HLTPhoId_LowHT = (TH1F*) file->Get("h_MET_HLTPhoId_LowHT");
  h_MET_HLTPhoId_LowHT->Rebin(3);
  TGraphAsymmErrors *Eff_MET_LowHT = new TGraphAsymmErrors;
  Eff_MET_LowHT->BayesDivide(h_MET_HLTPhoIdMet_LowHT, h_MET_HLTPhoId_LowHT, "");
  Eff_MET_LowHT->SetLineColor(kRed);
  Eff_MET_LowHT->Draw("SAME");

  TLegend *leg1 = new TLegend(0.614094,0.1538462,0.8993289,0.3321678,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->AddEntry(Eff_MET_HighHT,"HT>100 [GeV]","l");
  leg1->AddEntry(Eff_MET_LowHT,"HT<100 [GeV]","l");
  leg1->Draw();

  c1->SaveAs((infile+"_Efficiency_Trigger_MET_HT.pdf").c_str());
  c1->SaveAs((infile+"_Efficiency_Trigger_MET_HT.png").c_str());

}

void Fitting_Trigger_MET_PU(std::string infile){

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);

  gStyle->SetErrorX(0);

  TFile* file = TFile::Open((infile+".root").c_str());
  TH1F *h_MET_HLTPhoIdMet_HighPU = (TH1F*) file->Get("h_MET_HLTPhoIdMet_HighPU");
  h_MET_HLTPhoIdMet_HighPU->Rebin(3);
  TH1F *h_MET_HLTPhoId_HighPU = (TH1F*) file->Get("h_MET_HLTPhoId_HighPU");
  h_MET_HLTPhoId_HighPU->Rebin(3);
  TGraphAsymmErrors *Eff_MET_HighPU = new TGraphAsymmErrors;
  Eff_MET_HighPU->BayesDivide(h_MET_HLTPhoIdMet_HighPU, h_MET_HLTPhoId_HighPU, "");
  Eff_MET_HighPU->GetXaxis()->SetRangeUser(0, 800);
  Eff_MET_HighPU->GetYaxis()->SetTitle("Efficiency");
  Eff_MET_HighPU->GetXaxis()->SetTitle("Offline MET");
  Eff_MET_HighPU->SetLineColor(kBlue);
  Eff_MET_HighPU->SetTitle("Efficiency: HLT_Photon30_R9Id90_Met25 w.r.t HLT_Photon30_R9Id90");
  Eff_MET_HighPU->Draw("AP");

  TH1F *h_MET_HLTPhoIdMet_LowPU = (TH1F*) file->Get("h_MET_HLTPhoIdMet_LowPU");
  h_MET_HLTPhoIdMet_LowPU->Rebin(3);
  TH1F *h_MET_HLTPhoId_LowPU = (TH1F*) file->Get("h_MET_HLTPhoId_LowPU");
  h_MET_HLTPhoId_LowPU->Rebin(3);
  TGraphAsymmErrors *Eff_MET_LowPU = new TGraphAsymmErrors;
  Eff_MET_LowPU->BayesDivide(h_MET_HLTPhoIdMet_LowPU, h_MET_HLTPhoId_LowPU, "");
  Eff_MET_LowPU->SetLineColor(kRed);
  Eff_MET_LowPU->Draw("SAME");
  
  TLegend *leg1 = new TLegend(0.614094,0.1538462,0.8993289,0.3321678,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->AddEntry(Eff_MET_HighPU,"nVertices>15","l");
  leg1->AddEntry(Eff_MET_LowPU,"nVertices<15","l");
  leg1->Draw();

  c1->SaveAs((infile+"_Efficiency_Trigger_MET_PU.pdf").c_str());
  c1->SaveAs((infile+"_Efficiency_Trigger_MET_PU.png").c_str());

}

void Fitting_Trigger_PhotonID(std::string infile){

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);

  gStyle->SetErrorX(0);

  TFile* file = TFile::Open((infile+".root").c_str());
  TH1F *h_MET_HLTPho = (TH1F*) file->Get("h_PHO_HLTPho");
  h_MET_HLTPho->Rebin(3);
  TH1F *h_MET_HLTPhoId = (TH1F*) file->Get("h_PHO_HLTPhoId");
  h_MET_HLTPhoId->Rebin(3);
  TGraphAsymmErrors *Eff_PHO = new TGraphAsymmErrors;
  Eff_PHO->BayesDivide(h_PHO_HLTPhoId, h_PHO_HLTPho, "");
  Eff_PHO->SetLineColor(kBlue);
  Eff_PHO->GetXaxis()->SetRangeUser(0, 800);
  Eff_PHO->GetYaxis()->SetTitle("Efficiency");
  Eff_PHO->GetXaxis()->SetTitle("Offline Photon pT [GeV]");
  Eff_PHO->SetMinimum(0.0);
  Eff_PHO->SetTitle("Efficiency: HLT_Photon30_R9Id90 w.r.t HLT_Photon30");

  Eff_PHO->Draw("AP E");

  TF1 *fit = new TF1("fit","0.5*[2]*(1+TMath::Erf((x-[0])/[1]))",0.0, 800.0);
  fit->SetParameter(0, 20);
  fit->SetParameter(1, 40);
  fit->SetParameter(2, 0.8);
  fit->SetParLimits(2,0.0,1.0);
  Eff_PHO->Fit("fit", "", "", 0,800.0);

  fit->SetLineColor(kRed);
  fit->SetLineWidth(3);

  c1->SaveAs((infile+"_Efficiency_Trigger_PHO.pdf").c_str());
  c1->SaveAs((infile+"_Efficiency_Trigger_PHO.png").c_str());

}
