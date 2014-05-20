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
  his->Rebin(10);
  his->GetXaxis()->SetRangeUser(70, 110); 
  his->GetYaxis()->SetTitle("Events/0.1 GeV");
 
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
  his->GetYaxis()->SetTitle("Events/0.05 GeV");

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

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);
  c1->SetGridx();
  c1->SetGridy();
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderMode(-1);
  c1->GetFrame()->SetBorderSize(5);

  TFile* file = TFile::Open((infile+".root").c_str());
  TH1F *his = (TH1F*) file->Get("h_InvariantMass_Mu");
  his->GetXaxis()->SetRangeUser(8.0,10.0);
  his->GetYaxis()->SetTitle("Events/0.01 GeV");

  TF1 *fSignal = new TF1("fSignal","gaus",8.5,9.5);
  TF1 *fBackground = new TF1("fBackground","pol1", 8.0,10.0);
  fBackground->SetLineColor(kBlue);
  fBackground->SetLineWidth(2);


  TF1 *fSpectrum = new TF1("fSpectrum","gaus+pol1(3)",8.0,10.0);

  fSpectrum->SetParameters(2.0, his->GetMean(), his->GetRMS(), 10.0, 1.0);
  his->Fit("fSpectrum", "LL", "", 8.0,10.0); //log likelihood fit
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
  cout << "fBackground->Integral(8.0,10.0) = " << fBackground->Integral(8.0,10.0) << endl;
  cout << "fSignal->Integral(8.5, 9.5) = " << fSignal->Integral(8.5, 9.5) << endl;
  cout << "fSpectrum->Integral(8.0, 10.0) = " << fSpectrum->Integral(8.0, 10.0) << endl;


  TF1 *fitResult = his->GetFunction("fSpectrum");
  cout << "fitResult->GetChisquare() = " << fitResult->GetChisquare() << endl;
  cout << "fitResult->GetNDF() = " << fitResult->GetNDF() << endl;

  cout << "fitResult->GetChisquare() = " << fitResult->GetChisquare() << endl;
  cout << "fitResult->GetNDF() = " << fitResult->GetNDF() << endl;

  cout << "Sanity Check = " << fBackground->Integral(8.0, 10.0) + fSignal->Integral(8.5, 9.5) << " = Total = " << fSpectrum->Integral(8.0, 10.0) << endl;

  cout << "Fitted result = " << fSpectrum->Integral(8.0, 10.0)*10 << " Data = " << his->Integral() << endl;


  c1->SaveAs((infile+"_Upsilon_LL.pdf").c_str());
  c1->SaveAs((infile+"_Upsilon_LL.png").c_str());

}


void Fitting_Trigger_MET(std::string infile){

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);

  gStyle->SetErrorX(0);

  TFile* file = TFile::Open((infile+".root").c_str());
  TGraphAsymmErrors *gr1 = (TGraphAsymmErrors*) file->Get("Eff_MET");
  gr1->GetXaxis()->SetRangeUser(0, 800); 
  gr1->GetYaxis()->SetTitle("Efficiency");
  gr1->GetXaxis()->SetTitle("Offline MET");
  gr1->SetTitle("Efficiency: HLT_Photon30_R9Id90_Met25 w.r.t HLT_Photon30_R9Id90");
  
  gr1->Draw("AP");

  c1->SaveAs((infile+"_Efficiency_Trigger_MET.pdf").c_str());
  c1->SaveAs((infile+"_Efficiency_Trigger_MET.png").c_str());

}
            
