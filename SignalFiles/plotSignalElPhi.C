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


void plotSignalElPhi(){

  gROOT->SetStyle("Plain");

  TCanvas c1("c1","Signal Leading Electron pT",10,10,600,400);

  c1.SetLogy();

  TFile* file1 = TFile::Open("Output_LowPtSUSY_Tree_M1200_M2200_Delphes.root");
  TH1F *h1 = (TH1F*)file1->Get("h_el_phi_leading");
  if(h1->Integral()>0.0) h1->Scale(1/(h1->Integral()));
  h1->Rebin(2);
  //h1->GetXaxis()->SetRangeUser(0, 100);
  h1->GetYaxis()->SetTitleOffset(1.3);
  h1->SetLineColor(kBlack);
  h1->SetLineWidth(2);
  h1->SetMaximum(0.05);
  h1->SetStats(kFALSE);
  h1->Draw("HIST");

  TFile* file2 = TFile::Open("Output_LowPtSUSY_Tree_M1250_M2250_Delphes.root");
  TH1F *h2 = (TH1F*)file2->Get("h_el_phi_leading");
  if(h2->Integral()>0.0) h2->Scale(1/(h2->Integral()));
  h2->Rebin(2);
  h2->SetLineColor(kRed);
  h2->SetLineWidth(2);
  h2->Draw("hist same");

  TFile* file3 = TFile::Open("Output_LowPtSUSY_Tree_M1300_M2300_Delphes.root");
  TH1F *h3 = (TH1F*)file3->Get("h_el_phi_leading");
  if(h3->Integral()>0.0) h3->Scale(1/(h3->Integral()));
  h3->Rebin(2);
  h3->SetLineColor(kBlue);
  h3->SetLineWidth(2);
  h3->Draw("hist same");

  TFile* file4 = TFile::Open("Output_LowPtSUSY_Tree_M1350_M2350_Delphes.root");
  TH1F *h4 = (TH1F*)file4->Get("h_el_phi_leading");
  if(h4->Integral()>0.0) h4->Scale(1/(h4->Integral()));
  h4->Rebin(2);
  h4->SetLineColor(kGreen);
  h4->SetLineWidth(2);
  h4->Draw("hist same");

  TFile* file5 = TFile::Open("Output_LowPtSUSY_Tree_M1400_M2400_Delphes.root");
  TH1F *h5 = (TH1F*)file5->Get("h_el_phi_leading");
  if(h5->Integral()>0.0) h5->Scale(1/(h5->Integral()));
  h5->Rebin(2);
  h5->SetLineColor(kOrange+4);
  h5->SetLineWidth(2);
  h5->Draw("hist same");

  TFile* file6 = TFile::Open("Output_LowPtSUSY_Tree_M1450_M2450_Delphes.root");
  TH1F *h6 = (TH1F*)file6->Get("h_el_phi_leading");
  if(h6->Integral()>0.0) h6->Scale(1/(h6->Integral()));
  h6->Rebin(2);
  h6->SetLineColor(kCyan);
  h6->SetLineWidth(2);
  h6->Draw("hist same");

  TLegend *leg1 = new TLegend(0.5907383,0.6609651,0.8959732,0.8579088,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->AddEntry(h1,"M1_200 M2_200","l");
  leg1->AddEntry(h2,"M1_250 M2_250","l");
  leg1->AddEntry(h3,"M1_300 M2_300","l");
  leg1->AddEntry(h4,"M1_350 M2_350","l");
  leg1->AddEntry(h5,"M1_400 M2_400","l");
  leg1->AddEntry(h6,"M1_450 M2_450","l");
  leg1->Draw();
  
  c1.SaveAs("ElectronPhi_Delphes.pdf");
  c1.SaveAs("ElectronPhi_Delphes.png");

}
