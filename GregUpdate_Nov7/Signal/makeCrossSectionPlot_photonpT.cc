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

void makeCrossSectionPlot_photonpT(std::string outfile){

  TCanvas *c1 = new TCanvas("Cross section","The parametrized cross section",200,10,700,500);
  c1->SetLogy();
  //double cross_section[6] = {0.6773*1000, 0.2738*1000, 0.169*1000, 0.1264*1000, 0.06471*1000, 0.03566*1000};
  //double cross_section[6] = {335.9, 114.4, 65.41, 46.23, 20.95, 10.39};   
  double cross_section[6] = {1.032*1000, 0.688*1000, 0.5029*1000, 0.4*1000, 0.3171*1000, 0.268*1000};
  double phpT[6] = {5, 10, 15, 20, 25, 30};

  double x_error[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double y_error[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  gStyle->SetErrorX(0.5);

  TGraphErrors *Xsec_phpT = new TGraphErrors(6, phpT, cross_section, x_error, y_error);
  Xsec_phpT->SetTitle("");
  Xsec_phpT->SetLineColor(kRed);
  Xsec_phpT->SetMarkerColor(kRed);
  Xsec_phpT->SetMarkerStyle(20);
  Xsec_phpT->SetMarkerSize(1.0);
  Xsec_phpT->SetLineWidth(2);
  Xsec_phpT->SetLineStyle(1);
  Xsec_phpT->SetTitle("Photon pT versus cross section");
  Xsec_phpT->GetXaxis()->SetTitle("Photon pT [GeV]");
  Xsec_phpT->GetYaxis()->SetTitle("Cross section [pb]");
  Xsec_phpT->GetYaxis()->SetLabelSize(0.035);
  Xsec_phpT->GetXaxis()->SetLabelSize(0.035);
  Xsec_phpT->GetYaxis()->SetTitleSize(0.035);
  Xsec_phpT->GetXaxis()->SetTitleSize(0.035);
  Xsec_phpT->GetXaxis()->SetLabelFont(62);
  Xsec_phpT->GetYaxis()->SetLabelFont(62);
  Xsec_phpT->GetXaxis()->SetTitleFont(62);
  Xsec_phpT->GetYaxis()->SetTitleFont(62);
  Xsec_phpT->GetYaxis()->SetTitleOffset(1.4);
  Xsec_phpT->Draw("ALP");
  Xsec_phpT->SetName("Xsec_phpT");

  TLegend *leg1 = new TLegend(0.45,0.55,0.65,0.65,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0); 
  leg1->AddEntry("Xsec_phpT","Stop (m=200 GeV) Production with a photon", "lp"); 
  leg1->Draw();
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  Xsec_phpT->Write("Xsec_phpT");
  tFile->Close();

}
