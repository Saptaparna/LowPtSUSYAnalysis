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

void makeCrossSectionPlot_Nophoton(std::string outfile){

  TCanvas *c1 = new TCanvas("Cross section","The parametrized cross section",200,10,700,500);
  c1->SetLogy();
  //double cross_section[6] = {0.6773*1000, 0.2738*1000, 0.169*1000, 0.1264*1000, 0.06471*1000, 0.03566*1000};
  //double cross_section[6] = {335.9, 114.4, 65.41, 46.23, 20.95, 10.39};   
  double cross_section[12] = {335.9, 114.4, 65.41, 46.23, 20.95, 10.39, 1.086, 0.1908, 0.04378, 0.0121, 0.003805, 0.001284};
  double mass[12] = {100, 125, 140, 150, 175, 200, 300, 400, 500, 600, 700, 800};

  double x_error[12] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double y_error[12] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  gStyle->SetErrorX(0.5);

  TGraphErrors *Xsec_mass = new TGraphErrors(12, mass, cross_section, x_error, y_error);
  Xsec_mass->SetTitle("");
  Xsec_mass->SetLineColor(kRed);
  Xsec_mass->SetMarkerColor(kRed);
  Xsec_mass->SetMarkerStyle(20);
  Xsec_mass->SetMarkerSize(1.0);
  Xsec_mass->SetLineWidth(2);
  Xsec_mass->SetLineStyle(1);
  Xsec_mass->SetTitle("Mass versus cross section");
  Xsec_mass->GetXaxis()->SetTitle("Mass [GeV]");
  Xsec_mass->GetYaxis()->SetTitle("Cross section [pb]");
  Xsec_mass->GetYaxis()->SetLabelSize(0.035);
  Xsec_mass->GetXaxis()->SetLabelSize(0.035);
  Xsec_mass->GetYaxis()->SetTitleSize(0.035);
  Xsec_mass->GetXaxis()->SetTitleSize(0.035);
  Xsec_mass->GetXaxis()->SetLabelFont(62);
  Xsec_mass->GetYaxis()->SetLabelFont(62);
  Xsec_mass->GetXaxis()->SetTitleFont(62);
  Xsec_mass->GetYaxis()->SetTitleFont(62);
  Xsec_mass->GetYaxis()->SetTitleOffset(1.4);
  Xsec_mass->Draw("ALP");
  Xsec_mass->SetName("Xsec_mass");

  TLegend *leg1 = new TLegend(0.45,0.55,0.65,0.65,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0); 
  leg1->AddEntry("Xsec_mass","Stop Production without a photon", "lp"); 
  leg1->Draw();
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  Xsec_mass->Write("Xsec_mass");
  tFile->Close();

}
