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

using std::ifstream;
using std::string;
using std::cout;
using std::endl;
using std::istringstream;


void makeScatterPlot_Phi(std::string infile){


gROOT->SetStyle("Plain");

TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);

//c1->SetLogx();
//c1->SetLogy();

TFile* file1 = TFile::Open((infile+".root").c_str());
TH2F *h1 = (TH2F*)gDirectory->Get("h_ph_el1_Phi");
h1->SetMarkerStyle(kStar);
h1->SetMarkerSize(0.5);
h1->SetMarkerColor(kBlue);
h1->SetStats(kFALSE);
h1->GetXaxis()->SetLabelSize(0.03);
h1->GetYaxis()->SetLabelSize(0.03);
h1->GetYaxis()->SetTitleOffset(1.2);
h1->SetLineColor(kBlue);
h1->Draw();

TH2F *h2 = (TH2F*)gDirectory->Get("h_ph_el2_Phi");
h2->SetMarkerStyle(kCircle);
h2->SetMarkerSize(0.5);
h2->SetMarkerColor(kRed);
h2->Draw("same");

TLegend *leg1 = new TLegend(0.7080537,0.8181818,0.9832215,0.9965035,NULL,"brNDC");
leg1->SetBorderSize(0);
leg1->SetTextSize(0.03);
leg1->SetLineColor(1);
leg1->SetLineStyle(0);
leg1->SetLineWidth(1);
leg1->SetFillColor(10);
leg1->SetFillStyle(0);
leg1->AddEntry(h1,"Leading Electron","p");
leg1->AddEntry(h2,"Trailing Electron","p");
leg1->Draw();

c1->SaveAs((infile+"_ph_el_Phi.pdf").c_str());
c1->SaveAs((infile+"_ph_el_Phi.png").c_str());


   }
