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
#include <algorithm>
#include <TGraphAsymmErrors.h>
#include <TVector2.h>
#include <TF1.h>

int calculatedWeightedYields(std::string channel){

double ZGammaLL_SF = (1.2*132.6*7.3*1000)/6583032;
double ZGammaNU_SF = (1.2*123.9*7.3*1000)/3169096;
double WGamma_SF = (1.2*461.6*7.3*1000)/4802339;
double TTG_SF = (1.2*2.166*7.3*1000)/71328;
double DY_SF = (1.2*3532.8*7.3*1000)/27457743;
double TT_SF = (13.43*2.0*7.3*1000)/11843493;
double W_SF = (37509.2*1.2*7.3*1000)/18393019;
double WWG_SF = (0.528*1.2*7.3*1000)/214538;
double WGSE_SF = (5.873*1.2*7.3*1000)/314539;
double WGSMu_SF = (1.914*1.2*7.3*1000)/302312;
double WGSTau_SF = (0.336*1.2*7.3*1000)/49981;
double ZZ_SF = (0.1769*1.2*7.3*1000)/18393019;

if(channel=="ElMu")
{
  double yield_ZGammaLL = ZGammaLL_SF*70.5446;
  cout << "yield_ZGammaLL = " << yield_ZGammaLL << endl;
  double yield_WGamma = WGamma_SF*8.12154;
  cout << "yield_WGamma = " << yield_WGamma << endl;
  double yield_DY = DY_SF*91.5287;
  cout << "yield_DY = " << yield_DY << endl;
  double yield_TTbar = TT_SF*838.707;
  cout << "yield_TTbar = " << yield_TTbar << endl;
  double yield_TTG = TTG_SF*59.5822;
  cout << "yield_TTG = " << yield_TTG << endl;
  double yield_WWG = WWG_SF*97.0018;
  cout << "yield_WWG = " << yield_WWG << endl;
  double yield_WGSE = WGSE_SF*13.0749;
  cout << "yield_WGSE = " << yield_WGSE << endl;
  double yield_WGSMu = WGSMu_SF*1.50912;
  cout << "yield_WGSMu = " << yield_WGSMu << endl;
  double yield_WGSTau = WGSTau_SF*0.0;
  cout << "yield_WGSTau = " << yield_WGSTau << endl;  
  double yield_W = W_SF*1.40561;
  cout << "yield_W = " << yield_W << endl;
  double yield_ZZ = ZZ_SF*8691.91;
  cout << "yield_ZZ = " << yield_ZZ << endl;

}

if(channel=="MuMu")
{
  double yield_ZGammaLL = ZGammaLL_SF*4932.9;
  cout << "yield_ZGammaLL = " << yield_ZGammaLL << endl;
  double yield_WGamma = WGamma_SF*1.02557;
  cout << "yield_WGamma = " << yield_WGamma << endl;
  double yield_DY = DY_SF*897.51;
  cout << "yield_DY = " << yield_DY << endl;
  double yield_TTbar = TT_SF*1711.13;
  cout << "yield_TTbar = " << yield_TTbar << endl;
  double yield_TTG = TTG_SF*64.1012;
  cout << "yield_TTG = " << yield_TTG << endl;
  double yield_WWG = WWG_SF*97.41;
  cout << "yield_WWG = " << yield_WWG << endl;
  double yield_WGSE = WGSE_SF*0.0;
  cout << "yield_WGSE = " << yield_WGSE << endl;
  double yield_WGSMu = WGSMu_SF*211.018;
  cout << "yield_WGSMu = " << yield_WGSMu << endl;
  double yield_WGSTau = WGSTau_SF*2.37593;
  cout << "yield_WGSTau = " << yield_WGSTau << endl;
  double yield_W = W_SF*1.46941;
  cout << "yield_W = " << yield_W << endl;
  double yield_ZZ = ZZ_SF*37568.2;
  cout << "yield_ZZ = " << yield_ZZ << endl;

}




return 0;

} 
