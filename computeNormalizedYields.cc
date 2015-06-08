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
#include <TFractionFitter.h>


void computeNormalizedYields_HT_Mll()
{
  cout << "Scenario:  HT = 130 and  inv mass = 140 and HT_InvMass_Yield ElMu" << endl;
  //backgrounds
  double ZG_raw =  81.23; 
  double TTG_raw = 24.3962;
  double ZZ_raw = 362.853;
  double T_tW_raw = 0.38512;
  double Tbar_tW_raw = 2.46058;
  double WG_raw = 2.67064;
  double WZ_raw = 166.754;
  double WWG_raw = 99.0616; 
  double SS = 26;
  //signal
  double Stop80_raw = 1326.28;
  double Stop100_raw = 1570.34;
  double Stop120_raw = 1655.67;
  double Stop140_raw = 1791.06;
  double Stop160_raw = 1856.45;
  double Stop180_raw = 1831.15;
  double Stop200_raw = 1904.14;

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
  double signal_stop_SF_80 = (1.66*0.9301*1000*7.4*0.33*0.33*2)/99990;
  double signal_stop_SF_100 = (1.66*0.3983*1000*7.4*0.33*0.33*2)/99994;
  double signal_stop_SF_120 = (1.66*0.1943*1000*7.4*0.33*0.33*2)/99993;
  double signal_stop_SF_140 = (1.66*0.104*1000*7.4*0.33*0.33*2)/99996;
  double signal_stop_SF_160 = (1.66*0.05941*1000*7.4*0.33*0.33*2)/99990;
  double signal_stop_SF_180 = (1.66*0.03592*1000*7.4*0.33*0.33*2)/99987;
  double signal_stop_SF_200 = (1.66*0.02264*1000*7.4*0.33*0.33*2)/99991;

  double signal_SF = (0.0004318*(7.4)*1000*1000*1.66)/99990;

  cout << "Normalized ZGamma = " << ZG_raw*ZGammaLL_SF << endl;
  cout << "Normalized TTG = " << TTG_raw*TTG_SF << endl; 
  cout << "Normalized ZZ = " << ZZ_raw*ZZ_SF << endl;
  cout << "Normalized single top = " << (T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF) << endl;
  cout << "Normalized WGstar = " << WG_raw*WGSMu_SF << endl;
  cout << "Normalized WZ = " << WZ_raw*WZ_SF << endl;
  cout << "Normalized WWG = " << WWG_raw*WWG_SF << endl;
  cout << "Data driven = " << ZG_raw*ZGammaLL_SF*(116.688/2391.76) << endl;
  cout << "SS = " << SS << endl;

  cout << "Signal" << endl;
  cout << "Normalized Stop80 = " << Stop80_raw*signal_stop_SF_80 << " Significance = " << (Stop80_raw*signal_stop_SF_80)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop100 = " << Stop80_raw*signal_stop_SF_100 << " Significance = " << (Stop100_raw*signal_stop_SF_100)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop120 = " << Stop80_raw*signal_stop_SF_120  << " Significance = " << (Stop120_raw*signal_stop_SF_120)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop140 = " << Stop80_raw*signal_stop_SF_140  << " Significance = " << (Stop140_raw*signal_stop_SF_140)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop160 = " << Stop80_raw*signal_stop_SF_160  << " Significance = " << (Stop160_raw*signal_stop_SF_160)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop180 = " << Stop80_raw*signal_stop_SF_180  << " Significance = " << (Stop180_raw*signal_stop_SF_180)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop200 = " << Stop80_raw*signal_stop_SF_200  << " Significance = " << (Stop200_raw*signal_stop_SF_200)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Total background = " << (ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;
}

void computeNormalizedYields_HTb_Mll()
{
  cout << "Scenario:  HTb = 90 and  inv mass = 140 and HTb_InvMass_Yield ElMu" << endl;
  //backgrounds
  double ZG_raw =  87.6871;
  double TTG_raw = 43.284;
  double ZZ_raw = 412.322;
  double T_tW_raw = 2.09491;
  double Tbar_tW_raw = 2.46058;
  double WG_raw = 2.67064;
  double WZ_raw = 196.623;
  double WWG_raw = 116.694;
  double SS = 26;
  //signal
  double Stop80_raw = 1489.41;
  double Stop100_raw = 1757.61;
  double Stop120_raw = 1905.34;
  double Stop140_raw = 2077.74;
  double Stop160_raw = 2198.99;
  double Stop180_raw = 2194.35;
  double Stop200_raw = 2303.27;

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
  double signal_stop_SF_80 = (1.66*0.9301*1000*7.4*0.33*0.33*2)/99990;
  
  double signal_stop_SF_100 = (1.66*0.3983*1000*7.4*0.33*0.33*2)/99994;
  double signal_stop_SF_120 = (1.66*0.1943*1000*7.4*0.33*0.33*2)/99993;
  double signal_stop_SF_140 = (1.66*0.104*1000*7.4*0.33*0.33*2)/99996;
  double signal_stop_SF_160 = (1.66*0.05941*1000*7.4*0.33*0.33*2)/99990;
  double signal_stop_SF_180 = (1.66*0.03592*1000*7.4*0.33*0.33*2)/99987;
  double signal_stop_SF_200 = (1.66*0.02264*1000*7.4*0.33*0.33*2)/99991;

  double signal_SF = (0.0004318*(7.4)*1000*1000*1.66)/99990;

  cout << "Normalized ZGamma = " << ZG_raw*ZGammaLL_SF << endl;
  cout << "Normalized TTG = " << TTG_raw*TTG_SF << endl;
  cout << "Normalized ZZ = " << ZZ_raw*ZZ_SF << endl;
  cout << "Normalized single top = " << (T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF) << endl;
  cout << "Normalized WGstar = " << WG_raw*WGSMu_SF << endl;
  cout << "Normalized WZ = " << WZ_raw*WZ_SF << endl;
  cout << "Normalized WWG = " << WWG_raw*WWG_SF << endl;
  cout << "Data driven = " << ZG_raw*ZGammaLL_SF*(116.688/2391.76) << endl;
  cout << "SS = " << SS << endl;

  cout << "Signal" << endl;
  cout << "Normalized Stop80 = " << Stop80_raw*signal_stop_SF_80 << " Significance = " << (Stop80_raw*signal_stop_SF_80)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop100 = " << Stop80_raw*signal_stop_SF_100 << " Significance = " << (Stop100_raw*signal_stop_SF_100)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop120 = " << Stop80_raw*signal_stop_SF_120  << " Significance = " << (Stop120_raw*signal_stop_SF_120)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop140 = " << Stop80_raw*signal_stop_SF_140  << " Significance = " << (Stop140_raw*signal_stop_SF_140)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop160 = " << Stop80_raw*signal_stop_SF_160  << " Significance = " << (Stop160_raw*signal_stop_SF_160)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop180 = " << Stop80_raw*signal_stop_SF_180  << " Significance = " << (Stop180_raw*signal_stop_SF_180)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop200 = " << Stop80_raw*signal_stop_SF_200  << " Significance = " << (Stop200_raw*signal_stop_SF_200)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

 cout << "Total background = " << (ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

}

void percentUnc(){

  double btagDown = 43.4;
  double btagUp = 42.4782;
  double jecDown = 44.1676;
  double jecUp = 42.4838;
  double jerUp = 43.4;
  double jerDown = 43.4;

  double def = 43.284;
  
  double perbtagDown = (43.4-43.284)/43.284;
  double perbtagUp = (43.284-42.4782)/43.284;  
  double perjecUp = (44.1676-43.284)/43.284;
  double perjecDown = (43.284-42.4838)/43.284;
  double perjerUp = (43.4-43.284)/43.284;
  double perjerDown = (43.4-43.284)/43.284;

  cout << "perbtagDown = " << perbtagDown << endl;
  cout << "perbtagUp = " << perbtagUp << endl;
  cout << "perjecUp = " << perjecUp << endl;
  cout << "perjecDown = " << perjecDown << endl;
  cout << "perjerUp = " << perjerUp << endl;
  cout << "perjerDown = " << perjerDown << endl;
}

void computeNormalizedYields_HTb_Mll_WithoutZero()
{
  cout << "Scenario:  HTb = 105 and  inv mass = 130 and HT_InvMass_Yield ElMu" << endl;
  //backgrounds
  double ZG_raw =  0.0;
  double TTG_raw = 25.3945;
  double ZZ_raw = 13.2462;
  double T_tW_raw = 2.09491;
  double Tbar_tW_raw = 0.0;
  double WG_raw = 0.0;
  double WZ_raw = 1.97341;
  double WWG_raw = 0.780696;
  double SS = 25;
  //signal
  double Stop80_raw = 641.858;
  double Stop100_raw = 720.405;
  double Stop120_raw = 750.067;
  double Stop140_raw = 847.525;
  double Stop160_raw = 872.202;
  double Stop180_raw = 906.301;
  double Stop200_raw = 906.848;

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
  double signal_stop_SF_80 = (1.66*0.9301*1000*7.4*0.33*0.33*2)/99990;
  double signal_stop_SF_100 = (1.66*0.3983*1000*7.4*0.33*0.33*2)/99994;
  double signal_stop_SF_120 = (1.66*0.1943*1000*7.4*0.33*0.33*2)/99993;
  double signal_stop_SF_140 = (1.66*0.104*1000*7.4*0.33*0.33*2)/99996;
  double signal_stop_SF_160 = (1.66*0.05941*1000*7.4*0.33*0.33*2)/99990;
  double signal_stop_SF_180 = (1.66*0.03592*1000*7.4*0.33*0.33*2)/99987;
  double signal_stop_SF_200 = (1.66*0.02264*1000*7.4*0.33*0.33*2)/99991;

  double signal_SF = (0.0004318*(7.4)*1000*1000*1.66)/99990;

  cout << "Normalized ZGamma = " << ZG_raw*ZGammaLL_SF << endl;
  cout << "Normalized TTG = " << TTG_raw*TTG_SF << endl;
  cout << "Normalized ZZ = " << ZZ_raw*ZZ_SF << endl;
  cout << "Normalized single top = " << (T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF) << endl;
  cout << "Normalized WGstar = " << WG_raw*WGSMu_SF << endl;
  cout << "Normalized WZ = " << WZ_raw*WZ_SF << endl;
  cout << "Normalized WWG = " << WWG_raw*WWG_SF << endl;
  cout << "Data driven = " << ZG_raw*ZGammaLL_SF*(116.688/2391.76) << endl;
  cout << "SS = " << SS << endl;

  cout << "Signal" << endl;
  cout << "Normalized Stop80 = " << Stop80_raw*signal_stop_SF_80 << " Significance = " << (Stop80_raw*signal_stop_SF_80)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop100 = " << Stop80_raw*signal_stop_SF_100 << " Significance = " << (Stop100_raw*signal_stop_SF_100)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop120 = " << Stop80_raw*signal_stop_SF_120  << " Significance = " << (Stop120_raw*signal_stop_SF_120)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop140 = " << Stop80_raw*signal_stop_SF_140  << " Significance = " << (Stop140_raw*signal_stop_SF_140)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop160 = " << Stop80_raw*signal_stop_SF_160  << " Significance = " << (Stop160_raw*signal_stop_SF_160)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop180 = " << Stop80_raw*signal_stop_SF_180  << " Significance = " << (Stop180_raw*signal_stop_SF_180)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop200 = " << Stop80_raw*signal_stop_SF_200  << " Significance = " << (Stop200_raw*signal_stop_SF_200)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;
}

void computeNormalizedYields_HTb_Zero()
{
  cout << "Scenario:  HTb = 0.0 and  inv mass = 130 and HT_InvMass_Yield ElMu" << endl;
  //backgrounds
  double ZG_raw =  90.1995;
  double TTG_raw = 30.8962;
  double ZZ_raw = 464.72;
  double T_tW_raw = 0.0;
  double Tbar_tW_raw = 2.46058;
  double WG_raw = 2.67064;
  double WZ_raw = 242.201;
  double WWG_raw = 160.147;
  double SS = 13;
  //signal
  double Stop80_raw = 887.745;
  double Stop100_raw = 1135.18;
  double Stop120_raw =  1267.29;
  double Stop140_raw = 1332.03;
  double Stop160_raw = 1430.22;
  double Stop180_raw = 1380.06;
  double Stop200_raw = 1495.52;
          
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
  double signal_stop_SF_80 = (1.66*0.9301*1000*7.4*0.33*0.33*2)/99990;
  double signal_stop_SF_100 = (1.66*0.3983*1000*7.4*0.33*0.33*2)/99994;
  double signal_stop_SF_120 = (1.66*0.1943*1000*7.4*0.33*0.33*2)/99993;
  double signal_stop_SF_140 = (1.66*0.104*1000*7.4*0.33*0.33*2)/99996;
  double signal_stop_SF_160 = (1.66*0.05941*1000*7.4*0.33*0.33*2)/99990;
  double signal_stop_SF_180 = (1.66*0.03592*1000*7.4*0.33*0.33*2)/99987;
  double signal_stop_SF_200 = (1.66*0.02264*1000*7.4*0.33*0.33*2)/99991;

  double signal_SF = (0.0004318*(7.4)*1000*1000*1.66)/99990;

  cout << "Normalized ZGamma = " << ZG_raw*ZGammaLL_SF << endl;
  cout << "Normalized TTG = " << TTG_raw*TTG_SF << endl;
  cout << "Normalized ZZ = " << ZZ_raw*ZZ_SF << endl;
  cout << "Normalized single top = " << (T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF) << endl;
  cout << "Normalized WGstar = " << WG_raw*WGSMu_SF << endl;
  cout << "Normalized WZ = " << WZ_raw*WZ_SF << endl;
  cout << "Normalized WWG = " << WWG_raw*WWG_SF << endl;
  cout << "Data driven = " << ZG_raw*ZGammaLL_SF*(116.688/2391.76) << endl;
  cout << "SS = " << SS << endl;

  cout << "Signal" << endl;
  cout << "Normalized Stop80 = " << Stop80_raw*signal_stop_SF_80 << " Significance = " << (Stop80_raw*signal_stop_SF_80)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop100 = " << Stop80_raw*signal_stop_SF_100 << " Significance = " << (Stop100_raw*signal_stop_SF_100)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop120 = " << Stop80_raw*signal_stop_SF_120  << " Significance = " << (Stop120_raw*signal_stop_SF_120)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop140 = " << Stop80_raw*signal_stop_SF_140  << " Significance = " << (Stop140_raw*signal_stop_SF_140)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop160 = " << Stop80_raw*signal_stop_SF_160  << " Significance = " << (Stop160_raw*signal_stop_SF_160)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop180 = " << Stop80_raw*signal_stop_SF_180  << " Significance = " << (Stop180_raw*signal_stop_SF_180)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;

  cout << "Normalized Stop200 = " << Stop80_raw*signal_stop_SF_200  << " Significance = " << (Stop200_raw*signal_stop_SF_200)/TMath::Sqrt(ZG_raw*ZGammaLL_SF+TTG_raw*TTG_SF+ZZ_raw*ZZ_SF+(T_tW_raw*T_tW_SF + Tbar_tW_raw*Tbar_tW_SF)+ WG_raw*WGSMu_SF + WZ_raw*WZ_SF +  WWG_raw*WWG_SF + ZG_raw*ZGammaLL_SF*(116.688/2391.76) + SS) << endl;
}
 
