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

void makeCosPlot(){

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);
  TH1F *h_cosFunc = new TH1F("h_cosFunc", "h_cosFunc", 10, 0.0, 1.0);

  double array[32] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 53, 57, 59, 61, 63, 67, 69, 71};
  double array_den[32] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 53, 57, 59, 61, 63, 67, 69, 71, 73}; 
  
  for (int i=0; i<33; i++){
    //h_cosFunc->Fill(TMath::Cos(TMath::Sqrt(array[i])/TMath::Sqrt(array_den[i])));
    h_cosFunc->Fill(TMath::Cos(array[i]/array_den[i]));
  }

  h_cosFunc->Draw("");

}
