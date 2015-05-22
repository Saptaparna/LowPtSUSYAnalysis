#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


bool JecUnc(double jetPt, double jetEta, double &jecUp, double &jecDown)
{
  if (jetPt<10.0 || abs(jetEta)>5.4) return false;

  float jetEtalow;
  float jetEtahigh;
  std::ifstream file("Summer13_V5_DATA_UncertaintySources_AK5PF.txt");
  
  std::string s;

  while (!file.eof())
  {
    if (!getline(file, s)) break;

    std::stringstream ss(s);
    std::vector<std::string> line;
    while (ss)
    {
      std::string number;
      if (!getline(ss, number, ' ')) break; 
      line.push_back(number);
    }
    jetEtalow = atof(line.at(0).c_str());
    jetEtahigh = atof(line.at(1).c_str());
    if(jetEta>=jetEtalow && jetEta<jetEtahigh)
    {  
      cout << "jetEta = " << jetEtalow << " , " << jetEtahigh << jetPt << endl; 
      for(int i=0; i<43; i++)
      {
        if(jetPt>=atof(line.at(3*i+3).c_str()) && jetPt<=atof(line.at(3*i+6).c_str()))
        {
          jecUp = atof(line.at(3*i+4).c_str());
          jecDown = atof(line.at(3*i+5).c_str());
          return true;
        }
      }
      jecUp = atof(line.at(133).c_str()); 
      jecDown = atof(line.at(134).c_str()); 
      return true;
    }
  }
}

void JecUnc_Sapta()
{
  double thisJetpT = 4000;
  double thisJeteta = 1.3;
  double jecUp = 0;
  double jecDown = 0;

  if(JecUnc(thisJetpT, thisJeteta, jecUp, jecDown)){
   
    cout << "jecUp = " << jecUp << endl;
    cout << "jecDown = " << jecDown << endl;

  } 

}
