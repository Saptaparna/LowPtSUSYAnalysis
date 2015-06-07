#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


int JecUnc()
{
  //if (jetPt<10.0 || abs(jetEta)>5.4) return 1.0;

  //double number;
  double jetEtalow;
  double jetEtahigh;
  double jetPtlow;
  double jetPthigh;
  int jetPt_index;
  double up;
  double down;

  std::ifstream file("Summer13_V5_DATA_UncertaintySources_AK5PF.txt");
  
  std::string s;
  getline(file,s);

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
   jetEtalow = double(line.at(0));
   jetEtahigh = double(line.at(1));
   jetPtlow = double(line.at(3));
   up = double(line.at(4));
   down = double(line.at(5));
   jetPthigh = double(line.at(6));
   cout << "variables = " << line.at(0) << " , " << line.at(1) << " , " << line.at(3)  << " , " << line.at(4) << " , " <<  line.at(5) << " , " << line.at(6) << endl; 
  }

}
