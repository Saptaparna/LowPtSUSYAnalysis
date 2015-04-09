#!/bin/csh                                                                                                                                                                         
cp FlatTreeCreator_Template.C FlatTreeCreator.C
sed -i "s/SUFFIX/$1/g" FlatTreeCreator.C

cat > Executable.C << +EOF

void Executable(){

gSystem->Load("libFWCoreFWLite.so");
AutoLibraryLoader::enable();
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCPhysObject_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCMuon_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCMET_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCEGamma_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCElectron_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCTau_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCPhoton_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCJet_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCGenParticle_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCGenJet_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCPrimaryVtx_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCTriggerObject_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/plugins/HistManager_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/plugins/TriggerSelector_cc.so");
gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/lib/slc5_amd64_gcc462/libTauAnalysisSVfitStandalone.so");
gROOT->ProcessLine(".L FlatTreeCreator.C++");

}


+EOF

root -l -b -q Executable.C
rm -f Executable.C
rm -f FlatTreeCreator.C
mv FlatTreeCreator_C.so FlatTreeCreator_${1}_C.so

echo "FINISHED COMPILING"

cat > run_${1}.C <<EOF

    #include <iostream>
    #include <fstream>
    #include <string> 
    #include <vector> 
    #include <cstdlib>
    #include <stdio>
                      
    using namespace std;
        
    void run_${1}() {

    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCPhysObject_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCMuon_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCMET_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCEGamma_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCElectron_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCTau_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCPhoton_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCJet_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCGenParticle_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCGenJet_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCPrimaryVtx_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/Classes/src/TCTriggerObject_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/plugins/HistManager_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/plugins/TriggerSelector_cc.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/lib/slc5_amd64_gcc462/libTauAnalysisSVfitStandalone.so");
    gROOT->LoadMacro("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/FlatTreeCreator_${1}_C.so");
  
    TChain* fChain = new TChain("ntupleProducer/eventTree");
    ifstream sourceFiles("$1.txt");
    char line[128];
    int  count = 0;
    cout<< "Adding files from $1 to chain..."<< endl;
     while (sourceFiles >> line) {
        fChain->Add(line);
        ++count;
     }
    cout << count<<" files added!"<<endl;
    sourceFiles.close();
    TStopwatch timer;
    timer.Start();    
    fChain->Process("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/FlatTreeCreator");

    cout << "\n\nDone!" << endl;
    cout << "CPU Time : " << timer.CpuTime() <<endl;
    cout << "RealTime : " << timer.RealTime() <<endl;                             
    cout <<"\n";
}
EOF
root -l -b -q run_$1.C

echo "WROTE A NEW RUN_{SOURCE}.C FILE"
