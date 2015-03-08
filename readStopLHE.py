import sys
import ROOT as rt
import math
from LHEevent import *
from LHEfile import *
import plotTools

if __name__ == '__main__':

    #Stop histograms
    MChargino = rt.TH1D("MChargino", "MChargino", 50, 0., 150.0)
    MLSP = rt.TH1D("MLSP", "MLSP", 90, 0, 30)
    MStop = rt.TH1D("MStop", "MStop", 50, 100., 250.)
    Mb = rt.TH1F("Mb", "Mb", 18, 2.0, 8.0)
    PElectronPt = rt.TH1D("PElectronPt", "PElectronPt", 60, 0., 300) 
    MElectronPt = rt.TH1D("MElectronPt", "MElectronPt", 60, 0., 300)
    PMuonPt = rt.TH1D("PMuonPt", "PMuonPt", 60, 0., 300)
    MMuonPt = rt.TH1D("MMuonPt", "MMuonPt", 60, 0., 300)
    PTauPt = rt.TH1D("PTauPt", "PTauPt", 60, 0., 300)
    MTauPt = rt.TH1D("MTauPt", "MTauPt", 60, 0., 300)
    # find events in file
    myLHEfile = LHEfile(sys.argv[1])
    myLHEfile.setMax(100000)
    eventsReadIn = myLHEfile.readEvents()
    for oneEvent in eventsReadIn:
        sum_el = 0
        sum_mu = 0
        sum_ta = 0
        # read the event content
        myLHEevent = LHEevent()
        myLHEevent.fillEvent(oneEvent)
        # fill topology-specific histograms (this goes in a model loop)
        for i in range(0,len(myLHEevent.Particles)):
            p = myLHEevent.Particles[i]
            #print p['ID']
            if abs(p['ID']) == 1000024: MChargino.Fill(p['M'])
            if abs(p['ID']) == 1000006: MStop.Fill(p['M'])
            if (abs(p['ID']) ==5 and p['M'] > 0.0): Mb.Fill(p['M'])
            if abs(p['ID']) == 1000022: MLSP.Fill(p['M'])
            if (abs(p['ID']) == 11):
              sum_el += abs(p['ID'])
            if (abs(p['ID']) == 13):
              sum_mu += abs(p['ID'])
            if (abs(p['ID']) == 15):
              sum_ta += abs(p['ID'])
        for i in range(0,len(myLHEevent.Particles)):
            p = myLHEevent.Particles[i]
            if (sum_el == 22 and p['ID'] == 11): PElectronPt.Fill(math.sqrt(p['Px']*p['Px']+p['Py']*p['Py'])) 
            if (sum_el == 22 and p['ID'] == -11): MElectronPt.Fill(math.sqrt(p['Px']*p['Px']+p['Py']*p['Py'])) 
            if (sum_mu == 26 and p['ID'] == 13): PMuonPt.Fill(math.sqrt(p['Px']*p['Px']+p['Py']*p['Py']))
            if (sum_mu == 26 and p['ID'] == -13): MMuonPt.Fill(math.sqrt(p['Px']*p['Px']+p['Py']*p['Py']))
            if (sum_ta == 30 and p['ID'] == 15): PTauPt.Fill(math.sqrt(p['Px']*p['Px']+p['Py']*p['Py'])) 
            if (sum_ta == 30 and p['ID'] == -15): MTauPt.Fill(math.sqrt(p['Px']*p['Px']+p['Py']*p['Py']))
        del oneEvent, myLHEevent
        
    # write the histograms
    histoFILE = rt.TFile(sys.argv[2],"RECREATE")
    MChargino.Write()
    MLSP.Write()
    Mb.Write()
    MStop.Write()
    PElectronPt.Write()
    MElectronPt.Write()
    PMuonPt.Write()
    MMuonPt.Write()
    PTauPt.Write()
    MTauPt.Write()
    histoFILE.Close()
