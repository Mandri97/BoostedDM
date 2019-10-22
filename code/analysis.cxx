#include <iostream>
#include <fstream>
#include <string>


#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>


#include "event.hh"


#define NBIN_ENERGY 400
#define MIN_ENERGY 0
#define MAX_ENERGY 10

#define NBIN_SEGMENT 154
#define MIN_SEGMENT 0
#define MAX_SEGMENT 154


struct Statistics_t {
    float mean,
          std;
};


struct PSD_t {
    Statistics_t gammaBand,
                 neutronBand,
                 farany;
};


// Single pulse + Height (55 cm) + PileUp 2 us
PSD_t PSD_per_energy[10] = {
    { {0.1444, 0.01403}, {0.2429, 0.02336}, {0.2638, 0.01653} },       // Energy: 0.7 - 1 MeV
    { {0.1435, 0.01143}, {0.2202, 0.01755}, {0.2573, 0.01438} },       // Energy: 1 - 2 MeV 
    { {0.1423, 0.00857}, {0.2070, 0.01344}, {0.2508, 0.01312} },       // Energy: 2 - 3 Mev
    { {0.1397, 0.00786}, {0.2005, 0.01538}, {0.2498, 0.01102} },       // Energy: 3 - 4 MeV
    { {0.1395, 0.00728}, {0.1945, 0.01402}, {0.2429, 0.01266} },       // Energy: 4 - 5 MeV
    { {0.1393, 0.00684}, {0.1896, 0.01278}, {0.2382, 0.01382} },       // Energy: 5 - 6 MeV
    { {0.1391, 0.00652}, {0.1857, 0.01197}, {0.2299, 0.01749} },       // Energy: 6 - 7 MeV
    { {0.1388, 0.00628}, {0.1833, 0.01158}, {0.2280, 0.01599} },       // Energy: 7 - 8 MeV
    { {0.1384, 0.00611}, {0.1802, 0.01021}, {0.2180, 0.02036} },       // Energy: 8 - 9 MeV
    { {0.1383, 0.00593}, {0.1780, 0.01008}, {0.2116, 0.02357} }        // Energy: 9 - 10 MeV
};

using namespace std;

/* Prototypes {{{ */
 
/* Argument requirement to launch the program */
inline void helper (char *prgramName);


/* Consider each root file separately */
void analysis (char* filename, char* outname);


/* Analyze root file */
void analyzeRootFile (string rootFile);


void CutEvents (Events *events);

float timeWindow (Pulse_t *a, Pulse_t *b);
float heightDifference (Pulse_t *a, Pulse_t *b);


/* Set of Cuts */
bool neutronAdjacent (int iCurrentEvent, Events *allEvents, int PID, float time);
bool muonAdjacent (int iCurrentEvent, Events *allEvents, float time);
bool pileUp (int iCurrentEvent, Events *allEvents, float time);
bool correlatedDecayRnPo(int iCurrentEvent, Events *allEvents, float time, float height);
bool correlatedDecayBiPo(int iCurrentEvent, Events *allEvents, float time, float height);
bool doubleFiducialCut (int segment);

/* }}} */


/* Histograms {{{ */

auto hist_EnergyPerEvent   = new TH1F("hist_EnergyPerEvent",   "All events",   NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hist_SinglePulseEvent = new TH1F("hist_SinglePulseEvent", "Single Pulse", NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);

auto hist_Signal                             = new TH1F("hist_Signal",                             "Signal", 
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hist_Signal_Segmet_z_DoubleFV           = new TH1F("hist_Signal_Segmet_z_DoubleFV",           "Segment-z double fiducial cut", 
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hist_Signal_MuonAdjacent_time5          = new TH1F("hist_Signal_MuonAdjacent_time5",          "Muon Adjacent veto within #pm 5 #mus", 
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hist_Signal_NeutronRecoilAdjacent_time5 = new TH1F("hist_Signal_NeutronRecoilAdjacent_time5", "Neutron recoil veto within #pm 5 #mus", 
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hist_Signal_NLiCaptureAdjacent_time400  = new TH1F("hist_Signal_NLiCaptureAdjacent_time400",  "nLi capture veto within #pm 400 #mus", 
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hist_Signal_PileUp_time2                = new TH1F("hist_Signal_PileUp_time2",                "PileUp veto within #pm 2 #mus", 
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hist_Signal_RnPoCorrelatedDecay         = new TH1F("hist_Signal_RnPoCorrelatedDecay",         "Rn-Po correlated decay (#pm 25 cm & #pm 15000 #mus)", 
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hist_Signal_BiPoCorrelatedDecay         = new TH1F("hist_Signal_BiPoCorrelatedDecay",         "Bi-Po correlated decay (#pm 25 cm & - 1200 #mus)", 
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hist_BetaDecay                          = new TH1F("hist_BetaDecay", "", 
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hist_Energy_vs_PSD_noPIDCut             = new TH2F("hist_Energy_vs_PSD_noPIDCut",             "Energy vs PSD",
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY, 200, 0, 1);
auto hist_Signal_PSD                         = new TH1F("hist_Signal_PSD", "PSD value", 
                                                        200, 0, 0.5);
auto hist_default_PSD                        = new TH1F("hist_default_PSD", "Potential signals using default PSD", 
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hist_custom_PSD                         = new TH1F("hist_custom_PSD", "Potential signals using custom PSD", 
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
TH1F* hist_PSD_Energy[10];

// count
auto hist_liveSegment                    = new TH1F("hist_liveSegment",                    "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);
auto hist_liveSegment_Segment_z_DoubleFV = new TH1F("hist_liveSegment_Segment_z_DoubleFV", "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);

/* }}} */


/* Count variables {{{ */

long double totalRunTime = 0.0;
 
 /* }}} */

int main(int argc, char* argv[]){
    // Check for argument 
    if (argc < 3){
        cout << "Error: Not enough argument\n";
        helper(argv[0]);

        return 1;
    }else if (argc > 3){
        cout << "Too many argument \n";
        helper(argv[0]);

        return 1;
    }

    // Initialize histograms
    for (int i = 0; i < 10; i++)
        hist_PSD_Energy[i] = new TH1F(Form("hist_PSD_Energy_%i", i), Form("PSD in %i to %i;PSD", i, i + 1), 200, 0, 0.5);

    analysis(argv[1], argv[2]);

    return 0;
}


void analyzeRootFile(string rootFile){
    auto _file = new TFile(rootFile.c_str());

    cout << "Processing " << rootFile << " ... \n";

    Events *events = new Events();

    auto tree = (TTree*)_file->Get("PhysPulse");

    /* Declare variables for branches */
    Long64_t t_event;           // Keep track of events
    int      t_PID,             // PID
             t_segment;         // segment 
    float    t_PSD,             // Pulse Shape Discrimination parameter
             t_energy,          // Energy in [MeV]
             t_height,          // width in the segment in [mm]
             t_dtime;           // dt between two pulse
    double   t_time;            // absolute time [ns]

    /* Set branch address */
    tree->SetBranchAddress("evt",  &t_event);
    tree->SetBranchAddress("PID",  &t_PID);
    tree->SetBranchAddress("seg",  &t_segment);
    tree->SetBranchAddress("PSD",  &t_PSD);
    tree->SetBranchAddress("E",    &t_energy);
    tree->SetBranchAddress("z",    &t_height);
    tree->SetBranchAddress("t",    &t_time);
    tree->SetBranchAddress("dt",   &t_dtime);

    long int nentries = tree->GetEntries();

    // memory for event
    Long64_t lastEventID = 0;

    Event oneEvent;

    // Loop over all events -- Store
    for (long int ientry = 0; ientry < nentries; ientry++){
        tree->GetEntry(ientry);

        // Skip (pretended) dead segments
        switch(t_segment){
            case 2: case 4: case 6: case 11: case 13: case 18: case 21: case 32: case 44: case 79:
                continue;

                break;
        }

        // determine if a pulse belong to an event
        bool sameEvent = false;

        if (!ientry){
            lastEventID = t_event;
            sameEvent = true;

            // substract total running time
            totalRunTime -= t_time;
        }else{
            sameEvent = lastEventID == t_event;
            lastEventID = t_event;
        }

        // Store the previous event
        if (!sameEvent){
            events->addEvent(oneEvent);
            
            // reset
            oneEvent.Initialize();
        }

        Pulse_t pulse;
        
        pulse.PID     = t_PID;
        pulse.PSD     = t_PSD;
        pulse.time    = t_time;
        pulse.dtime   = t_dtime;
        pulse.height  = t_height;
        pulse.energy  = t_energy;
        pulse.segment = t_segment;

        oneEvent.addPulse(pulse);

        hist_liveSegment->Fill(t_segment);
    }

    events->addEvent(oneEvent);

    totalRunTime += t_time;

    CutEvents(events);

    delete events; 

    _file->Close();
}


inline void helper(char *programName){
    cout << "\tUsage: " << programName << " filename output\n" << endl;
}


void analysis(char* filename, char* outname){ // {{{

    ifstream textFile(filename);

    string oneRootFile;

    // Consider root file separately
    while (textFile >> oneRootFile) analyzeRootFile(oneRootFile);

    auto outFile = new TFile(Form("%s.root", outname), "recreate");

    // Save histogram
    hist_Signal->Write();
    hist_BetaDecay->Write();
    hist_EnergyPerEvent->Write();
    hist_SinglePulseEvent->Write();
    hist_Signal_PileUp_time2->Write();
    hist_Signal_Segmet_z_DoubleFV->Write();
    hist_Signal_MuonAdjacent_time5->Write();
    hist_Signal_RnPoCorrelatedDecay->Write();
    hist_Signal_BiPoCorrelatedDecay->Write();
    hist_Signal_NLiCaptureAdjacent_time400->Write();
    hist_Signal_NeutronRecoilAdjacent_time5->Write();

    hist_Energy_vs_PSD_noPIDCut->Write();
    hist_Signal_PSD->Write();
    
    hist_default_PSD->Write();
    hist_custom_PSD->Write();

    for (int i = 0; i < 10; i++) hist_PSD_Energy[i]->Write();

    hist_liveSegment->Write();
    hist_liveSegment_Segment_z_DoubleFV->Write();

    outFile->Close();

    totalRunTime = (totalRunTime - 18419.3) * 1e-9;

    cout << "Analysis finished.\n\n";
    cout << "**********************\n";
    cout << "Total running time: " << totalRunTime << " s\n";

} //}}}


void CutEvents (Events *events){ // {{{
    
        for (long int iEvent = 0, iMax = events->getNumberOfEvents(); iEvent < iMax; iEvent++){

        Event *event = events->getEvent(iEvent);

        float energyEvent = event->getEnergyEvent();

        hist_EnergyPerEvent->Fill(energyEvent);

        if (event->isBetaDecayEvent()) hist_BetaDecay->Fill(energyEvent);

        if (event->isSinglePulse() != 0) hist_SinglePulseEvent->Fill(energyEvent);

        if (event->isSinglePulse() == 4 || event->isSinglePulse() == 6){
            hist_default_PSD->Fill(energyEvent);
        }

        int iHist = -1;

        // get energy range
        if (energyEvent < 10 && energyEvent > 0 && event->isSinglePulse() != 0) iHist = (int) energyEvent;

        if (iHist == -1) continue;
        
        PSD_t psdEnergy = PSD_per_energy[iHist];
        float neutronBandMin = psdEnergy.neutronBand.mean - psdEnergy.neutronBand.std;
        float neutronBandMax = psdEnergy.neutronBand.mean + psdEnergy.neutronBand.std;

        float faranyMin = psdEnergy.farany.mean - psdEnergy.farany.std;
        float faranyMax = psdEnergy.farany.mean + psdEnergy.farany.std;
    
        if ((event->getPulse(0)->PSD >= neutronBandMin && event->getPulse(0)->PSD <= neutronBandMax) ||
            (event->getPulse(0)->PSD >= faranyMin  && event->getPulse(0)->PSD <= faranyMax          )){
            hist_custom_PSD->Fill(energyEvent);
        }

        // Select potential signals
        if (event->isSinglePulse() == 4 || event->isSinglePulse() == 6){ // No PID requirement for potential signals
            hist_Signal->Fill(energyEvent);

            /*
            if (abs(event->getPulse(0)->height) < 550 && pileUp(iEvent, events, 2)){
                int iHist = -1;

                // get energy range
                if (energyEvent < 10) iHist = (int) energyEvent;

                if (energyEvent >= 0.7 && energyEvent < 1) hist_PSD_Energy[0]->Fill(event->getPulse(0)->PSD);
                else if (iHist > 0)                        hist_PSD_Energy[iHist]->Fill(event->getPulse(0)->PSD);

                if (energyEvent >= 3.5 && energyEvent < 10) hist_Signal_PSD->Fill(event->getPulse(0)->PSD);
            }
            */
            
            // Segment and z double fiducial cuts
            if (doubleFiducialCut(event->getPulse(0)->segment) && abs(event->getPulse(0)->height) < 200){
                hist_Signal_Segmet_z_DoubleFV->Fill(energyEvent);

                hist_liveSegment_Segment_z_DoubleFV->Fill(event->getPulse(0)->segment);

                // Muon adjacent veto within 5 us
                if (muonAdjacent(iEvent, events, 5)){
                    hist_Signal_MuonAdjacent_time5->Fill(energyEvent);

                    // Neutron recoil adjacent veto within 5 us
                    if (neutronAdjacent(iEvent, events, 4, 5)){
                        hist_Signal_NeutronRecoilAdjacent_time5->Fill(energyEvent);

                        // nLi capture adjacent veto within 400 us
                        if (neutronAdjacent(iEvent, events, 6, 400)){
                            hist_Signal_NLiCaptureAdjacent_time400->Fill(energyEvent);

                            // PileUp veto within 2 us
                            if (pileUp(iEvent, events, 2)){
                                hist_Signal_PileUp_time2->Fill(energyEvent);

                                // Rn-Po Correlated decay +-250 mm +- 15000 us
                                if (correlatedDecayRnPo(iEvent, events, 15000, 250)){
                                    hist_Signal_RnPoCorrelatedDecay->Fill(energyEvent);

                                    // Bi-Po Correlated Decay +- 250 mm - 1200 us
                                    if (correlatedDecayBiPo(iEvent, events, 1200, 250)){
                                        hist_Signal_BiPoCorrelatedDecay->Fill(energyEvent);

                                        hist_Energy_vs_PSD_noPIDCut->Fill(energyEvent, event->getPulse(0)->PSD);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
} // }}}


float timeWindow (Pulse_t *a, Pulse_t *b){ // {{{
    // return the time difference in microseconds between 2 pulses
    
    return abs(a->time - b->time) * 1e-3;
} // }}}
 

float heightDifference (Pulse_t *a, Pulse_t *b){ // {{{
    // return the height difference in mm between 2 pulses

    return abs(a->height - b->height);
} // }}}


bool neutronAdjacent (int iCurrentEvent, Events *allEvents, int PID, float time){ // {{{

    /* Check if there is a neutron event adjacent to the pulse within certain time interval 
     *
     * @param:
     *        - iCurrentEvent: index of the current event in @allEvents
     *        - allEvents: all events in the root file
     *        - PID: PID of the neutron signature (4: neutron recoil, 6 nLi capture)
     *        - time: time interval
     *
     * @return:
     *        - true: if no neutron adjacent (PID requirement) within @time us
     *        - false; if there is a neutron adjancet (PID requirement) within @time us
     */

    // Current pulses
    Event *event = allEvents->getEvent(iCurrentEvent);

    bool prevNeutronAdjacent = true;
    
    long int iPrev = iCurrentEvent - 1;

    for ( ; iPrev >= 0; iPrev--){
        
        Event *temp = allEvents->getEvent(iPrev);

        // PID requirement
        if (PID == 4){

            // Looking for event containing neutron recoil
            if (temp->isContainingNeutronRecoil()){
                prevNeutronAdjacent = timeWindow(event->getPulse(0),
                    temp->getPulse(temp->getNumberOfPulses() - 1)) > time;

                break;
            }

        } else if (PID == 6){

            // Looking for event containing nLi capture
            if (temp->isContainingNLiCapture()){
                prevNeutronAdjacent = timeWindow(event->getPulse(0),
                    temp->getPulse(temp->getNumberOfPulses() - 1)) > time;

                break;
            }

        } else {
            prevNeutronAdjacent = timeWindow(event->getPulse(0),
                temp->getPulse(temp->getNumberOfPulses() - 1)) > time;

            break;
        }
    }

    bool nextNeutronAdjacent = true;

    long int iNext = iCurrentEvent + 1,
             iMax  = allEvents->getNumberOfEvents();

    for ( ; iNext < iMax; iNext++){
        
        Event *temp = allEvents->getEvent(iNext);

        // PID requirement
        if (PID == 4){

            // Looking for event containing neutron recoil
            if (temp->isContainingNeutronRecoil()){
                nextNeutronAdjacent = timeWindow(event->getPulse(0), temp->getPulse(0)) > time;

                break;
            }

        }else if (PID == 6){

            // Looking for event containing nLi capture
            if (temp->isContainingNLiCapture()){
                nextNeutronAdjacent = timeWindow(event->getPulse(0), temp->getPulse(0)) > time;

                break;
            }

        } else {
            nextNeutronAdjacent = timeWindow(event->getPulse(0),
                temp->getPulse(0)) > time;

            break;
        }
    }
    
    return prevNeutronAdjacent && nextNeutronAdjacent;
} // }}}


bool muonAdjacent (int iCurrentEvent, Events *allEvents, float time){ // {{{

    /* Check if here is a muon event adjacent to the pulse within a certain time interval
     * @param:
     *        - iCurrentEvent: index of the current event
     *        - allEvents: all events in the root file
     *        - time: time interval
     *
     * @return:
     *        - true if no muon adjacent within @time us
     *        - false if there is a muon adjacent within @time us
     */

    // Current pulses
    Event *event = allEvents->getEvent(iCurrentEvent);

    bool prevMuonAdjacent = true;
    
    long int iPrev = iCurrentEvent - 1;

    for ( ; iPrev >= 0; iPrev--){
        
        Event *temp = allEvents->getEvent(iPrev);

        // Energy requirement
        if (temp->getEnergyEvent() > 15){

            prevMuonAdjacent =
                timeWindow(event->getPulse(0), temp->getPulse(temp->getNumberOfPulses() - 1)) > time;

            break;
        }
    }

    bool nextMuonAdjacent = true;

    long int iNext = iCurrentEvent + 1,
             iMax  = allEvents->getNumberOfEvents();

    for ( ; iNext < iMax; iNext++){

        Event *temp = allEvents->getEvent(iNext);

        // Energy requirement
        if (temp->getEnergyEvent() > 15){

            nextMuonAdjacent =
                timeWindow(event->getPulse(0), temp->getPulse(0)) > time;

            break;
        }
    }

    return prevMuonAdjacent && nextMuonAdjacent;
} // }}} 


bool pileUp (int iCurrentEvent, Events *allEvents, float time){
    /* Check if there is any event within @time microseconds
     * @param:
     *        - iCurrentEvent: index of the current event in @allEvents
     *        - allEvents: all events in the root file
     *        - time: time interval 
     * Return:
     *        - true: if no event adjacent to the current event within @time microseconds
     *        - false: otherwise
     */

    // Current pulses
    Event *event = allEvents->getEvent(iCurrentEvent);
    
    long int iPrev = iCurrentEvent - 1;

    bool prevPileUp = true;

    if (iPrev > -1){
        Event *temp = allEvents->getEvent(iPrev);

        prevPileUp = timeWindow(event->getPulse(0), temp->getPulse(temp->getNumberOfPulses() - 1)) > time;
    }

    long int iNext = iCurrentEvent + 1;

    bool nextPileUp = true;

    if (iNext < allEvents->getNumberOfEvents()){
        Event *temp = allEvents->getEvent(iNext);

        nextPileUp = timeWindow(event->getPulse(0), temp->getPulse(0)) > time;
    }

    return prevPileUp && nextPileUp;
}

bool correlatedDecayRnPo(int iCurrentEvent, Events *allEvents, float time, float height){

    // Current pulses
    Pulse_t *signal = allEvents->getEvent(iCurrentEvent)->getPulse(0);
    
    long int iPrev = iCurrentEvent - 1;

    bool prevCorrelated = true;

    for (; iPrev >= 0; iPrev--){
        Event *prevEvent = allEvents->getEvent(iPrev);

        // single pulse neutron recoil
        if (prevEvent->isSinglePulse() == 4){
            Pulse_t* prevPulse = prevEvent->getPulse(0);

            if (timeWindow(signal, prevPulse) < time){

                // same height and same segment
                if (signal->segment == prevPulse->segment){
                      prevCorrelated = heightDifference(signal, prevPulse) > height;

                      break;
                } 
            } else break;
        }
    }

    long int iNext = iCurrentEvent + 1,
             iMax  = allEvents->getNumberOfEvents();

    bool nextCorrelated = true;

    for (; iNext < iMax; iNext++){
        Event *nextEvent = allEvents->getEvent(iNext);

        // single pulse neutron recoil
        if (nextEvent->isSinglePulse() == 4){
            Pulse_t* nextPulse = nextEvent->getPulse(0);

            if (timeWindow(signal, nextPulse) < time){

                // same height and same segment
                if (signal->segment == nextPulse->segment){
                    nextCorrelated = heightDifference(signal, nextPulse) > height;

                    break;
                }
            } else break;
        }
    }

    return prevCorrelated && nextCorrelated;
}

bool correlatedDecayBiPo(int iCurrentEvent, Events *allEvents, float time, float height){

    // Current pulses
    Pulse_t* signal = allEvents->getEvent(iCurrentEvent)->getPulse(0);
    
    long int iPrev = iCurrentEvent - 1;

    bool prevCorrelated = true,
         foundRequiredPulse = false;

    for (; iPrev >= 0; iPrev--){
        if (foundRequiredPulse) break;

        Event *prevEvent = allEvents->getEvent(iPrev);

        if ( prevEvent->getEnergyEvent() >= 0.25 && prevEvent->getEnergyEvent() < 3.25 &&  // Energy restriction
            prevEvent->isBetaDecayEvent() ){                                               // Beta decay only

            for (int iPulse = 0, nbPulses = prevEvent->getNumberOfPulses(); iPulse < nbPulses; iPulse++){
                Pulse_t *pulse = prevEvent->getPulse(iPulse);
                
                if (timeWindow(signal, pulse) < time){

                    if ( signal->segment == pulse->segment     ||
                         signal->segment == pulse->segment - 1 ||
                         signal->segment == pulse->segment + 1 ||
                         signal->segment == pulse->segment - 1 - 14 ||
                         signal->segment == pulse->segment     - 14 ||
                         signal->segment == pulse->segment + 1 - 14 ||
                         signal->segment == pulse->segment - 1 + 14 ||
                         signal->segment == pulse->segment     + 14 ||
                         signal->segment == pulse->segment + 1 + 14 ){

                        prevCorrelated = heightDifference(signal, pulse) > height;

                        foundRequiredPulse = true;

                        break;
                    }
                } else {
                    foundRequiredPulse = true;

                    break;
                }
            }
        } 
    }

    // TODO: Remove this
    return prevCorrelated;

    long int iNext = iCurrentEvent + 1,
             iMax = allEvents->getNumberOfEvents();

    bool nextCorrelated = true;
         foundRequiredPulse = false;

    for (; iNext < iMax; iNext++){
        if (foundRequiredPulse) break;

        Event *prevEvent = allEvents->getEvent(iNext);

        // n-pulses
        if (prevEvent->isSinglePulse() == 0){
            for (int iPulse = 0, nbPulses = prevEvent->getNumberOfPulses(); iPulse < nbPulses; iPulse++){
                Pulse_t *pulse = prevEvent->getPulse(iPulse);
                
                if (timeWindow(signal, pulse) < time){

                    // same height and same segment
                    if (signal->segment == pulse->segment){

                        nextCorrelated = heightDifference(signal, pulse) > height;

                        foundRequiredPulse = true;

                        break;
                    }
                } else {
                    foundRequiredPulse = true;

                    break;
                }
            }
        }
    }

    return prevCorrelated && nextCorrelated;
}

/*
* 140 141 142 143 144 145 146 147 148 149 150 151 152 153 
* 126 127 128 129 130 131 132 133 134 135 136 137 138 139
* 112 113 114 115 116 117 118 119 120 121 122 123 124 125
* 98  99  100 101 102 103 104 105 106 107 108 109 110 111
* 84  85  86  87  88  89  90  91  92  93  94  95  96  97
* 70  71  72  73  74  75  76  77  78  79  80  81  82  83
* 56  57  58  59  60  61  62  63  64  65  66  67  68  69
* 42  43  44  45  46  47  48  49  50  51  52  53  54  55
* 28  29  30  31  32  33  34  35  36  37  38  39  40  41
* 14  15  16  17  18  19  20  21  22  23  24  25  26  27
* 0   1   2   3   4   5   6   7   8   9   10  11  12  13 
*/

bool doubleFiducialCut (int segment){
    return ((segment >= 30  && segment <= 39)  ||
            (segment >= 44  && segment <= 53)  ||
            (segment >= 58  && segment <= 67)  ||
            (segment >= 72  && segment <= 81)  ||
            (segment >= 86  && segment <= 95)  ||
            (segment >= 100 && segment <= 109) ||
            (segment >= 114 && segment <= 123) );
}
