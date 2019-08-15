#include <iostream>
#include <fstream>
#include <string>


#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>


#include "event.hh"


#define NBIN_ENERGY 200
#define MIN_ENERGY 0
#define MAX_ENERGY 10

#define NBIN_SEGMENT 154
#define MIN_SEGMENT 0
#define MAX_SEGMENT 154

using namespace std;

/* Prototypes {{{ */
 
/* Argument requirement to launch the program */
inline void helper (char *prgramName);


/* Consider each root file separately */
void analysis (char* filename);


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
auto hist_Signal_RnPoCorrelatedDecay         = new TH1F("hist_Signal_RnPoCorrelatedDecay",         "Rn-Po correlated decay (#pm 25 cm & #pm 7500 #mus)", 
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hist_Signal_BiPoCorrelatedDecay         = new TH1F("hist_Signal_BiPoCorrelatedDecay",         "Bi-Po correlated decay (#pm 30 cm & - 750 #mus)", 
                                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);

// count
auto hist_liveSegment                    = new TH1F("hist_liveSegment",                    "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);
auto hist_liveSegment_Segment_z_DoubleFV = new TH1F("hist_liveSegment_Segment_z_DoubleFV", "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);

/* }}} */


/* Count variables {{{ */

long double totalRunTime = 0.0;
 
 /* }}} */
int main(int argc, char* argv[]){

    // Check for argument 
    if (argc < 2){
        cout << "Error: Not enough argument\n";
        helper(argv[0]);

        return 1;
    }else if (argc > 2){
        cout << "Too many argument \n";
        helper(argv[0]);

        return 1;
    }

    analysis(argv[1]);

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
             t_height;          // width in the segment in [mm]
    double   t_time;            // absolute time [ns]

    /* Set branch address */
    tree->SetBranchAddress("evt",  &t_event);
    tree->SetBranchAddress("PID",  &t_PID);
    tree->SetBranchAddress("seg",  &t_segment);
    tree->SetBranchAddress("PSD",  &t_PSD);
    tree->SetBranchAddress("E",    &t_energy);
    tree->SetBranchAddress("z",    &t_height);
    tree->SetBranchAddress("t",    &t_time);

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
    cout << "\tUsage: " << programName << " filename\n" << endl;
}


void analysis(char* filename){ // {{{

    ifstream textFile(filename);

    string oneRootFile;

    // Consider root file separately
    while (textFile >> oneRootFile) analyzeRootFile(oneRootFile);

    auto outFile = new TFile("analysis.root", "recreate");

    // Save histogram
    hist_Signal->Write();
    hist_EnergyPerEvent->Write();
    hist_SinglePulseEvent->Write();
    hist_Signal_PileUp_time2->Write();
    hist_Signal_Segmet_z_DoubleFV->Write();
    hist_Signal_MuonAdjacent_time5->Write();
    hist_Signal_RnPoCorrelatedDecay->Write();
    hist_Signal_BiPoCorrelatedDecay->Write();
    hist_Signal_NLiCaptureAdjacent_time400->Write();
    hist_Signal_NeutronRecoilAdjacent_time5->Write();

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

        if (event->isSinglePulse() != 0) hist_SinglePulseEvent->Fill(energyEvent);

        // Select signals
        if (event->isSinglePulse() == 4 || event->isSinglePulse() == 6){
            hist_Signal->Fill(energyEvent);

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

                        // nLi capture adjacent veto within 4000 us
                        if (neutronAdjacent(iEvent, events, 6, 400)){
                            hist_Signal_NLiCaptureAdjacent_time400->Fill(energyEvent);

                            // PileUp veto within 2 us
                            if (pileUp(iEvent, events, 2)){
                                hist_Signal_PileUp_time2->Fill(energyEvent);

                                // Rn-Po Correlated decay +-250 mm +- 7500 us
                                if (correlatedDecayRnPo(iEvent, events, 7500, 250)){
                                    hist_Signal_RnPoCorrelatedDecay->Fill(energyEvent);

                                    // Bi-Po Correlated Decay +- 300 mm - 750
                                    if (correlatedDecayBiPo(iEvent, events, 750, 300)){
                                        hist_Signal_BiPoCorrelatedDecay->Fill(energyEvent);
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
    Event *event = allEvents->getEvent(iCurrentEvent);
    
    long int iPrev = iCurrentEvent - 1;

    bool prevCorrelated = true;

    for (; iPrev >= 0; iPrev--){
        Event *temp = allEvents->getEvent(iPrev);

        // single pulse neutron recoil
        if (temp->isSinglePulse() == 4){
            if (timeWindow(event->getPulse(0), temp->getPulse(0)) < time){

                // same height and same segment
                if (event->getPulse(0)->segment == temp->getPulse(0)->segment){
                      prevCorrelated = heightDifference(event->getPulse(0), temp->getPulse(0)) > height;

                      break;
                } 
            } else break;
        }
    }

    long int iNext = iCurrentEvent + 1,
             iMax  = allEvents->getNumberOfEvents();

    bool nextCorrelated = true;

    for (; iNext < iMax; iNext++){
        Event *temp = allEvents->getEvent(iNext);

        // single pulse neutron recoil
        if (temp->isSinglePulse() == 4){
            if (timeWindow(event->getPulse(0), temp->getPulse(0)) < time){

                // same height and same segment
                if (event->getPulse(0)->segment == temp->getPulse(0)->segment){
                      nextCorrelated = heightDifference(event->getPulse(0), temp->getPulse(0)) > height;

                      break;
                } 
            } else break;
        }
    }

    return prevCorrelated && nextCorrelated;
}

bool correlatedDecayBiPo(int iCurrentEvent, Events *allEvents, float time, float height){

    // Current pulses
    Event *event = allEvents->getEvent(iCurrentEvent);
    
    long int iPrev = iCurrentEvent - 1;

    bool prevCorrelated = true,
         foundRequiredPulse = false;

    for (; iPrev >= 0; iPrev--){
        if (foundRequiredPulse) break;

        Event *temp = allEvents->getEvent(iPrev);

        // n-pulses
        if (temp->isSinglePulse() == 0){
            for (int iPulse = 0, nbPulses = temp->getNumberOfPulses(); iPulse < nbPulses; iPulse++){
                Pulse_t *pulse = temp->getPulse(iPulse);
                
                if (timeWindow(event->getPulse(0), pulse) < time){

                    // same height and same segment
		    // TODO: Change 15 August not same segment but +-1 segment
                    if (event->getPulse(0)->segment == pulse->segment ||
			event->getPulse(0)->segment == pulse->segment - 1 ||
			event->getPulse(0)->segment == pulse->segment + 1){
                        prevCorrelated = heightDifference(event->getPulse(0), pulse) > height;

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

    return prevCorrelated;

    // TODO: Need to modify this if you want to include pulse before signal

    long int iNext = iCurrentEvent + 1,
             iMax  = allEvents->getNumberOfEvents();

    bool nextCorrelated = true;

    if (iNext == iMax){
        prevCorrelated = true;
    } else {
        for (; iNext < iMax; iNext++){
            Event *temp = allEvents->getEvent(iNext);

            // n-pulses
            if (temp->isSinglePulse() > 0){
                
                for (int iPulse = 0, j = temp->getNumberOfPulses(); iPulse < j; iPulse++){

                    Pulse_t *pulse = temp->getPulse(iPulse);
                    
                    if (timeWindow(event->getPulse(0), pulse) < time){

                        // same height and same segment
                        if (event->getPulse(0)->segment == pulse->segment){
                            nextCorrelated = heightDifference(event->getPulse(0), pulse) > height;

                            break;
                        }
                    } else {
                         nextCorrelated = true;

                        break;
                    }
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
