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


using namespace std;

/* Prototypes {{{ */
 
/* Argument requirement to launch the program */
inline void helper (char *prgramName);


/* Consider each root file separately */
void analysis (char* filename, char* outname);


/* Analyze root file */
void analyzeRootFile (string rootFile);

void CutEvents (vector<Event> *events);

/* Histograms {{{ */
auto hEnergyPerEvent    = new TH1F("hEnergyPerEvent",   "All events",                                          NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hSinglePulseEvent  = new TH1F("hSinglePulseEvent", "Single Pulse",                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hSignalCandidate   = new TH1F("hSignalCandidate",  "Signal candidate",                                    NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hFiducialization   = new TH1F("hFiducialization",  "Segment-z double fiducial cut",                       NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hMuonAdjacent      = new TH1F("hMuonAdjacent",     "Muon Adjacent veto, #pm 5 #mus",                      NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hNeutronRecoil     = new TH1F("hNeutronRecoil",    "Neutron recoil veto, #pm 5 #mus",                     NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hNeutronCapture    = new TH1F("hNeutronCapture",   "NLi capture veto, #pm 1000 #mus",                     NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hRnPoDecay         = new TH1F("hRnPoDecay",        "Rn-Po correlated decay (#pm 25 cm & #pm 15000 #mus)", NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hBiPoDecay         = new TH1F("hBiPoDecay",        "Bi-Po correlated decay (#pm 25 cm & - 1200 #mus)",    NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hPileUp            = new TH1F("hPileUp",           "Pile Up veto, #pm 2 ns",                              NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hPulseCandidatePSD = new TH1F("hPulseCandidatePSD","PSD value", 200, 0, 0.5);

auto hPulseCandidateDefaultPSD = new TH1F("hPulseCandidateDefaultPSD", "Potential signals using default PSD", NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hPulseCandidateCustomPSD  = new TH1F("hPulseCandidateCustomPSD", "Potential signals using custom PSD", NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);

TH1F* hist_PSD_Energy[10];

// count
auto hLiveSegment         = new TH1F("hLiveSegment",         "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);
auto hLiveSegmentFiducial = new TH1F("hLiveSegmentFiducial", "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);

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

    vector<Event> *events;

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
            events->push_back(oneEvent);
            
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

        oneEvent.AddPulse(pulse);

        hLiveSegment->Fill(t_segment);
    }

    events->push_back(oneEvent);

    totalRunTime += t_time;

    CutEvents( events );

    events->clear( );

    delete events; 

    _file->Close( );
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
    hPileUp->Write();
    hRnPoDecay->Write();
    hBiPoDecay->Write();
    hLiveSegment->Write();
    hMuonAdjacent->Write();
    hNeutronRecoil->Write();
    hEnergyPerEvent->Write();
    hNeutronCapture->Write();
    hFiducialization->Write();
    hSignalCandidate->Write();
    hSinglePulseEvent->Write();
    hPulseCandidatePSD->Write();
    hLiveSegmentFiducial->Write();

    hPulseCandidateCustomPSD->Write();
    hPulseCandidateDefaultPSD->Write();

    for (int i = 0; i < 10; i++) hist_PSD_Energy[i]->Write();

    outFile->Close();

    totalRunTime = (totalRunTime - 18419.3) * 1e-9;

    cout << "Analysis finished.\n\n";
    cout << "**********************\n";
    cout << "Total running time: " << totalRunTime << " s\n";

} //}}}


void CutEvents (vector<Event> *events){ // {{{
    for (long int iEvent = 0, iMax = events->size(); iEvent < iMax; iEvent++){
        Event *event = &events->at(iEvent);

        float energyEvent = event->GetEnergyEvent();

        hEnergyPerEvent->Fill( energyEvent );

        if ( event->SinglePulseCut( ) )                                           { hSinglePulseEvent->Fill(energyEvent);
        if ( event->NeutronPulseCut( 4 ) || event->NeutronPulseCut( 6 ) )         { hSignalCandidate->Fill( energyEvent );
        if ( event->FiducialCut() && abs ( event->GetPulse( 0 )->height )  < 200 ){ hFiducialization->Fill( energyEvent );
                                                                                    hLiveSegmentFiducial->Fill( event->GetPulse( 0 )->segment );
        if ( event->MuonAdjacentCut( iEvent, events, 5 ) )                        { hMuonAdjacent->Fill( energyEvent );
        if ( event->NeutronAdjacentCut( iEvent, events, 5, 4 ) )                  { hNeutronRecoil->Fill( energyEvent );
        if ( event->NeutronAdjacentCut( iEvent, events, 1000, 6 ) )               { hNeutronCapture->Fill( energyEvent );
        if ( event->PileUpCut( iEvent, events, 2 ) )                              { hPileUp->Fill( energyEvent );
        if ( event->RnPoDecayCut( iEvent, events, 15000, 250 ) )                  { hRnPoDecay->Fill( energyEvent );
        if ( event->BiPoDecayCut( iEvent, events, 1200, 250 ) )                   { hBiPoDecay->Fill( energyEvent );

            }}}}}}}}}}
} // }}}




