#include <iostream>
#include <fstream>
#include <string>


#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>


#include "event.hh"
#include "cut.hh"


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

void CutEvents (Events *events);

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
auto hist_SegmentRateSignal                  = new TH1F("hist_SegmentRateSignal", "Signal segment rate",
                                                        154, 0, 154);

TH1F* hist_PSD_Energy[10];

// count
auto hist_liveSegment                       = new TH1F("hist_liveSegment",                    "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);

/* }}} */


/* Count variables {{{ */

long double totalRunTime = 0.0,
            deadTimeRnPo_d = 0.0,
            deadTimeBiPo_d = 0.0,
            deadTimeNeutronRecoil = 0.0,
            deadTimeNLi = 0.0,
            deadTimeMuon = 0.0,
            deadTimePile = 0.0;
 
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
    hist_SegmentRateSignal->Write();
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

    outFile->Close();

    totalRunTime = (totalRunTime - 18419.3) * 1e-9;
    
    // conversion
    deadTimeRnPo_d = deadTimeRnPo_d * 1e-4 / totalRunTime; 
    deadTimeBiPo_d = deadTimeBiPo_d * 1e-4 / totalRunTime; 
    deadTimeNeutronRecoil = deadTimeNeutronRecoil  * 1e-4 / totalRunTime;
    deadTimeNLi = deadTimeNLi  * 1e-4 / totalRunTime;
    deadTimeMuon = deadTimeMuon * 1e-4 / totalRunTime;
    deadTimePile = deadTimePile * 1e-4 / totalRunTime;


    cout << "Analysis finished.\n\n";
    cout << "**********************\n";
    cout << "Total running time: " << totalRunTime << " s\n";
    cout << "Dead time \n";
    cout << "\tRnPo: " << deadTimeRnPo_d << " %" << endl;
    cout << "\tBiPo: " << deadTimeBiPo_d << " %" << endl;
    cout << "\tNeutron Recoil: " << deadTimeNeutronRecoil << " %" << endl;
    cout << "\tNLi Capture: " << deadTimeNLi << " %" << endl;
    cout << "\tMuon: " << deadTimeMuon << " %" << endl;
    cout << "\tPileUp: " << deadTimePile << " %" << endl;

    cout << "\nTotal inefficiency: " << deadTimeRnPo_d + deadTimeBiPo_d + deadTimeNeutronRecoil +
        deadTimeNLi + deadTimeMuon + deadTimePile << endl;
} //}}}


void CutEvents (Events *events){

    Cut *cutEvent = new Cut(events, hist_EnergyPerEvent);

    cutEvent->addCut("singlePulseEvent", 4, 6, hist_Signal);
    cutEvent->addCut("fiducialAndHeight", 200, hist_Signal_Segmet_z_DoubleFV);
    cutEvent->addCut("muonAdjacent", 5, hist_Signal_MuonAdjacent_time5);
    cutEvent->addCut("neutronAdjacent", 4, 5, hist_Signal_NeutronRecoilAdjacent_time5);
    cutEvent->addCut("neutronAdjacent", 6, 400, hist_Signal_NLiCaptureAdjacent_time400);
    cutEvent->addCut("pileUp", 2, hist_Signal_PileUp_time2);
    cutEvent->addCut("RnPoDecay", 15000, 250, hist_Signal_RnPoCorrelatedDecay);
    cutEvent->addCut("BiPoDecay", 1200, 250, hist_Signal_BiPoCorrelatedDecay);
    hist_SegmentRateSignal = cutEvent->GetLiveSegment();

    cutEvent->Run();

    deadTimeRnPo_d = cutEvent->RnPoDeadTime();
    deadTimeBiPo_d = cutEvent->BiPoDeadTime();
    deadTimeNeutronRecoil = cutEvent->NeutronAdjacentDeadTime(4);
    deadTimeNLi = cutEvent->NeutronAdjacentDeadTime(6);
    deadTimeMuon = cutEvent->MuonAdjacentDeadTime();
    deadTimePile = cutEvent->PileUpDeadTime();
}
