#include <iostream>
#include <fstream>
#include <string>


#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TVectorT.h>
#include <TH1.h>
#include <TH2.h>
#include <TVectorT.h>

#include "cluster.hh"

#define NBIN_ENERGY 1970
#define MIN_ENERGY 15
#define MAX_ENERGY 1000

#define NBIN_SEGMENT 154
#define MIN_SEGMENT 0
#define MAX_SEGMENT 154

#define MUON_ANALYSIS 1
#define BDM_ANALYSIS 0

using namespace std;

/* Consider each root file separately */
void analysis (char* filename, char* outname);


/* Analyze root file */
void analyzeRootFile (string rootFile);

void CutEvents (vector<Cluster> *events);

/* Histograms */
auto hEnergyPerEvent    = new TH1F("hEnergyPerEvent",   "All events",                                          NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hMuonEvent         = new TH1F("hMuonEvent",        "Muon event;Energy (MeV);Count",                    NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hSinglePulseEvent  = new TH1F("hSinglePulseEvent", "Single Pulse",                                        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hSignalCandidate   = new TH1F("hSignalCandidate",  "Signal candidate",                                    NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hFiducialization   = new TH1F("hFiducialization",  "Segment-z double fiducial cut",                       NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hMuonAdjacent      = new TH1F("hMuonAdjacent",     "Muon Adjacent veto, #pm 5 #mus",                      NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hNeutronRecoil     = new TH1F("hNeutronRecoil",    "Neutron recoil veto, #pm 5 #mus",                     NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hNeutronCapture    = new TH1F("hNeutronCapture",   "NLi capture veto, #pm 1000 #mus",                     NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hRnPoDecay         = new TH1F("hRnPoDecay",        "Rn-Po correlated decay (#pm 25 cm & #pm 15000 #mus)", NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hBiPoDecay         = new TH1F("hBiPoDecay",        "Bi-Po correlated decay (#pm 25 cm & - 1200 #mus)",    NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hPileUp            = new TH1F("hPileUp",           "Pile Up veto, #pm 4 ns",                              NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hPulseCandidatePSD = new TH1F("hPulseCandidatePSD","PSD value", 200, 0, 0.5);

auto hPulseCandidateDefaultPSD = new TH1F("hPulseCandidateDefaultPSD", "Potential signals using default PSD", NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
auto hPulseCandidateCustomPSD  = new TH1F("hPulseCandidateCustomPSD",  "Potential signals using custom PSD",  NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);

// TODO: Change cuts before filling this histogram
TH1F* hist_PSD_Energy[10];

// count
auto hLiveSegment         = new TH1F("hLiveSegment",         "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);
auto hLiveSegmentFiducial = new TH1F("hLiveSegmentFiducial", "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);
auto hLiveSegmentSignal   = new TH1F("hLiveSegmentSignal",   "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);
auto hLiveSegmentSignal0  = new TH1F("hLiveSegmentSignal0",  "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);
auto hLiveSegmentSignal1  = new TH1F("hLiveSegmentSignal1",  "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);
auto hLiveSegmentNeutron  = new TH1F("hLiveSegmentNeutron",  "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);
auto hLiveSegmentPileUp   = new TH1F("hLiveSegmentPileUp",   "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);

 TVectorT<double> *runtime;
 
// Argument requirement to launch the program
inline void helper (char *prgramName);

// Combining pulses into cluster 
void ParsePhysPulse(TTree* p_tree, vector<Cluster>* p_cluster);

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
    //for (int i = 0; i < 10; i++)
    //    hist_PSD_Energy[i] = new TH1F(Form("hist_PSD_Energy_%i", i), Form("PSD in %i to %i;PSD", i, i + 1), 200, 0, 0.5);


    // Open root file
    auto _inFile = new TFile(argv[1]);

    if (_inFile->IsZombie()){
        cout << argv[1] << " could not be found ...\n";
        cout << "Exiting ...\n";

        exit(-1);
    }
    cout << "Processing " << argv[1] << " ...";

    // Create variable to contains all events
    vector<Cluster> *events = new vector<Cluster>;

    // Retrieve raw data
    auto tree = (TTree*)_inFile->Get("PhysPulse");

    // Combine pulse into cluster
    ParsePhyPulse(tree, events);

    if (tree == NULL){
	    _inFile->Close();
	    return -1;
    }


    // Prevent memory issues
    auto _outFile = new TFile(Form("%s.root", argv[2]), "recreate");
    _outFile->cd();

    // Create a tree to store selected events
    TTree *selectedEvts = new TTree("SelectedEvents", "");

    // Variables to store
    double c_energy;
    double c_segment;
    double c_time;

    selectedEvts->Branch("E",   &c_energy);
    selectedEvts->Branch("seg", &c_segment);
    selectedEvts->Branch("t",   &c_time);

    // Remove background
    for (long int iEvent = 0, iMax = events->size(); iEvent < iMax; iEvent++){
        Cluster *event = &events->at(iEvent);

        float energyEvent = event->GetClusterEnergy();

#ifdef MUON_ANALYSIS
        // Select Muon events
        if (energyEvent > 15){
            c_energy = energyEvent;
            // TODO: Assign the correct segment and time of the event
            c_segment = 0;
            c_time = 0;
        }
#endif

#ifdef BDM_ANALYSIS
        // avoid using nested if statement by using flags and continue
        if (event->SinglePulseCut()){ 
            // DO SOMETHING
        } else continue;
        
        if (event->NeutronPulseCut(4) || event->NeutronPulseCut(6)){ 
            // DO SOMETHING
        } else continue;

        if (event->FiducialCut() && abs(event->GetPulse(0)->height) < 200){ 
            // DO SOMETHING
        } else continue;

        if (event->MuonAdjacentCut(iEvent, events, 5)){
            // DO SOMETHING
        } else continue;
                        
        if (event->NeutronAdjacentCut(iEvent, events, 5, 4)){
            // DO SOMETHING
        } else continue;
                            
        if (event->NeutronAdjacentCut( iEvent, events, 1000, 6)){
            // DO SOMETHING
        } else continue;
                                
        if (event->PileUpCut( iEvent, events, 4)){
            // DO SOMETHING
        } else continue;

        if (event->RnPoDecayCut( iEvent, events, 15000, 250)){
            // DO SOMETHING
        } else continue;
                                        
        if (event->BiPoDecayCut(iEvent, events, 1200, 250)){
            // DO SOMETHING
        } else continue;
#endif

        selectedEvts->Fill();
    }

    cout << " done.\n";

    // Write tree into a file
    selectedEvts->Write();

    // Save runtime
    auto runtime = (TVectorT<double>*) _inFile->Get("runtime");
    runtime->Write("runtime");

    // Free memory
    events->clear();
    _inFile->Close();
    _outFile->Close();

    return 0;
}

void ParsePhysPulse(TTree* p_tree, vector<Cluster>* p_events){
    // Declare variables for branches 
    Long64_t t_event;           // Keep track of events
    int      t_PID,             // PID
             t_segment;         // segment 
    float    t_PSD,             // Pulse Shape Discrimination parameter
             t_energy,          // Energy in [MeV]
             t_height,          // width in the segment in [mm]
             t_dtime;           // dt between two pulse
    double   t_time;            // absolute time [ns]

    // Set branch address 
    p_tree->SetBranchAddress("evt",  &t_event);
    p_tree->SetBranchAddress("PID",  &t_PID);
    p_tree->SetBranchAddress("seg",  &t_segment);
    p_tree->SetBranchAddress("PSD",  &t_PSD);
    p_tree->SetBranchAddress("E",    &t_energy);
    p_tree->SetBranchAddress("z",    &t_height);
    p_tree->SetBranchAddress("t",    &t_time);
    p_tree->SetBranchAddress("dt",   &t_dtime);

    long int nentries = p_tree->GetEntries();

    // memory for event
    Long64_t lastEventID = 0;

    Cluster oneEvent;

    // Loop over all events -- Store
    for (long int ientry = 0; ientry < nentries; ientry++){
        p_tree->GetEntry(ientry);

        // Exclude data from these segments
        switch(t_segment){
	        case 0: case 1: case 2: case 3: case 4: case 5: case 6: case 9: case 11: case 13: case 18: case 21:
            case 23: case 24: case 27: case 32: case 40: case 44: case 68: case 73: case 79: case 102: case 107:
            case 122: case 127: case 130: case 139:
                continue;

                break;
        }

        // determine if a pulse belong to an event
        bool sameEvent = false;

        if (!ientry){
            lastEventID = t_event;
            sameEvent = true;
        }else{
            if ( t_dtime > 20 ) sameEvent = false;
            else sameEvent = lastEventID == t_event;

            lastEventID = t_event;
        }

        // Store the previous event
        if ( !sameEvent && oneEvent.GetNumberOfPulses() ){
            p_events->push_back(oneEvent);
            
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
    }

    p_events->push_back(oneEvent);

    CutEvents(p_events);

    p_events->clear( );

 // Save runtime
    auto runtime = (TVectorT<double>*) _inFile->Get("runtime");     
    auto abstime = (TVectorT<double>*) _inFile->Get("abstime");

    _inFile->Close( );

}

inline void helper(char *programName){
    cout << "\tUsage: " << programName << " filename output\n" << endl;
}


void analysis(char* filename, char* outname){ // {{{

    // Consider root file separately
    //while (textFile >> oneRootFile) 
    analyzeRootFile( string( filename ) );
  
    auto outFile = new TFile(Form("%s.root", outname), "recreate");

    // Save histogram
    hMuonEvent->Write();
    runtime->Write("runtime");

    //for (int i = 0; i < 10; i++) hist_PSD_Energy[i]->Write();

    outFile->Close();

    cout << "Analysis finished.\n\n";
} //}}}


void CutEvents (vector<Cluster> *events){ 
    for (long int iEvent = 0, iMax = events->size(); iEvent < iMax; iEvent++){
        Cluster *event = &events->at(iEvent);

        float energyEvent = event->GetClusterEnergy();

        hEnergyPerEvent->Fill( energyEvent );
	
	if ( event->MuonEvent( ) )						  { hMuonEvent->Fill(energyEvent);
}
} 
}
