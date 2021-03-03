#include <iostream>
#include <fstream>
#include <string>


#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>


#include "event.hh"


#define NBIN_ENERGY 560
#define MIN_ENERGY 0
#define MAX_ENERGY 14

#define NBIN_SEGMENT 154
#define MIN_SEGMENT 0
#define MAX_SEGMENT 154


using namespace std;

/* Prototypes {{{ */
 
/* Argument requirement to launch the program */
inline void helper (char *prgramName);

/* Analyze root file */
void analyzeRootFile (string rootFile, char *outname);

void CutAndSaveEvents (vector<Event> *events, char *outname);

auto hLiveSegment = new TH1F("hLiveSegment", "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);

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

    analyzeRootFile( string( argv[1] ), argv[2] );

    return 0;
}


void analyzeRootFile(string rootFile, char *outname){
    auto _file = new TFile(rootFile.c_str());

    cout << "Processing " << rootFile << " ... \n";

    vector<Event> *events = new vector<Event>;

    auto tree = (TTree*)_file->Get("PhysPulse");

    if (tree == NULL){
	_file->Close();
	return ;
    }


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

    CutAndSaveEvents( events, outname );

    events->clear( );

    _file->Close( );


    cout << "Analysis finished.\n\n";
}


inline void helper(char *programName){
    cout << "\tUsage: " << programName << " filename output\n" << endl;
}


void CutAndSaveEvents (vector<Event> *events, char *outname){     

	auto outFile = new TFile(Form("%s.root", outname), "recreate");

	/* Histograms {{{ */
	auto hEnergyPerEvent    = new TH1F("hEnergyPerEvent",   "All events",          NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
	auto hSinglePulseEvent  = new TH1F("hSinglePulseEvent", "Single Pulse",        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
	auto hSignalCandidate   = new TH1F("hSignalCandidate",  "Signal candidate",    NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
	auto hFiducialization   = new TH1F("hFiducialization",  "Fiducial cut",        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
	auto hMuonAdjacent      = new TH1F("hMuonAdjacent",     "Muon Adjacent veto",  NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
	auto hNeutronRecoil     = new TH1F("hNeutronRecoil",    "Neutron recoil veto", NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
	auto hNeutronCapture    = new TH1F("hNeutronCapture",   "NLi capture veto",    NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
	auto hPileUp            = new TH1F("hPileUp",           "Pile Up Veto",        NBIN_ENERGY, MIN_ENERGY, MAX_ENERGY);
	auto hPulseCandidatePSD = new TH1F("hPulseCandidatePSD","PSD value",	       200, 0, 0.5);

	// count
	auto hLiveSegmentFiducial = new TH1F("hLiveSegmentFiducial", "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);
	auto hLiveSegmentSignal   = new TH1F("hLiveSegmentSignal",   "", NBIN_SEGMENT, MIN_SEGMENT, MAX_SEGMENT);

	// TTrees
	auto tSelectedEvents =  new TTree("tSelectedEvents", "");
	auto tCountEvents    =  new TTree("tCountEvents", "");

	int tSeg = 0,
		tMuonCounts = 0,
		tRecoilCounts = 0,
		tCaptureCounts = 0,
		tClusterCounts = 0;

	double tEnergy = 0.0,
		   tPSD = 0.0,
		   tZ = 0.0,
		   tT = 0.0;

	// Initializing trees
    tSelectedEvents->Branch("E", &tEnergy);
    tSelectedEvents->Branch("PSD", &tPSD);
    tSelectedEvents->Branch("z", &tZ);
    tSelectedEvents->Branch("t", &tT);
    tSelectedEvents->Branch("seg", &tSeg);

    tCountEvents->Branch("muon", &tMuonCounts);
    tCountEvents->Branch("recoil", &tRecoilCounts);
    tCountEvents->Branch("capture", &tCaptureCounts);
    tCountEvents->Branch("cluster", &tClusterCounts);

	/* }}} */

    for (long int iEvent = 0, iMax = events->size(); iEvent < iMax; iEvent++){
        Event *event = &events->at(iEvent);

		// Count events
		if (event->IsMuonEvent()) tMuonCounts++;
		else if (event->IsRecoilEvent()) tRecoilCounts++;
		else if (event->IsCaptureEvent()) tCaptureCounts++;

        float energyEvent = event->GetEnergyEvent();

        hEnergyPerEvent->Fill( energyEvent );

		if ( event->SinglePulseCut() ){
			hSinglePulseEvent->Fill(energyEvent);
			hPulseCandidatePSD->Fill(event->GetPulse(0)->PSD);

			if (event->NeutronCut()){
				hSignalCandidate->Fill(energyEvent);

				if (event->FiducialCut() && event->HeightCut(200)) {
					hFiducialization->Fill(energyEvent);
					hLiveSegmentFiducial->Fill(event->GetPulse(0)->segment);

					if (event->MuonVeto(iEvent, events, 5)){
						hMuonAdjacent->Fill(energyEvent);

						if (event->RecoilVeto(iEvent, events, 5)){
							hNeutronRecoil->Fill(energyEvent);

							if (event->CaptureVeto(iEvent, events, 500)){
								hNeutronCapture->Fill(energyEvent);

								if (event->PileUpVeto(iEvent, events, 2)){
									hPileUp->Fill(energyEvent);
									hLiveSegmentSignal->Fill(event->GetPulse(0)->segment);

									//tZ = event->GetPulse(0)->height;
									//tT = event->GetPulse(0)->time;
									tSeg = event->GetPulse(0)->segment;
									tPSD = event->GetPulse(0)->PSD;
									tEnergy = event->GetPulse(0)->energy;

									tSelectedEvents->Fill();
								}
							}
						}
					}
				}
			}
		}
    }

    tClusterCounts = events->size();

    tCountEvents->Fill();

    // Save histogram
    hPileUp->Write();
    hLiveSegment->Write();
    hMuonAdjacent->Write();
    hNeutronRecoil->Write();
    hEnergyPerEvent->Write();
    hNeutronCapture->Write();
    hFiducialization->Write();
    hSignalCandidate->Write();
    hSinglePulseEvent->Write();
    hPulseCandidatePSD->Write();
    hLiveSegmentSignal->Write();
    hLiveSegmentFiducial->Write();

    tSelectedEvents->Write();
    tCountEvents->Write();

    outFile->Close();
} 
