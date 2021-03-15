#include "event.hh"

using namespace std;

struct Statistics_t {
    double mean,
           std;
};


struct PSD_t {
    Statistics_t gammaBand,
                 protonBand,
                 nucleusBand;
};

// Single pulse + Height (55 cm) + PileUp 2 us
/* Old parameters 1 MeV interval
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
*/

// New PSD parameters 0.5 MeV interval
Statistics_t protonRecoil[27] = {{2.610e-1, 2.057e-2}, {2.359e-1, 2.677e-2}, {2.159e-1, 1.635e-2}, {2.088e-1, 1.583e-2},
								 {2.048e-1, 1.519e-2}, {2.002e-1, 1.405e-2}, {1.975e-1, 1.417e-2}, {1.947e-1, 1.327e-2},
								 {1.923e-1, 1.297e-2}, {1.899e-1, 1.341e-2}, {1.877e-1, 1.184e-2}, {1.869e-1, 1.259e-2},
								 {1.850e-1, 1.207e-2}, {1.845e-1, 1.116e-2}, {1.822e-1, 1.116e-2}, {1.803e-1, 1.071e-2},
								 {1.801e-1, 1.066e-2}, {1.796e-1, 1.093e-2}, {1.778e-1, 1.085e-2}, {1.772e-1, 1.039e-2},
								 {1.778e-1, 9.558e-3}, {1.755e-1, 9.605e-3}, {1.753e-1, 1.152e-2}, {1.741e-1, 1.034e-2},
								 {1.743e-1, 1.090e-2}, {1.729e-1, 1.166e-2}, {1.707e-1, 1.137e-2}};

/* class Pulse_t {{{ */
Pulse_t::Pulse_t(int seg,  int PID,  double h,  
		 double E, double t, double dt, double PSD){

    PSD = PSD;  PID = PID; height = h;
    energy = E; time = t;  dtime = dt;
    segment = seg;
}

Pulse_t::Pulse_t( ){}

Pulse_t::~Pulse_t( ){}

double Pulse_t::TimeWindow (Pulse_t *a){
    // return the time difference in microseconds between 2 pulses
    return abs(a->time - this->time) * 1e-3;
} 

double Pulse_t::HeightDifference (Pulse_t *a){
    // return the height difference in mm between 2 pulses
    return abs(a->height - this->height);
}

/* }}} */


/* class Event {{{ */

Event::Event(){
    Initialize();
}

Event::~Event(){
    Initialize();
}

void Event::Initialize(){
    pulses.clear();
    energyEvent = 0.0;

    hasNeutronRecoil = false;
    hasNeutronCapture = false;
    isGammaEvent = true;
    isMuonEvent = false;
}


bool Event::SearchEventInTime( int iEvent, std::vector<Event> *allEvents, double time , double minE, int PID = -1){
    auto currentEvent = &allEvents->at(iEvent);

    assert (currentEvent->IsSinglePulse());

    double currentPulseTime = currentEvent->GetPulse(0)->time;

    int offset = time < 0 ? -1 : 1;

    time = abs(time);

    bool timeRequirement = false,
	 	 energyRequirement = false,
	 	 PIDRequirement = false;

    for ( long int i = iEvent + offset, j = allEvents->size();; i += offset ){
		// Limit reached
		if (i == -1 || i == j) return true;

		auto temp = &allEvents->at(i);	
		
		// time is the median time
		int size = temp->GetNumberOfPulses();

		double tempEventTime;

		if (size % 2 == 0) 
			tempEventTime = (temp->GetPulse(size / 2 - 1)->time + temp->GetPulse(size / 2)->time) / 2;
		else 
			tempEventTime = temp->GetPulse(size / 2)->time;

		// Time requirement
		timeRequirement =  abs(currentPulseTime - tempEventTime) * 1e-3 > time;

		if (timeRequirement) return true;

		// Energy requirement
		energyRequirement = temp->GetEnergyEvent() > minE ? true : false;

		if (!energyRequirement) continue;

		// PID and time requirements
		if ( PID == -1 ) PIDRequirement = true;
		else if (PID == 4) PIDRequirement = temp->IsRecoilEvent();
		else if (PID == 6) PIDRequirement = temp->IsCaptureEvent();

		if (!PIDRequirement) continue;
		
		return false;
    }
}

int Event::GetNumberOfPulses(){
    return pulses.size();
}

double Event::GetEnergyEvent(){
    return energyEvent;
}

Pulse_t* Event::GetPulse(int iPulse){
    // Make sure that there will no segmentation fault
    assert (iPulse >= 0 && iPulse < this->GetNumberOfPulses());

    return &pulses.at(iPulse);
}

void Event::AddPulse(Pulse_t pulse){
    pulses.push_back(pulse);

    // Update energyEvent
    energyEvent += pulse.energy;

    if (energyEvent > 15) isMuonEvent = true;

    // Check if it is a beta decay
    isGammaEvent = isGammaEvent && pulse.PID == 1;

    if (pulse.PID == 4) hasNeutronRecoil = true;
    if (pulse.PID == 6) hasNeutronCapture = true;
}

int Event::IsSinglePulse(){
    /* Determine if the event is a single pulse
     * 
     * Return:
     *      - 0 -- if not a single pulse
     *      - PID -- otherwise
     */

    assert (this->GetNumberOfPulses() > 0);

    if (this->GetNumberOfPulses() > 1) return 0;
    else return pulses.at(0).PID;
}

bool Event::SinglePulseCut( ){ return this->IsSinglePulse() != 0; }

bool Event::NeutronCut(TH1D *hMin, TH1D *hMax){

    if ( !this->SinglePulseCut( ) ) return false;

    double energyEvent = this->GetEnergyEvent( );

    int iHist = -1;
    if ( energyEvent < 14 && energyEvent >= 0.5 ) iHist = hMin->GetXaxis()->FindBin(energyEvent);
	else return false;

	double neutronMin = hMin->GetBinContent(iHist);
	double neutronMax = hMax->GetBinContent(iHist);

    double pulsePSD = this->GetPulse( 0 )->PSD;

   return pulsePSD >= neutronMin && pulsePSD <= neutronMax;
}

bool Event::IsBetaDecayEvent(){
    return isGammaEvent && this->IsSinglePulse() == 0;
}

bool Event::IsMuonEvent(){
    return isMuonEvent;
}

bool Event::IsCaptureEvent(){
    return hasNeutronCapture && this->IsSinglePulse() == 6;
}

bool Event::IsRecoilEvent(){
    return hasNeutronRecoil;
}

bool Event::FiducialCut (){
    int segment = this->GetPulse(0)->segment;
    
    int segX = segment % 14,
        segY = segment / 14;

    return ((segX >= 2 && segX <= 11) && (segY >= 3 && segY <= 8));
}

bool Event::HeightCut (double height){
    return abs(this->GetPulse(0)->height) < height;
}

bool Event::MuonVeto (int iEvent, std::vector<Event> *allEvents, double time){
    /* Check if here is a muon event adjacent to the pulse within a certain time interval
     * @param:
     *        - iEvent: index of the current event
     *        - allEvents: all events in the root file
     *        - time: time interval
     *
     * @return:
     *        - true if no muon adjacent within @time us
     *        - false if there is a muon adjacent within @time us
     */

    // Look for muon events before signal
    // The neutron single pulse event is likely 
    // to be create by muon events
    return SearchEventInTime(iEvent, allEvents, -time, 15);
}

bool Event::PileUpVeto(int iEvent, std::vector<Event> *allEvents, double time){
    /* Check if there is any event within @time microseconds
     * @param:
     *        - iEvent: index of the current event in @allEvents
     *        - allEvents: all events in the root file
     *        - time: time interval 
     * Return:
     *        - true: if no event adjacent to the current event within @time microseconds
     *        - false: otherwise
     */
    bool prevEvent = SearchEventInTime(iEvent, allEvents, -time, 0);
    bool nextEvent = SearchEventInTime(iEvent, allEvents, +time, 0);

    return prevEvent && nextEvent;
}

bool Event::RecoilVeto(int iEvent, std::vector<Event> *allEvents, double time){
    bool prevEvent = SearchEventInTime(iEvent, allEvents, -time, 0, 4);
    bool nextEvent = SearchEventInTime(iEvent, allEvents, +time, 0, 4);

    return prevEvent && nextEvent;
}

bool Event::CaptureVeto(int iEvent, std::vector<Event> *allEvents, double time){
    return SearchEventInTime(iEvent, allEvents, -time, 0, 6);
}

// TODO: Broken function - DO NOT USE
bool Event::RnPoDecayCut(int iEvent, std::vector<Event> *allEvents, double time, double height){
    return false;
    /*
    Pulse_t *signalCandidate = this->GetPulse(0);
    
    bool prevEvent = true;
    for (long int iPrev = iEvent - 1; iPrev >= 0; iPrev--){
        Event *temp = &allEvents->at(iPrev);

        // single pulse neutron recoil
        if ( temp->IsSinglePulse( ) == 4 ){
            Pulse_t* prevPulse = temp->GetPulse ( 0 );

            if ( signalCandidate->TimeWindow( prevPulse ) < time){
                // same height and same segment
                if ( signalCandidate->segment == prevPulse->segment ){
                      prevEvent = signalCandidate->HeightDifference( prevPulse ) > height;
                      break;
                } 
            } else break;
        }
    }


    bool nextEvent = true;
    for (long int iNext = iEvent + 1, iMax = allEvents->size(); iNext < iMax; iNext++){
        Event *temp = &allEvents->at(iNext);

        // single pulse neutron recoil
        if (temp->IsSinglePulse ( ) == 4){
            Pulse_t* nextPulse = temp->GetPulse( 0 );

            if (signalCandidate->TimeWindow( nextPulse ) < time){
                // same height and same segment
                if ( signalCandidate->segment == nextPulse->segment ){
                    nextEvent = signalCandidate->HeightDifference( nextPulse ) > height;
                    break;
                }
            } else break;
        }
    }

    bool returnValue = prevEvent && nextEvent;
    
    return returnValue;
    */
}


// TODO: Broken function - DO NOT USE
bool Event::BiPoDecayCut( int iEvent, std::vector<Event> *allEvents, double time, double height ){
    return false;

    /*
    assert(false);

    Pulse_t* signalCandidate = this->GetPulse(0);
    
    bool prevEvent = true,
         foundRequiredPulse = false;

    for ( long int iPrev = iEvent - 1; iPrev >= 0; iPrev-- ){
        if ( foundRequiredPulse ) break;

        Event *temp = &allEvents->at( iPrev );

        if ( temp->GetEnergyEvent() >= 0.25 && temp->GetEnergyEvent() < 3.25 && // Energy restriction
            temp->isBetaDecayEvent() ){                                         // Beta decay only

            for ( int iPulse = 0, nbPulses = temp->GetNumberOfPulses(); iPulse < nbPulses; iPulse++ ){
                Pulse_t *pulse = temp->GetPulse( iPulse );
                
                if ( signalCandidate->TimeWindow( pulse ) < time ){

                    if ( signalCandidate->segment == pulse->segment          ||
                         signalCandidate->segment == pulse->segment - 1      ||
                         signalCandidate->segment == pulse->segment + 1      ||
                         signalCandidate->segment == pulse->segment - 1 - 14 ||
                         signalCandidate->segment == pulse->segment     - 14 ||
                         signalCandidate->segment == pulse->segment + 1 - 14 ||
                         signalCandidate->segment == pulse->segment - 1 + 14 ||
                         signalCandidate->segment == pulse->segment     + 14 ||
                         signalCandidate->segment == pulse->segment + 1 + 14 ){

                        prevEvent = signalCandidate->HeightDifference( pulse ) > height;

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

    return prevEvent;
    */
}

/* }}} */
