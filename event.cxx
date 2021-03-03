#include "event.hh"
#include <cmath>

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
Statistics_t protonRecoil[27] = {{2.608e-1, 2.094e-2}, {2.156e-1, 1.541e-2}, {2.143e-1, 1.450e-2}, {2.087e-1, 1.325e-2},
				    			 {2.040e-1, 1.291e-2}, {2.011e-1, 1.510e-2}, {1.978e-1, 1.393e-2}, {1.952e-1, 1.340e-2},
				    			 {1.921e-1, 1.283e-2}, {1.900e-1, 1.268e-2}, {1.867e-1, 1.077e-2}, {1.853e-1, 1.087e-2},
				    			 {1.839e-1, 1.024e-2}, {1.838e-1, 1.172e-2}, {1.812e-1, 9.326e-3}, {1.796e-1, 9.306e-3},
				    			 {1.794e-1, 9.251e-3}, {1.787e-1, 9.301e-3}, {1.785e-1, 1.016e-2}, {1.766e-1, 8.876e-3},
				    			 {1.774e-1, 8.946e-3}, {1.759e-1, 9.287e-3}, {1.753e-1, 9.768e-3}, {1.741e-1, 1.104e-2},
				    			 {1.736e-1, 1.101e-2}, {1.722e-1, 1.272e-2}, {1.679e-1, 1.486e-2}};

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


// time of event is the median pulse
// 
//
bool Event::SearchEventInTime( int iEvent, std::vector<Event> *allEvents, double time , double minE, int PID = -1){
    auto currentEvent = &allEvents->at(iEvent);

    assert (currentEvent->IsSinglePulse());

    double currentPulseTime = currentEvent->GetPulse(0)->time;

    int offset = time < 0 ? -1 : 1;

    time = abs(time);

    long int i = iEvent + offset;

    bool foundPulse = false,
	 timeRequirement = false,
	 energyRequirement = false,
	 PIDRequirement = false;

    
    for ( long int i = iEvent + offset, j = allEvents->size();; i += offset ){
		// Limit reached
		if (i == 0 || i == j) return true;

		auto temp = &allEvents->at(i);	
		
		// Energy requirement
		energyRequirement = temp->GetEnergyEvent() > minE ? true : false;

		if (!energyRequirement) continue;

		double tempEventTime;

		// PID and time requirements
		if ( PID == -1 ) PIDRequirement = true;
		else if (PID == 4) PIDRequirement = this->HasNeutronRecoil();
		else if (PID == 6) PIDRequirement = this->HasNeutronCapture();

		if (!PIDRequirement) continue;
		
		// time is the median time
		int size = temp->GetNumberOfPulses();

		if (size % 2 == 0) 
			tempEventTime = (temp->GetPulse(size / 2 - 1)->time + temp->GetPulse(size / 2)->time) / 2;
		else 
			tempEventTime = temp->GetPulse(size / 2)->time;

		// Time requirement
		timeRequirement =  abs(currentPulseTime - tempEventTime) * 1e-3 > time;

		if (timeRequirement) return true;
		else {
		   foundPulse = PIDRequirement && energyRequirement;

		   if (foundPulse) return false;
		}
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

bool Event::NeutronCut(){

    if ( !this->SinglePulseCut( ) ) return false;

    double energyEvent = this->GetEnergyEvent( );

    int iHist = -1;
    if ( energyEvent < 14 && energyEvent >= 0.5 ) iHist = floor(energyEvent / 0.5) - 1;

    if (iHist == -1 ) return false;
    
	Statistics_t p = protonRecoil[iHist];

    double neutronBandMin = p.mean - 2 * p.std;
    double neutronBandMax = p.mean + 2 * p.std;

    double pulsePSD = this->GetPulse( 0 )->PSD;

   return pulsePSD >= neutronBandMin && pulsePSD <= neutronBandMax;
}

bool Event::HasNeutronRecoil(){
    return hasNeutronRecoil;
}

bool Event::HasNeutronCapture(){
    return hasNeutronCapture;
}

bool Event::IsBetaDecayEvent(){
    return isGammaEvent && this->IsSinglePulse() == 0;
}

bool Event::IsMuonEvent(){
    return isMuonEvent;
}

bool Event::IsCaptureEvent(){
    return hasNeutronCapture;
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
    return SearchEventInTime(iEvent, allEvents, +time, 0, 6);
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
