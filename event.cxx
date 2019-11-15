#include "event.hh"
#include <cmath>

using namespace std;

struct Statistics_t {
    float mean,
          std;
};


struct PSD_t {
    Statistics_t gammaBand,
                 neutronBand,
                 nLiBand;
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


double Event::_RnPoDeadTime = 0.0;
double Event::_BiPoDeadTime = 0.0;
double Event::_PileUpDeadTime = 0.0;
double Event::_MuonAdjacentDeadTime = 0.0;
double Event::_NeutronRecoilDeadTime = 0.0;
double Event::_NeutronCaptureDeadTime = 0.0;


/* class Pulse_t {{{ */
Pulse_t::Pulse_t(int seg,  int PID,
                float h,  float E, float t,
                float dt, float PSD){

    PSD = PSD;  PID = PID; height = h;
    energy = E; time = t;  dtime = dt;
    segment = seg;
}

Pulse_t::Pulse_t( ){}

Pulse_t::~Pulse_t( ){}

float Pulse_t::TimeWindow (Pulse_t *a){
    // return the time difference in microseconds between 2 pulses
    
    return abs(a->time - this->time) * 1e-3;
} 

float Pulse_t::HeightDifference (Pulse_t *a){
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
    if (pulses.size()) pulses.clear();
    energyEvent = 0.0;

    _hasNeutronRecoil = false;
    _hasNeutronCapture = false;
    _isGammaEvent = true;
}

int Event::GetNumberOfPulses(){
    return pulses.size();
}

float Event::GetEnergyEvent(){
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

    // Check if it is a beta decay
    _isGammaEvent = _isGammaEvent && pulse.PID == 1;

    if (pulse.PID == 4) _hasNeutronRecoil = true;
    if (pulse.PID == 6) _hasNeutronCapture = true;
}

int Event::isSinglePulse(){
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

bool Event::hasNeutronRecoil(){
    return _hasNeutronRecoil;
}

bool Event::hasNeutronCapture(){
    return _hasNeutronCapture;
}

bool Event::isBetaDecayEvent(){
    return _isGammaEvent && this->isSinglePulse() == 0;
}

bool Event::FiducialCut (){
    int segment = this->GetPulse(0)->segment;

    return ((segment >= 44  && segment <= 53)  ||
            (segment >= 58  && segment <= 67)  ||
            (segment >= 72  && segment <= 81)  ||
            (segment >= 86  && segment <= 95)  ||
            (segment >= 100 && segment <= 109) ||
            (segment >= 114 && segment <= 123) );
}

bool Event::HeightCut (float height){
    return abs(this->GetPulse(0)->height) < height;
}

bool Event::MuonAdjacentCut (int iEvent, std::vector<Event> *allEvents, float time){
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

    bool prevEvent = true;
    for (long int iPrev = iEvent - 1; iPrev >= 0; iPrev--){
        
        Event *temp = &allEvents->at(iPrev);

        // Energy requirement
        if (temp->GetEnergyEvent() > 15){

            prevEvent =
                this->GetPulse(0)->TimeWindow(temp->GetPulse(temp->GetNumberOfPulses() - 1)) > time;

            break;
        }
    }


    bool nextEvent = true;
    for (long int iNext = iEvent + 1, iMax = allEvents->size(); iNext < iMax; iNext++){

        Event *temp = &allEvents->at(iNext);

        // Energy requirement
        if (temp->GetEnergyEvent() > 15){

            nextEvent = this->GetPulse(0)->TimeWindow(temp->GetPulse(0)) > time;

            break;
        }
    }


    bool returnValue = prevEvent && nextEvent;
    if (!returnValue) _MuonAdjacentDeadTime += (2 * time);

    return returnValue;
}

bool Event::PileUpCut(int iEvent, std::vector<Event> *allEvents, float time){
    /* Check if there is any event within @time microseconds
     * @param:
     *        - iEvent: index of the current event in @allEvents
     *        - allEvents: all events in the root file
     *        - time: time interval 
     * Return:
     *        - true: if no event adjacent to the current event within @time microseconds
     *        - false: otherwise
     */
    long int iPrev = iEvent - 1;
    bool prevEvent = true;

    if (iPrev > -1){
        Event *temp = &allEvents->at(iPrev);

        prevEvent = this->GetPulse(0)->TimeWindow(temp->GetPulse(temp->GetNumberOfPulses() - 1)) > time;
    }


    long int iNext = iEvent + 1;
    bool nextEvent = true;

    if ( iNext < allEvents->size( ) ){
        Event *temp = &allEvents->at(iNext);

        nextEvent = this->GetPulse(0)->TimeWindow(temp->GetPulse(0)) > time;
    }


    bool returnValue = prevEvent && nextEvent;
    if (!returnValue) _PileUpDeadTime += ( 2 * time );

    return returnValue;
}


bool Event::NeutronAdjacentCut (int iEvent, std::vector<Event> *allEvents, float time, int PID){
    /* Check if there is a neutron event adjacent to the pulse within certain time interval 
     *
     * @param:
     *        - iEvent: index of the current event in @allEvents
     *        - allEvents: all events in the root file
     *        - PID: PID of the neutron signature (4: neutron recoil, 6 nLi capture)
     *        - time: time interval
     *
     * @return:
     *        - true: if no neutron adjacent (PID requirement) within @time us
     *        - false; if there is a neutron adjancet (PID requirement) within @time us
     */
    bool prevEvent = true;
    for (long int iPrev = iEvent - 1 ; iPrev >= 0; iPrev--){
        
        Event *temp = &allEvents->at(iPrev);

        if (PID == 4){
            // Looking for event containing neutron recoil
            if (temp->hasNeutronRecoil()){
                prevEvent =
                    ( this->GetPulse( 0 ) )->TimeWindow( temp->GetPulse(temp->GetNumberOfPulses() - 1) )
                        > time;
                break;
            }
        } else if (PID == 6){
            // Looking for event containing nLi capture
            if (temp->hasNeutronCapture()){
                prevEvent =
                    this->GetPulse(0)->TimeWindow( temp->GetPulse(temp->GetNumberOfPulses() - 1)) > time;
                break;
            }
        } else {
            prevEvent =
                this->GetPulse(0)->TimeWindow( temp->GetPulse(temp->GetNumberOfPulses() - 1)) > time;
            break;
        }
    }


    bool nextEvent = true;
    for (long int iNext = iEvent + 1, iMax = allEvents->size(); iNext < iMax; iNext++){

        Event *temp = &allEvents->at(iNext);

        if (PID == 4){
            // Looking for event containing neutron recoil
            if (temp->hasNeutronRecoil()){
                nextEvent = this->GetPulse(0)->TimeWindow( temp->GetPulse(0)) > time;
                break;
            }
        }else if (PID == 6){
            // Looking for event containing nLi capture
            if (temp->hasNeutronCapture()){
                nextEvent = this->GetPulse(0)->TimeWindow( temp->GetPulse(0)) > time;
                break;
            }
        } else {
            nextEvent = this->GetPulse(0)->TimeWindow( temp->GetPulse(0)) > time;
            break;
        }
    }
    
    bool returnValue =  prevEvent && nextEvent;

    if(!returnValue) {
        if (PID == 4) _NeutronRecoilDeadTime  += ( 2 * time );
        else          _NeutronCaptureDeadTime += ( 2 * time );
    }

    return returnValue;
}

bool Event::RnPoDecayCut(int iEvent, std::vector<Event> *allEvents, float time, float height){

    Pulse_t *signalCandidate = this->GetPulse(0);
    
    bool prevEvent = true;
    for (long int iPrev = iEvent - 1; iPrev >= 0; iPrev--){
        Event *temp = &allEvents->at(iPrev);

        // single pulse neutron recoil
        if ( temp->isSinglePulse( ) == 4 ){
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
        if (temp->isSinglePulse ( ) == 4){
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
    if ( !returnValue ) _RnPoDeadTime += ( 2 * time );
    
    return returnValue;
}


bool Event::BiPoDecayCut( int iEvent, std::vector<Event> *allEvents, float time, float height ){
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

    if (!prevEvent) _BiPoDeadTime += time;

    return prevEvent;
}

bool Event::SinglePulseCut( ){ return this->isSinglePulse() != 0; }

bool Event::NeutronPulseCut( int PID ){

    assert ( PID == 4 || PID == 6 );

    if ( !this->SinglePulseCut( ) ) return false;

    float energyEvent = this->GetEnergyEvent( );

    int iHist = -1;
    if ( energyEvent < 10 && energyEvent >= 0.5 ) iHist = (int) energyEvent;

    if (iHist == -1 ) return false;
    
    PSD_t psdEnergy = PSD_per_energy[iHist];

    float nLiBandMin     = psdEnergy.nLiBand.mean     - 2 * psdEnergy.nLiBand.std;
    float nLiBandMax     = psdEnergy.nLiBand.mean     + 2 * psdEnergy.nLiBand.std;
    float neutronBandMin = psdEnergy.neutronBand.mean - 2 * psdEnergy.neutronBand.std;
    float neutronBandMax = psdEnergy.neutronBand.mean + 2 * psdEnergy.neutronBand.std;

    float pulsePSD = this->GetPulse( 0 )->PSD;

    if (PID == 4) return pulsePSD >= neutronBandMin && pulsePSD <= neutronBandMax;
    else          return pulsePSD >= nLiBandMin  && pulsePSD <= nLiBandMax;
}

double Event::RnPoCutDeadTime()           { return _RnPoDeadTime; }
double Event::BiPCutoDeadTime()           { return _BiPoDeadTime; }
double Event::PileUpCutDeadTime()         { return _PileUpDeadTime; }
double Event::MuonAdjacentCutDeadTime()   { return _MuonAdjacentDeadTime; }
double Event::NeutronRecoilCutDeadTime()  { return _NeutronRecoilDeadTime; }
double Event::NeutronCaptureCutDeadTime() { return _NeutronCaptureDeadTime; }

/* }}} */
