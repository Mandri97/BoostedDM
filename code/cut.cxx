#include "cut.hh"
#include <cmath>

using namespace std;

Cut_t::Cut_t(Event *event, Events *data){
    events = data;
    currentEvent = event;
}

Cut_t::~Cut_t(){
    delete events, currentEvent;
}

bool Cut_t::correlatedDecayRnPo(float time, float height){

    // Get pulse
    Pulse_t *pulseCandidate = currentEvent->getPulse(0);

    // Previous events
    bool prevCorrelation = true;

    for (long int iPrev = currentEvent->index - 1; iPrev >= 0; iPrev--){
        Event *prevEvent = events->getEvent(iPrev);

        // single pulse neutron recoil
        if (prevEvent->isSinglePulse() == 4){
            Pulse_t* prevPulse = prevEvent->getPulse(0);

            if (timeWindow(pulseCandidate, prevPulse) < time){

                // same height and same segment
                if (pulseCandidate->segment == prevPulse->segment){
                      prevCorrelation = heightDifference(pulseCandidate, prevPulse) > height;

                      break;
                } 
            } else break;
        }
    }

    bool nextCorrelation = true;

    for (long int iNext = currentEvent->index + 1; iNext < events->getNumberOfEvents(); iNext++){
        Event *nextEvent = events->getEvent(iNext);

        // single pulse neutron recoil
        if (nextEvent->isSinglePulse() == 4){
            Pulse_t* nextPulse = nextEvent->getPulse(0);

            if (timeWindow(pulseCandidate, nextPulse) < time){

                // same height and same segment
                if (pulseCandidate->segment == nextPulse->segment){
                    nextCorrelation = heightDifference(pulseCandidate, nextPulse) > height;

                    break;
                }
            } else break;
        }
    }

    if (prevCorrelation && nextCorrelation == false) RnPo_d += 2 * time;
    return prevCorrelation && nextCorrelation;
}

bool Cut_t::correlatedDecayBiPo(float time, float height){

    // Current pulses
    Pulse_t* pulseCandidate = currentEvent->getPulse(0);
    
    bool prevCorrelation = true,
         foundRequiredPulse = false;

    for (long int iPrev = currentEvent->index - 1; iPrev >= 0; iPrev--){
        if (foundRequiredPulse) break;

        Event *prevEvent = events->getEvent(iPrev);

        if (prevEvent->getEnergyEvent() >= 0.25 && prevEvent->getEnergyEvent() < 3.25 &&  // Gamma-rays
            prevEvent->isBetaDecayEvent() ){                                               // Beta decay only

            for (int iPulse = 0, nbPulses = prevEvent->getNumberOfPulses(); iPulse < nbPulses; iPulse++){
                Pulse_t *pulse = prevEvent->getPulse(iPulse);
                
                if (timeWindow(pulseCandidate, pulse) < time){

                    if ( pulseCandidate->segment == pulse->segment     ||
                         pulseCandidate->segment == pulse->segment - 1 ||
                         pulseCandidate->segment == pulse->segment + 1 ||
                         pulseCandidate->segment == pulse->segment - 1 - 14 ||
                         pulseCandidate->segment == pulse->segment     - 14 ||
                         pulseCandidate->segment == pulse->segment + 1 - 14 ||
                         pulseCandidate->segment == pulse->segment - 1 + 14 ||
                         pulseCandidate->segment == pulse->segment     + 14 ||
                         pulseCandidate->segment == pulse->segment + 1 + 14 ){

                        prevCorrelation = heightDifference(pulseCandidate, pulse) > height;

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

    if (prevCorrelation == false) BiPo_d += time;

    return prevCorrelation;
}

bool Cut_t::doubleFiducialCut (){
    int segment = currentEvent->getPulse(0)->segment;

    return ((segment >= 30  && segment <= 39)  ||
            (segment >= 44  && segment <= 53)  ||
            (segment >= 58  && segment <= 67)  ||
            (segment >= 72  && segment <= 81)  ||
            (segment >= 86  && segment <= 95)  ||
            (segment >= 100 && segment <= 109) ||
            (segment >= 114 && segment <= 123) );
}

bool Cut_t::neutronAdjacent (int PID, float time){

    /* Check if there is a neutron event adjacent to the pulse within certain time interval 
     *
     * @param:
     *        - PID: PID of the neutron signature (4: neutron recoil, 6 nLi capture)
     *        - time: time interval
     *
     * @return:
     *        - true: if no neutron adjacent (PID requirement) within @time us
     *        - false; if there is a neutron adjancet (PID requirement) within @time us
     */
    bool prevNeutronAdjacent = true;

    for (long int iPrev = currentEvent->index - 1; iPrev >= 0; iPrev--){
        Event *temp = events->getEvent(iPrev);

        // PID requirement
        if (PID == 4){

            // Looking for event containing neutron recoil
            if (temp->isContainingNeutronRecoil()){
                prevNeutronAdjacent = timeWindow(currentEvent->getPulse(0),
                    temp->getPulse(temp->getNumberOfPulses() - 1)) > time;

                break;
            }

        } else if (PID == 6){

            // Looking for event containing nLi capture
            if (temp->isContainingNLiCapture()){
                prevNeutronAdjacent = timeWindow(currentEvent->getPulse(0),
                    temp->getPulse(temp->getNumberOfPulses() - 1)) > time;

                break;
            }

        } else {
            prevNeutronAdjacent = timeWindow(currentEvent->getPulse(0),
                temp->getPulse(temp->getNumberOfPulses() - 1)) > time;

            break;
        }
    }

    bool nextNeutronAdjacent = true;

    long int iMax  = events->getNumberOfEvents();

    for (long int iNext = currentEvent->index + 1; iNext < iMax; iNext++){
        
        Event *temp = events->getEvent(iNext);

        // PID requirement
        if (PID == 4){

            // Looking for event containing neutron recoil
            if (temp->isContainingNeutronRecoil()){
                nextNeutronAdjacent = timeWindow(currentEvent->getPulse(0), temp->getPulse(0)) > time;

                break;
            }

        }else if (PID == 6){

            // Looking for event containing nLi capture
            if (temp->isContainingNLiCapture()){
                nextNeutronAdjacent = timeWindow(currentEvent->getPulse(0), temp->getPulse(0)) > time;

                break;
            }

        } else {
            nextNeutronAdjacent = timeWindow(currentEvent->getPulse(0), temp->getPulse(0)) > time;

            break;
        }
    }
    
    // neutron recoil
    if (PID == 4){
        if (prevNeutronAdjacent && nextNeutronAdjacent == false) nRecoilAdjacent_d += 2 * time;
    }else if (PID == 6){
        if (prevNeutronAdjacent && nextNeutronAdjacent == false) nLiAdjacent_d += 2 * time;
    }else
        if (prevNeutronAdjacent && nextNeutronAdjacent == false) neutronAdjacent_d += 2 * time;

    return prevNeutronAdjacent && nextNeutronAdjacent;
} // }}}

bool Cut_t::muonAdjacent (float time){ // {{{

    /* Check if here is a muon event adjacent to the pulse within a certain time interval
     * @param:
     *        - time: time interval
     *
     * @return:
     *        - true if no muon adjacent within @time us
     *        - false if there is a muon adjacent within @time us
     */

    bool prevMuonAdjacent = true;
    
    for ( long int iPrev = currentEvent->index - 1; iPrev >= 0; iPrev--){
        
        Event *temp = events->getEvent(iPrev);

        // Energy requirement
        if (temp->getEnergyEvent() > 15){

            prevMuonAdjacent =
                timeWindow(currentEvent->getPulse(0), temp->getPulse(temp->getNumberOfPulses() - 1)) > time;

            break;
        }
    }

    bool nextMuonAdjacent = true;

    long int iMax  = events->getNumberOfEvents();

    for (long int iNext = currentEvent->index + 1; iNext < iMax; iNext++){

        Event *temp = events->getEvent(iNext);

        // Energy requirement
        if (temp->getEnergyEvent() > 15){

            nextMuonAdjacent =
                timeWindow(currentEvent->getPulse(0), temp->getPulse(0)) > time;

            break;
        }
    }

    if (prevMuonAdjacent && nextMuonAdjacent == false) muonAdjacent_d += 2 * time;
    return prevMuonAdjacent && nextMuonAdjacent;
} // }}} 

bool Cut_t::pileUp (float time){
    /* Check if there is any event within @time microseconds
     * @param:
     *        - time: time interval 
     * Return:
     *        - true: if no event adjacent to the current event within @time microseconds
     *        - false: otherwise
     */

    bool prevPileUp = true;

    long int iPrev = currentEvent->index - 1;

    if (iPrev > -1){
        Event *temp = events->getEvent(iPrev);

        prevPileUp = timeWindow(currentEvent->getPulse(0), temp->getPulse(temp->getNumberOfPulses() - 1)) > time;
    }

    long int iNext = currentEvent->index + 1;

    bool nextPileUp = true;

    if (iNext < events->getNumberOfEvents()){
        Event *temp = events->getEvent(iNext);

        nextPileUp = timeWindow(currentEvent->getPulse(0), temp->getPulse(0)) > time;
    }

    if (prevPileUp && nextPileUp == false) pileUp_d += 2 * time;
    return prevPileUp && nextPileUp;
}

bool Cut_t::height(int height){
    return currentEvent->getPulse(0)->height < 200;
}

float Cut_t::deadTimeRnPo(){
    return RnPo_d;
}

float Cut_t::deadTimeBiPo(){
    return BiPo_d;
}

float Cut_t::deadTimeNeutronAdjacent(int PID){
    switch(PID){
        case 4:
            return nRecoilAdjacent_d;

            break;
        case 6:
            return nLiAdjacent_d;

            break;
        default:
            return neutronAdjacent_d;

            break;
    }
}

float Cut_t::deadTimeMuonAdjacent(){
    return muonAdjacent_d;
}

float Cut_t::deadTimePileUp(){
    return pileUp_d;
}

float Cut_t::heightDifference (Pulse_t *a, Pulse_t *b){ // {{{
    // return the height difference in mm between 2 pulses

    return abs(a->height - b->height);
} // }}}

float Cut_t::timeWindow (Pulse_t *a, Pulse_t *b){ // {{{
    // return the time difference in microseconds between 2 pulses
    
    return abs(a->time - b->time) * 1e-3;
} // }}}
