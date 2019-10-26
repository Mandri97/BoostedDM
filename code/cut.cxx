#include "cut.hh"

#include <cmath>
#include <cstddef>
#include <assert.h>

using namespace std;

map<string, int> cutNameValue = {
    {"singlePulseEvent", 0},
    {"RnPoDecay", 1},
    {"BiPoDecay", 2},
    {"fiducialAndHeight", 3},
    {"neutronAdjacent", 4},
    {"muonAdjacent", 5},
    {"pileUp", 6}
};

Cut::Cut(Events *data, TH1F* hist){
    events = data;
    histogramEvent = hist;

    currentEvent = nullptr;
    pulseCandidate = nullptr;

    RnPo_d = .0; BiPo_d = .0; neutronAdjacent_d = .0; nRecoilAdjacent_d = .0;
    nLiAdjacent_d = .0; muonAdjacent_d = .0; pileUp_d = .0;

    rankCut = 0;

    cutsToBeApplied.clear();
}

Cut::~Cut(){
    cutsToBeApplied.clear();
    delete events; 
}

bool Cut::__SinglePulseCut__(float PID1, float PID2){

    if (currentEvent->isSinglePulse() == (int) PID1 || currentEvent->isSinglePulse() == (int) PID2){
        pulseCandidate = currentEvent->getPulse(0);

        return true;
    }

    return false;
}

bool Cut::__RnPoDecayCut__(float time, float height){
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

bool Cut::__BiPoDecayCut__(float time, float height){
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

bool Cut::__FiducializationAndHeightCut__ (float height){
    int segment = pulseCandidate->segment;

    bool F = ((segment >= 30  && segment <= 39)  ||
              (segment >= 44  && segment <= 53)  ||
              (segment >= 58  && segment <= 67)  ||
              (segment >= 72  && segment <= 81)  ||
              (segment >= 86  && segment <= 95)  ||
              (segment >= 100 && segment <= 109) ||
              (segment >= 114 && segment <= 123) );

    return F && this->__HeightCut__(height);
}

bool Cut::__NeutronAdjacentCut__ (float PID, float time){

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

    PID = (int) PID;

    bool prevNeutronAdjacent = true;

    for (long int iPrev = currentEvent->index - 1; iPrev >= 0; iPrev--){
        Event *temp = events->getEvent(iPrev);

        // PID requirement
        if (PID == 4){

            // Looking for event containing neutron recoil
            if (temp->isContainingNeutronRecoil()){
                prevNeutronAdjacent = timeWindow(pulseCandidate,
                    temp->getPulse(temp->getNumberOfPulses() - 1)) > time;

                break;
            }

        } else if (PID == 6){

            // Looking for event containing nLi capture
            if (temp->isContainingNLiCapture()){
                prevNeutronAdjacent = timeWindow(pulseCandidate,
                    temp->getPulse(temp->getNumberOfPulses() - 1)) > time;

                break;
            }

        } else {
            prevNeutronAdjacent = timeWindow(pulseCandidate,
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
                nextNeutronAdjacent = timeWindow(pulseCandidate, temp->getPulse(0)) > time;

                break;
            }

        }else if (PID == 6){

            // Looking for event containing nLi capture
            if (temp->isContainingNLiCapture()){
                nextNeutronAdjacent = timeWindow(pulseCandidate, temp->getPulse(0)) > time;

                break;
            }

        } else {
            nextNeutronAdjacent = timeWindow(pulseCandidate, temp->getPulse(0)) > time;

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

bool Cut::__MuonAdjacentCut__ (float time){ // {{{

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
                timeWindow(pulseCandidate, temp->getPulse(temp->getNumberOfPulses() - 1)) > time;

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
                timeWindow(pulseCandidate, temp->getPulse(0)) > time;

            break;
        }
    }

    if (prevMuonAdjacent && nextMuonAdjacent == false) muonAdjacent_d += 2 * time;
    return prevMuonAdjacent && nextMuonAdjacent;
} // }}} 

bool Cut::__PileUpCut__ (float time){
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

        prevPileUp = timeWindow(pulseCandidate, temp->getPulse(temp->getNumberOfPulses() - 1)) > time;
    }

    long int iNext = currentEvent->index + 1;

    bool nextPileUp = true;

    if (iNext < events->getNumberOfEvents()){
        Event *temp = events->getEvent(iNext);

        nextPileUp = timeWindow(pulseCandidate, temp->getPulse(0)) > time;
    }

    if (prevPileUp && nextPileUp == false) pileUp_d += 2 * time;

    return prevPileUp && nextPileUp;
}

void Cut::Run(){
    for (long int iEntry = 0, nEntries = events->getNumberOfEvents(); iEntry < nEntries; ++iEntry){
        currentEvent = events->getEvent(iEntry); 

        float energyEvent = currentEvent->getEnergyEvent();
        histogramEvent->Fill(energyEvent);

        if (energyEvent > 10 || energyEvent < 0) continue;

        bool satisfyCut = true;

        for (auto i = cutsToBeApplied.begin(), j = cutsToBeApplied.end(); i != j; ++i){
            auto cutParameters = i->second;

            int cutValue = cutParameters.cutValue;
            float arg1 = cutParameters.arg1;
            float arg2 = cutParameters.arg2;
            TH1F *histogram = cutParameters.hist;

            switch(cutValue){
                case 0:
                    // SinglePulseEvent
                    assert(arg1 != -1 && arg2 != -1);
                    satisfyCut *= this->__SinglePulseCut__(arg1, arg2);

                    break;
                case 1:
                    // RnPoDecay
                    assert(arg1 != -1 && arg2 != -1);
                    satisfyCut *= this->__RnPoDecayCut__(arg1, arg2);

                    break;
                case 2:
                    // BiPoDecay
                    assert(arg1 != -1 && arg2 != -1);
                    satisfyCut *= this->__BiPoDecayCut__(arg1, arg2);

                    break;
                case 3:
                    // fiducialAndHeight
                    assert(arg1 != -1);
                    satisfyCut *= this->__FiducializationAndHeightCut__(arg1);

                    break;
                case 4:
                    // neutronAdjacent
                    assert(arg1 != -1 && arg2 != -1);
                    satisfyCut *= this->__NeutronAdjacentCut__(arg1, arg2);

                    break;
                case 5:
                    // muonAdjacent
                    assert(arg1 != -1);
                    satisfyCut *= this->__MuonAdjacentCut__(arg1);

                    break;
                case 6:
                    // PileUp
                    assert(arg1 != -1);
                    satisfyCut *= this->__PileUpCut__(arg1);

                    break;
            }

            if (satisfyCut) histogram->Fill(currentEvent->getEnergyEvent());
            else break;
        }
         
    }
}

bool Cut::__HeightCut__(float height){
    return abs(pulseCandidate->height) < (int) height;
}

float Cut::RnPoDeadTime(){
    return RnPo_d;
}

float Cut::BiPoDeadTime(){
    return BiPo_d;
}

float Cut::NeutronAdjacentDeadTime(int PID){
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

float Cut::MuonAdjacentDeadTime(){
    return muonAdjacent_d;
}

float Cut::PileUpDeadTime(){
    return pileUp_d;
}

float Cut::heightDifference (Pulse_t *a, Pulse_t *b){ // {{{
    // return the height difference in mm between 2 pulses

    return abs(a->height - b->height);
} // }}}

float Cut::timeWindow (Pulse_t *a, Pulse_t *b){ // {{{
    // return the time difference in microseconds between 2 pulses
    
    return abs(a->time - b->time) * 1e-3;
} // }}}


void Cut::addCut(const char* cutName, TH1F* histogram){
    // Verify that the cut exist
    bool nameExist = false;

    string cutName_s = string(cutName);

    for (auto name : cutNameValue){
        if (name.first == cutName_s){
            nameExist = true;
            break;
        }
    }

    assert (nameExist);

    cutsToBeApplied[rankCut] = {cutNameValue[cutName_s], -1, -1, histogram};

    rankCut++;
}

void Cut::addCut(const char* cutName, float arg1, TH1F* histogram){
    this->addCut(cutName, histogram);

    cutsToBeApplied[rankCut - 1].arg1 = arg1;
}

 void Cut::addCut(const char* cutName, float arg1, float arg2, TH1F* histogram){
    this->addCut(cutName, arg1, histogram);

    cutsToBeApplied[rankCut - 1].arg2 = arg2;
}
