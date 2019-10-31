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

    if (currentEvent->isSinglePulse() == 0) return false;

    // Using default PID
    /*
    if (currentEvent->isSinglePulse() == (int) PID1 || currentEvent->isSinglePulse() == (int) PID2){
        pulseCandidate = currentEvent->getPulse(0);

        return true;
    }
    */

    // Using new PID -> 2 STD
    int iHist = (int) currentEvent->getEnergyEvent();

    PSD_t psdEnergy = PSD_per_energy[iHist];

    auto pulse = currentEvent->getPulse(0);

    float neutronBandMin = psdEnergy.neutronBand.mean - 2 * psdEnergy.neutronBand.std;
    float neutronBandMax = psdEnergy.neutronBand.mean + 2 * psdEnergy.neutronBand.std;

    float nLiBandMin = psdEnergy.nLiBand.mean - 2 * psdEnergy.nLiBand.std;
    float nLiBandMax = psdEnergy.nLiBand.mean + 2 * psdEnergy.nLiBand.std;

    if ((pulse->PSD >= neutronBandMin && pulse->PSD <= neutronBandMax) ||
        (pulse->PSD >= nLiBandMin     && pulse->PSD <= nLiBandMax    ) ){

        pulseCandidate = pulse;
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

    bool returnValue = prevCorrelation && nextCorrelation;

    if (!returnValue) RnPo_d += (2 * time);

    return returnValue;
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

    if (!prevCorrelation) BiPo_d += time;

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
    
    bool returnValue = prevNeutronAdjacent && nextNeutronAdjacent;

    if (!returnValue){
        if (PID == 4) nRecoilAdjacent_d += (2 * time);
        else if (PID == 6) nLiAdjacent_d += (2 * time);
        else  neutronAdjacent_d += (2 * time);
    }

    return returnValue;
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

    bool returnValue = prevMuonAdjacent && nextMuonAdjacent;

    if (!returnValue) muonAdjacent_d += (2 * time);

    return returnValue;
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

    bool returnValue = prevPileUp && nextPileUp;
    
    if (!returnValue) pileUp_d += (2 * time);

    return returnValue;
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

        if (satisfyCut) liveSegmentSignal->Fill(pulseCandidate->segment);
         
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

TH1F* Cut::GetLiveSegment(){
    return liveSegmentSignal;
}
