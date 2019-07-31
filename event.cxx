#include "event.hh"

using namespace std;

/* class Event {{{ */

Event::Event(){
    Initialize();
}

Event::~Event(){
    Initialize();
}

void Event::Initialize(){
    pulses.clear();
    energyEvent = 0;

    containNeutronRecoil = false;
    containNLiCapture = false;
}

int Event::getNumberOfPulses(){
    return pulses.size();
}

float Event::getEnergyEvent(){
    return energyEvent;
}

Pulse_t* Event::getPulse(int iPulse){
    // Make sure that there will no segmentation fault
    assert (iPulse >= 0 && iPulse < this->getNumberOfPulses());

    return &pulses.at(iPulse);
}

void Event::addPulse(Pulse_t pulse){
    pulses.push_back(pulse);

    // Update energyEvent
    energyEvent += pulse.energy;

    if (pulse.PID == 4) containNeutronRecoil = true;
    if (pulse.PID == 6) containNLiCapture = true;
}

int Event::isSinglePulse(){
    /* Determine if the event is a single pulse
     * 
     * Return:
     *      - 0 -- if not a single pulse
     *      - PID -- otherwise
     */

    assert (this->getNumberOfPulses() > 0);

    if (this->getNumberOfPulses() > 1) return 0;
    else return pulses.at(0).PID;
}

bool Event::isContainingNeutronRecoil(){
    return containNeutronRecoil;
}

bool Event::isContainingNLiCapture(){
    return containNLiCapture;
}

/* }}} */




/* class Events {{{ */
Events::Events(){
    Initialize();
}

Events::~Events(){
    Initialize();
}

void Events::Initialize(){
    events.clear();
}

int Events::getNumberOfEvents(){
    return events.size();
}

Event* Events::getEvent(int iEvent){
    // Make sure that there will no segmentation fault
    assert (iEvent >= 0 && this->getNumberOfEvents());

    return &events.at(iEvent);
}

void Events::addEvent(Event event){
    events.push_back(event);
}

/* }}} */
