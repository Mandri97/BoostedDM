#include <vector>
#include <assert.h>

#ifndef __EVENT_H_
#define __EVENT_H_

/* Structure to describe a pulse */

struct Pulse_t {
    int segment,
        PID;

    float height,
          energy,
          time,
          dtime,
          PSD;
};

/* class for one event */
class Event {
    private:
        std::vector<Pulse_t> pulses;     // vector containing all pulses in one event
        float energyEvent;               // energy of the event

        bool containNeutronRecoil,
             containNLiCapture;

    public:
        Event();
        ~Event();

        void Initialize();

        // Getter
        int      getNumberOfPulses();
        float    getEnergyEvent();
        Pulse_t* getPulse(int iPulse);

        // Setter
        void addPulse(Pulse_t pulse);

        int isSinglePulse();
        bool isContainingNeutronRecoil();
        bool isContainingNLiCapture();
};

class Events {
    private:
        std::vector<Event> events;

    public:
        Events();
        ~Events();

        void Initialize();

        // Getter
        int    getNumberOfEvents();
        Event* getEvent(int iEvent);


        // Setter
        void addEvent(Event event);
};

#endif 
