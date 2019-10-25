#include <vector>
#include <cstdint>
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
             containNLiCapture,
             allGammaRayPulses;

    protected:
        int64_t index;

    public:
        Event();
        ~Event();

        void Initialize();

        /******* Getter ********/
        float    getEnergyEvent();

            // Pulse realted
        int      getNumberOfPulses();
        Pulse_t* getPulse(int iPulse);

        /******** Setter ********/
        void addPulse(Pulse_t pulse);

        int isSinglePulse();
        bool isContainingNeutronRecoil();
        bool isContainingNLiCapture();
        bool isBetaDecayEvent();

        // friend class
        friend class Events;
        friend class Cut_t;
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
