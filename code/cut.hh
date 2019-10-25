#include "event.hh"

#ifndef __CUT_H__
#define __CUT_H__

class Cut_t {
    private:
        Events* events;
        Event*  currentEvent;

        // dead times
        float RnPo_d,
              BiPo_d,
              neutronAdjacent_d,
              nRecoilAdjacent_d,
              nLiAdjacent_d,
              muonAdjacent_d,
              pileUp_d;

    public:
        Cut_t(Event *event, Events *events);
        ~Cut_t();

        // IMPLEMENT AN INTERATOR FOR EVENTS
        /******** Cuts **********/
        bool correlatedDecayRnPo(float time, float height);
        bool correlatedDecayBiPo(float time, float height);
        bool doubleFiducialCut ();
        bool neutronAdjacent (int PID, float time);
        bool muonAdjacent (float time);
        bool pileUp (float time);
        bool height(int height);

        /******** Getter ********/
        float deadTimeRnPo();
        float deadTimeBiPo();
        float deadTimeNeutronAdjacent(int PID);
        float deadTimeMuonAdjacent();
        float deadTimePileUp();


        // Useful function
        float timeWindow (Pulse_t *a, Pulse_t *b);
        float heightDifference (Pulse_t *a, Pulse_t *b);
};

#endif
