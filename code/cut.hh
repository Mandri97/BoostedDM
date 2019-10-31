#include "event.hh"

#include <string>
#include <map>

#include <TH1.h>

#ifndef __CUT_H__
#define __CUT_H__

#define MAX_CUT 20

class Cut {
    private:
        Event* currentEvent;
        Events* events;
        Pulse_t* pulseCandidate;

        TH1F* histogramEvent;
        TH1F* liveSegmentSignal;

        struct __args__ {
            int cutValue;
            float arg1, arg2;
            TH1F* hist;
        };

        std::map<int, __args__> cutsToBeApplied;

        int rankCut;

        // dead times
        float RnPo_d,
              BiPo_d,
              neutronAdjacent_d,
              nRecoilAdjacent_d,
              nLiAdjacent_d,
              muonAdjacent_d,
              pileUp_d;

        bool __SinglePulseCut__ (float PID1, float PID2);
        bool __RnPoDecayCut__ (float time, float height);
        bool __BiPoDecayCut__ (float time, float height);
        bool __FiducializationAndHeightCut__ (float height);
        bool __NeutronAdjacentCut__ (float PID, float time);
        bool __MuonAdjacentCut__ (float time);
        bool __PileUpCut__ (float time);
        bool __HeightCut__ (float height);

    public:
        Cut(Events *events, TH1F* hist);
        ~Cut();

        // Add cut function
        void addCut(const char* cutName, TH1F* histogram);
        void addCut(const char* cutName, float arg1, TH1F* histogram);
        void addCut(const char* cutName, float arg1, float arg2, TH1F* histogram);

        /******** Getter ********/
        float RnPoDeadTime();
        float BiPoDeadTime();
        float NeutronAdjacentDeadTime(int PID);
        float MuonAdjacentDeadTime();
        float PileUpDeadTime();
        TH1F *GetLiveSegment(const char *name);

        void Run();

        // Useful function
        float timeWindow (Pulse_t *a, Pulse_t *b);
        float heightDifference (Pulse_t *a, Pulse_t *b);
};

#endif
