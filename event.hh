#include <vector>
#include <cstdint>
#include <assert.h>

#ifndef __EVENT_H_
#define __EVENT_H_

/* Structure to describe a pulse */

class Pulse_t {
    public:
        Pulse_t(int seg,  int PID,
                float h,  float E, float t,
                float dt, float PSD);
        Pulse_t();
        ~Pulse_t();

        float TimeWindow(Pulse_t* a);
        float HeightDifference(Pulse_t* a);

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

        bool _hasNeutronRecoil,
             _hasNeutronCapture,
             _isGammaEvent;

        // Dead Time
        double static _RnPoDeadTime,
                      _BiPoDeadTime,
                      _PileUpDeadTime,
                      _MuonAdjacentDeadTime,
                      _NeutronRecoilDeadTime,
                      _NeutronCaptureDeadTime;

    public:
        Event();
        ~Event();

        void Initialize();

        int      GetNumberOfPulses( );
        float    GetEnergyEvent( );
        double   RnPoCutDeadTime( );
        double   BiPCutoDeadTime( );
        double   PileUpCutDeadTime( );
        double   MuonAdjacentCutDeadTime( );
        double   NeutronRecoilCutDeadTime( );
        double   NeutronCaptureCutDeadTime( );


        Pulse_t* GetPulse(int iPulse);


        void AddPulse(Pulse_t pulse);


        int  isSinglePulse();
        bool hasNeutronRecoil();
        bool hasNeutronCapture();
        bool isBetaDecayEvent();


        // Cuts
        bool FiducialCut        ( );
        bool SinglePulseCut     ( );
        bool NeutronPulseCut    ( int PID );
        bool HeightCut          ( float height );
        bool MuonAdjacentCut    ( int iEvent, std::vector<Event> *allEvents, float time );
        bool PileUpCut          ( int iEvent, std::vector<Event> *allEvents, float time );
        bool NeutronAdjacentCut ( int iEvent, std::vector<Event> *allEvents, float time, int PID );
        bool RnPoDecayCut       ( int iEvent, std::vector<Event> *allEvents, float time, float height );
        bool BiPoDecayCut       ( int iEvent, std::vector<Event> *allEvents, float time, float height );
};

#endif 
