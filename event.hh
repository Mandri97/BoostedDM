#include <vector>
#include <cstdint>
#include <assert.h>
#include <cmath>

#include <TH1.h>

#ifndef __EVENT_H_
#define __EVENT_H_

/* Structure to describe a pulse */

class Pulse_t {
    public:
        Pulse_t(int seg,  int PID, double h,  
		double E, double t, double dt, double PSD);
        Pulse_t();
        ~Pulse_t();

        double TimeWindow(Pulse_t* a);
        double HeightDifference(Pulse_t* a);

        int segment,
            PID;

        double height,
               energy,
               time,
               dtime,
                PSD;
};

/* class for one event */
class Event {
    private:
        std::vector<Pulse_t> pulses;     // vector containing all pulses in one event
        double energyEvent;               // energy of the event

        bool hasNeutronRecoil,
             hasNeutronCapture;

		bool isMuonEvent,
             isGammaEvent;

		bool SearchEventInTime( int iEvent, std::vector<Event> *allEvents, double time, double minE, int PID);

    public:
        Event();
        ~Event();

        void Initialize();

        int      GetNumberOfPulses( );
        double   GetEnergyEvent( );

        Pulse_t* GetPulse(int iPulse);

        void AddPulse(Pulse_t pulse);

        int  IsSinglePulse();
        bool IsBetaDecayEvent();
		bool IsMuonEvent();
		bool IsGammaEvent();
		bool IsRecoilEvent();
		bool IsCaptureEvent();

        // Cuts
        bool FiducialCut    ( );
        bool SinglePulseCut ( );
        bool NeutronCut     ( TH1D *hMin, TH1D *hMax );
        bool HeightCut      ( double height );
        bool MuonVeto       ( int iEvent, std::vector<Event> *allEvents, double time );
        bool PileUpVeto     ( int iEvent, std::vector<Event> *allEvents, double time );
		bool RecoilVeto     ( int iEvent, std::vector<Event> *allEvents, double time);
		bool CaptureVeto    ( int iEvent, std::vector<Event> *allEvents, double time);

		// Broken function
        bool RnPoDecayCut   ( int iEvent, std::vector<Event> *allEvents, double time, double height );
        bool BiPoDecayCut   ( int iEvent, std::vector<Event> *allEvents, double time, double height );
};

#endif 
