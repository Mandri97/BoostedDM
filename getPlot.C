void Draw(vector<TH1F*> *hists);

TH1F* setupHistogram(TFile* file, char *name, Color_t color, float scaling);

void getPlot(char* rootFile){
    // open rootfile
    auto _file = new TFile(rootFile, "read");

    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptTitle(kFALSE);

    /* Scaling coefficient {{{ */
        
        auto hist_liveSegment = (TH1F*)_file->Get("hLiveSegment");
        auto hist_liveSegment_Segment_z_DoubleFV = (TH1F*)_file->Get("hLiveSegmentFiducial");

        int segmentAllEvent  = 0,
            segmentDoubleCut = 0;

        for (int i = 1; i < 155; i++){
            if (hist_liveSegment->GetBinContent(i)) segmentAllEvent++;
            if (hist_liveSegment_Segment_z_DoubleFV->GetBinContent(i)) segmentDoubleCut++;
        }

        // scale = Hz / MeV / kg
        float totalRunTime   = (329.406 + 22.014) * 3600, //202394 + 130341;
              binWidth       = 0.025,
              scaleAllEvent  = 1 / (totalRunTime * binWidth * segmentAllEvent * 24.2),
              scaleDoubleCut = 1 / (totalRunTime * binWidth * segmentDoubleCut * (40/117.8) * 24.2);
    /* }}} */


    /* All events {{{ */
    // Black
        auto hist_EnergyPerEvent = (TH1F*)_file->Get("hEnergyPerEvent");
        hist_EnergyPerEvent->SetLineColor(kBlack);
        hist_EnergyPerEvent->Scale(scaleAllEvent);
    /* }}} */


    /* Single Pulse Event {{{ */
    // Blue
        auto hist_SinglePulseEvent = (TH1F*)_file->Get("hSinglePulseEvent");
        hist_SinglePulseEvent->SetLineColor(kOrange + 3);
        hist_SinglePulseEvent->Scale(scaleAllEvent);
    /* }}} */
    

    /* Signal {{{ */
    // Red
        auto hist_Signal = (TH1F*)_file->Get("hSignalCandidate");
        hist_Signal->SetLineColor(kOrange - 2);
        hist_Signal->Scale(scaleAllEvent);
    /* }}} */
    

    /* Segment and height double fiducialization {{{ */
    // Purple
        auto hist_Signal_Segmet_z_DoubleFV = (TH1F*)_file->Get("hFiducialization");
        hist_Signal_Segmet_z_DoubleFV->SetLineColor(kYellow - 3);
        hist_Signal_Segmet_z_DoubleFV->Scale(scaleDoubleCut);
    /* }}} */


    /* Muon Adjacent veto {{{ */
    // kTeal +4
        auto hist_Signal_MuonAdjacent_time5 = (TH1F*)_file->Get("hMuonAdjacent");
        hist_Signal_MuonAdjacent_time5->SetLineColor(kMagenta + 2);
        hist_Signal_MuonAdjacent_time5->Scale(scaleDoubleCut);
    /* }}} */


    /* Neutron recoil Adjacent veto {{{ */
    // kCyan - 3
        auto hist_Signal_NeutronRecoilAdjacent_time5 = (TH1F*)_file->Get("hNeutronRecoil");
        hist_Signal_NeutronRecoilAdjacent_time5->SetLineColor(kTeal + 4);
        hist_Signal_NeutronRecoilAdjacent_time5->Scale(scaleDoubleCut);
    /* }}} */


    /* NLi Capture Adjacent veto {{{ */
    // kYellow - 3
        auto hist_Signal_NLiCaptureAdjacent_time400 = (TH1F*)_file->Get("hNeutronCapture");
        hist_Signal_NLiCaptureAdjacent_time400->SetLineColor(kCyan - 6);
        hist_Signal_NLiCaptureAdjacent_time400->Scale(scaleDoubleCut);
    /* }}} */


    /* Pile Up {{{ */
    // kOrange - 3
        auto hist_Signal_PileUp_time2 = (TH1F*)_file->Get("hPileUp");
	hist_Signal_PileUp_time2->SetTitle("Pile Up veto 4 #mus");
        hist_Signal_PileUp_time2->SetLineColor(kRed - 7);
        hist_Signal_PileUp_time2->Scale(scaleDoubleCut);
    /* }}} */


    /* Rn - Po Correlated Decay {{{ */
       auto hist_Signal_RnPoCorrelatedDecay = (TH1F*)_file->Get("hRnPoDecay") ;
       hist_Signal_RnPoCorrelatedDecay->Scale(scaleDoubleCut);
    /* }}} */


    /* Bi - Po Correlated Decay {{{ */
        auto hist_Signal_BiPoCorrelatedDecay = (TH1F*)_file->Get("hBiPoDecay");
        hist_Signal_BiPoCorrelatedDecay->Scale(scaleDoubleCut);
        hist_Signal_BiPoCorrelatedDecay->SetLineColor(kPink - 2);
    /* }}} */

    hist_EnergyPerEvent->GetYaxis()->SetRangeUser(1e-7, 0.8);
    hist_EnergyPerEvent->GetYaxis()->SetTitle("Rate [Hz/MeV/kg]");
    hist_EnergyPerEvent->GetYaxis()->SetTitleSize(0.03);
    hist_EnergyPerEvent->GetYaxis()->SetTitleOffset(1.55);
    hist_EnergyPerEvent->GetYaxis()->CenterTitle(kTRUE);
    hist_EnergyPerEvent->GetYaxis()->SetLabelSize(0.03);

    hist_EnergyPerEvent->GetXaxis()->SetRangeUser(0.5, 10);
    hist_EnergyPerEvent->GetXaxis()->SetTitle("Visible Energy [MeV]");
    hist_EnergyPerEvent->GetXaxis()->SetTitleSize(0.03);
    hist_EnergyPerEvent->GetXaxis()->CenterTitle(kTRUE);
    hist_EnergyPerEvent->GetXaxis()->SetLabelSize(0.03);

    hist_Signal->GetYaxis()->SetRangeUser(1e-7, 0.8);
    hist_Signal->GetYaxis()->SetTitle("Rate [Hz/MeV/kg]");
    hist_Signal->GetYaxis()->SetTitleSize(0.03);
    hist_Signal->GetYaxis()->SetTitleOffset(1.55);
    hist_Signal->GetYaxis()->CenterTitle(kTRUE);
    hist_Signal->GetYaxis()->SetLabelSize(0.03);

    hist_Signal->GetXaxis()->SetRangeUser(0.5, 10);
    hist_Signal->GetXaxis()->SetTitle("Visible Energy [MeV]");
    hist_Signal->GetXaxis()->SetTitleSize(0.03);
    hist_Signal->GetXaxis()->CenterTitle(kTRUE);
    hist_Signal->GetXaxis()->SetLabelSize(0.03);
    vector<TH1F*> toDraw;

    toDraw.push_back(hist_EnergyPerEvent);
    toDraw.push_back(hist_SinglePulseEvent);
    toDraw.push_back(hist_Signal);
    toDraw.push_back(hist_Signal_Segmet_z_DoubleFV);
    toDraw.push_back(hist_Signal_MuonAdjacent_time5);
    toDraw.push_back(hist_Signal_NeutronRecoilAdjacent_time5);
    toDraw.push_back(hist_Signal_NLiCaptureAdjacent_time400);
    toDraw.push_back(hist_Signal_PileUp_time2);
    toDraw.push_back(hist_Signal_RnPoCorrelatedDecay);
    toDraw.push_back(hist_Signal_BiPoCorrelatedDecay);

    Draw(&toDraw);
}


void Draw(vector<TH1F*> *hists){
    auto c = new TCanvas("c", "c", 1100, 800);
    c->SetLogy();

    // dx = 0.28, dy = 0.19
    auto l = new TLegend(0.4, 0.625, 0.6, 0.895);
    l->SetTextSize(0.025);
    l->SetLineWidth(0);

    for (int i = 0; i < hists->size(); i++){
        hists->at(i)->Draw("same hist");
        l->AddEntry(hists->at(i), hists->at(i)->GetTitle());
    }

    l->Draw();
}


TH1F* setupHistogram(TFile* file, char *name, Color_t color, float scaling){
    auto temp = (TH1F*) file->Get(name);
    temp->SetLineColor(color);
    temp->Scale(scaling);

    return temp;
}
