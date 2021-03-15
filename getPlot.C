void Draw(vector<TH1F*> *hists);

TH1F* setupHistogram(TFile* file, char *name, Color_t color, float scaling);

void getPlot(char* rootFile){
    // open rootfile
    auto _file = new TFile(rootFile, "read");

    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptTitle(kFALSE);

    /* Scaling coefficient {{{ */
        
        auto hLiveSegment = (TH1F*)_file->Get("hLiveSegment");
        auto hLiveSegmentFD = (TH1F*)_file->Get("hLiveSegmentFiducial");

        int segmentAllEvent  = 0,
            segmentDoubleCut = 0;

        for (int i = 1; i < 155; i++){
            if (hLiveSegment->GetBinContent(i)) segmentAllEvent++;
            if (hLiveSegmentFD->GetBinContent(i)) segmentDoubleCut++;
        }

        // scale = Hz / MeV / kg
        float totalRunTime   = (329.406 + 22.014) * 3600, //202394 + 130341;
              binWidth       = 0.025,
              scaleAllEvent  = 1 / (totalRunTime * binWidth * segmentAllEvent * 24.2),
              scaleDoubleCut = 1 / (totalRunTime * binWidth * segmentDoubleCut * (40/117.8) * 24.2);
    /* }}} */


    /* All events {{{ */
    // Black
        auto hAllEvent = (TH1F*)_file->Get("hEnergyPerEvent");
        hAllEvent->SetLineColor(kBlack);
        hAllEvent->Scale(scaleAllEvent);
    /* }}} */


    /* Single Pulse Event {{{ */
    // Blue
        auto hSinglePulse = (TH1F*)_file->Get("hSinglePulseEvent");
        hSinglePulse->SetLineColor(kOrange + 3);
        hSinglePulse->Scale(scaleAllEvent);
    /* }}} */
    

    /* Signal {{{ */
    // Red
        auto hSignalCandidate = (TH1F*)_file->Get("hSignalCandidate");
        hSignalCandidate->SetLineColor(kOrange - 2);
        hSignalCandidate->Scale(scaleAllEvent);
    /* }}} */
    

    /* Segment and height double fiducialization {{{ */
    // Purple
        auto hFiducialization = (TH1F*)_file->Get("hFiducialization");
        hFiducialization->SetLineColor(kYellow - 3);
        hFiducialization->Scale(scaleDoubleCut);
    /* }}} */


    /* Muon Adjacent veto {{{ */
    // kTeal +4
        auto hMuonVeto = (TH1F*)_file->Get("hMuonAdjacent");
        hMuonVeto->SetLineColor(kMagenta + 2);
        hMuonVeto->Scale(scaleDoubleCut);
    /* }}} */


    /* Neutron recoil Adjacent veto {{{ */
    // kCyan - 3
        auto hRecoiVeto = (TH1F*)_file->Get("hNeutronRecoil");
        hRecoiVeto->SetLineColor(kTeal + 4);
        hRecoiVeto->Scale(scaleDoubleCut);
    /* }}} */


    /* NLi Capture Adjacent veto {{{ */
    // kYellow - 3
        auto hCaptureVeto = (TH1F*)_file->Get("hNeutronCapture");
        hCaptureVeto->SetLineColor(kCyan - 6);
        hCaptureVeto->Scale(scaleDoubleCut);
    /* }}} */


    /* Pile Up {{{ */
    // kOrange - 3
        auto hPileUp = (TH1F*)_file->Get("hPileUp");
		hPileUp->SetTitle("Pile Up veto 4 #mus");
        hPileUp->SetLineColor(kRed - 7);
        hPileUp->Scale(scaleDoubleCut);
    /* }}} */


    /* Rn - Po Correlated Decay {{{ */
       auto hRnPoVeto = (TH1F*)_file->Get("hRnPoDecay") ;
       hRnPoVeto->Scale(scaleDoubleCut);
    /* }}} */


    /* Bi - Po Correlated Decay {{{ */
        auto hBiPoVeto = (TH1F*)_file->Get("hBiPoDecay");
        hBiPoVeto->Scale(scaleDoubleCut);
        hBiPoVeto->SetLineColor(kPink - 2);
    /* }}} */

    hAllEvent->GetYaxis()->SetRangeUser(1e-7, 0.8);
    hAllEvent->GetYaxis()->SetTitle("Rate [Hz/MeV/kg]");
    hAllEvent->GetYaxis()->SetTitleSize(0.03);
    hAllEvent->GetYaxis()->SetTitleOffset(1.55);
    hAllEvent->GetYaxis()->CenterTitle(kTRUE);
    hAllEvent->GetYaxis()->SetLabelSize(0.03);

    hAllEvent->GetXaxis()->SetRangeUser(0.5, 10);
    hAllEvent->GetXaxis()->SetTitle("Visible Energy [MeV]");
    hAllEvent->GetXaxis()->SetTitleSize(0.03);
    hAllEvent->GetXaxis()->CenterTitle(kTRUE);
    hAllEvent->GetXaxis()->SetLabelSize(0.03);

    hSignalCandidate->GetYaxis()->SetRangeUser(1e-7, 0.8);
    hSignalCandidate->GetYaxis()->SetTitle("Rate [Hz/MeV/kg]");
    hSignalCandidate->GetYaxis()->SetTitleSize(0.03);
    hSignalCandidate->GetYaxis()->SetTitleOffset(1.55);
    hSignalCandidate->GetYaxis()->CenterTitle(kTRUE);
    hSignalCandidate->GetYaxis()->SetLabelSize(0.03);

    hSignalCandidate->GetXaxis()->SetRangeUser(0.5, 10);
    hSignalCandidate->GetXaxis()->SetTitle("Visible Energy [MeV]");
    hSignalCandidate->GetXaxis()->SetTitleSize(0.03);
    hSignalCandidate->GetXaxis()->CenterTitle(kTRUE);
    hSignalCandidate->GetXaxis()->SetLabelSize(0.03);
    vector<TH1F*> toDraw;

    toDraw.push_back(hAllEvent);
    toDraw.push_back(hSinglePulse);
    toDraw.push_back(hSignalCandidate);
    toDraw.push_back(hFiducialization);
    toDraw.push_back(hMuonVeto);
    toDraw.push_back(hRecoiVeto);
    toDraw.push_back(hCaptureVeto);
    toDraw.push_back(hPileUp);
    toDraw.push_back(hRnPoVeto);
    toDraw.push_back(hBiPoVeto);

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
