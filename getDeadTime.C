void getDeadTime (const char* file){
    TFile *_file = new TFile(file);

    auto fiducialAndHeight = ((TH1F*) _file->Get("hFiducialization"))->GetEntries();
    auto muon = ((TH1F*) _file->Get("hMuonAdjacent"))->GetEntries();
    auto neutronRecoil = ((TH1F*) _file->Get("hNeutronRecoil"))->GetEntries();
    auto nLi = ((TH1F*) _file->Get("hNeutronCapture"))->GetEntries();
    auto pileUp = ((TH1F*) _file->Get("hPileUp"))->GetEntries();
    auto RnPo = ((TH1F*) _file->Get("hRnPoDecay"))->GetEntries();
    auto BiPo = ((TH1F*) _file->Get("hBiPoDecay"))->GetEntries();
    
    auto liveSegment__ = ((TH1F*) _file->Get("hLiveSegmentFiducial"));
    
    int liveSegment= 0;

    for (int i = 0; i < liveSegment_->GetNbinsX();i++)
        if(liveSegment__->GetBinContent(i)) liveSegment++;

    cout << "Live segment after double fiducialization cut: " << liveSegment << endl;

    float RnfractionVolume = (0.145 * 0.145 * 0.5) / (liveSegment * 0.145 * 0.145 * 1.176);
    float BifractionVolume =  9 * RnfractionVolume;

    auto MuonDeadTime = (fiducialAndHeight - muon) * 2 * 5 * 1e-6;
    auto NeutronRecoilDeadTime = (muon - neutronRecoil) * 2 * 5 * 1e-6;
    auto NLiDeadTime = (neutronRecoil - nLi) * 2 * 400 * 1e-6;
    auto PileUpDeadTime = (nLi - pileUp) * 2 * 2 * 1e-6;
    auto RnPoDeadTime = (pileUp - RnPo) * 2 * 15000 * 1e-6 * RnfractionVolume;
    auto BiPoDeadTime = (RnPo - BiPo) * 1200 * 1e-6 * BifractionVolume;

    cout << "Dead time" << endl;
    cout << "\tMuon: " << MuonDeadTime << " s, " << MuonDeadTime * 100 / 79253.3 << endl;
    cout << "\tNeutron Recoil: " << NeutronRecoilDeadTime << " s, " << NeutronRecoilDeadTime * 100/79253.3 << endl;
    cout << "\tNLi: " << NLiDeadTime << " s, " << NLiDeadTime * 100 / 79253.3 << endl;
    cout << "\tPile up: " << PileUpDeadTime << " s, " << PileUpDeadTime * 100 / 79253.3 << endl;
    cout << "\tRnPo: " << RnPoDeadTime << " s, " << RnPoDeadTime * 100 / 79253.3 << endl;
    cout << "\tBiPo: " << BiPoDeadTime << " s, " << BiPoDeadTime * 100 / 79253.3 << endl;

    cout << endl;

    auto total = (MuonDeadTime + NeutronRecoilDeadTime + NLiDeadTime + PileUpDeadTime + RnPoDeadTime +
                 BiPoDeadTime) * 100 / 79253.3;

    cout << "Total inefficiency: " << total << endl;
}
