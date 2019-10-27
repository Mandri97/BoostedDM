void getDeadTime.C(const char* file){
    TFile *_file = new TFile(file);

    auto fiducialAndHeight = ((TH1F*) _file->Get("hist_Signal_Segmet_z_DoubleFV"))->GetEntries();
    auto muon = ((TH1F*) _file->Get("hist_Signal_MuonAdjacent_time5"))->GetEntries();
    auto neutronRecoil = ((TH1F*) _file->Get("hist_Signal_NeutronRecoilAdjacent_time5"))->GetEntries();
    auto nLi = ((TH1F*) _file->Get("hist_Signal_NLiCaptureAdjacent_time400"))->GetEntries();
    auto pileUp = ((TH1F*) _file->Get("hist_Signal_PileUp_time2"))->GetEntries();
    auto RnPo = ((TH1F*) _file->Get("hist_Signal_RnPoCorrelatedDecay"))->GetEntries();
    auto BiPo = ((TH1F*) _file->Get("hist_Signal_BiPoCorrelatedDecay"))->GetEntries();

    auto MuonDeadTime = (fiducialAndHeight - muon) * 2 * 5 * 1e-6;
    auto NeutronRecoilDeadTime = (muon - neutronRecoil) * 2 * 5 * 1e-6;
    auto NLiDeadTime = (neutronRecoil - nLi) * 2 * 400 * 1e-6;
    auto PileUpDeadTime = (nLi - pileUp) * 2 * 2 * 1e-6;
    auto RnPoDeadTime = (pileUp - RnPo) * 2 * 15000 * 1e-6;
    auto BiPoDeadTime = (RnPo - BiPo) * 1200 * 1e-6;

    cout << "Dead time" << endl;
    cout << "\tMuon: " << MuonDeadTime << " s, " << MuonDeadTim * 100 / 79253.3 << endl;
    cout << "\tNeutron Recoil: " << NeutronRecoilDeadTime << " s, " << NeutronRecoilDeadTime * 100/79253.3 << endl;
    cout << "\tNLi: " << NLiDeadTime << " s, " << NLiDeadTime * 100 / 79253.3 << endl;
    cout << "\tPile up: " << PileUpDeadTime << " s, " << PileUpDeadTime * 100 / 79253.3 << endl;
    cout << "\tRnPo: " << RnPoDeadTime << " s, " << RnPoDeadTime * 100 / 79253.3 << endl;
    cout << "\tBiPo: " << BiPoDeadTime << " s, " << BiPoDeadTime * 100 / 79253.3 << endl;
}
