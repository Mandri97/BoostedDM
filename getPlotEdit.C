void Draw(vector<TH1F*> *hists){
    auto c = new TCanvas("c", "c", 1100, 800);

    for (int i = 0; i < hists->size(); i++){
        hists->at(i)->Draw("same hist");
        c->AddEntry(hists->at(i), hists->at(i)->GetTitle());
    }

    c->Draw();
    c->Write();
}
