void segmentVisualization(char *file){
    auto _file = new TFile(file);
    
    auto hist = (TH1F*)_file->Get("hLiveSegmentSignal");
    auto h2d  = new TH2F("h2d", ";;;", 14, 0, 14, 11, 0, 11);

    auto hNumber = new TH2F("hNummber", "", 14, 0, 14, 11, 0, 11);

    for (int i = 1; i < 155; i++){
   	int count = hist->GetBinContent(i);
	
	int segmentNumber = i - 1;
	
        int x = (segmentNumber % 14) + 1;
	int y = (segmentNumber / 14) + 1;
	
	hNumber->SetBinContent( x, y, segmentNumber );
	
        float totalRunTime   = (329.406 + 22.014) * 3600;

	h2d->SetBinContent( x, y, count/totalRunTime );
    }

    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(kViridis);

    auto c = new TCanvas("c", "c", 800, 600);
    c->SetGrid();
    c->SetLogz();

    h2d->GetZaxis()->SetTitle("Rate [Hz]");

    h2d->GetXaxis()->SetNdivisions(14);
    h2d->GetYaxis()->SetNdivisions(11);
    h2d->GetXaxis()->SetLabelSize(0);
    h2d->GetYaxis()->SetLabelSize(0);
    h2d->Draw("colz");
/*
    h2d->GetXaxis()->CenterLabels(kTRUE);
    h2d->GetYaxis()->CenterLabels(kTRUE);
    hNumber->GetXaxis()->CenterLabels(kTRUE);
    hNumber->GetYaxis()->CenterLabels(kTRUE);
*/
    hNumber->GetXaxis()->SetLabelSize(0);
    hNumber->GetYaxis()->SetLabelSize(0);
    hNumber->Draw("same Text");
}
