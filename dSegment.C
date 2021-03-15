void dSegment(){
    auto _fileBackground = new TFile("Result/background14MeV.root");
    auto _fileSignal = new TFile("Result/resultContinuousPSD.root");

	auto _TBackground = (TTree*)_fileBackground->Get("tSelectedEvents");
	auto _TSignal = (TTree*)_fileSignal->Get("tSelectedEvents");

	auto hBackground = new TH1F("hBackground", "", 154, 0, 154);
	auto hSignal = new TH1F("hSignal", "", 154, 0, 154);

	_TBackground->Draw("seg>>hBackground", "E>=3 && E<10", "goff");
	_TSignal->Draw("seg>>hSignal", "E>=3 && E<10", "goff");

    auto h2d  = new TH2F("h2d", ";;;", 14, 0, 14, 11, 0, 11);

    auto hNumber = new TH2F("hNummber", "", 14, 0, 14, 11, 0, 11);

	hBackground->Scale(1.0792);
	hSignal->Scale(1.0792);

    for (int i = 1; i < 155; i++){
		int count = hBackground->GetBinContent(i);
		
		int segmentNumber = i - 1;

		int x = (segmentNumber % 14) + 1;
		int y = (segmentNumber / 14) + 1;
		
		hNumber->SetBinContent( x, y, segmentNumber );
		
		//if (segmentNumber == 136) continue;

		float totalRunTime   = 1261513.745632;

		h2d->SetBinContent( x, y, count/totalRunTime );
    }
	
	for (int i = 1; i < 155; i++){
		int count = hSignal->GetBinContent(i);
		
		if (count < 1) continue;

		int segmentNumber = i - 1;
		
		int x = (segmentNumber % 14) + 1;
		int y = (segmentNumber / 14) + 1;
	
		hNumber->SetBinContent( x, y, segmentNumber );
		
		float totalRunTime   = 1261513.745632;

		h2d->SetBinContent( x, y, count/totalRunTime );
    }

	//printf("aaaaaaaaaa %f\n", hBackground->GetBinContent(153));

	gStyle->SetOptStat(kFALSE);
	gStyle->SetNumberContours(255);
    gStyle->SetPalette(kViridis);

	//h2d->Scale(1e3); // convert to mHz

    auto c = new TCanvas("c", "c", 800, 600);
    c->SetGrid();
    c->SetLogz();
	
	TGaxis::SetMaxDigits(3);

    h2d->GetZaxis()->SetTitle("Rate [Hz]");
	//h2d->GetZaxis()->SetRangeUser(4e-4, 5.5e-2);

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

	TLine *l1 = new TLine(2, 3, 2, 9);
	TLine *l2 = new TLine(2, 3, 12, 3);
	TLine *l3 = new TLine(12, 3, 12, 9);
	TLine *l4 = new TLine(2, 9, 12, 9);

	l1->SetLineWidth(3);
	l2->SetLineWidth(3);
	l3->SetLineWidth(3);
	l4->SetLineWidth(3);

	l1->SetLineColor(kRed);
	l2->SetLineColor(kRed);
	l3->SetLineColor(kRed);
	l4->SetLineColor(kRed);

	l1->Draw("same");
	l2->Draw("same");
	l3->Draw("same");
	l4->Draw("same");

	c->SaveAs("SignalWFiducializationZ.pdf");
}
