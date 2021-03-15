void bkgComparison(){
	TFile *_noSigFile = new TFile("noSignalBKGNoZ/noSignalBkg.root");
	TFile *_fullSigFile = new TFile("fullSignalBKGNoZ/fullSignalBkg.root");

	if (_noSigFile == nullptr || _fullSigFile == nullptr) exit(1);

	TH1D *hNoSigSpectrum = (TH1D*)_noSigFile->Get("hPileUp");
	TH1D *hFullSigSpectrum = (TH1D*)_fullSigFile->Get("hPileUp");

	auto hSegmentAlive = (TH1D*) _noSigFile->Get("hLiveSegmentFiducial");

	int aliveSegs = 0;

	for (int i = 1, j = hSegmentAlive->GetNbinsX(); i <= j; i++){
		if (hSegmentAlive->GetBinContent(i)) aliveSegs++;
	}

	
	double scaleFactor = hNoSigSpectrum->GetBinWidth(1) * aliveSegs * 24.2;

	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptTitle(kFALSE);

	auto c = new TCanvas("c", "", 950, 800);


	auto tPad = new TPad("tPad", "", 0, 0.35, 1, 1);
	tPad->SetBottomMargin(0);
	tPad->SetLogy();
	tPad->Draw();
	tPad->cd();

	auto tempNoSig = (TH1D*)hNoSigSpectrum->Clone("tempNoSig");
	auto tempFullSig = (TH1D*)hFullSigSpectrum->Clone("tempFullSig");

	tempNoSig->Scale(1/scaleFactor);
	tempFullSig->Scale(1/scaleFactor);

	tempNoSig->Scale(1.0792);
	tempFullSig->Scale(1.0792);

	double runtimeFullSignal = 255254.3598,
		   runtimeNoSignal = 244765.101;

	tempNoSig->Scale(1/runtimeNoSignal);
	tempFullSig->Scale(1/runtimeFullSignal);

	tempFullSig->SetLineColor(kRed);
	tempNoSig->SetLineColor(kBlack);

	tempFullSig->SetTitle(";;Rate [Hz/MeV/kg]");

	tempFullSig->GetXaxis()->SetRangeUser(1.5, 10);
	tempFullSig->GetYaxis()->SetRangeUser(9e-6, 3e-4);
	tempFullSig->GetYaxis()->SetTitleFont(43);
	tempFullSig->GetYaxis()->SetLabelFont(43);
	tempFullSig->GetYaxis()->CenterTitle(kTRUE);
	tempFullSig->GetYaxis()->SetTitleSize(22);
	tempFullSig->GetYaxis()->SetLabelSize(22);
	tempFullSig->GetYaxis()->SetTitleOffset(1.4);

	tempNoSig->GetXaxis()->SetRangeUser(1.5, 10);

	TLegend *l = new TLegend(0.67, 0.68, 0.88, 0.88);
	l->SetLineWidth(0);
	l->SetTextFont(43);
	l->SetTextSize(22);
	l->SetHeader("Non Fiducial segments", "C");
	l->AddEntry(tempFullSig, "00:00 - 04:00", "l");
	l->AddEntry(tempNoSig, "12:00 - 16:00", "l");

	tempFullSig->Draw("hist");
	tempNoSig->Draw("hist same");
	l->Draw();

	c->cd();

	auto bPad = new TPad("bPad", "", 0, 0, 1, 0.35);
	bPad->SetTopMargin(0);
	bPad->SetBottomMargin(0.2);
	bPad->Draw();
	bPad->cd();

	//hFullSigSpectrum->Scale(1/runtimeFullSignal);
	//hNoSigSpectrum->Scale(1/runtimeNoSignal);

	// Substract Backgrounds	
	auto hRatio = (TH1D*)hFullSigSpectrum->Clone();
	//hRatio->Add(hNoSigSpectrum, -1);
	hRatio->Divide(hNoSigSpectrum);
	hRatio->Scale(runtimeNoSignal/runtimeFullSignal);


	// Fit 1 is ax
	// Fit 2 is ax + b
	auto hRatioFit1 = (TH1D*)hRatio->Clone();
	auto hRatioFit2 = (TH1D*)hRatio->Clone();

	auto f1 = new TF1("f1", "pol0", 1.5, 10);
	auto f2 = new TF1("f2", "pol1", 1.5, 10);

	f2->SetLineColor(kBlue - 3);
	f2->SetLineWidth(2);

	// Error
	for (int i = 1, j = hNoSigSpectrum->GetNbinsX(); i <= j; i++){
		hRatioFit1->SetBinError(i, sqrt(1/hFullSigSpectrum->GetBinContent(i) + 1/hNoSigSpectrum->GetBinContent(i)) * hRatioFit1->GetBinContent(i));
		hRatioFit2->SetBinError(i, sqrt(1/hFullSigSpectrum->GetBinContent(i) + 1/hNoSigSpectrum->GetBinContent(i)) * hRatioFit2->GetBinContent(i));
	}

	//hRatioFit1->Scale(1/scaleFactor);
	//hRatioFit2->Scale(1/scaleFactor);

	hRatioFit1->Fit("f1", "R");
	hRatioFit2->Fit("f2", "R");

	printf("Fit 1:\n");
	printf("\t Chi-Square = %f\n", f1->GetChisquare());
	printf("\t Slope = %f\n\n", f1->GetParameter(0));

	printf("Fit 2:\n");
	printf("\t Chi-Square = %f\n", f2->GetChisquare());
	printf("\t Slope = %f\n\n", f2->GetParameter(0));
	printf("\t Intercept = %f\n\n", f2->GetParameter(1));

	hRatioFit2->SetTitle(";Reconstructed Pulse Energy [MeV]; Ratio");
	hRatioFit2->SetLineColor(kGray + 1);

	hRatioFit2->GetXaxis()->SetRangeUser(1.5, 10);
	hRatioFit2->GetXaxis()->CenterTitle(kTRUE);
	hRatioFit2->GetXaxis()->SetTitleFont(43);
	hRatioFit2->GetXaxis()->SetLabelFont(43);
	hRatioFit2->GetXaxis()->SetTitleSize(22);
	hRatioFit2->GetXaxis()->SetLabelSize(22);
	hRatioFit2->GetXaxis()->SetTickLength(0.06);
	hRatioFit2->GetXaxis()->SetTitleOffset(3);

	//hRatioFit2->GetYaxis()->SetRangeUser(-1.5e-5, 1.5e-5);
	hRatioFit2->GetYaxis()->CenterTitle(kTRUE);
	hRatioFit2->GetYaxis()->SetTitleFont(43);
	hRatioFit2->GetYaxis()->SetLabelFont(43);
	hRatioFit2->GetYaxis()->SetTitleSize(22);
	hRatioFit2->GetYaxis()->SetLabelSize(22);
	hRatioFit2->GetYaxis()->SetTitleOffset(1.4);
	hRatioFit2->GetYaxis()->SetNdivisions(505);

	TGaxis::SetMaxDigits(2);
	TGaxis::SetExponentOffset(0.02, -0.1, "y");

	hRatioFit2->SetMarkerStyle(8);
	hRatioFit2->SetMarkerSize(0.8);
	hRatioFit2->SetMarkerColor(kOrange + 8);

	hRatioFit2->Draw("hist P E1");
	f2->Draw("same");
}
