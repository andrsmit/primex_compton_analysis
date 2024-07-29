void plot_elas_tail() {
	
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	
	TString fname = "/work/halld/home/andrsmit/primex_compton_analysis/analyze_trees/phase1/analyze_data/rootFiles/Be_200nA_FIELDOFF.root";
	
	TFile *fIn = new TFile(fname.Data(), "READ");
	
	//--------------------------------------------------------------------//
	// DeltaE vs DeltaK:
	
	
	TH2F *h_dEdK = (TH2F*)fIn->Get("deltaE_vs_deltaK");
	h_dEdK->GetXaxis()->SetTitle("#DeltaK [GeV]");
	h_dEdK->GetXaxis()->SetTitleSize(0.05);
	h_dEdK->GetXaxis()->SetTitleOffset(0.9);
	h_dEdK->GetXaxis()->CenterTitle(true);
	h_dEdK->GetYaxis()->SetTitle("#DeltaE [GeV]");
	h_dEdK->GetYaxis()->SetTitleSize(0.05);
	h_dEdK->GetYaxis()->SetTitleOffset(0.85);
	h_dEdK->GetYaxis()->CenterTitle(true);
	h_dEdK->SetTitle("");
	h_dEdK->SetMinimum(0.);
	h_dEdK->GetXaxis()->SetRangeUser(-8.0, 8.0);
	h_dEdK->GetYaxis()->SetRangeUser(-8.0, 8.0);
	
	for(int xbin=1; xbin<=h_dEdK->GetXaxis()->GetNbins(); xbin++) {
		for(int ybin=1; ybin<=h_dEdK->GetXaxis()->GetNbins(); ybin++) {
			double loc_deltaE = h_dEdK->GetYaxis()->GetBinCenter(ybin);
			if(fabs(loc_deltaE)<1.0) h_dEdK->SetBinContent(xbin,ybin,0.);
		}
	}
	
	/*
	TH2F *h_dEdK = (TH2F*)fIn->Get("deltaK_vs_deltaE/deltaK_vs_deltaE_8");
	h_dEdK->GetXaxis()->SetTitle("#DeltaE [GeV]");
	h_dEdK->GetXaxis()->SetTitleSize(0.05);
	h_dEdK->GetXaxis()->SetTitleOffset(0.9);
	h_dEdK->GetXaxis()->CenterTitle(true);
	h_dEdK->GetYaxis()->SetTitle("#DeltaK [GeV]");
	h_dEdK->GetYaxis()->SetTitleSize(0.05);
	h_dEdK->GetYaxis()->SetTitleOffset(0.85);
	h_dEdK->GetYaxis()->CenterTitle(true);
	h_dEdK->SetTitle("");
	h_dEdK->SetMinimum(0.);
	h_dEdK->GetXaxis()->SetRangeUser(-8.0, 8.0);
	h_dEdK->GetYaxis()->SetRangeUser(-8.0, 8.0);
	
	for(int xbin=1; xbin<=h_dEdK->GetXaxis()->GetNbins(); xbin++) {
		for(int ybin=1; ybin<=h_dEdK->GetXaxis()->GetNbins(); ybin++) {
			double loc_deltaE = h_dEdK->GetXaxis()->GetBinCenter(xbin);
			if(fabs(loc_deltaE)<1.0) h_dEdK->SetBinContent(xbin,ybin,0.);
		}
	}
	*/
	TCanvas *c_dEdK = new TCanvas("c_dEdK", "dE_vs_dK", 650, 500);
	c_dEdK->SetTopMargin(0.07);
	c_dEdK->SetBottomMargin(0.13);
	c_dEdK->SetTickx(); c_dEdK->SetTicky();
	c_dEdK->SetLogz();
	c_dEdK->cd();
	h_dEdK->Draw("COLZ");
	
	return;
}
