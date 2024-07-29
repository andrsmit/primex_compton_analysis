void plot_rf_dt() {
	
	gStyle->SetOptStat(0);
	
	TString fname = "/work/halld/home/andrsmit/primex_compton_analysis/analyze_trees/phase1/analyze_data/rootFiles/Be_200nA_FIELDOFF.root";
	
	TFile *fIn = new TFile(fname.Data(), "READ");
	
	//--------------------------------------------------------------------//
	// CCAL:
	
	TH1F *h_ccal = (TH1F*)fIn->Get("ccal_rf_dt");
	h_ccal->SetLineWidth(2);
	h_ccal->SetLineColor(kBlue+2);
	h_ccal->SetFillColor(kYellow);
	h_ccal->GetXaxis()->SetTitle("t_{CCAL} - t_{RF} [ns]");
	h_ccal->GetXaxis()->SetTitleSize(0.05);
	h_ccal->GetXaxis()->SetTitleOffset(0.9);
	h_ccal->GetXaxis()->CenterTitle(true);
	h_ccal->GetYaxis()->SetTitle("counts / 0.02 ns");
	h_ccal->GetYaxis()->SetTitleSize(0.05);
	h_ccal->GetYaxis()->SetTitleOffset(0.85);
	h_ccal->GetXaxis()->SetRangeUser(-10.0,10.0);
	h_ccal->SetTitle("");
	h_ccal->GetYaxis()->SetMaxDigits(2);
	
	TCanvas *c_ccal = new TCanvas("c_ccal", "ccal", 600, 500);
	c_ccal->SetTopMargin(0.07);
	c_ccal->SetBottomMargin(0.13);
	c_ccal->SetTickx(); c_ccal->SetTicky();
	c_ccal->cd();
	h_ccal->Draw("hist");
	c_ccal->Update();
	
	TLine *lc[2];
	lc[0] = new TLine(-2.0, gPad->GetUymin(), -2.0, gPad->GetUymax());
	lc[1] = new TLine( 2.0, gPad->GetUymin(),  2.0, gPad->GetUymax());
	for(int i=0; i<2; i++) {
		lc[i]->SetLineColor(kBlack);
		lc[i]->SetLineStyle(2);
		lc[i]->SetLineWidth(2);
		lc[i]->Draw("same");
	}
	
	//--------------------------------------------------------------------//
	// FCAL:
	
	TH1F *h_fcal = (TH1F*)fIn->Get("fcal_rf_dt");
	h_fcal->SetLineWidth(2);
	h_fcal->SetLineColor(kBlue+2);
	h_fcal->SetFillColor(kYellow);
	h_fcal->GetXaxis()->SetTitle("t_{FCAL} - t_{RF} [ns]");
	h_fcal->GetXaxis()->SetTitleSize(0.05);
	h_fcal->GetXaxis()->SetTitleOffset(0.9);
	h_fcal->GetXaxis()->CenterTitle(true);
	h_fcal->GetYaxis()->SetTitle("counts / 0.02 ns");
	h_fcal->GetYaxis()->SetTitleSize(0.05);
	h_fcal->GetYaxis()->SetTitleOffset(0.85);
	h_fcal->GetXaxis()->SetRangeUser(-10.0,10.0);
	h_fcal->SetTitle("");
	h_fcal->GetYaxis()->SetMaxDigits(2);
	
	TCanvas *c_fcal = new TCanvas("c_fcal", "fcal", 600, 500);
	c_fcal->SetTopMargin(0.07);
	c_fcal->SetBottomMargin(0.13);
	c_fcal->SetTickx(); c_fcal->SetTicky();
	c_fcal->cd();
	h_fcal->Draw("hist");
	c_fcal->Update();
	
	TLine *lf[2];
	lf[0] = new TLine(-2.0, gPad->GetUymin(), -2.0, gPad->GetUymax());
	lf[1] = new TLine( 2.0, gPad->GetUymin(),  2.0, gPad->GetUymax());
	for(int i=0; i<2; i++) {
		lf[i]->SetLineColor(kBlack);
		lf[i]->SetLineStyle(2);
		lf[i]->SetLineWidth(2);
		lf[i]->Draw("same");
	}
	
	//--------------------------------------------------------------------//
	// Beam:
	
	TH1F *h_beam = (TH1F*)fIn->Get("beam_rf_dt");
	h_beam->SetLineWidth(2);
	h_beam->SetLineColor(kBlue+2);
	h_beam->SetFillColor(kYellow);
	h_beam->GetXaxis()->SetTitle("t_{#gamma} - t_{RF} [ns]");
	h_beam->GetXaxis()->SetTitleSize(0.05);
	h_beam->GetXaxis()->SetTitleOffset(0.9);
	h_beam->GetXaxis()->CenterTitle(true);
	h_beam->GetYaxis()->SetTitle("counts / 0.02 ns");
	h_beam->GetYaxis()->SetTitleSize(0.05);
	h_beam->GetYaxis()->SetTitleOffset(0.85);
	h_beam->GetXaxis()->SetRangeUser(-10.0,10.0);
	h_beam->SetTitle("");
	h_beam->GetYaxis()->SetMaxDigits(2);
	
	TCanvas *c_beam = new TCanvas("c_beam", "beam", 600, 500);
	c_beam->SetTopMargin(0.07);
	c_beam->SetBottomMargin(0.13);
	c_beam->SetTickx(); c_beam->SetTicky();
	c_beam->cd();
	h_beam->Draw("hist");
	c_beam->Update();
	
	TLine *lb[2];
	lb[0] = new TLine(-2.004, gPad->GetUymin(), -2.004, gPad->GetUymax());
	lb[1] = new TLine( 2.004, gPad->GetUymin(),  2.004, gPad->GetUymax());
	for(int i=0; i<2; i++) {
		lb[i]->SetLineColor(kBlack);
		lb[i]->SetLineStyle(2);
		lb[i]->SetLineWidth(2);
		lb[i]->Draw("same");
	}
	
	TH1F *h_beam_cut = (TH1F*)fIn->Get("beam_rf_dt_cut");
	h_beam_cut->SetLineWidth(2);
	h_beam_cut->SetLineColor(kBlue+2);
	h_beam_cut->SetFillColor(kYellow);
	h_beam_cut->GetXaxis()->SetTitle("t_{#gamma} - t_{RF} [ns]");
	h_beam_cut->GetXaxis()->SetTitleSize(0.05);
	h_beam_cut->GetXaxis()->SetTitleOffset(0.9);
	h_beam_cut->GetXaxis()->CenterTitle(true);
	h_beam_cut->GetYaxis()->SetTitle("counts / 0.02 ns");
	h_beam_cut->GetYaxis()->SetTitleSize(0.05);
	h_beam_cut->GetYaxis()->SetTitleOffset(0.85);
	h_beam_cut->GetXaxis()->SetRangeUser(-10.0,10.0);
	h_beam_cut->SetTitle("");
	h_beam_cut->GetYaxis()->SetMaxDigits(2);
	
	TCanvas *c_beam_cut = new TCanvas("c_beam_cut", "beam (cut)", 600, 500);
	c_beam_cut->SetTopMargin(0.07);
	c_beam_cut->SetBottomMargin(0.13);
	c_beam_cut->SetTickx(); c_beam_cut->SetTicky();
	c_beam_cut->cd();
	h_beam_cut->Draw("hist");
	c_beam_cut->Update();
	
	TLine *lb_cut[2];
	lb_cut[0] = new TLine(-2.004, gPad->GetUymin(), -2.004, gPad->GetUymax());
	lb_cut[1] = new TLine( 2.004, gPad->GetUymin(),  2.004, gPad->GetUymax());
	for(int i=0; i<2; i++) {
		lb_cut[i]->SetLineColor(kBlack);
		lb_cut[i]->SetLineStyle(2);
		lb_cut[i]->SetLineWidth(2);
		lb_cut[i]->Draw("same");
	}
	
	TLatex lat;
	lat.SetTextColor(kRed);
	lat.SetTextFont(52);
	lat.DrawLatexNDC(0.130, 0.875, "#scale[1.0]{Compton cuts applied}");
	
	return;
}
