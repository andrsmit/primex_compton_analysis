/*
Alignment procedure to determine relative offset of the CCAL from the beam.

For each event, consider the line constructed by the FCAL shower position: (x1, y1) and 
the interaction vertex: (0, 0). In this scenario, the FCAL shower position is already
corrected from the FCAL-Compton alignment procedure and the coordinates are given in a system
where the vertex is at the origin.

Then consider the distance between this line and the point given by the CCAL shower position:
(x2+x0, y2+y0). In this scenario, x2 and y2 are the coordinates reconstructed by the 
DCCALShower_factory (remove any alignment correction applied in the reconstruction), and x0, y0 
is the relative offset between the CCAL and the Beam that needs to be applied to the shower
position. We minimize the sum of squares of this distance for every event (d/dx0(D^2)=0 &
d/dy0(D^2)=0) to get a system of linear equations to solve for x0, y0.

Then the CCAL position is given by: x_ccal = x0 + x_beam & y_ccal = y0 + y_beam.
*/


Double_t crys_ball_fit(Double_t *x, Double_t *par);

void alignCCAL() {
	
	char pathName[256];
	sprintf(pathName, "/work/halld/home/andrsmit/ccal_calib/phaseIII/alignment/ccalTrees");
	
	char infile[256];
	sprintf(infile, "%s/110600.root", pathName);
	
	//---------------------------------------------//
	
	const double m_e = 0.510998928e-3;
	
	const double m_fcalX_ccdb =    0.529;
	const double m_fcalY_ccdb =   -0.002;
	const double m_fcalZ      =  624.906;
	
	const double m_beamX_ccdb =    0.;
	const double m_beamY_ccdb =    0.;
	const double m_beamZ      =   64.935;
	
	const double m_ccalX_ccdb =   -0.0225;
	const double m_ccalY_ccdb =    0.0073;
	const double m_ccalZ      = 1274.274;
	
	// results of FCAL-Compton alignment:
	
	double m_fcalX_aligned = 0.408;
	double m_fcalY_aligned = 0.027;
	double m_beamX_aligned = 0.151;
	double m_beamY_aligned = 0.012;
	
	//---------------------------------------------//
	
	double min_deltaE_cut = -0.2;
	double max_deltaE_cut =  0.4;
	double min_deltaK_cut = -1.0;
	double max_deltaK_cut =  0.6;
	
	double MIN_ENERGY_CUT =  3.0;
	double MAX_ENERGY_CUT = 12.0;
	
	double fcal_block_size = 4.0157;
	double ccal_block_size = 2.09;
	
	//---------------------------------------------//
	
	if(gSystem->AccessPathName(infile)) {
		cout << "Specified input file does not exist." << endl;
		return;
	}
	
	//---------  Histograms  ---------//
	
	TH2F *hPhi_old = new TH2F("phi_old", "Before Correction", 500, 50., 130., 1000, -90., 90.);
	hPhi_old->GetXaxis()->SetTitle("|#phi_{1}-#phi_{2}| / 2 [deg.]");
	hPhi_old->GetYaxis()->SetTitle("(#phi_{1}+#phi_{2}) / 2 [deg.]");
	
	TH2F *hPhi_new = new TH2F("phi_new", "After Correction",  500, 50., 130., 1000, -90., 90.);
	hPhi_new->GetXaxis()->SetTitle( "|#phi_{1}-#phi_{2}| / 2 [deg.]");
	hPhi_new->GetYaxis()->SetTitle( "(#phi_{1}+#phi_{2}) / 2 [deg.]");
	
	TH1F *h_deltaE     = new TH1F("deltaE",     "E_{1} + E_{2} - E_{#gamma}; [GeV]", 
		2000, -2.0, 2.0);
	TH1F *h_deltaK_old = new TH1F("deltaK_old", "E_{Comp} - E_{#gamma}; [GeV]", 
		2000, -6.0, 6.0);
	TH1F *h_deltaK_new = new TH1F("deltaK_new", "E_{Comp} - E_{#gamma}; [GeV]", 
		2000, -6.0, 6.0);
	
	TLatex lat;
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetPalette(1);
	
	//----------  Begin Analysis  ----------//
	
	TChain *t = new TChain("ccal_compton");
	
	double eb;
	double e1, x1, y1, z1;
	double e2, x2, y2, z2;
	int tof_match;
	
	t->SetBranchAddress("eb", &eb);
	t->SetBranchAddress("e1", &e1);
	t->SetBranchAddress("x1", &x1);
	t->SetBranchAddress("y1", &y1);
	t->SetBranchAddress("z1", &z1);
	t->SetBranchAddress("e2", &e2);
	t->SetBranchAddress("x2", &x2);
	t->SetBranchAddress("y2", &y2);
	t->SetBranchAddress("z2", &z2);
	t->SetBranchAddress("tof_match", &tof_match);
	
	t->Add(infile);
	
	int n_entries = (int)t->GetEntries();
	cout << "Total number of entries: " << n_entries << endl;
	
	double A = 0.;
	double B = 0.;
	double C = 0.;
	double F = 0.;
	double G = 0.;
	
	int counter = 0;
	vector<double> countVec, x0Vec, y0Vec;
	countVec.clear();
	x0Vec.clear();
	y0Vec.clear();
	
	double x0 = 0.; // x0 = x_ccal - x_beam
	double y0 = 0.; // y0 = y_ccal - y_beam
	
	int num_events_used = 0;
	
	for(int ievt = 1; ievt < n_entries; ievt++) {
		
		t->GetEvent(ievt);
		if(ievt%1000000 == 0) cout << "  Processing event " << ievt << endl; 
		
		if(!(MIN_ENERGY_CUT < eb && eb < MAX_ENERGY_CUT)) continue;
		if(e1 < 1.0 || e2 < 3.0) continue;
		
		/*
		The position stored in the tree comes directly from the DF(C)CALShower_factory.
		We want to apply a fiducial cut to remove the inner layer of the calorimeter. 
		In order to do so, we have to project both shower positions to the face of the 
		detectors.
		*/
		
		//----------------------------------------------//
		// Correct FCAL shower coordinates from the result of its alignment procedure:
		
		double x1_cor = x1 - m_fcalX_ccdb + m_fcalX_aligned - m_beamX_aligned;
		double y1_cor = y1 - m_fcalY_ccdb + m_fcalY_aligned - m_beamY_aligned;
		double z1_cor = z1 - m_beamZ;
		double x2_cor = x2 - m_beamX_aligned;
		double y2_cor = y2 - m_beamY_aligned;
		double z2_cor = z2 - m_beamZ;
		
		//----------------------------------------------//
		// Apply fiducial cuts:
		
		double fcal_fid_cut = 3.5*fcal_block_size;
		double ccal_fid_cut = 2.0*ccal_block_size;
		
		double x1_fcal_face = m_beamX_aligned + (x1_cor * (m_fcalZ-m_beamZ)/z1_cor);
		double y1_fcal_face = m_beamY_aligned + (y1_cor * (m_fcalZ-m_beamZ)/z1_cor);
		x1_fcal_face       -= m_fcalX_aligned;
		y1_fcal_face       -= m_fcalY_aligned;
		
		// project CCAL shower coordinates to the face of calorimeter (for this step, assume 
		// that the alignment data is correct):
		
		double x2_face = m_beamX_aligned + (x2_cor * (m_ccalZ-m_beamZ)/z2_cor);
		double y2_face = m_beamY_aligned + (y2_cor * (m_ccalZ-m_beamZ)/z2_cor);
		x2_face       -= m_ccalX_ccdb;
		y2_face       -= m_ccalY_ccdb;
		
		if(((-1.*fcal_fid_cut)<x1_fcal_face&&x1_fcal_face<fcal_fid_cut) && 
			((-1.*fcal_fid_cut)<y1_fcal_face&&y1_fcal_face<fcal_fid_cut)) continue;
		
		if(((-1.*ccal_fid_cut)<x2_face&&x2_face<ccal_fid_cut) && 
			((-1.*ccal_fid_cut)<y2_face&&y2_face<ccal_fid_cut)) continue;
		
		//----------------------------------------------//
		// Apply elasticity cut:
		
		double deltaE = (e1+e2) - (eb+m_e);
		h_deltaE->Fill(deltaE);
		
		if(deltaE<min_deltaE_cut || deltaE>max_deltaE_cut) continue;
		
		//----------------------------------------------//
		// Now project FCAL shower coordinates to the z-position of the CCAL shower:
		
		double loc_x1 = x1_cor * (z2 - m_beamZ)/z1_cor;
		double loc_y1 = y1_cor * (z2 - m_beamZ)/z1_cor;
		double loc_z1 = z2 - m_beamZ;
		
		// remove alignment correction applied to x2,y2 during reconstruction:
		
		double loc_x2 = x2 - m_ccalX_ccdb;
		double loc_y2 = y2 - m_ccalY_ccdb;
		double loc_z2 = z2 - m_beamZ;
		
		/*
		Now 'loc_x1' and 'loc_y1' give the position of the FCAL shower at the z-position of 
		the CCAL shower in a coordinate system where the vertex is at the origin.
		'loc_x2' and 'loc_y2' give the position of the CCAL shower relative to the center of the CCAL.
		*/
		
		double theta_e   = atan2(sqrt(loc_x2*loc_x2 + loc_y2*loc_y2), loc_z1);
		double theta_gam = atan2(sqrt(loc_x1*loc_x1 + loc_y1*loc_y1), loc_z2);
		
		double Ecompton  = (m_e * sin(theta_gam) 
			/ (2. * pow(sin(theta_gam/2.),2.0) * tan(theta_e))) - m_e;
		double deltaK    = Ecompton - eb;
		
		if(deltaK<min_deltaK_cut || deltaK>max_deltaK_cut) continue;
		
		num_events_used++;
		
		//----------------------------------------------//
		
		double r1sq = loc_x1*loc_x1 + loc_y1*loc_y1;
		
		A += (loc_y1*loc_y1) / r1sq;
		B += -1. * (loc_x1*loc_y1) / r1sq;
		C += loc_y1 * (loc_x1*loc_y2 - loc_x2*loc_y1) / r1sq;
		F += (loc_x1*loc_x1) / r1sq;
		G += loc_x1 * (loc_x2*loc_y1 - loc_x1*loc_y2) / r1sq;
		
		double loc_x0 = (C*F - B*G) / (A*F - B*B);
		double loc_y0 = (B*C - A*G) / (B*B - A*F);
		
		if((counter+1)%100 == 0 ) {
			countVec.push_back((double)counter);
			x0Vec.push_back(loc_x0);
			y0Vec.push_back(loc_y0);
		}
		counter++;
	}
	
	double m_ccal_beam_dX = x0Vec[(int)x0Vec.size()-1];
	double m_ccal_beam_dY = y0Vec[(int)y0Vec.size()-1];
	
	for(int ievt = 1; ievt < n_entries; ievt++) {
		
		t->GetEvent(ievt);
		if(ievt%1000000 == 0) cout << "  Processing event " << ievt << endl; 
		
		if(!(MIN_ENERGY_CUT < eb && eb < MAX_ENERGY_CUT)) continue;
		if(e1 < 1.0 || e2 < 3.0) continue;
		
		//----------------------------------------------//
		// Correct FCAL shower coordinates from the result of its alignment procedure:
		
		double x1_cor = x1 - m_fcalX_ccdb + m_fcalX_aligned - m_beamX_aligned;
		double y1_cor = y1 - m_fcalY_ccdb + m_fcalY_aligned - m_beamY_aligned;
		double z1_cor = z1 - m_beamZ;
		double x2_cor = x2 - m_beamX_aligned;
		double y2_cor = y2 - m_beamY_aligned;
		double z2_cor = z2 - m_beamZ;
		
		//----------------------------------------------//
		// Apply fiducial cuts:
		
		double fcal_fid_cut = 3.5*fcal_block_size;
		double ccal_fid_cut = 2.0*ccal_block_size;
		
		double x1_fcal_face = m_beamX_aligned + (x1_cor * (m_fcalZ-m_beamZ)/z1_cor);
		double y1_fcal_face = m_beamY_aligned + (y1_cor * (m_fcalZ-m_beamZ)/z1_cor);
		x1_fcal_face       -= m_fcalX_aligned;
		y1_fcal_face       -= m_fcalY_aligned;
		
		// project CCAL shower coordinates to the face of calorimeter (for this step, assume 
		// that the alignment data is correct):
		
		double x2_face = m_beamX_aligned + (x2_cor * (m_ccalZ-m_beamZ)/z2_cor);
		double y2_face = m_beamY_aligned + (y2_cor * (m_ccalZ-m_beamZ)/z2_cor);
		x2_face       -= m_ccalX_ccdb;
		y2_face       -= m_ccalY_ccdb;
		
		if(((-1.*fcal_fid_cut)<x1_fcal_face&&x1_fcal_face<fcal_fid_cut) && 
			((-1.*fcal_fid_cut)<y1_fcal_face&&y1_fcal_face<fcal_fid_cut)) continue;
		
		if(((-1.*ccal_fid_cut)<x2_face&&x2_face<ccal_fid_cut) && 
			((-1.*ccal_fid_cut)<y2_face&&y2_face<ccal_fid_cut)) continue;
		
		//----------------------------------------------//
		// Apply elasticity cut:
		
		double deltaE = (e1+e2) - (eb+m_e);
		h_deltaE->Fill(deltaE);
		
		if(deltaE<min_deltaE_cut || deltaE>max_deltaE_cut) continue;
		
		//----------------------------------------------//
		// Now project FCAL shower coordinates to the z-position of the CCAL shower:
		
		double loc_x1 = x1_cor * (z2 - m_beamZ)/z1_cor;
		double loc_y1 = y1_cor * (z2 - m_beamZ)/z1_cor;
		double loc_z1 = z2 - m_beamZ;
		
		// remove alignment correction applied to x2,y2 during reconstruction:
		
		double loc_x2 = x2 - m_ccalX_ccdb;
		double loc_y2 = y2 - m_ccalY_ccdb;
		double loc_z2 = z2 - m_beamZ;
		
		double loc_x2_cor = loc_x2 + m_ccal_beam_dX;
		double loc_y2_cor = loc_y2 + m_ccal_beam_dY;
		
		double phi1     = atan2(loc_y1,     loc_x1)     * (180./TMath::Pi());
		double phi2_old = atan2(loc_y2,     loc_x2)     * (180./TMath::Pi());
		double phi2_new = atan2(loc_y2_cor, loc_x2_cor) * (180./TMath::Pi());
		
		double sumPhi_old = 0.5 *     (phi2_old + phi1);
		double delPhi_old = 0.5 * fabs(phi2_old - phi1);
		double sumPhi_new = 0.5 *     (phi2_new + phi1);
		double delPhi_new = 0.5 * fabs(phi2_new - phi1);
		
		hPhi_old->Fill(delPhi_old, sumPhi_old);
		hPhi_new->Fill(delPhi_new, sumPhi_new);
		
		double theta1     = atan2(sqrt(pow(loc_x1,2.0)     + pow(loc_y1,2.0)),     loc_z1);
		double theta2     = atan2(sqrt(pow(loc_x2,2.0)     + pow(loc_y2,2.0)),     loc_z2);
		double theta2_new = atan2(sqrt(pow(loc_x2_cor,2.0) + pow(loc_y2_cor,2.0)), loc_z2);
		
		double ecomp1      =  1. / ((1./eb) + (1./m_e)*(1.-cos(theta1)));
		double ecomp2      =  1. / ((1./eb) + (1./m_e)*(1.-cos(theta2)));
		double ecomp2_new  =  1. / ((1./eb) + (1./m_e)*(1.-cos(theta2_new)));
		
		h_deltaK_old->Fill(ecomp1 + ecomp2     - eb);
		h_deltaK_new->Fill(ecomp1 + ecomp2_new - eb);
	}
	
	int nevts = (int)countVec.size();
	
	double *counts = new double[nevts];
	double *x0s    = new double[nevts];
	double *y0s    = new double[nevts];
	
	for(int i=0; i<nevts; i++) {
		counts[i] = countVec[i];
		x0s[i]    =  x0Vec[i];
		y0s[i]    =  y0Vec[i];
	}
	
	TGraph *gX = new TGraph(nevts, counts, x0s);
	gX->SetMarkerColor(kBlack);
	gX->SetMarkerStyle(8);
	gX->SetTitle("CompCal Center X Position");
	gX->GetXaxis()->SetTitle("Number Events");
	gX->GetYaxis()->SetTitle("x_{0} [cm]");
	
	TCanvas *cX = new TCanvas("cX", "cX", 1200, 600);
	cX->SetTickx(); cX->SetTicky(); cX->SetGrid();
	gX->Draw("AP");
	cX->Update();
	TLine *lX = new TLine(gPad->GetUxmin(), -0.0225, gPad->GetUxmax(), -0.0225);
	lX->SetLineStyle(2);
	lX->SetLineWidth(2);
	lX->SetLineColor(kRed);
	lX->Draw("same");
	
	TGraph *gY = new TGraph(nevts, counts, y0s);
	gY->SetMarkerColor(kBlack);
	gY->SetMarkerStyle(8);
	gY->SetTitle("CompCal Center Y Position");
	gY->GetXaxis()->SetTitle("Number Events");
	gY->GetYaxis()->SetTitle("y_{0} [cm]");
	
	TCanvas *cY = new TCanvas("cY", "cY", 1200, 600);
	cY->SetTickx(); cY->SetTicky(); cY->SetGrid();
	gY->Draw("AP");
	cY->Update();
	TLine *lY = new TLine(gPad->GetUxmin(), 0.0073, gPad->GetUxmax(), 0.0073);
	lY->SetLineStyle(2);
	lY->SetLineWidth(2);
	lY->SetLineColor(kRed);
	lY->Draw("same");
	
	//----------   DeltaE   ----------//
	
	h_deltaE->SetTitle("#DeltaE");
	h_deltaE->GetXaxis()->SetTitle("E_{1} + E_{2} - E_{#gamma} [GeV]");
	h_deltaE->GetXaxis()->SetTitleOffset(1.2);
	//h_deltaE->GetXaxis()->SetRangeUser(-1.0, 1.0);
	h_deltaE->SetLineColor(kBlack);
	
	TF1 *f_deltaE = new TF1("f_deltaE", crys_ball_fit, -0.6, 0.4, 9);
	
	f_deltaE->SetParName(0, "N");
	f_deltaE->SetParName(1, "#mu");
	f_deltaE->SetParName(2, "#sigma");
	f_deltaE->SetParName(3, "#alpha");
	f_deltaE->SetParName(4, "n");
	
	f_deltaE->SetParameters(h_deltaE->GetMaximum(), 
		h_deltaE->GetBinCenter(h_deltaE->GetMaximumBin()), 0.15, 1.5, 1.5);
	f_deltaE->FixParameter(5, 0.);
	f_deltaE->FixParameter(6, 0.);
	f_deltaE->FixParameter(7, 0.);
	f_deltaE->FixParameter(8, 0.);
	
	h_deltaE->Fit("f_deltaE", "R0QL");
	f_deltaE->SetLineColor(kRed);
	
	TCanvas *ce = new TCanvas("ce", "ce", 600, 600);
	ce->SetTickx(); ce->SetTicky();
	h_deltaE->Draw();
	f_deltaE->Draw("same");
	
	ce->Update();
	TLine *lE1 = new TLine(min_deltaE_cut, gPad->GetUymin(), 
		min_deltaE_cut, gPad->GetUymax());
	lE1->SetLineColor(kRed); lE1->SetLineStyle(2); lE1->Draw("same");
	TLine *lE2 = new TLine(max_deltaE_cut, gPad->GetUymin(), 
		max_deltaE_cut, gPad->GetUymax());
	lE2->SetLineColor(kRed); lE2->SetLineStyle(2); lE2->Draw("same");
	
	//----------   DeltaK   ----------//
	
	h_deltaK_old->SetTitle("#DeltaK");
	h_deltaK_old->GetXaxis()->SetTitle("E_{Compton} - E_{#gamma} [GeV]");
	h_deltaK_old->GetXaxis()->SetTitleOffset(1.2);
	//h_deltaK_old->GetXaxis()->SetRangeUser(-1.5, 1.5);
	h_deltaK_old->SetLineColor( kBlack );
	
	TF1 *f_deltaK_old = new TF1("f_deltaK_old", "gaus", -0.3, 0.3);
	f_deltaK_old->SetParName(0, "N");
	f_deltaK_old->SetParName(1, "#mu");
	f_deltaK_old->SetParName(2, "#sigma");
	f_deltaK_old->SetParameters(h_deltaK_old->GetMaximum(), 
		h_deltaK_old->GetBinCenter(h_deltaK_old->GetMaximumBin()), 0.3);
	h_deltaK_old->Fit("f_deltaK_old", "R0QL");
	f_deltaK_old->SetLineColor(kRed);
	
	h_deltaK_new->SetTitle("#DeltaK");
	h_deltaK_new->GetXaxis()->SetTitle("E_{Compton} - E_{#gamma} [GeV]");
	h_deltaK_new->GetXaxis()->SetTitleOffset(1.2);
	//h_deltaK_new->GetXaxis()->SetRangeUser(-1.5, 1.5);
	h_deltaK_new->SetLineColor( kBlack );
	
	TF1 *f_deltaK_new = new TF1("f_deltaK_new", "gaus", -0.3, 0.3);
	f_deltaK_new->SetParName(0, "N");
	f_deltaK_new->SetParName(1, "#mu");
	f_deltaK_new->SetParName(2, "#sigma");
	f_deltaK_new->SetParameters(h_deltaK_new->GetMaximum(), 
		h_deltaK_new->GetBinCenter(h_deltaK_new->GetMaximumBin()), 0.3);
	h_deltaK_new->Fit("f_deltaK_new", "R0QL");
	f_deltaK_new->SetLineColor(kRed);
	
	TCanvas *ck = new TCanvas("ck", "ck", 1200, 600);
	ck->Divide(2,1);
	
	TPad *ck1 = (TPad*)ck->cd(1);
	ck1->SetTickx(); ck1->SetTicky();
	h_deltaK_old->Draw();
	f_deltaK_old->Draw("same");
	
	TPad *ck2 = (TPad*)ck->cd(2);
	ck2->SetTickx(); ck2->SetTicky();
	h_deltaK_new->Draw();
	f_deltaK_new->Draw("same");
	/*
	ck1->Update();
	TLine *lK1 = new TLine(min_deltaK_cut, gPad->GetUymin(), 
		min_deltaK_cut, gPad->GetUymax());
	lK1->SetLineColor(kRed); lK1->SetLineStyle(2); lK1->Draw("same");
	TLine *lK2 = new TLine(max_deltaK_cut, gPad->GetUymin(), 
		max_deltaK_cut, gPad->GetUymax());
	lK2->SetLineColor(kRed); lK2->SetLineStyle(2); lK2->Draw("same");
	*/
	//--------------------------------//
	
	TCanvas *cPhi = new TCanvas("cPhi", "cPhi", 1200, 600);
	cPhi->Divide(2,1);
	
	TPad *cPhi1 = (TPad*)cPhi->cd(1);
	cPhi1->SetTickx(); cPhi1->SetTicky(); cPhi1->SetGrid();
	hPhi_old->GetXaxis()->SetRangeUser(70., 110.);
	hPhi_old->Draw("COL");
	TLine *lv1 = new TLine(90., -90., 90., 90.);
	lv1->SetLineColor(kBlack);
	lv1->SetLineWidth(2);
	lv1->Draw("same");
	
	TPad *cPhi2 = (TPad*)cPhi->cd(2);
	cPhi2->SetTickx(); cPhi2->SetTicky(); cPhi2->SetGrid();
	hPhi_new->GetXaxis()->SetRangeUser(70., 110.);
	hPhi_new->Draw("COL");
	lv1->Draw("same");
	
	TH1I *h1_new = (TH1I*)hPhi_new->ProjectionX("h1_new");
	h1_new->SetLineColor(kBlack);
	h1_new->SetLineWidth(2);
	h1_new->SetTitle("#Delta#phi (after correction)");
	h1_new->GetXaxis()->SetRangeUser(70., 110.);
	h1_new->GetXaxis()->SetTitle("0.5 * | #phi_{1} - #phi_{2} |  [deg.]");
	h1_new->GetXaxis()->SetTitleOffset(1.3);
	
	TH1I *h1_old = (TH1I*)hPhi_old->ProjectionX("h1_old");
	h1_old->SetLineColor(kBlack);
	h1_old->SetLineWidth(2);
	h1_old->SetTitle("#Delta#phi");
	h1_old->GetXaxis()->SetRangeUser(70., 110.);
	h1_old->GetXaxis()->SetTitle("0.5 * | #phi_{1} - #phi_{2} |  [deg.]");
	h1_old->GetXaxis()->SetTitleOffset(1.3);
	
	TCanvas *cdPhi = new TCanvas("cdPhi", "cdPhi", 1200, 600);
	cdPhi->Divide(2,1);
	
	TPad *pdPhi1 = (TPad*)cdPhi->cd(1);
	pdPhi1->SetTickx(); pdPhi1->SetTicky();
	h1_old->Draw();
	
	TPad *pdPhi2 = (TPad*)cdPhi->cd(2);
	pdPhi2->SetTickx(); pdPhi2->SetTicky();
	h1_new->Draw();
	
	//----------     Fit Distributions    ----------//
	
	pdPhi1->cd();
	
	TF1 *f_gaus_old = new TF1("f_gaus_old", "gaus", 85., 95.);
	f_gaus_old->SetParameters((double)h1_old->GetMaximum(), 
		h1_old->GetBinCenter(h1_old->GetMaximumBin()), 3.0);
	h1_old->Fit("f_gaus_old", "R0QL");
	
	f_gaus_old->Draw("same");
	
	pdPhi2->cd();
	
	TF1 *f_gaus_new = new TF1("f_gaus_new", "gaus", 85., 95.);
	f_gaus_new->SetParameters((double)h1_new->GetMaximum(), 
		h1_new->GetBinCenter(h1_new->GetMaximumBin()), 3.0);
	h1_new->Fit("f_gaus_new", "R0QL");
	
	f_gaus_new->Draw("same");
	
	//--------------------------------//
	
	cout << "\n\n\n\n";
	cout << "x_CCAL - x_BEAM = " << m_ccal_beam_dX << endl;
	cout << "y_CCAL - y_BEAM = " << m_ccal_beam_dY << endl;
	cout << "   x_CCAL = " << (m_ccal_beam_dX + m_beamX_aligned) << endl;
	cout << "   y_CCAL = " << (m_ccal_beam_dY + m_beamY_aligned) << endl;
	cout << "\n\n\n";
	
	return;
}

Double_t crys_ball_fit(Double_t *x, Double_t *par) 
{
	Double_t xx  =   x[0];
	
	Double_t N   = par[0];
	Double_t mu  = par[1];
	Double_t sig = par[2];
	Double_t a   = par[3];
	Double_t n   = par[4];
	
	Double_t p0  = par[5];
	Double_t p1  = par[6];
	Double_t p2  = par[7];
	Double_t p3  = par[8];
	
	Double_t A   = 	pow(n/fabs(a), n) * exp(-0.5*pow(fabs(a),2.0));
	Double_t B   =  (n/fabs(a)) - fabs(a);
	
	Double_t loc_x = (xx - mu) / sig;
	Double_t f;
	
	if(loc_x > -a) {
		f = N * exp(-0.5*pow(loc_x,2.0));
	} else {
		f = N * A * pow(B - loc_x, -n);
	}
	
	f += p0 + p1*xx + p2*xx*xx + p3*xx*xx*xx;
	
	return f;
}
