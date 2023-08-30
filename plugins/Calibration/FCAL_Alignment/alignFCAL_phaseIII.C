
Double_t bkgd_fit(Double_t *x, Double_t *par);
Double_t crys_ball_fit(Double_t *x, Double_t *par);
Double_t deltaK_fit(Double_t *x, Double_t *par);
Double_t deltaK_fit_exp(Double_t *x, Double_t *par);

void alignFCAL() {
	
	gStyle->SetPalette(1);
	
	char pathName[256];
	sprintf(pathName, "/work/halld/home/andrsmit/ccal_calib/phaseIII/alignment/fcalTrees");
	
	char infile[256];
	sprintf(infile, "%s/110609.root", pathName);
	
	//----------   Initialize Constants and Cuts   ----------//
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	const int rebins = 1;
	
	double m_fcalZ      = 624.906;
	double m_fcalX_ccdb =   0.529;
	double m_fcalY_ccdb =  -0.002;
	
	double m_beamZ      =  64.939;
	double m_beamX_ccdb =   0.;
	double m_beamY_ccdb =   0.;
	
	double m_fcalX, m_fcalY;
	double m_beamX, m_beamY;
	
	// from phase-II alignment result:
	
	m_fcalX = 0.408;
	m_fcalY = 0.027;
	
	//-------------------------------------------------------//
	
	bool USE_FACE_POSITIONS = true;
	
	double MIN_ENERGY_CUT = 2.0;
	double MAX_ENERGY_CUT = 4.0;
	
	double deltaE_min_cut = -0.3;
	double deltaE_max_cut =  1.0;
	double deltaK_min_cut = -1.0;
	double deltaK_max_cut =  0.5;
	
	double fcal_block_size = 4.0157;
	double m_e = 0.510998928e-3;
	
	double fid_cut = 2.0*fcal_block_size;
	
	//-------------------------------------------------------//
	
	if(gSystem->AccessPathName(infile)) {
		cout << "Input file specified does not exist. Check pathname." << endl;
		exit(-1);
	}
	
	//---------  Histograms  ---------//
	
	TH1F *h_deltaPhi_old = new TH1F("deltaPhi_old", "#Delta#phi_{old}; [deg.]", 
		1000, 100., 260.);
	TH1F *h_deltaPhi_new = new TH1F("deltaPhi_new", "#Delta#phi_{new}; [deg.]", 
		1000, 100., 260.);
	
	TH1F *h_deltaK_old = new TH1F("deltaK_old", "#DeltaK_{old}; [GeV]", 1000, -4., 4.);
	TH1F *h_deltaK_new = new TH1F("deltaK_new", "#DeltaK_{new}; [GeV]", 1000, -4., 4.);
	
	TH1F *h_deltaE = new TH1F("deltaE", "#DeltaE; [GeV]", 1000, -4., 4.);
	
	TH2F *h_xy1  =  new TH2F( "xy1", "Shower 1 Position; x [cm]; y [cm]", 
		500, -30., 30., 500, -30., 30. );
	TH2F *h_xy2  =  new TH2F( "xy2", "Shower 2 Position; x [cm]; y [cm]", 
		500, -30., 30., 500, -30., 30. );
	
	TH1F *h_e1 = new TH1F("e1", "Energy of Shower 1", 1000, 0., 5.);
	TH1F *h_e2 = new TH1F("e2", "Energy of Shower 2", 1000, 0., 5.);
	
	TH2F *h_sumPhi_vs_delPhi_old = new TH2F("sumPhi_vs_delPhi_old", "Before Alignment", 
		500, 50., 130., 1000, -90., 90.);
	h_sumPhi_vs_delPhi_old->GetXaxis()->SetTitle( "|#phi_{1}-#phi_{2}| / 2 [deg.]" );
	h_sumPhi_vs_delPhi_old->GetYaxis()->SetTitle( "(#phi_{1}+#phi_{2}) / 2 [deg.]" );
	
	TH2F *h_sumPhi_vs_delPhi_new = new TH2F("sumPhi_vs_delPhi_new", "After Alignment", 
		500, 50., 130., 1000, -90., 90.);
	h_sumPhi_vs_delPhi_new->GetXaxis()->SetTitle( "|#phi_{1}-#phi_{2}| / 2 [deg.]" );
	h_sumPhi_vs_delPhi_new->GetYaxis()->SetTitle( "(#phi_{1}+#phi_{2}) / 2 [deg.]" );
	
	//----------  Begin Analysis  ----------//
	
	TChain *t = new TChain("fcal_compton");
	
	double eb, tb, rfTime;
	double e1, x1, y1, z1, t1;
	double e2, x2, y2, z2, t2;
	int tof_match1, tof_match2;
	
	t->SetBranchAddress("eb",     &eb);
	t->SetBranchAddress("tb",     &tb);
	t->SetBranchAddress("rfTime", &rfTime);
	t->SetBranchAddress("e1",         &e1);
	t->SetBranchAddress("x1",         &x1);
	t->SetBranchAddress("y1",         &y1);
	t->SetBranchAddress("z1",         &z1);
	t->SetBranchAddress("t1",         &t1);
	t->SetBranchAddress("tof_match1", &tof_match1);
	t->SetBranchAddress("e2",         &e2);
	t->SetBranchAddress("x2",         &x2);
	t->SetBranchAddress("y2",         &y2);
	t->SetBranchAddress("z2",         &z2);
	t->SetBranchAddress("t2",         &t2);
	t->SetBranchAddress("tof_match2", &tof_match2);
	
	t->Add(infile);
	
	int n_entries = t->GetEntries();
	cout << "Total number of entries: " << n_entries << endl;
	
	double A = 0., B = 0., C = 0., F = 0., G = 0.;
	
	int n_count = 0;
	
	vector<double> countVec, x0Vec, y0Vec;
	
	double x0 = 0.; // x0 = x_beam - x_fcal
	double y0 = 0.; // y0 = y_beam - y_fcal
	
	int num_events_used = 0;
	
	for(int ievt = 0; ievt < n_entries; ievt++) {
		
		t->GetEvent(ievt);
		if(ievt%100000 == 0) cout << "Processing event " << ievt << endl;
		
		if(!(MIN_ENERGY_CUT < eb && eb < MAX_ENERGY_CUT)) continue;
		if(e1 < 1.0 || e2 < 1.0) continue;
		
		/*
		First, we want to apply a fiducial cut to remove the inner part of the calorimeter. 
		In order to do so, we have to project both shower positions to the face of the 
		detectors. This requires knowledge of the beam position, and FCAL position, but that 
		obviously isn't yet known at this point. 
		It should be fine at this stage to just assume the beam is (0,0), and the FCAL position 
		is what's stored in the ccdb. If we're wrong about this, it will have a minimal effect
		on the final result, as it is only a cut to remove showers on the inner part, and these
		positions are not what's used to calculate the relative offset.
		*/
		
		double z_face = m_fcalZ - m_beamZ;
		
		double x1_cor = x1 - m_beamX_ccdb;
		double y1_cor = y1 - m_beamY_ccdb;
		double z1_cor = z1 - m_beamZ;
		double x2_cor = x2 - m_beamX_ccdb;
		double y2_cor = y2 - m_beamY_ccdb;
		double z2_cor = z2 - m_beamZ;
		
		double x1_face = m_beamX_ccdb + (x1_cor*z_face/z1_cor);
		double y1_face = m_beamY_ccdb + (y1_cor*z_face/z1_cor);
		double x2_face = m_beamX_ccdb + (x2_cor*z_face/z2_cor);
		double y2_face = m_beamY_ccdb + (y2_cor*z_face/z2_cor);
		
		x1_face -= m_fcalX_ccdb;
		y1_face -= m_fcalY_ccdb;
		x2_face -= m_fcalX_ccdb;
		y2_face -= m_fcalY_ccdb;
		
		if(((-1.*fid_cut)<x1_face&&x1_face<fid_cut) && 
			((-1.*fid_cut)<y1_face&&y1_face<fid_cut)) continue;
		if(((-1.*fid_cut)<x2_face&&x2_face<fid_cut) && 
			((-1.*fid_cut)<y2_face&&y2_face<fid_cut)) continue;
		
		/*
		The position stored in the tree comes directly from the DFCALShower_factory.
		This position already has the offset of the FCAL which is stored in the CCDB applied to it.
		Our goal in this analysis is to determine this offset (relative to the beam), so to start, 
		let's remove this initial offset correction:
		*/
		
		double loc_x1, loc_y1, loc_z1;
		double loc_x2, loc_y2, loc_z2;
		
		if(USE_FACE_POSITIONS) {
			
			// Option 1: Use shower positions projected to the surface of the FCAL:
			
			loc_x1 = x1_face;
			loc_y1 = y1_face;
			loc_z1 = m_fcalZ - m_beamZ;
			loc_x2 = x2_face;
			loc_y2 = y2_face;
			loc_z2 = m_fcalZ - m_beamZ;
			
		} else {
			
			// Option 2: Use shower positions as reconstructed inside the calorimeter:
			
			loc_x1 = x1 - m_fcalX_ccdb;
			loc_y1 = y1 - m_fcalY_ccdb;
			loc_z1 = z1 - m_beamZ;
			loc_x2 = x2 - m_fcalX_ccdb;
			loc_y2 = y2 - m_fcalY_ccdb;
			loc_z2 = z2 - m_beamZ;
		}
		
		//---------------------------------------------------//
		// Apply TOF Veto:
		
		int tof_veto = 0;
		if(tof_match1 && !tof_match2) tof_veto = 1;
		else if(tof_match2 && !tof_match1) tof_veto = 1;
		
		if(!tof_veto) continue;
		
		//---------------------------------------------------//
		// DeltaE and DeltaK Cuts:
		
		double deltaE = (e1+e2) - (eb+m_e);
		h_deltaE->Fill(deltaE);
		
		if(deltaE<deltaE_min_cut || deltaE>deltaE_max_cut) continue;
		
		double theta_g, theta_e;
		
		if(tof_match1) {
			theta_g = atan2(sqrt(loc_x2*loc_x2 + loc_y2*loc_y2), loc_z2);
			theta_e = atan2(sqrt(loc_x1*loc_x1 + loc_y1*loc_y1), loc_z1);
		} else {
			theta_g = atan2(sqrt(loc_x1*loc_x1 + loc_y1*loc_y1), loc_z1);
			theta_e = atan2(sqrt(loc_x2*loc_x2 + loc_y2*loc_y2), loc_z2);
		}
		
		double Ecompton = (m_e * sin(theta_g) 
			/ (2. * pow(sin(theta_g/2.),2.0) * tan(theta_e))) - m_e;
		double deltaK   = Ecompton - eb;
		
		h_deltaK_old->Fill(deltaK);
		
		if(deltaK<deltaK_min_cut || deltaK>deltaK_max_cut) continue;
		
		double phi1     = atan2(loc_y1, loc_x1) * (180. / TMath::Pi());
		double phi2     = atan2(loc_y2, loc_x2) * (180. / TMath::Pi());
		double deltaPhi = fabs(phi2 - phi1);
		
		h_deltaPhi_old->Fill(deltaPhi);
		h_sumPhi_vs_delPhi_old->Fill(0.5*fabs(phi2-phi1), 0.5*(phi2+phi1));
		
		h_xy1->Fill(x1, y1);
		h_xy2->Fill(x2, y2);
		
		num_events_used++;
		
		h_e1->Fill(e1);
		h_e2->Fill(e2);
		
		//---------------------------------------------------//
		
		double r2 = pow(loc_x2-loc_x1,2.0) + pow(loc_y2-loc_y1,2.0);
		
		A += pow(loc_y1-loc_y2,2.0) / r2;
		B += ((loc_x1-loc_x2) * (loc_y2-loc_y1)) / r2;
		C += (loc_x1*loc_y2*loc_y2 + loc_x2*loc_y1*loc_y1 - loc_y1*loc_y2*(loc_x1+loc_x2))/ r2;
		
		F += pow(loc_x1-loc_x2,2.0) / r2;
		G += (loc_x1*loc_x1*loc_y2 + loc_x2*loc_x2*loc_y1 - loc_x1*loc_x2*(loc_y1+loc_y2)) / r2;
		
		double loc_x0 = (F*C - B*G) / (A*F - B*B);
		double loc_y0 = (A*G - B*C) / (A*F - B*B);
		
		if(n_count%10 == 0 && n_count > 1) {
			x0Vec.push_back(loc_x0);
			y0Vec.push_back(loc_y0);
			countVec.push_back(n_count);
		}
		n_count++;
	}
	
	int n_bins = (int)countVec.size();
	
	double *counts = new double[n_bins];
	double *x0s    = new double[n_bins];
	double *y0s    = new double[n_bins];
	
	for(int i=0; i<n_bins; i++) {
		counts[i] = countVec[i];
		x0s[i]    = x0Vec[i];
		y0s[i]    = y0Vec[i];
	}
	
	TLatex lat;
	
	TGraph *gx0 = new TGraph(n_bins, counts, x0s);
	gx0->SetMarkerColor( kBlack );
	gx0->SetMarkerStyle(3);
	gx0->GetXaxis()->SetTitle("number of Compton events" );
	gx0->GetYaxis()->SetTitle("x_{0} [cm]");
	gx0->SetTitle("Center x_{0}");
	gx0->GetYaxis()->SetTitleSize(0.075);
	gx0->GetYaxis()->SetTitleOffset(0.5);
	gx0->GetYaxis()->SetLabelSize(0.05);
	gx0->GetXaxis()->SetTitleSize(0.06);
	gx0->GetXaxis()->SetTitleOffset(0.8);
	gx0->GetXaxis()->SetLabelSize(0.05);
	
	TGraph *gy0 = new TGraph( n_bins, counts, y0s );
	gy0->SetMarkerColor( kBlack );
	gy0->SetMarkerStyle(3);
	gy0->GetXaxis()->SetTitle("number of Compton events");
	gy0->GetYaxis()->SetTitle("y_{0} [cm]");
	gy0->SetTitle("Center y_{0}");
	gy0->GetYaxis()->SetTitleSize(0.075);
	gy0->GetYaxis()->SetTitleOffset(0.5);
	gy0->GetYaxis()->SetLabelSize(0.05);
	gy0->GetXaxis()->SetTitleSize(0.06);
	gy0->GetXaxis()->SetTitleOffset(0.8);
	gy0->GetXaxis()->SetLabelSize(0.05);
	
	TCanvas *c = new TCanvas("c", "c", 1000, 600);
	c->Divide(1,2);
	
	gx0->GetYaxis()->SetRangeUser(-1.,0.2);
	gy0->GetYaxis()->SetRangeUser(-0.6,0.6);
	
	TPad *p1 = (TPad*)c->cd(1);
	p1->SetTickx(); p1->SetTicky();
	gx0->Draw( "AP" );
	lat.DrawLatexNDC(0.55, 0.2, Form("#scale[1.5]{x_{Beam} - x_{FCAL} = %.3f cm}", x0s[n_bins-1]));
	
	TPad *p2 = (TPad*)c->cd(2);
	p2->SetTickx(); p2->SetTicky();
	gy0->Draw( "AP" );
	lat.DrawLatexNDC(0.55, 0.2, Form("#scale[1.5]{y_{Beam} - y_{FCAL} = %.3f cm}", y0s[n_bins-1]));
	
	x0 = x0Vec[n_bins-1];
	y0 = y0Vec[n_bins-1];
	
	for(int ievt = 0; ievt < n_entries; ievt++) {
		
		t->GetEvent(ievt);
		if(ievt%100000 == 0) cout << "Processing event " << ievt << endl;
		
		if(!(MIN_ENERGY_CUT < eb && eb < MAX_ENERGY_CUT)) continue;
		if(e1 < 1.0 || e2 < 1.0) continue;
		
		/*
		First, we want to apply a fiducial cut to remove the inner part of the calorimeter. 
		In order to do so, we have to project both shower positions to the face of the 
		detectors. This requires knowledge of the beam position, and FCAL position, but that 
		obviously isn't yet known at this point. 
		It should be fine at this stage to just assume the beam is (0,0), and the FCAL position 
		is what's stored in the ccdb. If we're wrong about this, it will have a minimal effect
		on the final result, as it is only a cut to remove showers on the inner part, and these
		positions are not what's used to calculate the relative offset.
		*/
		
		double z_face = m_fcalZ - m_beamZ;
		
		double x1_cor = x1 - m_beamX_ccdb;
		double y1_cor = y1 - m_beamY_ccdb;
		double z1_cor = z1 - m_beamZ;
		double x2_cor = x2 - m_beamX_ccdb;
		double y2_cor = y2 - m_beamY_ccdb;
		double z2_cor = z2 - m_beamZ;
		
		double x1_face = m_beamX_ccdb + (x1_cor*z_face/z1_cor);
		double y1_face = m_beamY_ccdb + (y1_cor*z_face/z1_cor);
		double x2_face = m_beamX_ccdb + (x2_cor*z_face/z2_cor);
		double y2_face = m_beamY_ccdb + (y2_cor*z_face/z2_cor);
		
		x1_face -= m_fcalX_ccdb;
		y1_face -= m_fcalY_ccdb;
		x2_face -= m_fcalX_ccdb;
		y2_face -= m_fcalY_ccdb;
		
		if(((-1.*fid_cut)<x1_face&&x1_face<fid_cut) && 
			((-1.*fid_cut)<y1_face&&y1_face<fid_cut)) continue;
		if(((-1.*fid_cut)<x2_face&&x2_face<fid_cut) && 
			((-1.*fid_cut)<y2_face&&y2_face<fid_cut)) continue;
		
		/*
		The position stored in the tree comes directly from the DFCALShower_factory.
		This position already has the offset of the FCAL which is stored in the CCDB applied to it.
		Our goal in this analysis is to determine this offset (relative to the beam), so to start, 
		let's remove this initial offset correction:
		*/
		
		double loc_x1, loc_y1, loc_z1;
		double loc_x2, loc_y2, loc_z2;
		
		if(USE_FACE_POSITIONS) {
			
			// Option 1: Use shower positions projected to the surface of the FCAL:
			
			loc_x1 = x1_face - x0;
			loc_y1 = y1_face - y0;
			loc_z1 = m_fcalZ - m_beamZ;
			loc_x2 = x2_face - x0;
			loc_y2 = y2_face - y0;
			loc_z2 = m_fcalZ - m_beamZ;
			
		} else {
			
			// Option 2: Use shower positions as reconstructed inside the calorimeter:
			
			loc_x1 = x1 - m_fcalX_ccdb - x0;
			loc_y1 = y1 - m_fcalY_ccdb - y0;
			loc_z1 = z1 - m_beamZ;
			loc_x2 = x2 - m_fcalX_ccdb - x0;
			loc_y2 = y2 - m_fcalY_ccdb - y0;
			loc_z2 = z2 - m_beamZ;
		}
		
		//---------------------------------------------------//
		// Apply TOF Veto:
		
		int tof_veto = 0;
		if(tof_match1 && !tof_match2) tof_veto = 1;
		else if(tof_match2 && !tof_match1) tof_veto = 1;
		
		if(!tof_veto) continue;
		
		//---------------------------------------------------//
		// DeltaE and DeltaK Cuts:
		
		double deltaE = (e1+e2) - (eb+m_e);
		if(deltaE<deltaE_min_cut || deltaE>deltaE_max_cut) continue;
		
		double theta_g, theta_e;
		
		if(tof_match1) {
			theta_g = atan2(sqrt(loc_x2*loc_x2 + loc_y2*loc_y2), loc_z2);
			theta_e = atan2(sqrt(loc_x1*loc_x1 + loc_y1*loc_y1), loc_z1);
		} else {
			theta_g = atan2(sqrt(loc_x1*loc_x1 + loc_y1*loc_y1), loc_z1);
			theta_e = atan2(sqrt(loc_x2*loc_x2 + loc_y2*loc_y2), loc_z2);
		}
		
		double Ecompton = (m_e * sin(theta_g) 
			/ (2. * pow(sin(theta_g/2.),2.0) * tan(theta_e))) - m_e;
		double deltaK   = Ecompton - eb;
		
		h_deltaK_new->Fill(deltaK);
		
		if(deltaK<deltaK_min_cut || deltaK>deltaK_max_cut) continue;
		
		double phi1     = atan2(loc_y1, loc_x1) * (180. / TMath::Pi());
		double phi2     = atan2(loc_y2, loc_x2) * (180. / TMath::Pi());
		double deltaPhi = fabs(phi2 - phi1);
		
		h_deltaPhi_new->Fill(deltaPhi);
		h_sumPhi_vs_delPhi_new->Fill(0.5*fabs(phi2-phi1), 0.5*(phi2+phi1));
	}
	
	
	//----------   DeltaE   ----------//
	
	h_deltaE->SetTitle("#DeltaE");
	h_deltaE->GetXaxis()->SetTitle("E_{1} + E_{2} - E_{#gamma} [GeV]");
	h_deltaE->GetXaxis()->SetTitleOffset(1.2);
	h_deltaE->GetXaxis()->SetRangeUser(-1.0, 1.0);
	h_deltaE->SetLineColor(kBlack);
	
	TCanvas *ce = new TCanvas("ce", "ce", 600, 600);
	ce->SetTickx(); ce->SetTicky();
	h_deltaE->Draw();
	
	ce->Update();
	TLine *lE1 = new TLine(deltaE_min_cut, gPad->GetUymin(), deltaE_min_cut, gPad->GetUymax());
	lE1->SetLineColor(kRed); lE1->SetLineStyle(2); lE1->Draw("same");
	TLine *lE2 = new TLine(deltaE_max_cut, gPad->GetUymin(), deltaE_max_cut, gPad->GetUymax());
	lE2->SetLineColor(kRed); lE2->SetLineStyle(2); lE2->Draw("same");
	
	
	//----------   DeltaK   ----------//
	
	h_deltaK_old->SetTitle("#DeltaK");
	h_deltaK_old->GetXaxis()->SetTitle("E_{1} + E_{2} - E_{#gamma} [GeV]");
	h_deltaK_old->GetXaxis()->SetTitleOffset(1.2);
	h_deltaK_old->GetXaxis()->SetRangeUser(-1.0, 1.0);
	h_deltaK_old->SetLineColor(kBlack);
	
	TCanvas *ck = new TCanvas("ck", "ck", 600, 600);
	ck->SetTickx(); ck->SetTicky();
	h_deltaK_old->Draw();
	
	ck->Update();
	TLine *lK1 = new TLine(deltaK_min_cut, gPad->GetUymin(), deltaK_min_cut, gPad->GetUymax());
	lK1->SetLineColor(kRed); lK1->SetLineStyle(2); lK1->Draw("same");
	TLine *lK2 = new TLine(deltaK_max_cut, gPad->GetUymin(), deltaK_max_cut, gPad->GetUymax());
	lK2->SetLineColor(kRed); lK2->SetLineStyle(2); lK2->Draw("same");
	
	
	//----------   Shower energy and position   ----------//
	
	TCanvas *cE12 = new TCanvas("cE12","cE12",1200,800);
	cE12->Divide(2,1);
	
	h_e1->GetXaxis()->SetTitle("E_{1} [GeV]");
	h_e1->GetXaxis()->SetRangeUser(0.,5.);
	h_e1->SetLineColor(kBlack);
	h_e2->GetXaxis()->SetTitle("E_{1} [GeV]");
	h_e2->GetXaxis()->SetRangeUser(0.,5.);
	h_e2->SetLineColor(kBlack);
	
	cE12->cd(1);
	h_e1->Draw();
	cE12->cd(2);
	h_e2->Draw();
	
	TCanvas *cx = new TCanvas( "cx", "cx", 1200, 600 );
	cx->Divide( 2,1 );
	TPad *pxy1 = (TPad*)cx->cd(1);
	pxy1->SetTickx(); pxy1->SetTicky(); pxy1->SetLogz();
	h_xy1->Draw( "colz" );
	TPad *pxy2 = (TPad*)cx->cd(2);
	pxy2->SetTickx(); pxy2->SetTicky(); pxy2->SetLogz();
	h_xy2->Draw( "colz" );
	
	
	//----- DeltaPhi Comparison (Before vs. After Alignment) -----//
	
	h_deltaPhi_old->Rebin(4);
	h_deltaPhi_new->Rebin(4);
	
	TF1 *f_old = new TF1( "f_old", "gaus", 170., 190. );
	f_old->SetParameters( h_deltaPhi_old->GetMaximum(), 180., 7.5 );
	h_deltaPhi_old->Fit( "f_old", "R0QL" );
	f_old->SetRange( f_old->GetParameter(1) - f_old->GetParameter(2), 
		f_old->GetParameter(1) + f_old->GetParameter(2) );
	h_deltaPhi_old->Fit( "f_old", "R0QL" );
	f_old->SetRange( 150., 210. );
	f_old->SetLineStyle( 2 );
	
	TF1 *f_new = new TF1( "f_new", "gaus", 170., 190. );
	f_new->SetParameters( h_deltaPhi_new->GetMaximum(), 180., 7.5 );
	h_deltaPhi_new->Fit( "f_new", "R0QL" );
	f_new->SetRange( f_new->GetParameter(1) - f_new->GetParameter(2), 
		f_new->GetParameter(1) + f_new->GetParameter(2) );
	h_deltaPhi_new->Fit( "f_new", "R0QL" );
	f_new->SetRange( 150., 210. );
	f_new->SetLineStyle( 2 );
	
	TCanvas *cPhi = new TCanvas( "cPhi", "cPhi", 1200, 800 );
	cPhi->Divide( 2,1 );
	
	h_deltaPhi_old->GetXaxis()->SetTitle( "|#phi_{1} - #phi_{2}| [deg.]" );
	h_deltaPhi_old->SetLineColor( kBlack );
	h_deltaPhi_old->GetXaxis()->SetRangeUser( 140., 220. );
	h_deltaPhi_new->GetXaxis()->SetTitle( "|#phi_{1} - #phi_{2}| [deg.]" );
	h_deltaPhi_new->SetLineColor( kBlack );
	h_deltaPhi_new->GetXaxis()->SetRangeUser( 140., 220. );
	
	TPad *pPhi1 = (TPad*)cPhi->cd(1);
	pPhi1->SetTickx(); pPhi1->SetTicky();
	h_deltaPhi_old->Draw();
	f_old->Draw("same");
	
	TPad *pPhi2 = (TPad*)cPhi->cd(2);
	pPhi2->SetTickx(); pPhi2->SetTicky();
	h_deltaPhi_new->Draw();
	f_new->Draw("same");
	
	//----- DeltaK Comparison (Before vs. After Alignment) -----//
	
	TCanvas *cK = new TCanvas("cK", "cK", 1200, 800);
	cK->Divide(2,1);
	
	h_deltaK_old->Rebin(4);
	h_deltaK_new->Rebin(4);
	
	h_deltaK_old->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
	h_deltaK_old->SetLineColor(kBlack);
	//h_deltaK_old->GetXaxis()->SetRangeUser(-2., 2.);
	h_deltaK_new->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
	h_deltaK_new->SetLineColor(kBlack);
	//h_deltaK_new->GetXaxis()->SetRangeUser(-2., 2.);
	
	TF1 *f_deltaK_old = new TF1("f_deltaK_old", deltaK_fit, -3., 0.5, 6);
	//f_deltaK_old->SetParameters(110., -0.5, 0.15, 1., -1., 0.);
	f_deltaK_old->SetParameters(110., -0.5, 0.15, 20., 0., 0.);
	h_deltaK_old->Fit(f_deltaK_old, "R0QL");
	
	TF1 *f_deltaK_new = new TF1("f_deltaK_new", deltaK_fit, -3., 0.5, 6);
	//f_deltaK_new->SetParameters(110., -0.5, 0.15, 1., -1., 0.);
	f_deltaK_new->SetParameters(110., -0.5, 0.15, 20., 0., 0.);
	h_deltaK_new->Fit(f_deltaK_new, "R0QL");
	
	TPad *pK1 = (TPad*)cK->cd(1);
	pK1->SetTickx(); pK1->SetTicky();
	h_deltaK_old->Draw();
	f_deltaK_old->Draw("same");
	
	TPad *pK2 = (TPad*)cK->cd(2);
	pK2->SetTickx(); pK2->SetTicky();
	h_deltaK_new->Draw();
	f_deltaK_new->Draw("same");
	
	
	//------ 2D Comparison of DeltaPhi ------//
	
	TCanvas *cPhi2D = new TCanvas("cPhi2D", "cPhi2D", 1200, 600);
	cPhi2D->Divide(2,1);
	
	TPad *cPhi2D1 = (TPad*)cPhi2D->cd(1);
	cPhi2D1->SetTickx(); cPhi2D1->SetTicky(); cPhi2D1->SetGrid();
	h_sumPhi_vs_delPhi_old->GetXaxis()->SetRangeUser(70., 110.);
	h_sumPhi_vs_delPhi_old->Draw("COL");
	
	TPad *cPhi2D2 = (TPad*)cPhi2D->cd(2);
	cPhi2D2->SetTickx(); cPhi2D2->SetTicky(); cPhi2D2->SetGrid();
	h_sumPhi_vs_delPhi_new->GetXaxis()->SetRangeUser(70., 110.);
	h_sumPhi_vs_delPhi_new->Draw("COL");
	
	
	cout << "\n\n";
	cout << "BEAM - FCAL X: " << x0 << endl;
	cout << "BEAM - FCAL Y: " << y0 << endl;
	cout << "    BEAM X (Lab Frame): " << (x0 + m_fcalX) << endl;
	cout << "    BEAM Y (Lab Frame): " << (y0 + m_fcalY) << endl;
	cout << "    FCAL X (Lab Frame): " << (m_fcalX) << endl;
	cout << "    FCAL Y (Lab Frame): " << (m_fcalY) << endl;
	//cout << "    FCAL X (Lab Frame): " << (m_beamX - x0) << endl;
	//cout << "    FCAL Y (Lab Frame): " << (m_beamY - y0) << endl;
	cout << "\n\n";
	
	
	
	return;
}



Double_t bkgd_fit( Double_t *x, Double_t *par ) 
{
	
	Double_t xx  =   x[0];
	
	if( !(xx < par[4] || xx > 0.5) ) {
		TF1::RejectPoint();
		return 0.;
	}
	
	Double_t p0  = par[0];
	Double_t p1  = par[1];
	Double_t p2  = par[2];
	Double_t p3  = par[3];
	
	Double_t f   = p0 + p1*xx + p2*xx*xx + p3*xx*xx*xx;
	
	
	return f;
}



Double_t crys_ball_fit( Double_t *x, Double_t *par ) 
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
	
	
	Double_t A   = 	pow( n/fabs(a), n ) * exp( -0.5*pow(fabs(a),2.0) );
	Double_t B   =  ( n/fabs(a) ) - fabs(a);
	
	Double_t loc_x = ( xx - mu ) / sig;
	Double_t f;
	
	if( loc_x > -a ) {
		f = N * exp( -0.5*pow(loc_x,2.0) );
	} else {
		f = N * A * pow( B - loc_x, -n );
	}
	
	f += p0 + p1*xx + p2*xx*xx + p3*xx*xx*xx;
	
	
	return f;
}



Double_t deltaK_fit_exp( Double_t *x, Double_t *par ) 
{
	Double_t xx  =   x[0];
	
	Double_t A   = par[0];
	Double_t mu  = par[1];
	Double_t sig = par[2];
	
	Double_t p0  = par[3];
	Double_t p1  = par[4];
	Double_t p2  = par[5];
	
	Double_t f_bkgd = p0*exp(p1*xx + p2);
	Double_t f_gaus = A*exp( -0.5*pow((xx-mu)/sig,2.0) );
	
	Double_t f = f_bkgd + f_gaus;
	
	return f;
}



Double_t deltaK_fit( Double_t *x, Double_t *par ) 
{
	Double_t xx  =   x[0];
	
	Double_t A   = par[0];
	Double_t mu  = par[1];
	Double_t sig = par[2];
	
	Double_t p0  = par[3];
	Double_t p1  = par[4];
	Double_t p2  = par[5];
	
	Double_t f_bkgd = p0 + p1*xx + p2*xx*xx;
	Double_t f_gaus = A*exp( -0.5*pow((xx-mu)/sig,2.0) );
	
	Double_t f = f_bkgd + f_gaus;
	
	return f;
}
