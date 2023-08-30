
void alignCCAL_runs() {
	
	char pathName[256];
	sprintf(pathName, "/work/halld/home/andrsmit/ccal_calib/phaseIII/alignment/ccalTrees");
	
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
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetPalette(1);
	
	//----------  Begin Analysis  ----------//
	
	vector<double> runVec, x0Vec, y0Vec;
	runVec.clear();
	x0Vec.clear();
	y0Vec.clear();
	
	for(int irun=110600; irun<=110621; irun++) {
		
		char infile[256];
		sprintf(infile, "%s/%d.root", pathName, irun);
		if(gSystem->AccessPathName(infile)) continue;
		
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
		//cout << "Total number of entries: " << n_entries << endl;
		
		double A = 0.;
		double B = 0.;
		double C = 0.;
		double F = 0.;
		double G = 0.;
		
		int num_events_used = 0;
		double loc_x0 = 0., loc_y0 = 0.;
		
		for(int ievt = 1; ievt < n_entries; ievt++) {
			
			t->GetEvent(ievt);
			//if(ievt%1000000 == 0) cout << "  Processing event " << ievt << endl; 
			
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
			'loc_x2' and 'loc_y2' give the position of the CCAL shower relative to the center of 
			the CCAL.
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
			
			loc_x0 = (C*F - B*G) / (A*F - B*B);
			loc_y0 = (B*C - A*G) / (B*B - A*F);
		}
		
		if(num_events_used > 10000) {
			runVec.push_back((double)irun);
			x0Vec.push_back(loc_x0);
			y0Vec.push_back(loc_y0);
		}
		
		cout << "  num_events_used = " << num_events_used << endl;
	}
	
	int n_runs = (int)runVec.size();
	
	double *runs = new double[n_runs];
	double *x0s  = new double[n_runs];
	double *y0s  = new double[n_runs];
	
	double avg_counter = 0.;
	double x0_mean = 0.;
	double y0_mean = 0.;
	
	for(int i=0; i<n_runs; i++) {
		runs[i] = runVec[i] - 110000;
		x0s[i]  = x0Vec[i] + m_beamX_aligned;
		y0s[i]  = y0Vec[i] + m_beamY_aligned;
		
		avg_counter += 1.;
		x0_mean += (x0Vec[i] + m_beamX_aligned);
		y0_mean += (y0Vec[i] + m_beamY_aligned);
	}
	
	x0_mean /= avg_counter;
	y0_mean /= avg_counter;
	
	TGraph *gx0 = new TGraph(n_runs, runs, x0s);
	gx0->GetXaxis()->SetTitle("Run Number (-110000)");
	gx0->GetXaxis()->SetTitleSize(0.05);
	gx0->GetXaxis()->SetTitleOffset(1.);
	gx0->GetYaxis()->SetTitle("x_{CCAL} [cm]");
	gx0->GetYaxis()->SetTitleSize(0.055);
	gx0->GetYaxis()->SetTitleOffset(0.75);
	gx0->SetTitle("CCAL-X Position from Compton Alignment");
	gx0->SetMarkerStyle(8);
	gx0->SetMarkerColor(kBlack);
	
	TF1 *fx0 = new TF1("fx0", "[0]", runs[0], runs[n_runs-1]);
	fx0->SetLineColor(kRed);
	fx0->SetParameter(0,x0_mean);
	
	TGraph *gy0 = new TGraph(n_runs, runs, y0s);
	gy0->GetXaxis()->SetTitle("Run Number (-110000)");
	gy0->GetXaxis()->SetTitleSize(0.05);
	gy0->GetXaxis()->SetTitleOffset(1.);
	gy0->GetYaxis()->SetTitle("y_{CCAL} [cm]");
	gy0->GetYaxis()->SetTitleSize(0.055);
	gy0->GetYaxis()->SetTitleOffset(0.75);
	gy0->SetTitle("CCAL-Y Position from Compton Alignment");
	gy0->SetMarkerStyle(8);
	gy0->SetMarkerColor(kBlack);
	
	TF1 *fy0 = new TF1("fy0", "[0]", runs[0], runs[n_runs-1]);
	fy0->SetLineColor(kRed);
	fy0->SetParameter(0,y0_mean);
	
	//gx0->GetYaxis()->SetRangeUser(-1.,1.);
	//gy0->GetYaxis()->SetRangeUser(-1.,1.);
	
	TLatex lat;
	
	TCanvas *cxy = new TCanvas("cxy", "cxy", 1000, 1000);
	cxy->Divide(1,2);
	cxy->cd(1);
	gx0->Draw("AP");
	fx0->Draw("same");
	
	cxy->cd(2);
	gy0->Draw("AP");
	fy0->Draw("same");
	
	cout << "\n\n\n";
	cout << "x0: " << x0_mean << "(x_CCAL = " << (x0_mean) << ")" << endl;
	cout << "y0: " << y0_mean << "(y_CCAL = " << (y0_mean) << ")" << endl;
	cout << "\n\n\n";
	
	return;
}
