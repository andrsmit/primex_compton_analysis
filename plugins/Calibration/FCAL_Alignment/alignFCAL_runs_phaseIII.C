
Double_t bkgd_fit(Double_t *x, Double_t *par);
Double_t crys_ball_fit(Double_t *x, Double_t *par);
Double_t deltaK_fit(Double_t *x, Double_t *par);
Double_t deltaK_fit_exp(Double_t *x, Double_t *par);

void alignFCAL_runs() {
	
	gStyle->SetPalette(1);
	
	char pathName[256];
	sprintf(pathName, "/work/halld/home/andrsmit/ccal_calib/phaseIII/alignment/fcalTrees");
	
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
	
	//----------  Begin Analysis  ----------//
	
	vector<double> runVec, x0Vec, y0Vec;
	
	for(int irun=110600; irun<=110621; irun++) {
		
		char infile[256];
		sprintf(infile, "%s/%d.root", pathName, irun);
		if(gSystem->AccessPathName(infile)) continue;
		
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
		//cout << "Total number of entries: " << n_entries << endl;
		
		double A = 0., B = 0., C = 0., F = 0., G = 0.;
		
		int num_events_used = 0;
		double loc_x0 = 0., loc_y0 = 0.;
		
		for(int ievt = 0; ievt < n_entries; ievt++) {
			
			t->GetEvent(ievt);
			//if(ievt%100000 == 0) cout << "Processing event " << ievt << endl;
			
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
			if(deltaK<deltaK_min_cut || deltaK>deltaK_max_cut) continue;
			
			num_events_used++;
			
			//---------------------------------------------------//
			
			double r2 = pow(loc_x2-loc_x1,2.0) + pow(loc_y2-loc_y1,2.0);
			
			A += pow(loc_y1-loc_y2,2.0) / r2;
			B += ((loc_x1-loc_x2) * (loc_y2-loc_y1)) / r2;
			C += (loc_x1*loc_y2*loc_y2 + loc_x2*loc_y1*loc_y1 - loc_y1*loc_y2*(loc_x1+loc_x2))/ r2;
			
			F += pow(loc_x1-loc_x2,2.0) / r2;
			G += (loc_x1*loc_x1*loc_y2 + loc_x2*loc_x2*loc_y1 - loc_x1*loc_x2*(loc_y1+loc_y2)) / r2;
			
			loc_x0 = (F*C - B*G) / (A*F - B*B);
			loc_y0 = (A*G - B*C) / (A*F - B*B);
		}
		
		if(num_events_used > 100) {
			runVec.push_back((double)irun);
			x0Vec.push_back(loc_x0 + m_fcalX);
			y0Vec.push_back(loc_y0 + m_fcalY);
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
		x0s[i]  = x0Vec[i];
		y0s[i]  = y0Vec[i];
		
		avg_counter += 1.;
		x0_mean += x0Vec[i];
		y0_mean += y0Vec[i];
	}
	
	x0_mean /= avg_counter;
	y0_mean /= avg_counter;
	
	
	TGraph *gx0 = new TGraph(n_runs, runs, x0s);
	gx0->GetXaxis()->SetTitle("Run Number (-110000)");
	gx0->GetXaxis()->SetTitleSize(0.05);
	gx0->GetXaxis()->SetTitleOffset(1.);
	gx0->GetYaxis()->SetTitle("x_{Beam} [cm]");
	gx0->GetYaxis()->SetTitleSize(0.055);
	gx0->GetYaxis()->SetTitleOffset(0.75);
	gx0->SetTitle("Beam-X Position from FCAL-Compton Alignment");
	gx0->SetMarkerStyle(8);
	gx0->SetMarkerColor(kBlack);
	
	TF1 *fx0 = new TF1("fx0", "[0]", runs[0], runs[n_runs-1]);
	fx0->SetLineColor(kRed);
	fx0->SetParameter(0,x0_mean);
	
	TGraph *gy0 = new TGraph(n_runs, runs, y0s);
	gy0->GetXaxis()->SetTitle("Run Number (-110000)");
	gy0->GetXaxis()->SetTitleSize(0.05);
	gy0->GetXaxis()->SetTitleOffset(1.);
	gy0->GetYaxis()->SetTitle("y_{Beam} - y_{CCAL} [cm]");
	gy0->GetYaxis()->SetTitleSize(0.055);
	gy0->GetYaxis()->SetTitleOffset(0.75);
	gy0->SetTitle("Beam-Y Position from FCAL-Compton Alignment");
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
	cout << "x0: " << x0_mean << "(x_Beam = " << (x0_mean) << ")" << endl;
	cout << "y0: " << y0_mean << "(y_Beam = " << (y0_mean) << ")" << endl;
	cout << "\n\n\n";
	
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
