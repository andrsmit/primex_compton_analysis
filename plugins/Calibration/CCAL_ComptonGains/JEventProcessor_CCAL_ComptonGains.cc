
/*-------------------------------------------------------------------------
	
	This plugin can be used to calibrate CCAL gain constants using 
	Compton Scattering events.
	
	If running on simulation, you will probably want to set the 
	'BYPASS_TRIG' parameter to 1 (default 0). This can be done at run time
	by doing (for example):
		
		hd_root -PPLUGINS=CCAL_ComptonGains -PCCAL_GAINS:BYPASS_TRIG=1 <input_file>.hddm
	
--------------------------------------------------------------------------*/

#include "JEventProcessor_CCAL_ComptonGains.h"

extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
    	app->AddProcessor(new JEventProcessor_CCAL_ComptonGains());
}
} // "C"


//==========================================================
//
//   Constructor
//
//==========================================================

JEventProcessor_CCAL_ComptonGains::JEventProcessor_CCAL_ComptonGains()
{
	// Set defaults:
	
	BYPASS_TRIG = false;
	
	MIN_FCAL_ENERGY_CUT  =  0.35;
	MIN_CCAL_ENERGY_CUT  =  3.00;
	MIN_BEAM_ENERGY_CUT  =  5.00;
	
	COPL_CUT             = 50.0;
	DELTA_K_CUT          =  0.5;
	
	MIN_CALIB_ENERGY     =  4.75;
	MAX_CALIB_ENERGY     =  5.25;
	
	gPARMS->SetDefaultParameter("CCAL_ComptonGains:BYPASS_TRIG", BYPASS_TRIG, 
			"Set to true to bypass trigger information (for simulation)");
	
	gPARMS->SetDefaultParameter("CCAL_ComptonGains:MIN_FCAL_ENERGY_CUT", MIN_FCAL_ENERGY_CUT,
			"Minimum energy of FCAL Shower");
	gPARMS->SetDefaultParameter("CCAL_ComptonGains:MIN_CCAL_ENERGY_CUT", MIN_CCAL_ENERGY_CUT,
			"Minimum energy of CCAL Shower");
	gPARMS->SetDefaultParameter("CCAL_ComptonGains:MIN_BEAM_ENERGY_CUT", MIN_BEAM_ENERGY_CUT,
			"Minimum beam photon energy");
	
	gPARMS->SetDefaultParameter("CCAL_ComptonGains:COPL_CUT", COPL_CUT,
			"Coplanarity Cut (degrees)");
	gPARMS->SetDefaultParameter("CCAL_ComptonGains:DELTA_K_CUT", DELTA_K_CUT);
	
	gPARMS->SetDefaultParameter("CCAL_ComptonGains:MIN_CALIB_ENERGY", MIN_CALIB_ENERGY);
	gPARMS->SetDefaultParameter("CCAL_ComptonGains:MAX_CALIB_ENERGY", MAX_CALIB_ENERGY);
}


//==========================================================
//
//   init
//
//==========================================================

jerror_t JEventProcessor_CCAL_ComptonGains::init(void)
{
  	
	TDirectory *dir_CCAL_ComptonGains = new TDirectoryFile("CCAL_ComptonGains", 
		"CCAL_ComptonGains");
	dir_CCAL_ComptonGains->cd();
	
	// Trigger histograms:
	
	h_trig            = new TH1F("trig",   "GTP Bits",        34, -0.5, 33.5);
	h_fptrig          = new TH1F("fptrig", "FP Trigger Bits", 34, -0.5, 33.5);
	
	// TOF timing histograms:
	
	h_tof_rf_dt       = new TH1F("tof_rf_dt",     "t_{TOF} - t_{RF}; [ns]", 5000, -100., 100.);
	h_tof_rf_dt_cut   = new TH1F("tof_rf_dt_cut", "t_{TOF} - t_{RF}; [ns]", 5000, -100., 100.);
	
	h_fcal_tof_dt     = new TH1F("fcal_tof_dt",     "t_{FCAL} - t_{TOF}; [ns]", 5000, -100., 100.);
	h_fcal_tof_dt_cut = new TH1F("fcal_tof_dt_cut", "t_{FCAL} - t_{TOF}; [ns]", 5000, -100., 100.);
	
	// RF Timing Histograms
	
	h_fcal_rf_dt   = new TH1F("fcal_rf_dt", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	h_ccal_rf_dt   = new TH1F("ccal_rf_dt", "t_{CCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	h_beam_rf_dt   = new TH1F("beam_rf_dt", "t_{BEAM} - t_{RF}; [ns]", 2000, -100., 100.);
	
	h_ccal_beam_dt = new TH1F("ccal_beam_dt", "t_{#gamma} - t_{CCAL}; [ns]", 2000, -100., 100.);
	h_fcal_beam_dt = new TH1F("fcal_beam_dt", "t_{#gamma} - t_{FCAL}; [ns]", 2000, -100., 100.);
	
	// Compton scattering distribtions:
	
	h_deltaT   = new TH1F("deltaT", "#DeltaT; t_{FCAL} - t_{CCAL} [ns]", 2000, -20., 20.);
	h_deltaPhi = new TH1F("deltaPhi", "#Delta#phi; #phi_{FCAL} - #phi_{CCAL} [deg.]", 
		3600, 0., 360.);
	h_deltaE   = new TH1F("deltaE", "#DeltaE; E_{FCAL} + E_{CCAL} - E_{#gamma} [GeV]", 
		4000, -4., 4.);
	h_deltaK   = new TH1F("deltaK", "#DeltaK; E_{Compton} - E_{#gamma} [GeV]", 
		4000, -4., 4.);
	
	// Histograms used for calibration:
	
	for(int ihist = 0; ihist < n_cuts; ihist++) {
		
		// Elasticity ratios vs channel number:
		
		h_CompRatio[ihist]   = new TH2F(Form("CompRatio_%d",ihist),
			"Compton Energy Ratio; CCAL Channel Number; E_{CCAL} / E_{Calc.}",
			144, -0.5, 143.5, 2000, 0., 2.);
		h_CompRatio_g[ihist] = new TH2F(Form("CompRatio_g_%d",ihist),
			"Compton Energy Ratio (Photons Only); CCAL Channel Number; E_{CCAL} / E_{Calc.}",
			144, -0.5, 143.5, 2000, 0., 2.);
		h_CompRatio_e[ihist] = new TH2F(Form("CompRatio_e_%d",ihist),
			"Compton Energy Ratio (Electrons Only); CCAL Channel Number; E_{CCAL} / E_{Calc.}",
			144, -0.5, 143.5, 2000, 0., 2.);
		
		h_ElasRatio[ihist]   = new TH2F(Form("ElasRatio_%d",ihist),
			"Elasticity Ratio; CCAL Channel Number; E_{CCAL} / (E_{#gamma} - E_{FCAL})",
			144, -0.5, 143.5, 2000, 0., 2.);
		h_ElasRatio_g[ihist] = new TH2F(Form("ElasRatio_g_%d",ihist),
			"Elasticity Ratio (Photons Only); CCAL Channel Number; E_{CCAL} / (E_{#gamma} - E_{FCAL})",
			144, -0.5, 143.5, 2000, 0., 2.);
		h_ElasRatio_e[ihist] = new TH2F(Form("ElasRatio_e_%d",ihist),
			"Elasticity Ratio (Electrons Only); CCAL Channel Number; E_{CCAL} / (E_{#gamma} - E_{FCAL})",
			144, -0.5, 143.5, 2000, 0., 2.);
		
		// Elasticity ratios vs channel number for selected energy range:
		
		h_CompRatio_cut[ihist]   = new TH2F(Form("CompRatio_cut_%d",ihist),
			"Compton Energy Ratio; CCAL Channel Number; E_{CCAL} / E_{Calc.}",
			144, -0.5, 143.5, 2000, 0., 2.);
		h_CompRatio_g_cut[ihist] = new TH2F(Form("CompRatio_g_cut_%d",ihist),
			"Compton Energy Ratio (Photons Only); CCAL Channel Number; E_{CCAL} / E_{Calc.}",
			144, -0.5, 143.5, 2000, 0., 2.);
		h_CompRatio_e_cut[ihist] = new TH2F(Form("CompRatio_e_cut_%d",ihist),
			"Compton Energy Ratio (Electrons Only); CCAL Channel Number; E_{CCAL} / E_{Calc.}",
			144, -0.5, 143.5, 2000, 0., 2.);
		
		h_ElasRatio_cut[ihist]   = new TH2F(Form("ElasRatio_cut_%d",ihist),
			"Elasticity Ratio; CCAL Channel Number; E_{CCAL} / (E_{#gamma} - E_{FCAL})",
			144, -0.5, 143.5, 2000, 0., 2.);
		h_ElasRatio_g_cut[ihist] = new TH2F(Form("ElasRatio_g_cut_%d",ihist),
			"Elasticity Ratio (Photons Only); CCAL Channel Number; E_{CCAL} / (E_{#gamma} - E_{FCAL})",
			144, -0.5, 143.5, 2000, 0., 2.);
		h_ElasRatio_e_cut[ihist] = new TH2F(Form("ElasRatio_e_cut_%d",ihist),
			"Elasticity Ratio (Electrons Only); CCAL Channel Number; E_{CCAL} / (E_{#gamma} - E_{FCAL})",
			144, -0.5, 143.5, 2000, 0., 2.);
		
		// Elasticity ratios vs energy (nonlinearity):
		
		h_CompRatio_vs_E[ihist] = new TH2F(Form("CompRatio_vs_E_%d",ihist),
			"CCAL Compton Energy Ratio; E_{Calc.} [GeV]; E_{CCAL} / E_{Calc.}",
			120, 0., 12., 2000, 0., 2.);
		h_CompRatio_vs_E_g[ihist] = new TH2F(Form("CompRatio_vs_E_g_%d",ihist),
			"CCAL Compton Energy Ratio (Photons Only); E_{Calc.} [GeV]; E_{CCAL} / E_{Calc.}",
			120, 0., 12., 2000, 0., 2.);
		h_CompRatio_vs_E_e[ihist] = new TH2F(Form("CompRatio_vs_E_e_%d",ihist),
			"CCAL Compton Energy Ratio (Electrons Only); E_{Calc.} [GeV]; E_{CCAL} / E_{Calc.}",
			120, 0., 12., 2000, 0., 2.);
		
		h_ElasRatio_vs_E[ihist] = new TH2F(Form("ElasRatio_vs_E_%d",ihist),
			"CCAL Elasticity Ratio; E_{#gamma} - E_{FCAL} [GeV]; E_{CCAL} / (E_{#gamma} - E_{FCAL})",
			120, 0., 12., 2000, 0., 2.);
		h_ElasRatio_vs_E_g[ihist] = new TH2F(Form("ElasRatio_vs_E_g_%d",ihist),
			"CCAL Elasticity Ratio (Photons Only); E_{#gamma} - E_{FCAL} [GeV]; E_{CCAL} / (E_{#gamma} - E_{FCAL})",
			120, 0., 12., 2000, 0., 2.);
		h_ElasRatio_vs_E_e[ihist] = new TH2F(Form("ElasRatio_vs_E_e_%d",ihist),
			"CCAL Elasticity Ratio (Electrons Only); E_{#gamma} - E_{FCAL} [GeV]; E_{CCAL} / (E_{#gamma} - E_{FCAL})",
			120, 0., 12., 2000, 0., 2.);
	}
	
	// Distribution of FCAL and CCAL showers:
	
	h_xy_fcal   = new TH2F("xy_fcal", 
		"Distribution of FCAL Showers; x_{FCAL} [cm]; y_{FCAL} [cm]", 
		500, -125., 125., 500, -125., 125.);
	h_xy_ccal   = new TH2F("xy_ccal", 
		"Distribution of CCAL Showers; x_{CCAL} [cm]; y_{CCAL} [cm]", 
		500, -12.5, 12.5, 500, -12.5, 12.5);
	
	dir_CCAL_ComptonGains->cd("../");
	
	return NOERROR;
}


///==========================================================
//
//   brun
//
//==========================================================

jerror_t JEventProcessor_CCAL_ComptonGains::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	DGeometry*   dgeom = NULL;
  DApplication* dapp = dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  if(dapp)     dgeom = dapp->GetDGeometry(runnumber);
	
	if(dgeom){
		dgeom->GetTargetZ(m_beamZ);
		dgeom->GetFCALPosition(m_fcalX, m_fcalY, m_fcalZ);
		dgeom->GetCCALPosition(m_ccalX, m_ccalY, m_ccalZ);
  } else{
		cerr << "No geometry accessbile to ccal_timing monitoring plugin." << endl;
		return RESOURCE_UNAVAILABLE;
  }
	
	jana::JCalibration *jcalib = japp->GetJCalibration(runnumber);
  std::map<string, float> beam_spot;
  jcalib->Get("PHOTON_BEAM/beam_spot", beam_spot);
  m_beamX = beam_spot.at("x");
  m_beamY = beam_spot.at("y");
	
	// apply correction to FCAL and CCAL positions from Compton-scattering alignment procedure:
	
	if(runnumber>60000 && runnumber<69999) {
		
		// PhaseI:
		
		// (2/4/2024): Correction to alignment after Simon updated beam spot with new CDC alignment:
		
		m_fcalX_new =  0.455;
		m_fcalY_new = -0.032;
		
		m_ccalX_new = -0.082;
		if(runnumber<61483) 
			m_ccalY_new =  0.061;
		else
			m_ccalY_new =  0.051;
		
		if(runnumber<61483) {
			m_beamX =  0.027;
			m_beamY = -0.128;
		} else if(runnumber<61774) {
			m_beamX =  0.001;
			m_beamY = -0.077;
		} else {
			m_beamX =  0.038;
			m_beamY = -0.095;
		}
		
		/*
		m_fcalX_new =  0.617;
		m_fcalY_new = -0.002;
		
		if(runnumber<61483) {
			m_ccalX_new = 0.083;
			m_ccalY_new = 0.148;
		} else {
			m_ccalX_new = 0.083;
			m_ccalY_new = 0.119;
		}
		
		if(runnumber<61483) {
			m_beamX =  0.190;
			m_beamY = -0.074;
		} else if(runnumber<61774) {
			m_beamX =  0.165;
			m_beamY = -0.024;
		} else {
			m_beamX =  0.202;
			m_beamY = -0.042;
		}
		*/
		
	} else if(runnumber>80000 && runnumber<89999) {
		
		// PhaseII:
		
		m_fcalX_new =  0.455;
		m_fcalY_new = -0.032;
		m_ccalX_new = m_ccalX;
		m_ccalY_new = m_ccalY;
		
	} else if(runnumber>110000) {
		
		// Phase III:
		
		m_fcalX_new =  0.455;
		m_fcalY_new = -0.032;
		m_ccalX_new = m_ccalX;
		m_ccalY_new = m_ccalY;
	}
	
	m_fcal_correction.SetXYZ(m_fcalX_new-m_fcalX, m_fcalY_new-m_fcalY, 0.);
	m_ccal_correction.SetXYZ(m_ccalX_new-m_ccalX, m_ccalY_new-m_ccalY, 0.);
	
	return NOERROR;
}


//==========================================================
//
//   evnt
//
//==========================================================

jerror_t JEventProcessor_CCAL_ComptonGains::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
	//--------------------------------------------------------------------//
	// Check trigger:
	
	bool is_ccal_trig = false;
	uint32_t trigmask, fp_trigmask;
	
	if(BYPASS_TRIG) {}
	else {
		const DL1Trigger *locTrig = NULL;
		try {
			eventLoop->GetSingle(locTrig);
		} catch (...) {}
		if (locTrig == NULL) { return NOERROR; }
		
		trigmask    = locTrig->trig_mask;
		fp_trigmask = locTrig->fp_trig_mask;
		
		// only select events where the CCAL trigger bit is set:
		if(trigmask & (1 << 0)) is_ccal_trig = true;
	}
	
	//--------------------------------------------------------------------//
	// Check RF time:
	
	const DEventRFBunch *locRFBunch = NULL;
	try { 
	  	eventLoop->GetSingle(locRFBunch, "CalorimeterOnly");
	} catch (...) { return NOERROR; }
	double locRFTime = locRFBunch->dTime;
	
	//--------------------------------------------------------------------//
	// Get releveant data objects:
	
	DVector3 locVertex;
	locVertex.SetXYZ(m_beamX, m_beamY, m_beamZ);
	
	vector<const DBeamPhoton*> locBeamPhotons;
	eventLoop->Get(locBeamPhotons);
	
	vector<const DCCALShower*> locCCALShowers;
	eventLoop->Get(locCCALShowers);
	
	vector<const DFCALShower*> locFCALShowers;
	eventLoop->Get(locFCALShowers);
	
	vector<const DTOFPoint*> locTOFPoints;
	eventLoop->Get(locTOFPoints);
	
	//--------------------------------------------------------------------//
	// Fill lock for multi-threaded running:
	
	japp->RootFillLock(this);  // Acquire root lock
	
	//--------------------------------------------------------------------//
	// Fill histograms with trigger information:
	
	for(int ibit = 0; ibit < 34; ibit++) {
		if(trigmask & (1 << ibit))    h_trig->Fill(ibit);
		if(fp_trigmask & (1 << ibit)) h_fptrig->Fill(ibit);
	}
	
	//--------------------------------------------------------------------//
	// Skip events that aren't triggered by the CCAL:
	
	if(fp_trigmask || !is_ccal_trig || locRFBunch->dNumParticleVotes < 2) {
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	
	//--------------------------------------------------------------------//
	// Select events with correct multiplicities:
	
	int locNFCALShowers = 0, locNCCALShowers = 0;
	
	vector<const DFCALShower*> locGoodFCALShowers;
	locGoodFCALShowers.clear();
	vector<const DCCALShower*> locGoodCCALShowers;
	locGoodCCALShowers.clear();
	
	for(vector<const DFCALShower*>::const_iterator show = locFCALShowers.begin(); 
		show != locFCALShowers.end(); show++) {
		
		DVector3 loc_pos = (*show)->getPosition_log() - locVertex + m_fcal_correction;
		double loc_rf_dt = (*show)->getTime() - (loc_pos.Mag()/c) - locRFTime;
		
		h_fcal_rf_dt->Fill(loc_rf_dt);
		
		double loc_e = (*show)->getEnergy();
		
		if(fabs(loc_rf_dt)<RF_TIME_CUT) {
			locNFCALShowers++;
			if(fabs(loc_rf_dt)<2.0 && loc_e>MIN_FCAL_ENERGY_CUT) {
				locGoodFCALShowers.push_back((*show));
			}
		}
	}
	
	for(vector<const DCCALShower*>::const_iterator show = locCCALShowers.begin(); 
		show != locCCALShowers.end(); show++) {
		
		DVector3 loc_pos((*show)->x1, (*show)->y1, (*show)->z);
		loc_pos += (m_ccal_correction - locVertex);
		double loc_rf_dt = (*show)->time - (loc_pos.Mag()/c) - locRFTime;
		
		h_ccal_rf_dt->Fill(loc_rf_dt);
		
		double loc_e = (*show)->E;
		
		if(fabs(loc_rf_dt)<RF_TIME_CUT) {
			locNCCALShowers++;
			if(fabs(loc_rf_dt)<2.0 && loc_e>MIN_CCAL_ENERGY_CUT) {
				locGoodCCALShowers.push_back((*show));
			}
		}
	}
	
	/*
	if(locNFCALShowers!=1 || locFCCALShowers!=1) {
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	*/
	
	//--------------------------------------------------------------------//
	// Evaluate 2-cluster pair in the FCAL:
	
	for(vector<const DCCALShower*>::const_iterator show1 = locGoodCCALShowers.begin(); 
		show1 != locGoodCCALShowers.end(); show1++) {
		
		DVector3 pos1((*show1)->x1, (*show1)->y1, (*show1)->z);
		pos1 += (m_ccal_correction - locVertex);
		
		double e1     = (*show1)->E;
		double t1     = (*show1)->time - (pos1.Mag()/c);
		double phi1   = pos1.Phi() * (180./TMath::Pi());
		double theta1 = pos1.Theta();
		
		int idmax = (*show1)->idmax;
		
		for(vector<const DFCALShower*>::const_iterator show2 = locGoodFCALShowers.begin(); 
			show2 != locGoodFCALShowers.end(); show2++) {
			
			DVector3 pos2 = (*show2)->getPosition_log() - locVertex + m_fcal_correction;
			
			double e2     = (*show2)->getEnergy();
			double t2     = (*show2)->getTime() - (pos2.Mag()/c);
			double phi2   = pos2.Phi() * (180./TMath::Pi());
			double theta2 = pos2.Theta();
			
			double tof_dx, tof_dy, tof_dt;
			check_TOF_match(pos2, locRFTime, locVertex, locTOFPoints, tof_dx, tof_dy, tof_dt, 1.0);
			
			int loc_tof_match = 0;
			if(sqrt(pow(tof_dx,2.0) + pow(tof_dy,2.0)) < 8.0) loc_tof_match = 1;
			
			//--------------------------------//
			
			double deltaT   = t1 - t2;
			double deltaPhi = fabs(phi1-phi2);
			
			h_deltaT->Fill(deltaT);
			h_deltaPhi->Fill(deltaPhi);
			
			//------------------------------------------------------------//
			// Apply a Coplanarity cut:
			
			if(fabs(deltaPhi-180.) > COPL_CUT) continue;
			
			//------------------------------------------------------------//
			// Loop over beam photons:
			
			for(vector<const DBeamPhoton*>::const_iterator gam = locBeamPhotons.begin();
				gam != locBeamPhotons.end(); gam++) {
				
				double eb = (*gam)->lorentzMomentum().E();
				if(eb < MIN_BEAM_ENERGY_CUT) continue;
				
				double tb    = (*gam)->time();
				double brfdt = tb - locRFTime;
				
				h_beam_rf_dt->Fill(brfdt);
				h_ccal_beam_dt->Fill(t1-tb);
				h_fcal_beam_dt->Fill(t2-tb);
				
				double fill_weight = 0.;
				
				if(fabs(brfdt) < 2.004)
					fill_weight =  1.0;
				else if(( -(6.012 + 5.*4.008) <= brfdt && brfdt <= -(6.012 + 3.*4.008))
					|| ((6.012 + 3.*4.008) <= brfdt && brfdt <=  (6.012 + 5.*4.008)))
					fill_weight = -0.25;
				else 
					continue;
				
				double ecomp1 = 1. / ((1./eb) + (1./m_e)*(1.-cos(theta1)));
				double ecomp2 = 1. / ((1./eb) + (1./m_e)*(1.-cos(theta2)));
				
				double deltaE = (e1+e2) - (eb+m_e);
				double deltaK = (ecomp1+ecomp2) - (eb+m_e);
				
				double ecomp1_ratio = e1/ecomp1;
				double elas_ratio   = e1/(eb-e2);
				
				h_deltaE->Fill(deltaE, fill_weight);
				h_deltaK->Fill(deltaK, fill_weight);
				
				if(fabs(deltaK) > DELTA_K_CUT) continue;
				
				//----------------------------------------------------------//
				
				int cut_vals[n_cuts];
				for(int icut = 0; icut < n_cuts; icut++) cut_vals[icut] = 0;
				
				// first cut is trivial:
				cut_vals[0] = 1;
				
				// second cut is a fiducial cut on FCAL:
				if(!fcal_fiducial_cut(pos2,locVertex,2.0)) {
					
					cut_vals[1] = 1;
					
					// third cut is requirement that there is only one shower in FCAL and one in CCAL:
					if(locNFCALShowers==1 && locNCCALShowers==1) {
						
						cut_vals[2] = 1;
					}
				}
				
				//-----------------------------------------------------------------//
				// Fill Histograms:
				
				for(int icut = 0; icut < n_cuts; icut++) {
					if(cut_vals[icut]) {
						
						h_CompRatio[icut]->Fill(idmax, ecomp1_ratio, fill_weight);
						if(MIN_CALIB_ENERGY<ecomp1 && ecomp1<MAX_CALIB_ENERGY) {
							h_CompRatio_cut[icut]->Fill(idmax, ecomp1_ratio, fill_weight);
						}
						
						h_ElasRatio[icut]->Fill(idmax, elas_ratio, fill_weight);
						if(MIN_CALIB_ENERGY<((eb+m_e)-e2) && ((eb+m_e)-e2)<MAX_CALIB_ENERGY) {
							h_ElasRatio_cut[icut]->Fill(idmax, elas_ratio, fill_weight);
						}
						
						h_CompRatio_vs_E[icut]->Fill(ecomp1, ecomp1_ratio, fill_weight);
						h_ElasRatio_vs_E[icut]->Fill((eb+m_e)-e2, elas_ratio, fill_weight);
						
						if(loc_tof_match) {
							
							// Photon in CCAL:
							
							h_CompRatio_g[icut]->Fill(idmax, ecomp1_ratio, fill_weight);
							if(MIN_CALIB_ENERGY<ecomp1 && ecomp1<MAX_CALIB_ENERGY) {
								h_CompRatio_g_cut[icut]->Fill(idmax, ecomp1_ratio, fill_weight);
							}
							
							h_ElasRatio_g[icut]->Fill(idmax, elas_ratio, fill_weight);
							if(MIN_CALIB_ENERGY<((eb+m_e)-e2) && ((eb+m_e)-e2)<MAX_CALIB_ENERGY) {
								h_ElasRatio_g_cut[icut]->Fill(idmax, elas_ratio, fill_weight);
							}
							
							h_CompRatio_vs_E_g[icut]->Fill(ecomp1, ecomp1_ratio, fill_weight);
							h_ElasRatio_vs_E_g[icut]->Fill((eb+m_e)-e2, elas_ratio, fill_weight);
							
						} else {
							
							// Electron in CCAL:
							
							h_CompRatio_e[icut]->Fill(idmax, ecomp1_ratio, fill_weight);
							if(MIN_CALIB_ENERGY<ecomp1 && ecomp1<MAX_CALIB_ENERGY) {
								h_CompRatio_e_cut[icut]->Fill(idmax, ecomp1_ratio, fill_weight);
							}
							
							h_ElasRatio_e[icut]->Fill(idmax, elas_ratio, fill_weight);
							if(MIN_CALIB_ENERGY<((eb+m_e)-e2) && ((eb+m_e)-e2)<MAX_CALIB_ENERGY) {
								h_ElasRatio_e_cut[icut]->Fill(idmax, elas_ratio, fill_weight);
							}
							
							h_CompRatio_vs_E_e[icut]->Fill(ecomp1, ecomp1_ratio, fill_weight);
							h_ElasRatio_vs_E_e[icut]->Fill((eb+m_e)-e2, elas_ratio, fill_weight);
						}
						
					}
				}
				
			} // end loop over beam photons
		} // end loop over fcal showers
	} // end loop over ccal showers
	
	japp->RootFillUnLock(this);  // Release root lock
	
	return NOERROR;
}

//------------------
// fcal_fiducial_cut
//------------------

int JEventProcessor_CCAL_ComptonGains::fcal_fiducial_cut(DVector3 pos, DVector3 vertex, 
	double cut_layers) {
	
  int fid_cut = 0;
  
  double fcal_inner_layer_cut = (1.5 + cut_layers) * 4.0157;
  
  double fcal_face_x = vertex.X() + (pos.X() * (m_fcalZ - vertex.Z())/pos.Z());
  double fcal_face_y = vertex.Y() + (pos.Y() * (m_fcalZ - vertex.Z())/pos.Z());
  
  fcal_face_x -= m_fcalX;
  fcal_face_y -= m_fcalY;
  
  if((-1.*fcal_inner_layer_cut < fcal_face_x && fcal_face_x < fcal_inner_layer_cut)
     && (-1.*fcal_inner_layer_cut < fcal_face_y 
	 && fcal_face_y < fcal_inner_layer_cut)) fid_cut = 1;
  
  return fid_cut;
}

//------------------
// check_TOF_match
//------------------

void JEventProcessor_CCAL_ComptonGains::check_TOF_match(DVector3 pos1, double rfTime, 
	DVector3 vertex, vector<const DTOFPoint*> tof_points, double &dx_min, double &dy_min, 
	double &dt_min, double rf_time_cut) {
	
	dx_min = 1000.;
	dy_min = 1000.;
	dt_min = 1000.;
	
	for(vector< const DTOFPoint* >::const_iterator tof = tof_points.begin();
		tof != tof_points.end(); tof++) {
		
		double xt = (*tof)->pos.X() - vertex.X();
		double yt = (*tof)->pos.Y() - vertex.Y();
		double zt = (*tof)->pos.Z() - vertex.Z();		
		double rt = sqrt(xt*xt + yt*yt + zt*zt);
		double tt = (*tof)->t - (rt/c);
		double dt = tt - rfTime;
		xt *= pos1.Z() / zt;
		yt *= pos1.Z() / zt;
		double dx = pos1.X() - xt;
		double dy = pos1.Y() - yt;
		
		if(fabs(dt) < rf_time_cut) {
			if((dx*dx + dy*dy) < (dx_min*dx_min + dy_min*dy_min)) {
				dx_min = dx;
				dy_min = dy;
				dt_min = dt;
			}
		}
	}
	
	return;
}
