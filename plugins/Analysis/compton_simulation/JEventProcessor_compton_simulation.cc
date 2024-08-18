// $Id$
//
//    File: JEventProcessor_compton_simulation.cc
// Created: Thu Feb 10 21:33:24 EST 2022
// Creator: andrsmit (on Linux ifarm1802.jlab.org 3.10.0-1160.11.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_compton_simulation.h"

// Routine used to create our JEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_compton_simulation());
}
} // "C"

JEventProcessor_compton_simulation::JEventProcessor_compton_simulation() {
	
	m_USE_REACTION_WEIGHT = 0;
	gPARMS->SetDefaultParameter("compton_simulation:USE_REACTION_WEIGHT", 
		m_USE_REACTION_WEIGHT);
	
	m_REACTION_CUT_WEIGHT = 1.e4;
	gPARMS->SetDefaultParameter("compton_simulation:REACTION_CUT_WEIGHT", 
		m_REACTION_CUT_WEIGHT);
}

//------------------
// init
//------------------
jerror_t JEventProcessor_compton_simulation::init(void)
{
	rndm_gen = new TRandom3(0);
	
	TDirectory *dir_compton = new TDirectoryFile("compton_simulation", 
		"compton_simulation");
	dir_compton->cd();
	
	h_beam           = new TH1F("beam", "Is there a beam photon?", 2, -0.5, 1.5);
	h_tagh_flux      = new TH1F("tagh_flux", "TAGH Flux", 274, 0.5, 274.5);
	h_tagm_flux      = new TH1F("tagm_flux", "TAGM Flux", 102, 0.5, 102.5);
	h_double_compton = new TH1F("double_compton", "Is e2 > 0?", 2, -0.5, 1.5);
	h_beam_energy    = new TH1F("beam_energy", "Photon Beam Energy", 12000, 0., 12.);
	
	h_vertex           = new TH1F("vertex",          
		"Vertex Z Position (unweighted)", 1000, 0., 100.);
	h_vertex_weight    = new TH1F("vertex_weight", 
		"Event Weight", 1000, 0., 2.);
	h_vertex_accepted  = new TH1F("vertex_accepted", 
		"Vertex Z Position (weighted)",   1000, 0., 100.);
	h_vertex_xy        = new TH2F("vertex_xy", "Vertex Y vs. X; x [cm]; y [cm]", 
		1000, -5., 5., 1000, -5., 5.);
	h_reaction_weight  = new TH1F("reaction_weight", 
		"Event Reaction Weight (All Events)", 1000, 0., 1.e7);
	
	//--------------------------------------------------------------------//
	
	h_fcal_rf_dt  = new TH1F("fcal_rf_dt", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	h_ccal_rf_dt  = new TH1F("ccal_rf_dt", "t_{CCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	h_beam_rf_dt  = new TH1F("beam_rf_dt", "t_{Beam} - t_{RF}; [ns]", 2000, -100., 100.);
	h_beam_rf_dt_cut = new TH1F("beam_rf_dt_cut", "t_{Beam} - t_{RF} (Compton Cuts); [ns]", 
		2000, -100., 100.);
	h_beam_rf_dt_tagh = new TH2F("beam_rf_dt_tagh", "t_{Beam} - t_{RF}; TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -100., 100.);
	h_beam_rf_dt_tagm = new TH2F("beam_rf_dt_tagm", "t_{Beam} - t_{RF}; TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -100., 100.);
	
	//--------------------------------------------------------------------//
	
	h_nccal = new TH1F("nccal", "number of CCAL Showers", 10, -0.5, 9.5);
	h_nfcal = new TH1F("nfcal", "number of FCAL Showers", 10, -0.5, 9.5);
	h_n_ccal_fcal     = new TH1F("n_ccal_fcal",     
		"Is there a shower in both CCAL+FCAL", 2, -0.5, 1.5);
	h_n_ccal_fcal_cut = new TH1F("n_ccal_fcal_cut", 
		"Is there a shower in both CCAL+FCAL", 2, -0.5, 1.5);
	
	h_ccal_nblocks = new TH1F("ccal_nblocks", "Number of hits in CCAL Cluster", 51, -0.5, 50.5);
	
	h_ccalE     = new TH2F("ccalE", "Energy of CCAL Shower; E_{Beam} [GeV]; E_{CCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	h_fcalE     = new TH2F("fcalE", "Energy of FCAL Shower; E_{Beam} [GeV]; E_{FCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	h_ccalE_cut = new TH2F("ccalE_cut", "Energy of CCAL Shower; E_{Beam} [GeV]; E_{CCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	h_fcalE_cut = new TH2F("fcalE_cut", "Energy of FCAL Shower; E_{Beam} [GeV]; E_{FCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	
	//--------------------------------------------------------------------//
	
	h_deltaE_vs_deltaK         = new TH2F("deltaE_vs_deltaK", 
		"#DeltaE vs. #DeltaK; E_{Compton} - E_{#gamma} [GeV]; E_{1} + E_{2} - E_{#gamma} [GeV]", 
		500, -8.0, 8.0, 500, -4.0, 4.0);
	h_deltaE_vs_deltaK_smeared = new TH2F("deltaE_vs_deltaK_smeared", 
		"#DeltaE vs. #DeltaK; E_{Compton} - E_{#gamma} [GeV]; E_{1} + E_{2} - E_{#gamma} [GeV]", 
		500, -8.0, 8.0, 500, -4.0, 4.0);
	
	/*
	h_deltaR_tagh = new TH2F("deltaR_tagh", 
		"Cluster Separation; TAGH Counter; #DeltaR [cm]", 274, 0.5, 274.5, 2000, 0., 100.);
	h_deltaR_tagm = new TH2F("deltaR_tagm", 
		"Cluster Separation; TAGM Counter; #DeltaR [cm]", 274, 0.5, 274.5, 2000, 0., 100.);
	h_deltaR = new TH2F("deltaR", 
		"Cluster Separation; E_{#gamma} [GeV]; #DeltaR [cm]", 120, 0., 12., 2000, 0., 100.);
	*/
	h_deltaE_tagh = new TH2F("deltaE_tagh", 
		"#DeltaE; TAGH Counter; E_{1}+E_{2} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0);
	h_deltaE_tagm = new TH2F("deltaE_tagm", 
		"#DeltaE; TAGM Counter; E_{1}+E_{2} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0);
	h_deltaE_smeared_tagh = new TH2F("deltaE_smeared_tagh", 
		"#DeltaE; TAGH Counter; E_{1}+E_{2} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0);
	h_deltaE_smeared_tagm = new TH2F("deltaE_smeared_tagm", 
		"#DeltaE; TAGM Counter; E_{1}+E_{2} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0);
	
	h_deltaPhi_tagh = new TH2F("deltaPhi_tagh", 
		"#Delta#phi; TAGH Counter; |#phi_{1}-#phi_{2}| [deg.]", 
		274, 0.5, 274.5, 3600, 0.0, 360.0);
	h_deltaPhi_tagm = new TH2F("deltaPhi_tagm", 
		"#Delta#phi; TAGM Counter; |#phi_{1}-#phi_{2}| [deg.]", 
		102, 0.5, 102.5, 3600, 0.0, 360.0);
	h_deltaPhi_smeared_tagh = new TH2F("deltaPhi_smeared_tagh", 
		"#Delta#phi; TAGH Counter; |#phi_{1}-#phi_{2}| [deg.]", 
		274, 0.5, 274.5, 3600, 0.0, 360.0);
	h_deltaPhi_smeared_tagm = new TH2F("deltaPhi_smeared_tagm", 
		"#Delta#phi; TAGM Counter; |#phi_{1}-#phi_{2}| [deg.]", 
		102, 0.5, 102.5, 3600, 0.0, 360.0);
	
	h_deltaK_tagh = new TH2F("deltaK_tagh", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm = new TH2F("deltaK_tagm", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	h_deltaK_shifted_tagh = new TH2F("deltaK_shifted_tagh", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_shifted_tagm = new TH2F("deltaK_shifted_tagm", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	h_deltaK_smeared_tagh = new TH2F("deltaK_smeared_tagh", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_smeared_tagm = new TH2F("deltaK_smeared_tagm", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	
	h_deltaK_tagh_cut = new TH2F("deltaK_tagh_cut", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm_cut = new TH2F("deltaK_tagm_cut", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	h_deltaK_shifted_tagh_cut = new TH2F("deltaK_shifted_tagh_cut", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_shifted_tagm_cut = new TH2F("deltaK_shifted_tagm_cut", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	h_deltaK_smeared_tagh_cut = new TH2F("deltaK_smeared_tagh_cut", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_smeared_tagm_cut = new TH2F("deltaK_smeared_tagm_cut", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	
	h_deltaK_tagh_cut_main = new TH2F("deltaK_tagh_cut_main", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm_cut_main = new TH2F("deltaK_tagm_cut_main", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	
	h_deltaK_tagh_cut_acc = new TH2F("deltaK_tagh_cut_acc", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm_cut_acc = new TH2F("deltaK_tagm_cut_acc", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	
	h_fcal_xy = new TH2F("fcal_xy", 
		"FCAL Y vs. X; x_{FCAL} [cm]; y_{FCAL} [cm]", 500, -100., 100., 500, -100., 100.);
	h_ccal_xy = new TH2F("ccal_xy", 
		"CCAL Y vs. X; x_{CCAL} [cm]; y_{CCAL} [cm]", 500,  -15.,  15., 500,  -15.,  15.);
	
	dir_compton->cd("../");
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_compton_simulation::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	
	DGeometry*   dgeom = NULL;
	DApplication* dapp = dynamic_cast< DApplication* >(eventLoop->GetJApplication());
	if(dapp)     dgeom = dapp->GetDGeometry(runnumber);
	
	if(dgeom){
		dgeom->GetTargetZ(m_beamZ);
		dgeom->GetFCALPosition(m_fcalX, m_fcalY, m_fcalZ);
		dgeom->GetCCALPosition(m_ccalX, m_ccalY, m_ccalZ);
	} else{
		cerr << "No geometry accessbile to PrimExComptonAnalysis plugin." << endl;
		return RESOURCE_UNAVAILABLE;
	}
	
	jana::JCalibration *jcalib = japp->GetJCalibration(runnumber);
	std::map<string, float> beam_spot;
	jcalib->Get("PHOTON_BEAM/beam_spot", beam_spot);
	m_beamX = beam_spot.at("x");
	m_beamY = beam_spot.at("y");
	
	if(runnumber>60000 && runnumber<69999) {
		
		phase_val = 1;
		
		if(runnumber<61355) {
			m_target_length  = 1.7755;
			m_target_density = 1.85;
			m_atten          = 0.01172;
		} else {
			m_target_length  = 29.535;
			m_target_density = 0.1217;
			m_atten          = 0.00821;
		}
		
		// (2/4/2024): Correction to alignment after Simon updated beam spot with new CDC alignment:
		
		m_fcalX_new =  0.455;
		m_fcalY_new = -0.032;
		
		m_ccalX_new = -0.082;
		if(runnumber<61483) m_ccalY_new = 0.061;
		else                m_ccalY_new = 0.051;
		
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
		
		phase_val = 2;
		
		if(runnumber<81384) {
			m_target_length  = 1.7755;
			m_target_density = 1.85;
			m_atten          = 0.01172;
		} else {
			m_target_length  = 29.535;
			m_target_density = 0.1217;
			m_atten          = 0.00821;
		}
		
		m_fcalX_new = 0.408;
		m_fcalY_new = 0.027;
		
		m_ccalX_new = 0.108;
		m_ccalY_new = 0.130;
		
		if(runnumber<81471) {
			m_beamX =  0.129;
			m_beamY = -0.038;
		} else {
			bfield_val = 1;
			m_beamX =  0.139874;
			m_beamY = -0.040895;
		}
	} else {
		
		phase_val = 3;
		
		if(runnumber<110622) {
			m_target_length  = 1.7755;
			m_target_density = 1.85;
			m_atten          = 0.01172;
		} else {
			m_target_length  = 29.535;
			m_target_density = 0.1217;
			m_atten          = 0.00821;
		}
		/*
		m_fcalX_new = 0.408;
		m_fcalY_new = 0.027;
		m_ccalX_new = 0.135;
		m_ccalY_new = 0.135;
		m_beamX     = 0.146;
		m_beamY     = 0.017;
		*/
		m_fcalX_new = 0.408;
		m_fcalY_new = 0.027;
		m_ccalX_new = 0.184;
		m_ccalY_new = 0.110;
		m_beamX     = 0.151;
		m_beamY     = 0.012;
	}
	
	fcal_correction.SetXYZ(m_fcalX_new-m_fcalX, m_fcalY_new-m_fcalY, 0.);
	ccal_correction.SetXYZ(m_ccalX_new-m_ccalX, m_ccalY_new-m_ccalY, 0.);
	
	set_cuts(runnumber);
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_compton_simulation::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
	//-----   Data Objects   -----//
	
	DVector3 vertex;
	vertex.SetXYZ(m_beamX, m_beamY, m_beamZ);
	
	vector<const DBeamPhoton*> beam_photons;
	eventLoop->Get(beam_photons);
	
	vector<const DCCALShower*> ccal_showers;
	eventLoop->Get(ccal_showers);
	
	vector<const DFCALShower*> fcal_showers;
	eventLoop->Get(fcal_showers);
	
	vector<const DMCThrown*> mc_thrown;
	eventLoop->Get(mc_thrown);
	
	vector<const DMCReaction*> mc_reaction;
	eventLoop->Get(mc_reaction);
	
	vector<const DFCALShower*> good_fcal_showers;
	vector<const DCCALShower*> good_ccal_showers;
	
	double thrown_energy   = mc_reaction[0]->beam.energy();
	
	double vertex_z        = mc_thrown[0]->position().Z();
	double vertex_weight   = get_vertex_weight(vertex_z);
	double reaction_weight = mc_reaction[0]->weight;
	
	japp->RootFillLock(this);  // Acquire root lock
	
	h_vertex->Fill(vertex_z);
	h_vertex_weight->Fill(vertex_weight);
	h_reaction_weight->Fill(reaction_weight);
	
	// apply accept-reject filter based on vertex-event weight:
	
	if(vertex_weight < rndm_gen->Uniform()) {
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	if(m_USE_REACTION_WEIGHT && reaction_weight > m_REACTION_CUT_WEIGHT) {
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	
	h_vertex_accepted->Fill(vertex_z);
	h_vertex_xy->Fill(mc_thrown[0]->position().X(), mc_thrown[0]->position().Y());
	
	//-----   RF Bunch   -----//
	
	const DEventRFBunch *locRFBunch = NULL;
	try { 
	  	eventLoop->GetSingle(locRFBunch, "CalorimeterOnly");
	} catch (...) { 
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	double rfTime = locRFBunch->dTime;
	
	
	int n_beam_photons = (int)beam_photons.size();
	if(n_beam_photons) 
		h_beam->Fill(1);
	else 
		h_beam->Fill(0);
	
	for( vector< const DBeamPhoton* >::const_iterator gam = beam_photons.begin();
		gam != beam_photons.end(); gam++ ) {
		
		int counter = (*gam)->dCounter;
		
		DetectorSystem_t sys = (*gam)->dSystem;
		if(sys==SYS_TAGH)      h_tagh_flux->Fill(counter);
		else if(sys==SYS_TAGM) h_tagm_flux->Fill(counter);
		
		h_beam_rf_dt->Fill((*gam)->time() - rfTime);
		h_beam_energy->Fill((*gam)->lorentzMomentum().E());
		if(sys==SYS_TAGH) h_beam_rf_dt_tagh->Fill(counter, (*gam)->time()-rfTime);
		if(sys==SYS_TAGM) h_beam_rf_dt_tagm->Fill(counter, (*gam)->time()-rfTime);
		
	}
	
	int n_fcal_showers = 0, n_fcal_showers_cut = 0;	
	for(vector<const DFCALShower*>::const_iterator show = fcal_showers.begin(); 
		show != fcal_showers.end(); show++) {
		
		DVector3 loc_pos = (*show)->getPosition_log() - vertex + fcal_correction;
		double loc_t     = (*show)->getTime() - (loc_pos.Mag()/c) - rfTime;
		
		int fid_cut = fcal_fiducial_cut(loc_pos, vertex);
		
		if((fabs(loc_t) < FCAL_RF_time_cut) && !fid_cut) {
		  	n_fcal_showers_cut++;
			if((*show)->getEnergy() > FCAL_min_energy_cut) {
				n_fcal_showers++;
		  		good_fcal_showers.push_back((*show));
			}
		}
		
		h_fcal_rf_dt->Fill(loc_t);
	}
	
	int n_ccal_showers = 0, n_ccal_showers_cut = 0;
	for(vector<const DCCALShower*>::const_iterator show = ccal_showers.begin();
		show != ccal_showers.end(); show++) {
		
		DVector3 loc_pos((*show)->x1, (*show)->y1, (*show)->z);
		loc_pos = loc_pos - vertex + ccal_correction;
		double loc_t = (*show)->time - (loc_pos.Mag()/c) - rfTime;
		
		int fid_cut = ccal_fiducial_cut(loc_pos, vertex);
		
		if((fabs(loc_t) < CCAL_RF_time_cut) && !fid_cut) {
			n_ccal_showers_cut++;
			if((*show)->E > CCAL_min_energy_cut) {
				n_ccal_showers++;
				good_ccal_showers.push_back((*show));
			}
		}
		
		h_ccal_rf_dt->Fill(loc_t);
	}
	
	h_nccal->Fill(n_ccal_showers);
	h_nfcal->Fill(n_fcal_showers);
	
	if(ccal_showers.size() && fcal_showers.size())
		h_n_ccal_fcal->Fill(1.);
	else 
		h_n_ccal_fcal->Fill(0.);
	
	if(n_ccal_showers && n_fcal_showers)
		h_n_ccal_fcal_cut->Fill(1.);
	else 
		h_n_ccal_fcal_cut->Fill(0.);
	
	if(locRFBunch->dNumParticleVotes < 2) {
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	
	// plot energy of showers when there are two good showers in the FCAL/CCAL:
	
	if(n_ccal_showers_cut==1 && n_fcal_showers_cut==1) {
		
		for(vector< const DFCALShower*>::const_iterator show = fcal_showers.begin(); 
			show != fcal_showers.end(); show++) {
			
			DVector3 loc_pos = (*show)->getPosition_log() - vertex + fcal_correction;
			double loc_t     = (*show)->getTime() - (loc_pos.Mag()/c) - rfTime;
			
			int fid_cut = fcal_fiducial_cut(loc_pos, vertex);
			
			if((fabs(loc_t) < FCAL_RF_time_cut) && !fid_cut)
				h_fcalE->Fill(thrown_energy, (*show)->getEnergy());
		}
		
		for(vector< const DCCALShower*>::const_iterator show = ccal_showers.begin(); 
			show != ccal_showers.end(); show++) {
			
			DVector3 loc_pos((*show)->x1, (*show)->y1, (*show)->z);
			loc_pos = loc_pos - vertex + ccal_correction;
			double loc_t = (*show)->time - (loc_pos.Mag()/c) - rfTime;
			
			int fid_cut = ccal_fiducial_cut(loc_pos, vertex);
			
			if((fabs(loc_t) < CCAL_RF_time_cut) && !fid_cut) 
				h_ccalE->Fill(thrown_energy, (*show)->E);
		}
	}
	
	
	//----------     Check FCAL-CCAL Pairs     ----------//
	
	vector<ComptonCandidate_t> candidates;
	
	for(vector<const DFCALShower*>::const_iterator show1 = good_fcal_showers.begin(); 
		show1 != good_fcal_showers.end(); show1++) {
		
		double e1     = (*show1)->getEnergy();
		DVector3 pos1 = (*show1)->getPosition_log() - vertex + fcal_correction;
		
		double t1     = (*show1)->getTime() - (pos1.Mag()/c);
		double phi1   = pos1.Phi() * (180./TMath::Pi());
		double theta1 = pos1.Theta();
		
		for(vector<const DCCALShower*>::const_iterator show2 = good_ccal_showers.begin(); 
			show2 != good_ccal_showers.end(); show2++) {
			
			double e2  = (*show2)->E;
			DVector3 pos2( (*show2)->x1, (*show2)->y1, (*show2)->z );
			pos2       = pos2 - vertex + ccal_correction;
			
			double t2     = (*show2)->time - (pos2.Mag()/c);
			double phi2   = pos2.Phi() * (180./TMath::Pi());
			double theta2 = pos2.Theta();
			
			// calculate deltaPhi and deltaT:
			
			double deltaPhi = fabs(phi2 - phi1);
			double deltaT   = t2 - t1;
			
			// loop over beam photons:
			
			for(vector<const DBeamPhoton*>::const_iterator gam = beam_photons.begin();
				gam != beam_photons.end(); gam++) {
				
				double eb = (*gam)->lorentzMomentum().E();
				double tb = (*gam)->time();
				
				double brfdt = tb - rfTime;
				
				int bunch_val;
				
				if(fabs(brfdt) < 2.004)
					bunch_val = 1;
				else if((-(2.004 + 10.*4.008)<=brfdt && brfdt<=-(2.004 + 5.*4.008))
					||((2.004 + 5.*4.008)<=brfdt && brfdt<=(2.004 + 10.*4.008)))
					bunch_val = 0;
				else 
					continue;
				
				if(eb < BEAM_min_energy_cut) continue;
				/*
				double ecomp1 =  1. / ((1./eb) + (1./m_e)*(1.-cos(theta1)));
				double ecomp2 =  1. / ((1./eb) + (1./m_e)*(1.-cos(theta2)));
				double deltaK = (ecomp1 + ecomp2) - (eb + m_e);
				*/
				double deltaE = (e1 + e2) - (eb + m_e);
				double deltaK = m_e * sin(theta1+theta2) 
					/ (sin(theta1) + sin(theta2) - sin(theta1+theta2));
				deltaK       -= eb;
				
				if(fabs(deltaPhi-180.) < 60. && fabs(deltaE)<1.0 && fabs(deltaK)<1.5) {
					h_ccalE_cut->Fill(thrown_energy, e2);
					h_fcalE_cut->Fill(thrown_energy, e1);
				}
				
				ComptonCandidate_t loc_Cand;
				
				loc_Cand.bunch_val    = bunch_val;
				
				loc_Cand.e1           = e1;
				loc_Cand.x1           = pos1.X();
				loc_Cand.y1           = pos1.Y();
				loc_Cand.z1           = pos1.Z();
				loc_Cand.e2           = e2;
				loc_Cand.x2           = pos2.X();
				loc_Cand.y2           = pos2.Y();
				loc_Cand.z2           = pos2.Z();
				
				loc_Cand.deltaPhi     = deltaPhi;
				loc_Cand.deltaT       = deltaT;
				loc_Cand.deltaE       = deltaE;
				loc_Cand.deltaK       = deltaK;
				loc_Cand.ccal_nblocks = (*show2)->dime;
				
				loc_Cand.vz           = vertex_z;
				
				loc_Cand.eb           = eb;
				loc_Cand.brfdt        = brfdt;
				loc_Cand.tag_counter  = (*gam)->dCounter;
				
				DetectorSystem_t sys  = (*gam)->dSystem;
				if(sys==SYS_TAGH)      loc_Cand.tag_sys = 0;
				else if(sys==SYS_TAGM) loc_Cand.tag_sys = 1;
				
				loc_Cand.event_weight = reaction_weight;
				
				candidates.push_back(loc_Cand);
				
			} // end DBeamPhoton loop
		} // end DCCALShower loop
	} // end DFCALShower loop
	
	
	fill_histograms(candidates, vertex);
	
	japp->RootFillUnLock(this);  // Release root lock
	
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_compton_simulation::erun(void)
{
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_compton_simulation::fini(void)
{
	
	return NOERROR;
}


int JEventProcessor_compton_simulation::fcal_fiducial_cut(DVector3 pos, DVector3 vertex)
{
	int fid_cut = 0;
	
	double fcal_inner_layer_cut = 2.5 * 4.0157;
	
	double fcal_face_x = vertex.X() + (pos.X() * (m_fcalZ - vertex.Z())/pos.Z());
	double fcal_face_y = vertex.Y() + (pos.Y() * (m_fcalZ - vertex.Z())/pos.Z());
	
	fcal_face_x -= m_fcalX_new;
	fcal_face_y -= m_fcalY_new;
	
	if((-1.*fcal_inner_layer_cut < fcal_face_x && fcal_face_x < fcal_inner_layer_cut)
		&& (-1.*fcal_inner_layer_cut < fcal_face_y 
		&& fcal_face_y < fcal_inner_layer_cut)) fid_cut = 1;
	
	// only apply the next fiducial cut for runs from phase-I:
	
	if(phase_val < 2) {
	
	if((-32.<fcal_face_y && fcal_face_y<-20.) && (-8.<fcal_face_x && fcal_face_x<4.))
			fid_cut = 1;
	
	}
	
	return fid_cut;
}


int JEventProcessor_compton_simulation::ccal_fiducial_cut(DVector3 pos, DVector3 vertex)
{
	int fid_cut = 0;
	
	double ccal_inner_layer_cut = 2.0 * 2.09;
	
	double ccal_face_x = vertex.X() + (pos.X() * (m_ccalZ - vertex.Z())/pos.Z());
	double ccal_face_y = vertex.Y() + (pos.Y() * (m_ccalZ - vertex.Z())/pos.Z());
	
	ccal_face_x -= m_ccalX_new;
	ccal_face_y -= m_ccalY_new;
	
	if((-1.*ccal_inner_layer_cut < ccal_face_x && ccal_face_x < ccal_inner_layer_cut)
		&& (-1.*ccal_inner_layer_cut < ccal_face_y 
		&& ccal_face_y < ccal_inner_layer_cut)) fid_cut = 1;
	
	if(ccal_face_x<-8.36 || ccal_face_x>10.45 
		|| ccal_face_y<-10.45 || ccal_face_y>10.45) fid_cut = 1;
	
	return fid_cut;
}



void JEventProcessor_compton_simulation::fill_histograms(
	vector<ComptonCandidate_t> Comp_Candidates, DVector3 vertex) 
{
	int n_candidates = static_cast<int>(Comp_Candidates.size());
	
	for(int ic = 0; ic < n_candidates; ic++) {
		
		ComptonCandidate_t loc_Cand = Comp_Candidates[ic];
		
		//-------------------------------------------//
		
		int bunch_val   = loc_Cand.bunch_val;
		double eb       = loc_Cand.eb;
		double brfdt    = loc_Cand.brfdt;
		int tag_sys     = loc_Cand.tag_sys;
		int tag_counter = loc_Cand.tag_counter;
		
		double deltaPhi = loc_Cand.deltaPhi;
		double deltaE   = loc_Cand.deltaE;
		double deltaK   = loc_Cand.deltaK;
		
		double x1 = loc_Cand.x1;
		double y1 = loc_Cand.y1;
		double x2 = loc_Cand.x2;
		double y2 = loc_Cand.y2;
		
		int ccal_nblocks = loc_Cand.ccal_nblocks;
		
		double event_weight = loc_Cand.event_weight;
		
		//--------------     Cuts      --------------//
		
		double deltaE_smeared, deltaPhi_smeared;
		double deltaK_smeared, deltaK_shifted;
		
		//--------------------------------------------------------//
		// DeltaE Cuts:
		
		double deltaE_mu_mc    = f_deltaE_mu_mc->Eval(eb);
		double deltaE_mu_data  = f_deltaE_mu_data->Eval(eb);
		double deltaE_sig_mc   = eb * f_deltaE_sig_mc->Eval(eb);
		double deltaE_sig_data = eb * f_deltaE_sig_data->Eval(eb);
		
		deltaE_smeared = deltaE + (deltaE_mu_data - deltaE_mu_mc);
		if(deltaE_sig_data > deltaE_sig_mc) {
			double loc_deltaE_smear = sqrt(pow(deltaE_sig_data,2.0)
				- pow(deltaE_sig_mc,2.0));
			deltaE_smeared += rndm_gen->Gaus(0.,loc_deltaE_smear);
		}
		
		//--------------------------------------------------------//
		// DeltaPhi Cuts:
		
		double deltaPhi_mu_mc    = f_deltaPhi_mu_mc->Eval(eb);
		double deltaPhi_mu_data  = f_deltaPhi_mu_data->Eval(eb);
		double deltaPhi_sig_mc   = f_deltaPhi_sig_mc->Eval(eb);
		double deltaPhi_sig_data = f_deltaPhi_sig_data->Eval(eb);
		
		deltaPhi_smeared = deltaPhi + (deltaPhi_mu_data - deltaPhi_mu_mc);
		if(deltaPhi_sig_data > deltaPhi_sig_mc) {
			double loc_deltaPhi_smear = sqrt(pow(deltaPhi_sig_data,2.0)
				- pow(deltaPhi_sig_mc,2.0));
			deltaPhi_smeared += rndm_gen->Gaus(0.,loc_deltaPhi_smear);
		}
		
		//--------------------------------------------------------//
		// DeltaK Cuts:
		
		double deltaK_mu_mc    = f_deltaK_mu_mc->Eval(eb);
		double deltaK_mu_data  = f_deltaK_mu_data->Eval(eb);
		double deltaK_sig_mc   = f_deltaK_sig_mc->Eval(eb);
		double deltaK_sig_data = f_deltaK_sig_data->Eval(eb);
		
		deltaK_shifted = deltaK + (deltaK_mu_data - deltaK_mu_mc);
		deltaK_smeared = deltaK_shifted;
		if(deltaK_sig_data > deltaK_sig_mc) {
			double loc_deltaK_smear = sqrt(pow(deltaK_sig_data,2.0)
				- pow(deltaK_sig_mc,2.0));
			deltaK_smeared += rndm_gen->Gaus(0.,loc_deltaK_smear);
		}
		
		int e_cut = 0;
		if(fabs(deltaE_smeared - deltaE_mu_data) < 5.0*deltaE_sig_data) e_cut = 1;
		
		int p_cut = 0;
		if(bfield_val) {
			if(fabs(deltaPhi - 180.) < 50.) {
				p_cut = 1;
			}
		} else {
			if(fabs(deltaPhi_smeared - deltaPhi_mu_data) < 5.0*deltaPhi_sig_data) p_cut = 1;
		}
		
		int k_cut = 0;
		if(fabs(deltaK - deltaK_mu_data) < 5.0*deltaK_sig_data) k_cut = 1;
		
		int k_cut_smeared = 0;
		if(fabs(deltaK_smeared - deltaK_mu_data) < 5.0*deltaK_sig_data) 
			k_cut_smeared = 1;
		
		int k_cut_shifted = 0;
		if(fabs(deltaK_shifted - deltaK_mu_data) < 5.0*deltaK_sig_data) 
			k_cut_shifted = 1;
		
		//-------------------------------------------//
		
		double fill_weight;
		if(bunch_val) fill_weight =  1.0;
		else          fill_weight = -0.1;
		
		if(m_USE_REACTION_WEIGHT) {
			fill_weight *= event_weight;
		}
		
		h_deltaE_vs_deltaK->Fill(deltaK, deltaE, fill_weight);
		h_deltaE_vs_deltaK_smeared->Fill(deltaK_smeared, deltaE_smeared, fill_weight);
		
		if(tag_sys==0) {
			h_deltaE_tagh->Fill(tag_counter, deltaE, fill_weight);
			h_deltaE_smeared_tagh->Fill(tag_counter, deltaE_smeared, fill_weight);
			
			if(e_cut) {
				h_deltaPhi_tagh->Fill(tag_counter, deltaPhi, fill_weight);
				h_deltaPhi_smeared_tagh->Fill(tag_counter, deltaPhi_smeared, fill_weight);
				
				if(p_cut) {
					h_deltaK_tagh->Fill(tag_counter, deltaK, fill_weight);
					h_deltaK_shifted_tagh->Fill(tag_counter, deltaK_shifted, fill_weight);
					h_deltaK_smeared_tagh->Fill(tag_counter, deltaK_smeared, fill_weight);
					if(k_cut) {
						h_deltaK_tagh_cut->Fill(tag_counter, deltaK, fill_weight);
						if(bunch_val) { 
							h_deltaK_tagh_cut_main->Fill(tag_counter, deltaK);
						} else {
							h_deltaK_tagh_cut_acc->Fill(tag_counter, deltaK, 0.5);
						}
						h_beam_rf_dt_cut->Fill(brfdt);
					}
					if(k_cut_shifted) {
						h_deltaK_shifted_tagh_cut->Fill(tag_counter, deltaK_shifted, fill_weight);
					}
					if(k_cut_smeared) {
						h_deltaK_smeared_tagh_cut->Fill(tag_counter, deltaK_smeared, fill_weight);
					}
				}
			}
			
		} else {
			h_deltaE_tagm->Fill(tag_counter, deltaE, fill_weight);
			h_deltaE_smeared_tagm->Fill(tag_counter, deltaE_smeared, fill_weight);
			
			if(e_cut) {
				h_deltaPhi_tagm->Fill(tag_counter, deltaPhi, fill_weight);
				h_deltaPhi_smeared_tagm->Fill(tag_counter, deltaPhi_smeared, 
					fill_weight);
				
				if(p_cut) {
					h_deltaK_tagm->Fill(tag_counter, deltaK, fill_weight);
					h_deltaK_shifted_tagm->Fill(tag_counter, deltaK_shifted, fill_weight);
					h_deltaK_smeared_tagm->Fill(tag_counter, deltaK_smeared, fill_weight);
					if(k_cut) {
						h_deltaK_tagm_cut->Fill(tag_counter, deltaK, fill_weight);
						if(bunch_val) { 
							h_deltaK_tagm_cut_main->Fill(tag_counter, deltaK);
						} else {
							h_deltaK_tagm_cut_acc->Fill(tag_counter, deltaK, 0.5);
						}
						h_beam_rf_dt_cut->Fill(brfdt);
					} 
					if(k_cut_shifted) {
						h_deltaK_shifted_tagm_cut->Fill(tag_counter, deltaK_shifted, fill_weight);
					}
					if(k_cut_smeared) {
						h_deltaK_smeared_tagm_cut->Fill(tag_counter, deltaK_smeared, fill_weight);
					}
				}
			}
		}
		
		if(e_cut && p_cut && k_cut) {
			h_fcal_xy->Fill(x1, y1);
			h_ccal_xy->Fill(x2, y2);
			
			h_ccal_nblocks->Fill(ccal_nblocks);
		}
	}
	
	return;
}



void JEventProcessor_compton_simulation::set_cuts(int32_t runnumber)
{
	if(runnumber>60000 && runnumber<61355) {
	  
		// Phase I, Be Target
		
		deltaE_mu_p0_mc      = -1.69536e-02;  deltaE_mu_p0_data    =  8.33517e-03;
		deltaE_mu_p1_mc      = -1.51864e-02;  deltaE_mu_p1_data    =  2.09025e-03;
		deltaE_mu_p2_mc      =  1.18343e-03;  deltaE_mu_p2_data    = -1.09342e-04;
		deltaE_mu_p3_mc      =  0.;           deltaE_mu_p3_data    =  0.;
		
		deltaE_sig_p0_mc     =  8.49000e-03;  deltaE_sig_p0_data   =  8.37004e-03;
		deltaE_sig_p1_mc     =  4.06306e-02;  deltaE_sig_p1_data   =  4.56259e-02;
		deltaE_sig_p2_mc     =  0.;           deltaE_sig_p2_data   =  0.;
		//----------------------------------------------------------------------//
		deltaPhi_mu_p0_mc    =  1.80231e+02;  deltaPhi_mu_p0_data  =  1.79943e+02;
		deltaPhi_mu_p1_mc    = -1.53158e-01;  deltaPhi_mu_p1_data  = -2.11766e-02;
		deltaPhi_mu_p2_mc    =  1.47909e-02;  deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_mc    = -5.69442e-04;  deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_mc   =  8.94978e+00;  deltaPhi_sig_p0_data =  1.20139e+01;
		deltaPhi_sig_p1_mc   = -8.13230e-01;  deltaPhi_sig_p1_data = -1.75486e+00;
		deltaPhi_sig_p2_mc   =  7.04521e-02;  deltaPhi_sig_p2_data =  1.67515e-01;
		deltaPhi_sig_p3_mc   = -2.13704e-03;  deltaPhi_sig_p3_data = -5.48316e-03;
		//----------------------------------------------------------------------//
		deltaK_mu_p0_mc      =  2.55352e-01;  deltaK_mu_p0_data    = -9.36095e-02;
		deltaK_mu_p1_mc      = -7.43055e-02;  deltaK_mu_p1_data    =  5.48923e-02;
		deltaK_mu_p2_mc      =  3.58192e-03;  deltaK_mu_p2_data    = -1.19844e-02;
		deltaK_mu_p3_mc      = -1.44681e-04;  deltaK_mu_p3_data    =  4.38188e-04;
		
		deltaK_sig_p0_mc     =  8.56721e-01;  deltaK_sig_p0_data   =  6.68283e-01;
		deltaK_sig_p1_mc     = -1.61018e-01;  deltaK_sig_p1_data   = -8.45642e-02;
		deltaK_sig_p2_mc     =  2.48317e-02;  deltaK_sig_p2_data   =  1.61255e-02;
		deltaK_sig_p3_mc     = -9.24893e-04;  deltaK_sig_p3_data   = -5.93363e-04;
		
	} else if(runnumber>60000 && runnumber<69999) {
		
		// Phase I, He Target
		
		deltaE_mu_p0_mc      = -4.70073e-02;  deltaE_mu_p0_data    =  4.07223e-02;
		deltaE_mu_p1_mc      = -9.47891e-03;  deltaE_mu_p1_data    = -1.78574e-02;
		deltaE_mu_p2_mc      =  9.16184e-04;  deltaE_mu_p2_data    =  1.71081e-03;
		deltaE_mu_p3_mc      =  0.;           deltaE_mu_p3_data    = -5.38583e-05;
		
		deltaE_sig_p0_mc     =  1.03019e-02;  deltaE_sig_p0_data   =  1.08608e-02;
		deltaE_sig_p1_mc     =  3.88701e-02;  deltaE_sig_p1_data   =  4.32721e-02;
		deltaE_sig_p2_mc     =  4.92845e-09;  deltaE_sig_p2_data   =  4.73705e-08;
		//----------------------------------------------------------------------//
		deltaPhi_mu_p0_mc    =  1.79371e+02;  deltaPhi_mu_p0_data  =  1.79797e+02;
		deltaPhi_mu_p1_mc    =  1.58594e-01;  deltaPhi_mu_p1_data  = -9.76515e-03;
		deltaPhi_mu_p2_mc    = -2.13751e-02;  deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_mc    =  8.09262e-04;  deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_mc   =  1.12974e+01;  deltaPhi_sig_p0_data =  1.21151e+01;
		deltaPhi_sig_p1_mc   = -1.68146e+00;  deltaPhi_sig_p1_data = -1.84430e+00;
		deltaPhi_sig_p2_mc   =  1.68430e-01;  deltaPhi_sig_p2_data =  1.75617e-01;
		deltaPhi_sig_p3_mc   = -5.74504e-03;  deltaPhi_sig_p3_data = -5.68551e-03;
		//----------------------------------------------------------------------//
		deltaK_mu_p0_mc      = -4.86844e-02;  deltaK_mu_p0_data    = -5.93615e-01;
		deltaK_mu_p1_mc      =  4.65378e-02;  deltaK_mu_p1_data    =  2.57995e-01;
		deltaK_mu_p2_mc      = -1.09603e-02;  deltaK_mu_p2_data    = -3.70844e-02;
		deltaK_mu_p3_mc      =  4.23499e-04;  deltaK_mu_p3_data    =  1.44465e-03;
		
		deltaK_sig_p0_mc     =  3.18796e-01;  deltaK_sig_p0_data   =  3.80560e-01;
		deltaK_sig_p1_mc     =  2.79811e-02;  deltaK_sig_p1_data   =  2.15649e-02;
		deltaK_sig_p2_mc     =  3.01475e-03;  deltaK_sig_p2_data   =  3.03526e-03;
		deltaK_sig_p3_mc     = -7.64643e-05;  deltaK_sig_p3_data   = -4.06248e-05;
		
	} else if(runnumber>80000 && runnumber<81400) {
		
		// Phase II, Be Target
		
		deltaE_mu_p0_mc      = -3.74142e-02;  deltaE_mu_p0_data    = -0.05;
		deltaE_mu_p1_mc      = -1.23213e-02;  deltaE_mu_p1_data    =  0.0;
		deltaE_mu_p2_mc      =  1.17682e-03;  deltaE_mu_p2_data    =  0.0;
		deltaE_mu_p3_mc      =  0.;           deltaE_mu_p3_data    =  0.0;
		
		deltaE_sig_p0_mc     =  8.19775e-03;  deltaE_sig_p0_data   =  1.40025e-02;
		deltaE_sig_p1_mc     =  3.90388e-02;  deltaE_sig_p1_data   =  3.55640e-02;
		deltaE_sig_p2_mc     =  0.;           deltaE_sig_p2_data   =  5.15193e-04;
		//----------------------------------------------------------------------//
		deltaPhi_mu_p0_mc    =  1.79400e+02;  deltaPhi_mu_p0_data  =  1.79816e+02;
		deltaPhi_mu_p1_mc    =  1.62370e-01;  deltaPhi_mu_p1_data  = -7.06279e-03;
		deltaPhi_mu_p2_mc    = -2.35714e-02;  deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_mc    =  9.77592e-04;  deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_mc   =  1.30598e+01;  deltaPhi_sig_p0_data =  1.21214e+01;
		deltaPhi_sig_p1_mc   = -2.19866e+00;  deltaPhi_sig_p1_data = -1.84405e+00;
		deltaPhi_sig_p2_mc   =  2.24090e-01;  deltaPhi_sig_p2_data =  1.81923e-01;
		deltaPhi_sig_p3_mc   = -7.85765e-03;  deltaPhi_sig_p3_data = -6.13418e-03;
		//----------------------------------------------------------------------//
		deltaK_mu_p0_mc      = -1.74619e-01;  deltaK_mu_p0_data    = -1.74619e-01;
		deltaK_mu_p1_mc      =  1.01328e-01;  deltaK_mu_p1_data    =  1.01328e-01;
		deltaK_mu_p2_mc      = -1.84593e-02;  deltaK_mu_p2_data    = -1.84593e-02;
		deltaK_mu_p3_mc      =  7.13390e-04;  deltaK_mu_p3_data    =  7.13390e-04;
		
		deltaK_sig_p0_mc     =  3.88777e-01;  deltaK_sig_p0_data   =  3.88777e-01;
		deltaK_sig_p1_mc     =  2.50588e-02;  deltaK_sig_p1_data   =  2.50588e-02;
		deltaK_sig_p2_mc     =  2.20573e-03;  deltaK_sig_p2_data   =  2.20573e-03;
		deltaK_sig_p3_mc     = -2.76677e-05;  deltaK_sig_p3_data   = -2.76677e-05;
		
	} else if(runnumber>80000 && runnumber<81473) {
		
		// Phase II, He Target, Field OFF (stand-in values, adjust later)
		
		deltaE_mu_p0_mc      = -3.74142e-02;  deltaE_mu_p0_data    = -0.05;
		deltaE_mu_p1_mc      = -1.23213e-02;  deltaE_mu_p1_data    =  0.0;
		deltaE_mu_p2_mc      =  1.17682e-03;  deltaE_mu_p2_data    =  0.0;
		deltaE_mu_p3_mc      =  0.;           deltaE_mu_p3_data    =  0.0;
		
		deltaE_sig_p0_mc     =  8.19775e-03;  deltaE_sig_p0_data   =  1.40025e-02;
		deltaE_sig_p1_mc     =  3.90388e-02;  deltaE_sig_p1_data   =  3.55640e-02;
		deltaE_sig_p2_mc     =  0.;           deltaE_sig_p2_data   =  5.15193e-04;
		//----------------------------------------------------------------------//
		deltaPhi_mu_p0_mc    =  1.79400e+02;  deltaPhi_mu_p0_data  =  1.79816e+02;
		deltaPhi_mu_p1_mc    =  1.62370e-01;  deltaPhi_mu_p1_data  = -7.06279e-03;
		deltaPhi_mu_p2_mc    = -2.35714e-02;  deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_mc    =  9.77592e-04;  deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_mc   =  1.30598e+01;  deltaPhi_sig_p0_data =  1.21214e+01;
		deltaPhi_sig_p1_mc   = -2.19866e+00;  deltaPhi_sig_p1_data = -1.84405e+00;
		deltaPhi_sig_p2_mc   =  2.24090e-01;  deltaPhi_sig_p2_data =  1.81923e-01;
		deltaPhi_sig_p3_mc   = -7.85765e-03;  deltaPhi_sig_p3_data = -6.13418e-03;
		//----------------------------------------------------------------------//
		deltaK_mu_p0_mc      = -1.74619e-01;  deltaK_mu_p0_data    = -1.74619e-01;
		deltaK_mu_p1_mc      =  1.01328e-01;  deltaK_mu_p1_data    =  1.01328e-01;
		deltaK_mu_p2_mc      = -1.84593e-02;  deltaK_mu_p2_data    = -1.84593e-02;
		deltaK_mu_p3_mc      =  7.13390e-04;  deltaK_mu_p3_data    =  7.13390e-04;
		
		deltaK_sig_p0_mc     =  3.88777e-01;  deltaK_sig_p0_data   =  3.88777e-01;
		deltaK_sig_p1_mc     =  2.50588e-02;  deltaK_sig_p1_data   =  2.50588e-02;
		deltaK_sig_p2_mc     =  2.20573e-03;  deltaK_sig_p2_data   =  2.20573e-03;
		deltaK_sig_p3_mc     = -2.76677e-05;  deltaK_sig_p3_data   = -2.76677e-05;
		
	} else if(runnumber>80000 && runnumber<89999) {
		
		// Phase II, He Target, Field ON (stand-in values, adjust later)
		
		deltaE_mu_p0_mc      = -3.74142e-02;  deltaE_mu_p0_data    = -0.05;
		deltaE_mu_p1_mc      = -1.23213e-02;  deltaE_mu_p1_data    =  0.0;
		deltaE_mu_p2_mc      =  1.17682e-03;  deltaE_mu_p2_data    =  0.0;
		deltaE_mu_p3_mc      =  0.;           deltaE_mu_p3_data    =  0.0;
		
		deltaE_sig_p0_mc     =  8.19775e-03;  deltaE_sig_p0_data   =  1.40025e-02;
		deltaE_sig_p1_mc     =  3.90388e-02;  deltaE_sig_p1_data   =  3.55640e-02;
		deltaE_sig_p2_mc     =  0.;           deltaE_sig_p2_data   =  5.15193e-04;
		//----------------------------------------------------------------------//
		deltaPhi_mu_p0_mc    =  1.79400e+02;  deltaPhi_mu_p0_data  =  1.79816e+02;
		deltaPhi_mu_p1_mc    =  1.62370e-01;  deltaPhi_mu_p1_data  = -7.06279e-03;
		deltaPhi_mu_p2_mc    = -2.35714e-02;  deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_mc    =  9.77592e-04;  deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_mc   =  1.30598e+01;  deltaPhi_sig_p0_data =  1.21214e+01;
		deltaPhi_sig_p1_mc   = -2.19866e+00;  deltaPhi_sig_p1_data = -1.84405e+00;
		deltaPhi_sig_p2_mc   =  2.24090e-01;  deltaPhi_sig_p2_data =  1.81923e-01;
		deltaPhi_sig_p3_mc   = -7.85765e-03;  deltaPhi_sig_p3_data = -6.13418e-03;
		//----------------------------------------------------------------------//
		deltaK_mu_p0_mc      = -1.74619e-01;  deltaK_mu_p0_data    = -1.74619e-01;
		deltaK_mu_p1_mc      =  1.01328e-01;  deltaK_mu_p1_data    =  1.01328e-01;
		deltaK_mu_p2_mc      = -1.84593e-02;  deltaK_mu_p2_data    = -1.84593e-02;
		deltaK_mu_p3_mc      =  7.13390e-04;  deltaK_mu_p3_data    =  7.13390e-04;
		
		deltaK_sig_p0_mc     =  3.88777e-01;  deltaK_sig_p0_data   =  3.88777e-01;
		deltaK_sig_p1_mc     =  2.50588e-02;  deltaK_sig_p1_data   =  2.50588e-02;
		deltaK_sig_p2_mc     =  2.20573e-03;  deltaK_sig_p2_data   =  2.20573e-03;
		deltaK_sig_p3_mc     = -2.76677e-05;  deltaK_sig_p3_data   = -2.76677e-05;
		
	} else if(runnumber>110000 && runnumber<110584) {
		
		// Phase III, Be Target, Field ON (stand-in values, adjust later)
		
		deltaE_mu_p0_mc      = -4.64132e-02;  deltaE_mu_p0_data    = -1.15974e-01;
		deltaE_mu_p1_mc      = -7.27578e-03;  deltaE_mu_p1_data    =  3.09559e-02;
		deltaE_mu_p2_mc      =  8.61276e-04;  deltaE_mu_p2_data    = -2.18315e-03;
		deltaE_mu_p3_mc      =  0.; 			    deltaE_mu_p3_data    =  0.;
		
		deltaE_sig_p0_mc     =  8.54695e-03;  deltaE_sig_p0_data   =  1.59557e-02;
		deltaE_sig_p1_mc     =  4.02650e-02;  deltaE_sig_p1_data   =  3.36354e-02;
		deltaE_sig_p2_mc     =  0.; 			    deltaE_sig_p2_data   =  0.;
		//----------------------------------------------------------------------//
		deltaPhi_mu_p0_mc    =  1.80029e+02;  deltaPhi_mu_p0_data  =  1.80032e+02;
		deltaPhi_mu_p1_mc    = -1.02623e-01;  deltaPhi_mu_p1_data  = -9.61299e-03;
		deltaPhi_mu_p2_mc    =  1.15760e-02;  deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_mc    = -5.22398e-04;  deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_mc   =  1.19447e+01;  deltaPhi_sig_p0_data =  1.26943e+01;
		deltaPhi_sig_p1_mc   = -1.75646e+00;  deltaPhi_sig_p1_data = -1.95900e+00;
		deltaPhi_sig_p2_mc   =  1.67672e-01;  deltaPhi_sig_p2_data =  1.87668e-01;
		deltaPhi_sig_p3_mc   = -5.50310e-03;  deltaPhi_sig_p3_data = -6.08281e-03;
		//----------------------------------------------------------------------//
		deltaK_mu_p0_mc      =  6.80186e-02;  deltaK_mu_p0_data    =  3.52591e-01;
		deltaK_mu_p1_mc      = -9.72781e-03;  deltaK_mu_p1_data    = -9.12544e-02;
		deltaK_mu_p2_mc      = -3.79565e-03;  deltaK_mu_p2_data    =  5.17767e-03;
		deltaK_mu_p3_mc      =  1.31657e-04;  deltaK_mu_p3_data    = -2.20872e-04;
		
		deltaK_sig_p0_mc     =  9.82692e-01;  deltaK_sig_p0_data   =  5.97017e-01;
		deltaK_sig_p1_mc     = -1.92373e-01;  deltaK_sig_p1_data   = -5.77323e-02;
		deltaK_sig_p2_mc     =  2.69333e-02;  deltaK_sig_p2_data   =  1.29216e-02;
		deltaK_sig_p3_mc     = -9.51459e-04;  deltaK_sig_p3_data   = -4.54356e-04;
		
	} else if(runnumber>110000 && runnumber<110622) {
		
		// Phase III, Be Target, Field OFF
		
		deltaE_mu_p0_mc      = -2.47457e-02;  deltaE_mu_p0_data    = -1.16581e-01;
		deltaE_mu_p1_mc      = -1.32566e-02;  deltaE_mu_p1_data    =  3.40885e-02;
		deltaE_mu_p2_mc      =  1.20307e-03;  deltaE_mu_p2_data    = -2.77050e-03;
		deltaE_mu_p3_mc      =  0.; 			    deltaE_mu_p3_data    =  2.09654e-05;
		
		deltaE_sig_p0_mc     =  9.09206e-03;  deltaE_sig_p0_data   =  1.65687e-02;
		deltaE_sig_p1_mc     =  3.99628e-02;  deltaE_sig_p1_data   =  2.65247e-02;
		deltaE_sig_p2_mc     =  2.94773e-08;  deltaE_sig_p2_data   =  4.38584e-02;
		//----------------------------------------------------------------------//
		deltaPhi_mu_p0_mc    =  1.79914e+02;  deltaPhi_mu_p0_data  =  1.79825e+02;
		deltaPhi_mu_p1_mc    = -5.75977e-02;  deltaPhi_mu_p1_data  = -1.10610e-02;
		deltaPhi_mu_p2_mc    =  6.32321e-03;  deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_mc    = -3.34789e-04;  deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_mc   =  1.11697e+01;  deltaPhi_sig_p0_data =  1.33099e+01;
		deltaPhi_sig_p1_mc   = -1.47298e+00;  deltaPhi_sig_p1_data = -2.14366e+00;
		deltaPhi_sig_p2_mc   =  1.33343e-01;  deltaPhi_sig_p2_data =  2.06214e-01;
		deltaPhi_sig_p3_mc   = -4.12731e-03;  deltaPhi_sig_p3_data = -6.71967e-03;
		//----------------------------------------------------------------------//
		deltaK_mu_p0_mc      =  2.90606e-01;  deltaK_mu_p0_data    =  2.09726e-01;
		deltaK_mu_p1_mc      = -8.47192e-02;  deltaK_mu_p1_data    = -4.71020e-02;
		deltaK_mu_p2_mc      =  4.36136e-03;  deltaK_mu_p2_data    =  3.81198e-04;
		deltaK_mu_p3_mc      = -1.58040e-04;  deltaK_mu_p3_data    = -4.82342e-05;
		
		deltaK_sig_p0_mc     =  8.05682e-01;  deltaK_sig_p0_data   =  5.19636e-01;
		deltaK_sig_p1_mc     = -1.37154e-01;  deltaK_sig_p1_data   = -3.35925e-02;
		deltaK_sig_p2_mc     =  2.14922e-02;  deltaK_sig_p2_data   =  1.03144e-02;
		deltaK_sig_p3_mc     = -7.82710e-04;  deltaK_sig_p3_data   = -3.63616e-04;
		
	} else if(runnumber>110000 && runnumber<111969) {
		
		// Phase III, He Target, Field ON (stand-in values, adjust later)
		
		deltaE_mu_p0_mc      = -4.64132e-02;  deltaE_mu_p0_data    = -1.15974e-01;
		deltaE_mu_p1_mc      = -7.27578e-03;  deltaE_mu_p1_data    =  3.09559e-02;
		deltaE_mu_p2_mc      =  8.61276e-04;  deltaE_mu_p2_data    = -2.18315e-03;
		deltaE_mu_p3_mc      =  0.; 			    deltaE_mu_p3_data    =  0.;
		
		deltaE_sig_p0_mc     =  8.54695e-03;  deltaE_sig_p0_data   =  1.59557e-02;
		deltaE_sig_p1_mc     =  4.02650e-02;  deltaE_sig_p1_data   =  3.36354e-02;
		deltaE_sig_p2_mc     =  0.; 			    deltaE_sig_p2_data   =  0.;
		//----------------------------------------------------------------------//
		deltaPhi_mu_p0_mc    =  1.80029e+02;  deltaPhi_mu_p0_data  =  1.80032e+02;
		deltaPhi_mu_p1_mc    = -1.02623e-01;  deltaPhi_mu_p1_data  = -9.61299e-03;
		deltaPhi_mu_p2_mc    =  1.15760e-02;  deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_mc    = -5.22398e-04;  deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_mc   =  1.19447e+01;  deltaPhi_sig_p0_data =  1.26943e+01;
		deltaPhi_sig_p1_mc   = -1.75646e+00;  deltaPhi_sig_p1_data = -1.95900e+00;
		deltaPhi_sig_p2_mc   =  1.67672e-01;  deltaPhi_sig_p2_data =  1.87668e-01;
		deltaPhi_sig_p3_mc   = -5.50310e-03;  deltaPhi_sig_p3_data = -6.08281e-03;
		//----------------------------------------------------------------------//
		deltaK_mu_p0_mc      =  6.80186e-02;  deltaK_mu_p0_data    =  3.52591e-01;
		deltaK_mu_p1_mc      = -9.72781e-03;  deltaK_mu_p1_data    = -9.12544e-02;
		deltaK_mu_p2_mc      = -3.79565e-03;  deltaK_mu_p2_data    =  5.17767e-03;
		deltaK_mu_p3_mc      =  1.31657e-04;  deltaK_mu_p3_data    = -2.20872e-04;
		
		deltaK_sig_p0_mc     =  9.82692e-01;  deltaK_sig_p0_data   =  5.97017e-01;
		deltaK_sig_p1_mc     = -1.92373e-01;  deltaK_sig_p1_data   = -5.77323e-02;
		deltaK_sig_p2_mc     =  2.69333e-02;  deltaK_sig_p2_data   =  1.29216e-02;
		deltaK_sig_p3_mc     = -9.51459e-04;  deltaK_sig_p3_data   = -4.54356e-04;
		
	} else if(runnumber>110000 && runnumber<119999) {
		
		// Phase III, He Target, Field OFF (stand-in values, adjust later)
		
		deltaE_mu_p0_mc      = -4.64132e-02;  deltaE_mu_p0_data    = -1.15974e-01;
		deltaE_mu_p1_mc      = -7.27578e-03;  deltaE_mu_p1_data    =  3.09559e-02;
		deltaE_mu_p2_mc      =  8.61276e-04;  deltaE_mu_p2_data    = -2.18315e-03;
		deltaE_mu_p3_mc      =  0.; 			    deltaE_mu_p3_data    =  0.;
		
		deltaE_sig_p0_mc     =  8.54695e-03;  deltaE_sig_p0_data   =  1.59557e-02;
		deltaE_sig_p1_mc     =  4.02650e-02;  deltaE_sig_p1_data   =  3.36354e-02;
		deltaE_sig_p2_mc     =  0.; 			    deltaE_sig_p2_data   =  0.;
		//----------------------------------------------------------------------//
		deltaPhi_mu_p0_mc    =  1.80029e+02;  deltaPhi_mu_p0_data  =  1.80032e+02;
		deltaPhi_mu_p1_mc    = -1.02623e-01;  deltaPhi_mu_p1_data  = -9.61299e-03;
		deltaPhi_mu_p2_mc    =  1.15760e-02;  deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_mc    = -5.22398e-04;  deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_mc   =  1.19447e+01;  deltaPhi_sig_p0_data =  1.26943e+01;
		deltaPhi_sig_p1_mc   = -1.75646e+00;  deltaPhi_sig_p1_data = -1.95900e+00;
		deltaPhi_sig_p2_mc   =  1.67672e-01;  deltaPhi_sig_p2_data =  1.87668e-01;
		deltaPhi_sig_p3_mc   = -5.50310e-03;  deltaPhi_sig_p3_data = -6.08281e-03;
		//----------------------------------------------------------------------//
		deltaK_mu_p0_mc      =  6.80186e-02;  deltaK_mu_p0_data    =  3.52591e-01;
		deltaK_mu_p1_mc      = -9.72781e-03;  deltaK_mu_p1_data    = -9.12544e-02;
		deltaK_mu_p2_mc      = -3.79565e-03;  deltaK_mu_p2_data    =  5.17767e-03;
		deltaK_mu_p3_mc      =  1.31657e-04;  deltaK_mu_p3_data    = -2.20872e-04;
		
		deltaK_sig_p0_mc     =  9.82692e-01;  deltaK_sig_p0_data   =  5.97017e-01;
		deltaK_sig_p1_mc     = -1.92373e-01;  deltaK_sig_p1_data   = -5.77323e-02;
		deltaK_sig_p2_mc     =  2.69333e-02;  deltaK_sig_p2_data   =  1.29216e-02;
		deltaK_sig_p3_mc     = -9.51459e-04;  deltaK_sig_p3_data   = -4.54356e-04;
		
	} else {
		
		// Placeholder for non-PrimEx runs:
		
		deltaE_mu_p0_mc      = -4.64132e-02;  deltaE_mu_p0_data    = -1.15974e-01;
		deltaE_mu_p1_mc      = -7.27578e-03;  deltaE_mu_p1_data    =  3.09559e-02;
		deltaE_mu_p2_mc      =  8.61276e-04;  deltaE_mu_p2_data    = -2.18315e-03;
		deltaE_mu_p3_mc      =  0.; 			    deltaE_mu_p3_data    =  0.;
		
		deltaE_sig_p0_mc     =  8.54695e-03;  deltaE_sig_p0_data   =  1.59557e-02;
		deltaE_sig_p1_mc     =  4.02650e-02;  deltaE_sig_p1_data   =  3.36354e-02;
		deltaE_sig_p2_mc     =  0.; 			    deltaE_sig_p2_data   =  0.;
		//----------------------------------------------------------------------//
		deltaPhi_mu_p0_mc    =  1.80029e+02;  deltaPhi_mu_p0_data  =  1.80032e+02;
		deltaPhi_mu_p1_mc    = -1.02623e-01;  deltaPhi_mu_p1_data  = -9.61299e-03;
		deltaPhi_mu_p2_mc    =  1.15760e-02;  deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_mc    = -5.22398e-04;  deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_mc   =  1.19447e+01;  deltaPhi_sig_p0_data =  1.26943e+01;
		deltaPhi_sig_p1_mc   = -1.75646e+00;  deltaPhi_sig_p1_data = -1.95900e+00;
		deltaPhi_sig_p2_mc   =  1.67672e-01;  deltaPhi_sig_p2_data =  1.87668e-01;
		deltaPhi_sig_p3_mc   = -5.50310e-03;  deltaPhi_sig_p3_data = -6.08281e-03;
		//----------------------------------------------------------------------//
		deltaK_mu_p0_mc      =  6.80186e-02;  deltaK_mu_p0_data    =  3.52591e-01;
		deltaK_mu_p1_mc      = -9.72781e-03;  deltaK_mu_p1_data    = -9.12544e-02;
		deltaK_mu_p2_mc      = -3.79565e-03;  deltaK_mu_p2_data    =  5.17767e-03;
		deltaK_mu_p3_mc      =  1.31657e-04;  deltaK_mu_p3_data    = -2.20872e-04;
		
		deltaK_sig_p0_mc     =  9.82692e-01;  deltaK_sig_p0_data   =  5.97017e-01;
		deltaK_sig_p1_mc     = -1.92373e-01;  deltaK_sig_p1_data   = -5.77323e-02;
		deltaK_sig_p2_mc     =  2.69333e-02;  deltaK_sig_p2_data   =  1.29216e-02;
		deltaK_sig_p3_mc     = -9.51459e-04;  deltaK_sig_p3_data   = -4.54356e-04;
	}
	
	
	f_deltaE_mu_mc = new TF1("f_deltaE_mu_mc", "pol3", 3.0, 12.0);
	f_deltaE_mu_mc->SetParameters(deltaE_mu_p0_mc, deltaE_mu_p1_mc, 
		deltaE_mu_p2_mc, deltaE_mu_p3_mc);
	
	f_deltaE_sig_mc = new TF1("f_deltaE_sig_mc", 
		"sqrt([0]*[0] + ([1]/sqrt(x))*([1]/sqrt(x)) + ([2]/x)*([2]/x))", 3.0, 12.0);
	f_deltaE_sig_mc->SetParameters(deltaE_sig_p0_mc, deltaE_sig_p1_mc, deltaE_sig_p2_mc);
	
	
	f_deltaE_mu_data = new TF1("f_deltaE_mu_data", "pol3", 3.0, 12.0);
	f_deltaE_mu_data->SetParameters(deltaE_mu_p0_data, deltaE_mu_p1_data, 
		deltaE_mu_p2_data, deltaE_mu_p3_data);
	
	f_deltaE_sig_data = new TF1("f_deltaE_sig_data", 
		"sqrt([0]*[0] + ([1]/sqrt(x))*([1]/sqrt(x)) + ([2]/x)*([2]/x))", 3.0, 12.0);
	f_deltaE_sig_data->SetParameters(deltaE_sig_p0_data, deltaE_sig_p1_data, 
		deltaE_sig_p2_data);
	
	//--------------------------------------------------------------------------------------//
	
	f_deltaPhi_mu_mc = new TF1("f_deltaPhi_mu_mc", "pol3", 3.0, 12.0);
	f_deltaPhi_mu_mc->SetParameters(deltaPhi_mu_p0_mc, deltaPhi_mu_p1_mc, 
		deltaPhi_mu_p2_mc, deltaPhi_mu_p3_mc);
	
	f_deltaPhi_sig_mc = new TF1("f_deltaPhi_sig_mc", "pol3", 3.0, 12.0);
	f_deltaPhi_sig_mc->SetParameters(deltaPhi_sig_p0_mc, deltaPhi_sig_p1_mc, 
		deltaPhi_sig_p2_mc, deltaPhi_sig_p3_mc);
	
	
	f_deltaPhi_mu_data = new TF1("f_deltaPhi_mu_data", "pol3", 3.0, 12.0);
	f_deltaPhi_mu_data->SetParameters(deltaPhi_mu_p0_data, deltaPhi_mu_p1_data, 
		deltaPhi_mu_p2_data, deltaPhi_mu_p3_data);
	
	f_deltaPhi_sig_data = new TF1("f_deltaPhi_sig_data", "pol3", 3.0, 12.0);
	f_deltaPhi_sig_data->SetParameters(deltaPhi_sig_p0_data, deltaPhi_sig_p1_data, 
		deltaPhi_sig_p2_data, deltaPhi_sig_p3_data);
	
	//--------------------------------------------------------------------------------------//
	
	f_deltaK_mu_mc = new TF1("f_deltaK_mu_mc", "pol3", 3.0, 12.0);
	f_deltaK_mu_mc->SetParameters(deltaK_mu_p0_mc, deltaK_mu_p1_mc, 
		deltaK_mu_p2_mc, deltaK_mu_p3_mc);
	
	f_deltaK_sig_mc = new TF1("f_deltaK_sig_mc", "pol3", 3.0, 12.0);
	f_deltaK_sig_mc->SetParameters(deltaK_sig_p0_mc, deltaK_sig_p1_mc, 
		deltaK_sig_p2_mc, deltaK_sig_p3_mc);
	
	
	f_deltaK_mu_data = new TF1("f_deltaK_mu_data", "pol3", 3.0, 12.0);
	f_deltaK_mu_data->SetParameters(deltaK_mu_p0_data, deltaK_mu_p1_data, 
		deltaK_mu_p2_data, deltaK_mu_p3_data);
	
	f_deltaK_sig_data = new TF1("f_deltaK_sig_data", "pol3", 3.0, 12.0);
	f_deltaK_sig_data->SetParameters(deltaK_sig_p0_data, deltaK_sig_p1_data, 
		deltaK_sig_p2_data, deltaK_sig_p3_data);
	
	
	return;
}


double JEventProcessor_compton_simulation::get_vertex_weight(double vertex_z) {
	
	double loc_weight = 1.;
	
	// shift coordinate system so that upstream entrance of target is at z=0:
	
	double loc_z = vertex_z - m_beamZ + (m_target_length/2);
	
	// use attenuation length from XCOM database to calculate probability of photon 
	// absorption.
	
	if(loc_z<0.) {
		return 1.;
	} else if(loc_z>m_target_length) {
		return TMath::Exp(-m_atten * m_target_density * m_target_length);
	}
	loc_weight = TMath::Exp(-m_atten * m_target_density * loc_z);
	
	return loc_weight;
}
