// $Id$
//
//    File: JEventProcessor_compton_analysis_TOF_new.cc
// Created: Wed May 22 03:17:27 PM EDT 2024
// Creator: andrsmit (on Linux ifarm180302.jlab.org 5.14.0-362.13.1.el9_3.x86_64 x86_64)
//

#include "JEventProcessor_compton_analysis_TOF_new.h"

// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_compton_analysis_TOF_new());
  }
} // "C"

JEventProcessor_compton_analysis_TOF_new::JEventProcessor_compton_analysis_TOF_new() {
	
	m_USE_REACTION_WEIGHT = 0;
	gPARMS->SetDefaultParameter("compton_analysis_TOF_new:USE_REACTION_WEIGHT", 
		m_USE_REACTION_WEIGHT);
	
	m_REACTION_CUT_WEIGHT = 1.e4;
	gPARMS->SetDefaultParameter("compton_analysis_TOF_new:REACTION_CUT_WEIGHT", 
		m_REACTION_CUT_WEIGHT);
}

//------------------
// init
//------------------
jerror_t JEventProcessor_compton_analysis_TOF_new::init(void)
{
	rndm_gen = new TRandom3(0);
	
	TDirectory *dir_compton = new TDirectoryFile("compton_analysis_TOF", "compton_analysis_TOF");
	dir_compton->cd();
	
	//---------------------------------------------//
	
	h_beam            = new TH1F("beam", 
		"Is there a beam photon?", 2, -0.5, 1.5);
	h_vertex          = new TH1F("vertex",          
		"Vertex Z Position (unweighted)", 1000, 0., 100.);
	h_vertex_weight   = new TH1F("vertex_weight", 
		"Vertex Weight", 1000, 0., 2.);
	h_vertex_accepted = new TH1F("vertex_accepted", 
		"Vertex Z Position (weighted)",   1000, 0., 100.);
	h_vertex_xy       = new TH2F("vertex_xy", "Vertex Y vs. X; x [cm]; y [cm]", 
		1000, -5., 5., 1000, -5., 5.);
	h_reaction_weight = new TH1F("reaction_weight", 
		"Event Reaction Weight (All Events)", 1000, 0., 1.e7);
	
	//---------------------------------------------//
	
	h_n_showers          = new TH1F("n_showers", 
		"Number of showers in FCAL and CCAL", 20, -0.5, 19.5);
	h_n_showers_ccal     = new TH1F("n_showers_fcal", 
		"Number of showers in FCAL", 20, -0.5, 19.5);
	h_n_showers_fcal     = new TH1F("n_showers_fcal", 
		"Number of showers in FCAL", 20, -0.5, 19.5);
	h_n_showers_cut      = new TH1F("n_showers", 
		"Number of showers in FCAL and CCAL", 20, -0.5, 19.5);
	h_n_showers_ccal_cut = new TH1F("n_showers_fcal_cut", 
		"Number of showers in FCAL", 20, -0.5, 19.5);
	h_n_showers_fcal_cut = new TH1F("n_showers_fcal_cut", 
		"Number of showers in FCAL", 20, -0.5, 19.5);
	
	h_n_good_showers          = new TH1F("n_good_showers", 
		"Number of showers in FCAL and CCAL", 20, -0.5, 19.5);
	h_n_good_showers_ccal     = new TH1F("n_good_showers_fcal", 
		"Number of showers in FCAL", 20, -0.5, 19.5);
	h_n_good_showers_fcal     = new TH1F("n_good_showers_fcal", 
		"Number of showers in FCAL", 20, -0.5, 19.5);
	h_n_good_showers_cut      = new TH1F("n_good_showers", 
		"Number of showers in FCAL and CCAL", 20, -0.5, 19.5);
	h_n_good_showers_ccal_cut = new TH1F("n_good_showers_fcal_cut", 
		"Number of showers in FCAL", 20, -0.5, 19.5);
	h_n_good_showers_fcal_cut = new TH1F("n_good_showers_fcal_cut", 
		"Number of showers in FCAL", 20, -0.5, 19.5);
	
	h_tof_match_case     = new TH1F("tof_match_case", 
		"TOF Match Case", 4, -0.5, 3.5);
	h_tof_match_case_cut = new TH1F("tof_match_case_cut", 
		"TOF Match Case", 4, -0.5, 3.5);
	
	//---------------------------------------------//
	
	TDirectory *dir_deltaE = new TDirectoryFile("deltaE", "deltaE");
	dir_deltaE->cd();
	
	for(int ihist=0; ihist<n_cuts; ihist++) {
		
		h_deltaE_tagh[ihist] = new TH2F(Form("deltaE_tagh_%d",ihist),
			"#DeltaE; TAGH Counter; E_{1} + E_{2} - E_{#gamma} [GeV]",
			274, 0.5, 274.5, 4000, -8.0, 8.0);
		h_deltaE_tagm[ihist] = new TH2F(Form("deltaE_tagm_%d",ihist),
			"#DeltaE; TAGM Counter; E_{1} + E_{2} - E_{#gamma} [GeV]",
			102, 0.5, 102.5, 4000, -8.0, 8.0);
	}
	dir_deltaE->cd("../");
	
	TDirectory *dir_deltaK = new TDirectoryFile("deltaK", "deltaK");
	dir_deltaK->cd();
	
	for(int ihist=0; ihist<n_cuts; ihist++) {
		
		h_deltaK_tagh[ihist] = new TH2F(Form("deltaK_tagh_%d",ihist),
			"#DeltaK; TAGH Counter; E_{Compton} - E_{#gamma} [GeV]",
			274, 0.5, 274.5, 4000, -8.0, 8.0);
		h_deltaK_tagm[ihist] = new TH2F(Form("deltaK_tagm_%d",ihist),
			"#DeltaK; TAGM Counter; E_{Compton} - E_{#gamma} [GeV]",
			102, 0.5, 102.5, 4000, -8.0, 8.0);
	}
	dir_deltaK->cd("../");
	
	TDirectory *dir_deltaK_cut = new TDirectoryFile("deltaK_cut", "deltaK_cut");
	dir_deltaK_cut->cd();
	
	for(int ihist=0; ihist<n_cuts; ihist++) {
		
		h_deltaK_tagh_cut[ihist] = new TH2F(Form("deltaK_tagh_cut_%d",ihist),
			"#DeltaK; TAGH Counter; E_{Compton} - E_{#gamma} [GeV]",
			274, 0.5, 274.5, 4000, -8.0, 8.0);
		h_deltaK_tagm_cut[ihist] = new TH2F(Form("deltaK_tagm_cut_%d",ihist),
			"#DeltaK; TAGM Counter; E_{Compton} - E_{#gamma} [GeV]",
			102, 0.5, 102.5, 4000, -8.0, 8.0);
	}
	dir_deltaK_cut->cd("../");
	
	TDirectory *dir_elas = new TDirectoryFile("elas_vs_deltaE", "elas_vs_deltaE");
	dir_elas->cd();
	
	for(int ihist=0; ihist<n_cuts; ihist++) {
		
		h_elas_vs_deltaE[ihist] = new TH2F(Form("deltaK_tagh_cut_%d",ihist),
			"2-D Elasticity; E_{1}+E_{2} - E_{#gamma} [GeV]; E_{1}+E_{2} - E_{Compton} [GeV]",
			500, -8.0, 8.0, 500, -8.0, 8.0);
	}
	dir_elas->cd("../");
	
	TDirectory *dir_mgg = new TDirectoryFile("mgg_vs_deltaE", "mgg_vs_deltaE");
	dir_mgg->cd();
	
	for(int ihist=0; ihist<n_cuts; ihist++) {
		
		h_mgg_vs_deltaE[ihist] = new TH2F(Form("deltaK_tagh_cut_%d",ihist),
			"Invariant Mass vs. #DeltaE; E_{1}+E_{2} - E_{#gamma} [GeV]; Mass [GeV/c^{2}",
			500, -8.0, 8.0, 250, 0., 0.5);
	}
	dir_mgg->cd("../");
	
	//---------------------------------------------//
	
	TDirectory *dir_extra_fcal = new TDirectoryFile("extra_fcal", "extra_fcal");
	dir_extra_fcal->cd();
	
	h_extra_fcal_shower_energy          = new TH1F("shower_energy", 
		"Energy of extra FCAL shower; E_{extra} [GeV]", 1000, 0., 10.);
	h_extra_fcal_shower_distance        = new TH1F("shower_distance",
		"Distrance from main FCAL shower; #Deltad [cm]", 100, 0., 10.);
	h_extra_fcal_shower_deltaPhi        = new TH1F("shower_deltaPhi",
		"#phi_{extra} - #phi_{FCAL}; [deg.]", 3600, 0., 360.);
	h_extra_fcal_shower_deltaPhi_ccal   = new TH1F("shower_deltaPhi_ccal",
		"#phi_{extra} - #phi_{CCAL}; [deg.]", 3600, 0., 360.);
	h_extra_fcal_shower_elasticity      = new TH1F("shower_elasticity",
		"#DeltaE; E_{FCAL} + E_{CCAL} - E_{#gamma} [GeV]", 4000, -8.0, 8.0);
	h_extra_fcal_shower_elasticity_new  = new TH1F("shower_elasticity_new",
		"#DeltaE new; E_{CCAL} + E_{extra} - E_{#gamma} [GeV]", 4000, -8.0, 8.0);
	h_extra_fcal_shower_elasticity_corr = new TH1F("shower_elasticity_corr",
		"Corrected #DeltaE; E_{FCAL} + E_{CCAL} + E_{extra} - E_{#gamma} [GeV]", 4000, -8.0, 8.0);
	h_extra_fcal_shower_xy              = new TH2F("shower_xy",
		"Position of extra FCAL shower", 500, -12., 12., 500, -12., 12.);
	
	h_extra_fcal_shower_energy_cut          = new TH1F("shower_energy_cut", 
		"Energy of extra FCAL shower; E_{extra} [GeV]", 1000, 0., 10.);
	h_extra_fcal_shower_distance_cut        = new TH1F("shower_distance_cut",
		"Distrance from main FCAL shower; #Deltad [cm]", 100, 0., 10.);
	h_extra_fcal_shower_deltaPhi_cut        = new TH1F("shower_deltaPhi_cut",
		"#phi_{extra} - #phi_{FCAL}; [deg.]", 3600, 0., 360.);
	h_extra_fcal_shower_deltaPhi_ccal_cut   = new TH1F("shower_deltaPhi_ccal_cut",
		"#phi_{extra} - #phi_{CCAL}; [deg.]", 3600, 0., 360.);
	h_extra_fcal_shower_elasticity_cut      = new TH1F("shower_elasticity_cut",
		"#DeltaE; E_{FCAL} + E_{CCAL} - E_{#gamma} [GeV]", 4000, -8.0, 8.0);
	h_extra_fcal_shower_elasticity_new_cut  = new TH1F("shower_elasticity_new_cut",
		"#DeltaE new; E_{CCAL} + E_{extra} - E_{#gamma} [GeV]", 4000, -8.0, 8.0);
	h_extra_fcal_shower_elasticity_corr_cut = new TH1F("shower_elasticity_corr_cut",
		"Corrected #DeltaE; E_{FCAL} + E_{CCAL} + E_{extra} - E_{#gamma} [GeV]", 4000, -8.0, 8.0);
	h_extra_fcal_shower_xy_cut              = new TH2F("shower_xy_cut",
		"Position of extra FCAL shower", 500, -12., 12., 500, -12., 12.);
	
	dir_extra_fcal->cd("../");
	
	
	TDirectory *dir_extra_ccal = new TDirectoryFile("extra_ccal", "extra_ccal");
	dir_extra_ccal->cd();
	
	h_extra_ccal_shower_energy          = new TH1F("shower_energy", 
		"Energy of extra CCAL shower; E_{extra} [GeV]", 1000, 0., 10.);
	h_extra_ccal_shower_distance        = new TH1F("shower_distance",
		"Distrance from main CCAL shower; #Deltad [cm]", 100, 0., 10.);
	h_extra_ccal_shower_deltaPhi        = new TH1F("shower_deltaPhi",
		"#phi_{extra} - #phi_{CCAL}; [deg.]", 3600, 0., 360.);
	h_extra_ccal_shower_deltaPhi_fcal   = new TH1F("shower_deltaPhi_fcal",
		"#phi_{extra} - #phi_{FCAL}; [deg.]", 3600, 0., 360.);
	h_extra_ccal_shower_elasticity      = new TH1F("shower_elasticity",
		"#DeltaE; E_{FCAL} + E_{CCAL} - E_{#gamma} [GeV]", 4000, -8.0, 8.0);
	h_extra_ccal_shower_elasticity_new  = new TH1F("shower_elasticity_new",
		"#DeltaE new; E_{FCAL} + E_{extra} - E_{#gamma} [GeV]", 4000, -8.0, 8.0);
	h_extra_ccal_shower_elasticity_corr = new TH1F("shower_elasticity_corr",
		"Corrected #DeltaE; E_{FCAL} + E_{CCAL} + E_{extra} - E_{#gamma} [GeV]", 4000, -8.0, 8.0);
	h_extra_ccal_shower_xy              = new TH2F("shower_xy",
		"Position of extra CCAL shower", 500, -12., 12., 500, -12., 12.);
	
	h_extra_ccal_shower_energy_cut          = new TH1F("shower_energy_cut", 
		"Energy of extra CCAL shower; E_{extra} [GeV]", 1000, 0., 10.);
	h_extra_ccal_shower_distance_cut        = new TH1F("shower_distance_cut",
		"Distrance from main CCAL shower; #Deltad [cm]", 100, 0., 10.);
	h_extra_ccal_shower_deltaPhi_cut        = new TH1F("shower_deltaPhi_cut",
		"#phi_{extra} - #phi_{CCAL}; [deg.]", 3600, 0., 360.);
	h_extra_ccal_shower_deltaPhi_fcal_cut   = new TH1F("shower_deltaPhi_fcal_cut",
		"#phi_{extra} - #phi_{FCAL}; [deg.]", 3600, 0., 360.);
	h_extra_ccal_shower_elasticity_cut      = new TH1F("shower_elasticity_cut",
		"#DeltaE; E_{FCAL} + E_{CCAL} - E_{#gamma} [GeV]", 4000, -8.0, 8.0);
	h_extra_ccal_shower_elasticity_new_cut  = new TH1F("shower_elasticity_new_cut",
		"#DeltaE new; E_{FCAL} + E_{extra} - E_{#gamma} [GeV]", 4000, -8.0, 8.0);
	h_extra_ccal_shower_elasticity_corr_cut = new TH1F("shower_elasticity_corr_cut",
		"Corrected #DeltaE; E_{FCAL} + E_{CCAL} + E_{extra} - E_{#gamma} [GeV]", 4000, -8.0, 8.0);
	h_extra_ccal_shower_xy_cut              = new TH2F("shower_xy_cut",
		"Position of extra CCAL shower", 500, -12., 12., 500, -12., 12.);
	
	dir_extra_ccal->cd("../");
	
	
	
	dir_compton->cd("../");
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_compton_analysis_TOF_new::brun(JEventLoop *eventLoop, int32_t runnumber)
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
		
		// (2/4/2024): Correction to alignment after Simon updated beam spot:
		
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
		if(runnumber<61483) {
			m_ccalX_new = 0.083;
			m_ccalY_new = 0.148;
		} else {
			m_ccalX_new = 0.083;
			m_ccalY_new = 0.119;
		}
		
		m_fcalX_new =  0.617;
		m_fcalY_new = -0.002;
		
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
	
	m_fcal_correction.SetXYZ(m_fcalX_new-m_fcalX, m_fcalY_new-m_fcalY, 0.);
	m_ccal_correction.SetXYZ(m_ccalX_new-m_ccalX, m_ccalY_new-m_ccalY, 0.);
	
	set_cuts(runnumber);
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_compton_analysis_TOF_new::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
	//--------------------------------------------------------------------------------------//
	// Get all necessary data objects:
	
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
	
	vector<const DMCThrown*> locMCThrown;
	eventLoop->Get(locMCThrown);
	
	vector<const DMCReaction*> locMCReaction;
	eventLoop->Get(locMCReaction);
  
	japp->RootFillLock(this);
	
	//--------------------------------------------------------------------------------------//
	// If MC, apply accept-reject filter to make realistic distribution of target z-position:
	
	int loc_is_mc = 0;
	double loc_reaction_weight = 1.0;
	
	if(locMCThrown.size()) {
		
		loc_is_mc = 1;
		
		double vertex_z        = locMCThrown[0]->position().Z();
		double vertex_weight   = get_vertex_weight(vertex_z);
		
		h_vertex->Fill(vertex_z);
		h_vertex_weight->Fill(vertex_weight);
		
		if(vertex_weight < rndm_gen->Uniform()) {
			japp->RootFillUnLock(this);
			return NOERROR;
		}
		
		loc_reaction_weight = locMCReaction[0]->weight;
		h_reaction_weight->Fill(loc_reaction_weight);
		
		if(m_USE_REACTION_WEIGHT && loc_reaction_weight > m_REACTION_CUT_WEIGHT) {
			japp->RootFillUnLock(this);
			return NOERROR;
		}
		
		h_vertex_accepted->Fill(vertex_z);
		h_vertex_xy->Fill(locMCThrown[0]->position().X(), locMCThrown[0]->position().Y());
	}
	
	//--------------------------------------------------------------------------------------//
	// Check Trigger (if data):
	
	const DL1Trigger *locTrig = NULL;
	try {
		eventLoop->GetSingle(locTrig);
	} catch (...) {}
	if (locTrig == NULL && loc_is_mc==0) { 
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	
	//--------------------------------------------------------------------------------------//
	// Get RF Bunch (return NOERROR if none):
	
	const DEventRFBunch *locRFBunch = NULL;
	try { 
		eventLoop->GetSingle(locRFBunch, "CalorimeterOnly");
	} catch (...) { 
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	double locRFTime = locRFBunch->dTime;
	
	//--------------------------------------------------------------------------------------//
	// Check if there is a beam photon reconstructed (for MC):
	
	int locNBeamPhotons = (int)locBeamPhotons.size();
	if(locNBeamPhotons) h_beam->Fill(1);
	else                h_beam->Fill(0);
	
	if(locRFBunch->dNumParticleVotes < 2) {
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	
	//--------------------------------------------------------------------------------------//
	// Check number of "good" FCAL and CCAL showers:
	
	int locNFCALShowers = 0, locNGoodFCALShowers = 0;
	vector<const DFCALShower*> locGoodFCALShowers;
	
	for(vector<const DFCALShower*>::const_iterator show = locFCALShowers.begin(); 
		show != locFCALShowers.end(); show++) {
		
    	DVector3 loc_pos = (*show)->getPosition_log() - locVertex + m_fcal_correction;
    	double loc_t     = (*show)->getTime() - (loc_pos.Mag()/c) - locRFTime;
		
		int fid_cut = fcal_fiducial_cut(loc_pos, locVertex, 1.0);
		if((fabs(loc_t) < FCAL_RF_time_cut)) {
			locNFCALShowers++;
			if((*show)->getEnergy() > FCAL_min_energy_cut) {
				locNGoodFCALShowers++;
				if(!fid_cut) {
					locGoodFCALShowers.push_back((*show));
				}
			}
		}
	}
	
	int locNCCALShowers = 0, locNGoodCCALShowers = 0;
	vector<const DCCALShower*> locGoodCCALShowers;
	
	for(vector<const DCCALShower*>::const_iterator show = locCCALShowers.begin();
		show != locCCALShowers.end(); show++) {
		
		DVector3 loc_pos((*show)->x1, (*show)->y1, (*show)->z);
		loc_pos = loc_pos - locVertex + m_ccal_correction;
		double loc_t = (*show)->time - (loc_pos.Mag()/c) - locRFTime;
		
		int fid_cut = ccal_fiducial_cut(loc_pos, locVertex);
		if((fabs(loc_t) < CCAL_RF_time_cut)) {
			locNCCALShowers++;
			if((*show)->E > CCAL_min_energy_cut) {
				locNGoodCCALShowers++;
				if(!fid_cut) {
					locGoodCCALShowers.push_back((*show));
				}
			}
		}
	}
	
	//--------------------------------------------------------------------------------------//
	// Check FCAL-CCAL Pairs
	
	for(vector< const DFCALShower* >::const_iterator show1 = locGoodFCALShowers.begin(); 
		show1 != locGoodFCALShowers.end(); show1++) {
		
		double e1     = (*show1)->getEnergy();
		DVector3 pos1 = (*show1)->getPosition_log() - locVertex + m_fcal_correction;
		
		double theta1 = pos1.Theta();
		
		//--------------------------------------------------------------------------------//
		// check for matches with TOF:
		
		double tof_dx, tof_dy, tof_dt;
		check_TOF_match(pos1, locRFTime, locVertex, locTOFPoints, tof_dx, tof_dy, tof_dt, 1.0);
		double tof_dr = sqrt(pow(tof_dx,2.0) + pow(tof_dy,2.0));
		
		int tof_match = 0;
		if(tof_dr < 8.0) tof_match = 1;
		
		//--------------------------------------------------------------------------------//
		
		for(vector<const DCCALShower*>::const_iterator show2 = locGoodCCALShowers.begin(); 
			show2 != locGoodCCALShowers.end(); show2++) {
			
			double e2  = (*show2)->E;
			DVector3 pos2( (*show2)->x1, (*show2)->y1, (*show2)->z );
			pos2       = pos2 - locVertex + m_ccal_correction;
			
			double theta2 = pos2.Theta();
			
			// calculate deltaPhi:
			
			double deltaPhi = fabs(pos2.Phi() - pos1.Phi()) * (180./TMath::Pi());
			
			// calculate invariant mass (assuming massless particles):
			
			double cos12   = (pos1.X()*pos2.X() + pos1.Y()*pos2.Y() + pos1.Z()*pos2.Z()) 
				/ (pos1.Mag()*pos2.Mag());
			double invmass = sqrt(2.0*e1*e2*(1.-cos12));
			
			// loop over beam photons:
			
			for(vector<const DBeamPhoton*>::const_iterator gam = locBeamPhotons.begin();
				gam != locBeamPhotons.end(); gam++) {
					
				double eb = (*gam)->lorentzMomentum().E();
				double tb = (*gam)->time();
				
				double brfdt = tb - locRFTime;
				
				double fill_weight = 0.0;
				
				if(fabs(brfdt) < 2.004)
					fill_weight =  1.0;
				else if((-(2.004 + 10.*4.008)<=brfdt && brfdt<=-(2.004 + 5.*4.008))
					||((2.004 + 5.*4.008)<=brfdt && brfdt<=(2.004 + 10.*4.008)))
					fill_weight = -0.1;
				else 
					continue;
				
				if(eb < BEAM_min_energy_cut) continue;
				
				DetectorSystem_t tag_sys = (*gam)->dSystem;
				int tag_counter = (*gam)->dCounter;
				
				double deltaE = (e1 + e2) - (eb + m_e);
				
				// calculate expected shower energies from beam energy and scattering angles:
				/*
				double ecomp1 = e1 / (1. - (e1/m_e)*(1.-cos(theta1)));
				double ecomp2 = e2 / (1. - (e2/m_e)*(1.-cos(theta2)));
				
				double deltaK = m_e * sin(theta1+theta2) / 
					(sin(theta1) + sin(theta2) - sin(theta1+theta2));
				deltaK       -= eb;
				*/
				double deltaK = m_e*(sin(theta1) / (2.*pow(sin(theta1/2.),2.)*tan(theta2))) - m_e;
				deltaK       -= eb;
				
				//-------------------------------------------------------------//
				// Cuts:
				
				double deltaPhi_smeared = deltaPhi;
				double deltaE_smeared   = deltaE;
				double deltaK_smeared   = deltaK;
				
				double deltaE_mu_data  = f_deltaE_mu_data->Eval(eb);
				double deltaE_sig_data = eb * f_deltaE_sig_data->Eval(eb);
				
				double deltaPhi_mu_data  = f_deltaPhi_mu_data->Eval(eb);
				double deltaPhi_sig_data = f_deltaPhi_sig_data->Eval(eb);
				
				double deltaK_mu_data  = f_deltaK_mu_data->Eval(eb);
				double deltaK_sig_data = f_deltaK_sig_data->Eval(eb);
				
				if(loc_is_mc) {
					
					double deltaE_mu_mc    = f_deltaE_mu_mc->Eval(eb);
					double deltaE_sig_mc   = eb * f_deltaE_sig_mc->Eval(eb);
					
					double deltaPhi_mu_mc  = f_deltaPhi_mu_mc->Eval(eb);
					double deltaPhi_sig_mc = f_deltaPhi_sig_mc->Eval(eb);
					
					double deltaK_mu_mc    = f_deltaK_mu_mc->Eval(eb);
					double deltaK_sig_mc   = f_deltaK_sig_mc->Eval(eb);
					
					deltaE_smeared = deltaE + (deltaE_mu_data - deltaE_mu_mc);
					if(deltaE_sig_data > deltaE_sig_mc) {
						double loc_deltaE_smear = sqrt(pow(deltaE_sig_data,2.0)
							- pow(deltaE_sig_mc,2.0));
						deltaE_smeared += rndm_gen->Gaus(0.,loc_deltaE_smear);
					}
					
					deltaPhi_smeared = deltaPhi + (deltaPhi_mu_data - deltaPhi_mu_mc);
					if(deltaPhi_sig_data > deltaPhi_sig_mc) {
						double loc_deltaPhi_smear = sqrt(pow(deltaPhi_sig_data,2.0)
							- pow(deltaPhi_sig_mc,2.0));
						deltaPhi_smeared += rndm_gen->Gaus(0.,loc_deltaPhi_smear);
					}
					
					deltaK_smeared = deltaK + (deltaK_mu_data - deltaK_mu_mc);
					if(deltaK_sig_data > deltaK_sig_mc) {
						double loc_deltaK_smear = sqrt(pow(deltaK_sig_data,2.0)
							- pow(deltaK_sig_mc,2.0));
						deltaK_smeared += rndm_gen->Gaus(0.,loc_deltaK_smear);
					}
				}
				
				int phi_cut = 0;
				if(bfield_val) {
					if(fabs(deltaPhi_smeared - 180.) < 50.) {
						phi_cut = 1;
					}
				} else {
					if(fabs(deltaPhi_smeared - deltaPhi_mu_data) < 5.0*deltaPhi_sig_data) {
						phi_cut = 1;
					}
				}
				int e_cut = 0;
				if(fabs(deltaE_smeared - deltaE_mu_data) < 5.0*deltaE_sig_data) e_cut = 1;
				
				int k_cut = 0;
				if(fabs(deltaK_smeared - deltaK_mu_data) < 5.0*deltaK_sig_data) k_cut = 1;
				
				
				/*
				We want to look at how the following set of cuts change the data:
					1. Multiplicity cut (3 options)
						a. No extra showers
						b. Allow extra showers if they are below some energy threshold
						c. Allow any ammount of extra showers
					2. TOF Veto
						a. No match betweeen FCAL and TOF
						b. Maybe a match, maybe not (no information)
					3. FCAL fiducial cut
						a. One layer cut
						b. Two layer cut
				*/
				
				int cut_vals[n_cuts];
				for(int icut = 0; icut < n_cuts; icut++) cut_vals[icut] = 0.;
				
				int loc_fid_cut = fcal_fiducial_cut(pos1, locVertex, 2.0);
				
				cut_vals[0] = 1;
				if(!tof_match) {
					cut_vals[1] = 1;
					if(!loc_fid_cut) {
						cut_vals[2] = 1;
					}
				} 
				if(!loc_fid_cut) cut_vals[3] = 1;
				
				if(locNGoodFCALShowers==1 && locNGoodCCALShowers==1) {
					cut_vals[4] = 1;
					if(!tof_match) {
						cut_vals[5] = 1;
						if(!loc_fid_cut) {
							cut_vals[6] = 1;
						}
					}
					if(!loc_fid_cut) cut_vals[7] = 1;
				}
				
				if(locNFCALShowers==1 && locNCCALShowers==1) {
					cut_vals[8] = 1;
					if(!tof_match) {
						cut_vals[9] = 1;
						if(!loc_fid_cut) {
							cut_vals[10] = 1;
						}
					}
					if(!loc_fid_cut) cut_vals[11] = 1;
				}
				
				for(int icut = 0; icut < n_cuts; icut++) {
					
					if(tag_sys==SYS_TAGH) {
						
						if(cut_vals[icut]) {
							h_deltaE_tagh[icut]->Fill(tag_counter, deltaE, fill_weight);
							if(e_cut && phi_cut) {
								h_deltaK_tagh[icut]->Fill(tag_counter, deltaK, fill_weight);
								if(k_cut) {
									h_deltaK_tagh_cut[icut]->Fill(tag_counter, deltaK, 
										fill_weight);
								}
							}
						}
					} else if(tag_sys==SYS_TAGM) {
						
						if(cut_vals[icut]) {
							h_deltaE_tagm[icut]->Fill(tag_counter, deltaE, fill_weight);
							if(e_cut && phi_cut) {
								h_deltaK_tagm[icut]->Fill(tag_counter, deltaK, fill_weight);
								if(k_cut) {
									h_deltaK_tagm_cut[icut]->Fill(tag_counter, deltaK, 
										fill_weight);
								}
							}
						}
					}
					
					h_elas_vs_deltaE[icut]->Fill(deltaE, (deltaE-deltaK), fill_weight);
					h_mgg_vs_deltaE[icut]->Fill(deltaE, invmass, fill_weight);
				}
				
				//if(!e_cut || !k_cut || !phi_cut) continue;
				
				// If there are more than 1 shower, how many and where:
				
				h_n_showers->Fill(locNFCALShowers+locNCCALShowers);
				h_n_showers_ccal->Fill(locNCCALShowers);
				h_n_showers_fcal->Fill(locNFCALShowers);
				h_n_good_showers->Fill(locNGoodFCALShowers+locNGoodCCALShowers);
				h_n_good_showers_ccal->Fill(locNGoodCCALShowers);
				h_n_good_showers_fcal->Fill(locNGoodFCALShowers);
				
				if(e_cut && k_cut && phi_cut) {
					h_n_showers_cut->Fill(locNFCALShowers+locNCCALShowers);
					h_n_showers_ccal_cut->Fill(locNCCALShowers);
					h_n_showers_fcal_cut->Fill(locNFCALShowers);
					h_n_good_showers_cut->Fill(locNGoodFCALShowers+locNGoodCCALShowers);
					h_n_good_showers_ccal_cut->Fill(locNGoodCCALShowers);
					h_n_good_showers_fcal_cut->Fill(locNGoodFCALShowers);
				}
				
				if(locNFCALShowers==2 && locNCCALShowers==1) {
					
					// loop over FCAL showers and find the extra ones:
					for(vector<const DFCALShower*>::const_iterator loc_show=locFCALShowers.begin();
						loc_show != locFCALShowers.end(); loc_show++) {
						
						DVector3 loc_pos = (*loc_show)->getPosition_log() - locVertex 
							+ m_fcal_correction;
						double loc_t = (*loc_show)->getTime() - (loc_pos.Mag()/c);
						if((fabs(loc_t-locRFTime) >= FCAL_RF_time_cut)) continue;
						
						double loc_e = (*loc_show)->getEnergy();
						if(loc_e==e1) continue;
						
						// check if this shower can be matched to a TOF hit:
						double loc_tof_dx, loc_tof_dy, loc_tof_dt;
						check_TOF_match(loc_pos, locRFTime, locVertex, locTOFPoints, 
							loc_tof_dx, loc_tof_dy, loc_tof_dt, 1.0);
						double loc_tof_dr = sqrt(pow(loc_tof_dx,2.0) + pow(loc_tof_dy,2.0));
						int loc_tof_match = 0;
						if(loc_tof_dr < 8.0) loc_tof_match = 1;
						
						int tof_match_fill_val;
						if(!tof_match && !loc_tof_match) {
							tof_match_fill_val = 0;
						} else if(!tof_match && loc_tof_match) {
							tof_match_fill_val = 1;
						} else if(tof_match && !loc_tof_match) {
							tof_match_fill_val = 2;
						} else {
							tof_match_fill_val = 3;
						}
						h_tof_match_case->Fill(tof_match_fill_val, fill_weight);
						if(e_cut && k_cut && phi_cut) {
							h_tof_match_case_cut->Fill(tof_match_fill_val, fill_weight);
						}
						
						// project shower to same z-position:
						
						loc_pos *= (pos1.Z()/loc_pos.Z());
						
						// distance from current shower:
						double loc_dr = (loc_pos-pos1).Mag();
						
						double loc_deltaPhi      = (loc_pos.Phi()-pos1.Phi())*(180./TMath::Pi());
						double loc_deltaPhi_ccal = (loc_pos.Phi()-pos2.Phi())*(180./TMath::Pi());
						
						h_extra_fcal_shower_energy->Fill(loc_e, fill_weight);
						h_extra_fcal_shower_distance->Fill(loc_dr, fill_weight);
						h_extra_fcal_shower_deltaPhi->Fill(loc_deltaPhi, fill_weight);
						h_extra_fcal_shower_deltaPhi_ccal->Fill(loc_deltaPhi_ccal, fill_weight);
						h_extra_fcal_shower_elasticity->Fill(deltaE, fill_weight);
						h_extra_fcal_shower_elasticity_new->Fill(deltaE+loc_e-e1, fill_weight);
						h_extra_fcal_shower_elasticity_corr->Fill(deltaE+loc_e, fill_weight);
						h_extra_fcal_shower_xy->Fill(loc_pos.X(), loc_pos.Y(), fill_weight);
						
						if(e_cut && k_cut && phi_cut) {
							h_extra_fcal_shower_energy_cut->Fill(loc_e, fill_weight);
							h_extra_fcal_shower_distance_cut->Fill(loc_dr, fill_weight);
							h_extra_fcal_shower_deltaPhi_cut->Fill(loc_deltaPhi, fill_weight);
							h_extra_fcal_shower_deltaPhi_ccal_cut->Fill(loc_deltaPhi_ccal, 
								fill_weight);
							h_extra_fcal_shower_elasticity_cut->Fill(deltaE, fill_weight);
							h_extra_fcal_shower_elasticity_new_cut->Fill(deltaE+loc_e-e1, 
								fill_weight);
							h_extra_fcal_shower_elasticity_corr_cut->Fill(deltaE+loc_e, 
								fill_weight);
							h_extra_fcal_shower_xy_cut->Fill(loc_pos.X(), loc_pos.Y(), 
								fill_weight);
						}
					}
				}
				
				if(locNFCALShowers==1 && locNCCALShowers==2) {
					
					// loop over CCAL showers and find the extra ones:
					for(vector<const DCCALShower*>::const_iterator loc_show=locCCALShowers.begin();
						loc_show != locCCALShowers.end(); loc_show++) {
						
						DVector3 loc_pos((*loc_show)->x1, (*loc_show)->y1, (*loc_show)->z);
						loc_pos += (m_ccal_correction - locVertex);
						
						double loc_t = (*loc_show)->time - (loc_pos.Mag()/c);
						if((fabs(loc_t-locRFTime) >= CCAL_RF_time_cut)) continue;
						
						double loc_e = (*loc_show)->E;
						if(loc_e==e1) continue;
						
						// project shower to same z-position:
						
						loc_pos *= (pos2.Z()/loc_pos.Z());
						
						// distance from current shower:
						double loc_dr = (loc_pos-pos2).Mag();
						
						double loc_deltaPhi      = (loc_pos.Phi()-pos2.Phi())*(180./TMath::Pi());
						double loc_deltaPhi_fcal = (loc_pos.Phi()-pos1.Phi())*(180./TMath::Pi());
						
						h_extra_ccal_shower_energy->Fill(loc_e, fill_weight);
						h_extra_ccal_shower_distance->Fill(loc_dr, fill_weight);
						h_extra_ccal_shower_deltaPhi->Fill(loc_deltaPhi, fill_weight);
						h_extra_ccal_shower_deltaPhi_fcal->Fill(loc_deltaPhi_fcal, fill_weight);
						h_extra_ccal_shower_elasticity->Fill(deltaE, fill_weight);
						h_extra_ccal_shower_elasticity_new->Fill(deltaE+loc_e-e2, fill_weight);
						h_extra_ccal_shower_elasticity_corr->Fill(deltaE+loc_e, fill_weight);
						h_extra_ccal_shower_xy->Fill(loc_pos.X(), loc_pos.Y(), fill_weight);
						
						if(e_cut && k_cut && phi_cut) {
							h_extra_ccal_shower_energy_cut->Fill(loc_e, fill_weight);
							h_extra_ccal_shower_distance_cut->Fill(loc_dr, fill_weight);
							h_extra_ccal_shower_deltaPhi_cut->Fill(loc_deltaPhi, fill_weight);
							h_extra_ccal_shower_deltaPhi_fcal_cut->Fill(loc_deltaPhi_fcal, 
								fill_weight);
							h_extra_ccal_shower_elasticity_cut->Fill(deltaE, fill_weight);
							h_extra_ccal_shower_elasticity_new_cut->Fill(deltaE+loc_e-e2, 
								fill_weight);
							h_extra_ccal_shower_elasticity_corr_cut->Fill(deltaE+loc_e, 
								fill_weight);
							h_extra_ccal_shower_xy_cut->Fill(loc_pos.X(), loc_pos.Y(), 
								fill_weight);
						}
					}
				}
				
			} // end DBeamPhoton loop
		} // end DCCALShower loop
	} // end DFCALShower loop
	
	japp->RootFillUnLock(this);  // Release root lock
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_compton_analysis_TOF_new::erun(void)
{
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_compton_analysis_TOF_new::fini(void)
{
  return NOERROR;
}


int JEventProcessor_compton_analysis_TOF_new::fcal_fiducial_cut(DVector3 pos, DVector3 vertex, 
	double layer_cut)
{
	int fid_cut = 0;
	
	double fcal_inner_layer_cut = (1.5 + layer_cut) * 4.0157;
	
	double fcal_face_x = vertex.X() + (pos.X() * (m_fcalZ - vertex.Z())/pos.Z());
	double fcal_face_y = vertex.Y() + (pos.Y() * (m_fcalZ - vertex.Z())/pos.Z());
	
	fcal_face_x -= m_fcalX_new;
	fcal_face_y -= m_fcalY_new;
	
	if((-1.*fcal_inner_layer_cut < fcal_face_x && fcal_face_x < fcal_inner_layer_cut)
		&& (-1.*fcal_inner_layer_cut < fcal_face_y 
		&& fcal_face_y < fcal_inner_layer_cut)) fid_cut = 1;
	
	// only apply the next fiducial cut for runs from phase-I:
	
	if(phase_val < 2) {
		if((-32.<fcal_face_y && fcal_face_y<-20.) && (-8.<fcal_face_x && fcal_face_x<4.)) {
			fid_cut = 1;
		}
	}
	
	return fid_cut;
}


int JEventProcessor_compton_analysis_TOF_new::ccal_fiducial_cut(DVector3 pos, DVector3 vertex)
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


void JEventProcessor_compton_analysis_TOF_new::set_cuts(int32_t runnumber)
{
	
	if(runnumber > 60000 && runnumber < 61355) {
		
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
		
	} else if(runnumber>80000 && runnumber<81473) {
		
		// He Target, Field OFF, Phase-II
		
		deltaE_mu_p0_mc    = -4.52241e-02;
		deltaE_mu_p1_mc    = -9.08425e-03;
		deltaE_mu_p2_mc    =  8.94119e-04;
		deltaE_mu_p3_mc    =  0.;
		
		deltaE_sig_p0_mc   =  8.69176e-03;
		deltaE_sig_p1_mc   =  3.88600e-02;
		deltaE_sig_p2_mc   =  2.15792e-09;
		//------------------------------//
		deltaE_mu_p0_data  =  1.73121e-03;
		deltaE_mu_p1_data  = -6.33500e-03;
		deltaE_mu_p2_data  =  3.77085e-05;
		deltaE_mu_p3_data  =  0.;
		
		deltaE_sig_p0_data =  1.45365e-02;
		deltaE_sig_p1_data =  3.68813e-02;
		deltaE_sig_p2_data =  1.07375e-08;
		
		
		deltaPhi_mu_p0_mc    =  1.80041e+02;
		deltaPhi_mu_p1_mc    = -1.00150e-01;
		deltaPhi_mu_p2_mc    =  9.06181e-03;
		deltaPhi_mu_p3_mc    = -4.27409e-04;
		
		deltaPhi_sig_p0_mc   =  1.16985e+01;
		deltaPhi_sig_p1_mc   = -1.61493e+00;
		deltaPhi_sig_p2_mc   =  1.69828e-01;
		deltaPhi_sig_p3_mc   = -5.87253e-03;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.78347e+02;
		deltaPhi_mu_p1_data  =  5.11289e-01;
		deltaPhi_mu_p2_data  = -6.13319e-02;
		deltaPhi_mu_p3_data  =  2.33165e-03;
		
		deltaPhi_sig_p0_data =  1.28404e+01;
		deltaPhi_sig_p1_data = -2.01750e+00;
		deltaPhi_sig_p2_data =  1.95811e-01;
		deltaPhi_sig_p3_data = -6.46208e-03;
		
		
		deltaK_mu_p0_mc    =  6.53381e-02;
		deltaK_mu_p1_mc    = -1.53039e-02;
		deltaK_mu_p2_mc    =  1.31349e-04;
		deltaK_mu_p3_mc    = -6.34173e-06;
		
		deltaK_sig_p0_mc   =  1.69023e-02;
		deltaK_sig_p1_mc   =  2.17646e-02;
		deltaK_sig_p2_mc   =  2.03976e-04;
		deltaK_sig_p3_mc   = -2.23966e-05;
		//--------------------------------//
		deltaK_mu_p0_data    =  2.10762e-01;
		deltaK_mu_p1_data    = -5.67320e-02;
		deltaK_mu_p2_data    =  4.09237e-03;
		deltaK_mu_p3_data    = -1.32960e-04;
		
		deltaK_sig_p0_data   =  1.99286e-01;
		deltaK_sig_p1_data   = -3.50894e-02;
		deltaK_sig_p2_data   =  6.33112e-03;
		deltaK_sig_p3_data   = -2.39232e-04;
		
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
	
	//--------------------------------------------------------------------------------------//
	
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
	
	//--------------------------------------------------------------------------------------//
	
	return;
}

double JEventProcessor_compton_analysis_TOF_new::get_vertex_weight(double vertex_z) {
	
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

void JEventProcessor_compton_analysis_TOF_new::check_TOF_match(DVector3 pos1, double rfTime, 
	DVector3 vertex, vector<const DTOFPoint*> tof_points, double &dx_min, double &dy_min, 
	double &dt_min, double rf_time_cut) {
	
	dx_min = 1000.;
	dy_min = 1000.;
	dt_min = 1000.;
	
	for(vector<const DTOFPoint*>::const_iterator tof = tof_points.begin(); 
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
