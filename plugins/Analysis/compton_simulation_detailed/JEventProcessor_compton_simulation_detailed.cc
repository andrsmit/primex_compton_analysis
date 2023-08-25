// $Id$
//
//    File: JEventProcessor_compton_simulation_detailed.cc
// Created: Thu Feb 10 21:33:24 EST 2022
// Creator: andrsmit (on Linux ifarm1802.jlab.org 3.10.0-1160.11.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_compton_simulation_detailed.h"

// Routine used to create our JEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_compton_simulation_detailed());
}
} // "C"


//------------------
// init
//------------------
jerror_t JEventProcessor_compton_simulation_detailed::init(void)
{
	
	rndm_gen = new TRandom3(0);
	
	TDirectory *dir_reweight = new TDirectoryFile("compton_simulation_detailed", 
		"compton_simulation_detailed");
	dir_reweight->cd();
	
	//-------------------------------------------------------------------------------------------//
	
	for(int i=0; i<2; i++) {
		h_theta_thrown_photon[i] = new TH1F(Form("theta_thrown_photon_%d",i), 
			"Scattering Angle of Primary Photon", 1000, 0., 6.);
		h_theta_thrown_electron[i] = new TH1F(Form("theta_thrown_electron_%d",i), 
			"Scattering Angle of Recoil Electron", 1000, 0., 6.);
		h_theta_thrown_secondary[i] = new TH1F(Form("theta_thrown_secondary_%d",i), 
			"Scattering Angle of Secondary Photon", 1000, 0., 6.);
		
		h_phi12_thrown[i] = new TH1F(Form("phi12_thrown_%d",i), "#phi_{12}; [#circ]", 3600, 0., 360.);
		h_phi13_thrown[i] = new TH1F(Form("phi13_thrown_%d",i), "#phi_{13}; [#circ]", 3600, 0., 360.);
		h_phi23_thrown[i] = new TH1F(Form("phi23_thrown_%d",i), "#phi_{23}; [#circ]", 3600, 0., 360.);
	}
	
	h_deltaE_vs_deltaK[0] = new TH2F("deltaE_vs_deltaK", 
		"#DeltaE vs. #DeltaK; E_{Compton} - E_{#gamma} [GeV]; E_{1} + E_{2} - E_{#gamma} [GeV]", 
		500, -4., 4., 500, -4., 4.);
	h_deltaE_vs_deltaK[1] = new TH2F("deltaE_vs_deltaK_born_sv", 
		"#DeltaE vs. #DeltaK; E_{Compton} - E_{#gamma} [GeV]; E_{1} + E_{2} - E_{#gamma} [GeV]", 
		500, -4., 4., 500, -4., 4.);
	h_deltaE_vs_deltaK[2] = new TH2F("deltaE_vs_deltaK_double", 
		"#DeltaE vs. #DeltaK; E_{Compton} - E_{#gamma} [GeV]; E_{1} + E_{2} - E_{#gamma} [GeV]", 
		500, -4., 4., 500, -4., 4.);
	
	h_deltaE[0] = new TH1F("deltaE", 
		"#DeltaE; E_{1} + E_{2} - E_{#gamma} [GeV]", 2000, -4.0, 4.0);
	h_deltaE[1] = new TH1F("deltaE_born_sv", 
		"#DeltaE; E_{1} + E_{2} - E_{#gamma} [GeV]", 2000, -4.0, 4.0);
	h_deltaE[2] = new TH1F("deltaE_double", 
		"#DeltaE; E_{1} + E_{2} - E_{#gamma} [GeV]", 2000, -4.0, 4.0);
	
	h_deltaPhi[0] = new TH1F("deltaPhi", 
		"#Delta#phi; | #phi_{1} - #phi_{2} | [#circ]", 3600, 0.0, 360.0);
	h_deltaPhi[1] = new TH1F("deltaPhi_born_sv", 
		"#Delta#phi; | #phi_{1} - #phi_{2} | [#circ]", 3600, 0.0, 360.0);
	h_deltaPhi[2] = new TH1F("deltaPhi_double", 
		"#Delta#phi; | #phi_{1} - #phi_{2} | [#circ]", 3600, 0.0, 360.0);
	
	h_deltaK[0] = new TH1F("deltaK", 
		"#DeltaK; E_{Compton} - E_{#gamma} [GeV]", 2000, -4.0, 4.0);
	h_deltaK[1] = new TH1F("deltaK_born_sv", 
		"#DeltaK; E_{Compton} - E_{#gamma} [GeV]", 2000, -4.0, 4.0);
	h_deltaK[2] = new TH1F("deltaK_double", 
		"#DeltaK; E_{Compton} - E_{#gamma} [GeV]", 2000, -4.0, 4.0);
	
	h_deltax_photon[0] = new TH1F("deltax_photon", 
		"x_{thrown} - x_{rec}; [cm]", 2000, -100., 100.);
	h_deltay_photon[0] = new TH1F("deltay_photon", 
		"y_{thrown} - y_{rec}; [cm]", 2000, -100., 100.);
	h_deltax_electron[0] = new TH1F("deltax_electron", 
		"x_{thrown} - x_{rec}; [cm]", 2000, -100., 100.);
	h_deltay_electron[0] = new TH1F("deltay_electron", 
		"y_{thrown} - y_{rec}; [cm]", 2000, -100., 100.);
	h_deltaE_photon[0] = new TH1F("deltaE_photon", 
		"E_{thrown} - E_{rec}; [GeV]", 1000, -5., 5.);
	h_deltaE_electron[0] = new TH1F("deltaE_electron", 
		"E_{thrown} - E_{rec}; [GeV]", 1000, -5., 5.);
	
	h_deltax_photon[1] = new TH1F("deltax_photon_bad", 
		"x_{thrown} - x_{rec}; [cm]", 2000, -100., 100.);
	h_deltay_photon[1] = new TH1F("deltay_photon_bad", 
		"y_{thrown} - y_{rec}; [cm]", 2000, -100., 100.);
	h_deltax_electron[1] = new TH1F("deltax_electron_bad", 
		"x_{thrown} - x_{rec}; [cm]", 2000, -100., 100.);
	h_deltay_electron[1] = new TH1F("deltay_electron_bad", 
		"y_{thrown} - y_{rec}; [cm]", 2000, -100., 100.);
	h_deltaE_photon[1] = new TH1F("deltaE_photon_bad", 
		"E_{thrown} - E_{rec}; [GeV]", 1000, -5., 5.);
	h_deltaE_electron[1] = new TH1F("deltaE_electron_bad", 
		"E_{thrown} - E_{rec}; [GeV]", 1000, -5., 5.);
	
	//-------------------------------------------------------------------------------------------//
	
	h_vertex     = new TH1F("vertex", 
		"Vertex Z Position (unweighted)", 1000, 0., 100.);
	h_vertex_cut = new TH1F("vertex_cut", 
		"Vertex Z Position (unweighted, accepted events)", 1000, 0., 100.);
	
	h_vertex_accepted     = new TH1F("vertex_accepted", 
		"Vertex Z Position (weighted)", 1000, 0., 100.);
	
	h_vertex_weight     = new TH1F("vertex_weight", 
		"Event Weight", 1000, 0., 2.);
	h_vertex_weight_cut = new TH1F("vertex_weight_cut", 
		"Event Weight (accepted events)", 1000, 0., 2.);
	
	h_fcal_rf_dt  = new TH1F("fcal_rf_dt", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	h_ccal_rf_dt  = new TH1F("ccal_rf_dt", "t_{CCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	h_beam_rf_dt  = new TH1F("beam_rf_dt", "t_{Beam} - t_{RF}; [ns]", 2000, -100., 100.);
	h_beam_rf_dt_cut = new TH1F("beam_rf_dt_cut", "t_{Beam} - t_{RF} (Compton Cuts); [ns]", 
		2000, -100., 100.);
	
	h_beam        = new TH1F("beam", "Is there a beam photon?", 2, -0.5, 1.5);
	h_tagh_flux   = new TH1F("tagh_flux", "TAGH Flux", 274, 0.5, 274.5);
	h_tagm_flux   = new TH1F("tagm_flux", "TAGM Flux", 102, 0.5, 102.5);
	h_beam_energy = new TH1F( "beam_energy", "Photon Beam Energy", 12000, 0., 12. );
	
	h_nccal = new TH1F("nccal", "number of CCAL Showers", 10, -0.5, 9.5);
	h_nfcal = new TH1F("nfcal", "number of FCAL Showers", 10, -0.5, 9.5);
	h_n_ccal_fcal     = new TH1F("n_ccal_fcal",     
		"Is there a shower in both CCAL+FCAL", 2, -0.5, 1.5);
	h_n_ccal_fcal_cut = new TH1F("n_ccal_fcal_cut", 
		"Is there a shower in both CCAL+FCAL", 2, -0.5, 1.5);
	
	h_ccalE = new TH2F( "ccalE", "Energy of CCAL Shower; E_{Beam} [GeV]; E_{CCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	h_fcalE = new TH2F( "fcalE", "Energy of FCAL Shower; E_{Beam} [GeV]; E_{FCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	h_ccalE_cut = new TH2F( "ccalE_cut", "Energy of CCAL Shower; E_{Beam} [GeV]; E_{CCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	h_fcalE_cut = new TH2F( "fcalE_cut", "Energy of FCAL Shower; E_{Beam} [GeV]; E_{FCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	
	//====================================================================================//
	
	TDirectory *dir_onelayercut = new TDirectoryFile("onelayercut", "onelayercut");
	dir_onelayercut->cd();
	
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
		274, 0.5, 274.5, 2000, -4.0, 4.0);
	h_deltaK_tagm = new TH2F("deltaK_tagm", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0);
	h_deltaK_shifted_tagh = new TH2F("deltaK_shifted_tagh", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0);
	h_deltaK_shifted_tagm = new TH2F("deltaK_shifted_tagm", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0);
	h_deltaK_smeared_tagh = new TH2F("deltaK_smeared_tagh", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0);
	h_deltaK_smeared_tagm = new TH2F("deltaK_smeared_tagm", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0);
	
	h_deltaK_tagh_cut = new TH2F("deltaK_tagh_cut", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0);
	h_deltaK_tagm_cut = new TH2F("deltaK_tagm_cut", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0);
	h_deltaK_shifted_tagh_cut = new TH2F("deltaK_shifted_tagh_cut", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0);
	h_deltaK_shifted_tagm_cut = new TH2F("deltaK_shifted_tagm_cut", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0);
	h_deltaK_smeared_tagh_cut = new TH2F("deltaK_smeared_tagh_cut", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0);
	h_deltaK_smeared_tagm_cut = new TH2F("deltaK_smeared_tagm_cut", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0);
	
	h_deltaK_tagh_cut_main = new TH2F("deltaK_tagh_cut_main", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0);
	h_deltaK_tagm_cut_main = new TH2F("deltaK_tagm_cut_main", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0);
	
	h_fcal_xy = new TH2F("fcal_xy", 
		"FCAL Y vs. X; x_{FCAL} [cm]; y_{FCAL} [cm]", 500, -100., 100., 500, -100., 100.);
	h_ccal_xy = new TH2F("ccal_xy", 
		"CCAL Y vs. X; x_{CCAL} [cm]; y_{CCAL} [cm]", 500,  -15.,  15., 500,  -15.,  15.);
	
	dir_onelayercut->cd("../");
	
	//====================================================================================//
	
	dir_reweight->cd("../");
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_compton_simulation_detailed::brun(JEventLoop *eventLoop, int32_t runnumber)
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
		
		m_target_length  = 29.535;
		m_target_density = 0.1217;
		m_atten          = 0.00821;
		
		m_fcalX_new = 0.408;
		m_fcalY_new = 0.027;
		
		m_ccalX_new = 0.108;
		m_ccalY_new = 0.130;
	}
	
	fcal_correction.SetXYZ(m_fcalX_new-m_fcalX, m_fcalY_new-m_fcalY, 0.);
	ccal_correction.SetXYZ(m_ccalX_new-m_ccalX, m_ccalY_new-m_ccalY, 0.);
	
	set_cuts(runnumber);
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_compton_simulation_detailed::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
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
	
	double vertex_z      = mc_thrown[0]->position().Z();
	double vertex_weight = get_event_weight(vertex_z);
	
	japp->RootFillLock(this);  // Acquire root lock
	
	h_vertex->Fill(vertex_z);
	h_vertex_weight->Fill(vertex_weight);
	
	// apply accept-reject filter based on vertex-event weight:
	
	if(vertex_weight < rndm_gen->Uniform()) {
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	
	h_vertex_accepted->Fill(vertex_z);
	
	
	//------------------------------------------------------------------------------//
	
	
	bool is_double_compton = false;
	if(mc_thrown[0]->energy() > 0.) is_double_compton = true;
	
	// plot angle of photon and electron:
	
	double theta_photon    = mc_thrown[1]->momentum().Theta() * (180./TMath::Pi());
	double theta_electron  = mc_thrown[2]->momentum().Theta() * (180./TMath::Pi());
	double theta_secondary = mc_thrown[0]->momentum().Theta() * (180./TMath::Pi());
	double phi_photon      = mc_thrown[1]->momentum().Phi()   * (180./TMath::Pi());
	double phi_electron    = mc_thrown[2]->momentum().Phi()   * (180./TMath::Pi());
	double phi_secondary   = mc_thrown[0]->momentum().Phi()   * (180./TMath::Pi());
	
	if(!is_double_compton) {
		h_theta_thrown_photon[0]->Fill(theta_photon);
		h_theta_thrown_electron[0]->Fill(theta_electron);
		h_theta_thrown_secondary[0]->Fill(theta_secondary);
		
		h_phi12_thrown[0]->Fill(fabs(phi_photon-phi_electron));
		h_phi13_thrown[0]->Fill(fabs(phi_photon-phi_secondary));
		h_phi23_thrown[0]->Fill(fabs(phi_electron-phi_secondary));
		
	} else {
		h_theta_thrown_photon[1]->Fill(theta_photon);
		h_theta_thrown_electron[1]->Fill(theta_electron);
		h_theta_thrown_secondary[1]->Fill(theta_secondary);
		
		h_phi12_thrown[1]->Fill(fabs(phi_photon-phi_electron));
		h_phi13_thrown[1]->Fill(fabs(phi_photon-phi_secondary));
		h_phi23_thrown[1]->Fill(fabs(phi_electron-phi_secondary));
	}
	
	
	
	
	
	
	
	//------------------------------------------------------------------------------//
	
	
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
	
	for(vector< const DFCALShower* >::const_iterator show1 = good_fcal_showers.begin(); 
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
				else if((-(2.004 + 3.*4.008)<=brfdt && brfdt<=-(2.004 + 2.*4.008))
					||((2.004 + 2.*4.008)<=brfdt && brfdt<=(2.004 + 3.*4.008)))
					bunch_val = 0;
				else 
					continue;
				
				if(eb < BEAM_min_energy_cut) continue;
				
				double ecomp1 =  1. / ((1./eb) + (1./m_e)*(1.-cos(theta1)));
				double ecomp2 =  1. / ((1./eb) + (1./m_e)*(1.-cos(theta2)));
				double deltaE = (e1     + e2    ) - (eb + m_e);
				double deltaK = (ecomp1 + ecomp2) - (eb + m_e);
				
				
				if(fabs(deltaPhi-180.) < 60. && fabs(deltaE)<1.0 && fabs(deltaK)<1.0) {
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
				
				loc_Cand.vz           = vertex_z;
				
				loc_Cand.eb           = eb;
				loc_Cand.brfdt        = brfdt;
				loc_Cand.tag_counter  = (*gam)->dCounter;
				
				DetectorSystem_t tag_sys = (*gam)->dSystem;
				if(tag_sys==SYS_TAGH)      loc_Cand.tag_sys = 0;
				else if(tag_sys==SYS_TAGM) loc_Cand.tag_sys = 1;
				
				candidates.push_back(loc_Cand);
				
			} // end DBeamPhoton loop
		} // end DCCALShower loop
	} // end DFCALShower loop
	
	
	fill_histograms(candidates, vertex, is_double_compton, mc_thrown);
	
	japp->RootFillUnLock(this);  // Release root lock
	
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_compton_simulation_detailed::erun(void)
{
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_compton_simulation_detailed::fini(void)
{
	
	return NOERROR;
}


int JEventProcessor_compton_simulation_detailed::fcal_fiducial_cut(DVector3 pos, DVector3 vertex)
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


int JEventProcessor_compton_simulation_detailed::ccal_fiducial_cut(DVector3 pos, DVector3 vertex)
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



void JEventProcessor_compton_simulation_detailed::fill_histograms(
	vector<ComptonCandidate_t> Comp_Candidates, DVector3 vertex, bool is_double, 
	vector<const DMCThrown*> mc_thrown) {
	
	int n_candidates = static_cast<int>(Comp_Candidates.size());
	
	for(int ic = 0; ic < n_candidates; ic++) {
		
		ComptonCandidate_t loc_Cand = Comp_Candidates[ic];
		
		//-------------------------------------------//
		
		int bunch_val       = loc_Cand.bunch_val;
		
		double eb           = loc_Cand.eb;
		double brfdt        = loc_Cand.brfdt;
		int tag_sys         = loc_Cand.tag_sys;
		int tag_counter     = loc_Cand.tag_counter;
		
		double deltaPhi     = loc_Cand.deltaPhi;
		double deltaE       = loc_Cand.deltaE;
		double deltaK       = loc_Cand.deltaK;
		
		double e1           = loc_Cand.e1;
		double x1           = loc_Cand.x1;
		double y1           = loc_Cand.y1;
		double z1           = loc_Cand.z1;
		double e2           = loc_Cand.e2;
		double x2           = loc_Cand.x2;
		double y2           = loc_Cand.y2;
		double z2           = loc_Cand.z2;
		
		//--------------     Cuts      --------------//
		
		double deltaE_smeared,   deltaPhi_smeared,   deltaK_smeared;
		double deltaK_shifted;
		
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
		if(fabs(deltaE_smeared   - deltaE_mu_data  ) < 5.0*deltaE_sig_data) e_cut = 1;
		
		int p_cut = 0;
		if(bfield_val) {
			if(fabs(deltaPhi - 180.) < 50.) {
				p_cut = 1;
			}
		} else {
			if(fabs(deltaPhi - deltaPhi_mu_mc) < 5.0*deltaPhi_sig_mc) p_cut = 1;
		}
		
		int k_cut = 0;
		if(fabs(deltaK - deltaK_mu_data) < 5.0*deltaK_sig_data) k_cut   = 1;
		
		int k_cut_smeared = 0;
		if(fabs(deltaK_smeared - deltaK_mu_data) < 5.0*deltaK_sig_data) 
			k_cut_smeared = 1;
		
		int k_cut_shifted = 0;
		if(fabs(deltaK_shifted - deltaK_mu_data) < 5.0*deltaK_sig_data) 
			k_cut_shifted = 1;
		
		//-------------------------------------------//
		
		double fill_weight;
		if(bunch_val) fill_weight =  1.0;
		else          fill_weight = -0.5;
		
		//----------------------------------------------------------------------------------//
		
		h_deltaE_vs_deltaK[0]->Fill(deltaK, deltaE, fill_weight);
		h_deltaE[0]->Fill(deltaE, fill_weight);
		h_deltaPhi[0]->Fill(deltaPhi, fill_weight);
		h_deltaK[0]->Fill(deltaK, fill_weight);
		if(!is_double) {
			h_deltaE_vs_deltaK[1]->Fill(deltaK, deltaE, fill_weight);
			h_deltaE[1]->Fill(deltaE, fill_weight);
			h_deltaPhi[1]->Fill(deltaPhi, fill_weight);
			h_deltaK[2]->Fill(deltaK, fill_weight);
		} else {
			h_deltaE_vs_deltaK[2]->Fill(deltaK, deltaE, fill_weight);
			h_deltaE[2]->Fill(deltaE, fill_weight);
			h_deltaPhi[2]->Fill(deltaPhi, fill_weight);
			h_deltaK[2]->Fill(deltaK, fill_weight);
		}
		
		// calculate position of thrown electron at fcal/ccal:
		
		DVector3 thrown_pos = mc_thrown[0]->position();
		
		double x_photon_1 = thrown_pos.X() + 
			(z1-thrown_pos.Z()) * tan(mc_thrown[1]->momentum().Theta()) * 
			cos(mc_thrown[1]->momentum().Phi());
		double y_photon_1 = thrown_pos.Y() + 
			(z1-thrown_pos.Z()) * tan(mc_thrown[1]->momentum().Theta()) * 
			sin(mc_thrown[1]->momentum().Phi());
		double x_electron_1 = thrown_pos.X() + 
			(z1-thrown_pos.Z()) * tan(mc_thrown[2]->momentum().Theta()) * 
			cos(mc_thrown[1]->momentum().Phi());
		double y_electron_1 = thrown_pos.Y() + 
			(z1-thrown_pos.Z()) * tan(mc_thrown[2]->momentum().Theta()) * 
			sin(mc_thrown[1]->momentum().Phi());
		
		double x_photon_2 = thrown_pos.X() + 
			(z2-thrown_pos.Z()) * tan(mc_thrown[1]->momentum().Theta()) * 
			cos(mc_thrown[1]->momentum().Phi());
		double y_photon_2 = thrown_pos.Y() + 
			(z2-thrown_pos.Z()) * tan(mc_thrown[1]->momentum().Theta()) * 
			sin(mc_thrown[1]->momentum().Phi());
		double x_electron_2 = thrown_pos.X() + 
			(z2-thrown_pos.Z()) * tan(mc_thrown[2]->momentum().Theta()) * 
			cos(mc_thrown[1]->momentum().Phi());
		double y_electron_2 = thrown_pos.Y() + 
			(z2-thrown_pos.Z()) * tan(mc_thrown[2]->momentum().Theta()) * 
			sin(mc_thrown[1]->momentum().Phi());
		
		double deltax_photon = 0., deltay_photon = 0., deltax_electron = 0., deltay_electron;
		double deltaE_photon = 0., deltaE_electron = 0.;
		if(mc_thrown[1]->momentum().Theta() > mc_thrown[2]->momentum().Theta()) {
			deltax_photon   = (x_photon_1-vertex.X()) - x1;
			deltay_photon   = (y_photon_1-vertex.Y()) - y1;
			deltax_electron = (x_electron_2-vertex.X()) - x2;
			deltay_electron = (y_electron_2-vertex.Y()) - y2;
			deltaE_photon   = mc_thrown[1]->energy() - e1;
			deltaE_electron = mc_thrown[2]->energy() - e2;
		} else {
			deltax_photon   = (x_photon_2-vertex.X()) - x1;
			deltay_photon   = (y_photon_2-vertex.Y()) - y1;
			deltax_electron = (x_electron_1-vertex.X()) - x2;
			deltay_electron = (y_electron_1-vertex.Y()) - y2;
			deltaE_photon   = mc_thrown[1]->energy() - e2;
			deltaE_electron = mc_thrown[2]->energy() - e1;
		}
		
		h_deltax_photon[0]->Fill(deltax_photon);
		h_deltay_photon[0]->Fill(deltay_photon);
		h_deltax_electron[0]->Fill(deltax_electron);
		h_deltay_electron[0]->Fill(deltay_electron);
		h_deltaE_photon[0]->Fill(deltaE_photon);
		h_deltaE_electron[0]->Fill(deltaE_electron);
		
		if(deltaE<-0.5) {
			h_deltax_photon[1]->Fill(deltax_photon);
			h_deltay_photon[1]->Fill(deltay_photon);
			h_deltax_electron[1]->Fill(deltax_electron);
			h_deltay_electron[1]->Fill(deltay_electron);
			h_deltaE_photon[1]->Fill(deltaE_photon);
			h_deltaE_electron[1]->Fill(deltaE_electron);
		}
		
		//----------------------------------------------------------------------------------//
		
		if(tag_sys==0) {
			h_deltaE_tagh->Fill(tag_counter, deltaE, fill_weight);
			h_deltaE_smeared_tagh->Fill(tag_counter, deltaE_smeared, fill_weight);
			
			if(e_cut) {
				h_deltaPhi_tagh->Fill(tag_counter, deltaPhi, fill_weight);
				h_deltaPhi_smeared_tagh->Fill(tag_counter, deltaPhi_smeared, 
					fill_weight);
				
				if(p_cut) {
					h_deltaK_tagh->Fill(tag_counter, deltaK, fill_weight);
					h_deltaK_shifted_tagh->Fill(tag_counter, deltaK_shifted, 
						fill_weight);
					h_deltaK_smeared_tagh->Fill(tag_counter, deltaK_smeared, 
						fill_weight);
					
					if(k_cut) {
						h_deltaK_tagh_cut->Fill(tag_counter, 
							deltaK, fill_weight);
						if(bunch_val) { 
							h_deltaK_tagh_cut_main->Fill(tag_counter,
								deltaK);
						}
						h_beam_rf_dt_cut->Fill(brfdt);
					} 
					if(k_cut_shifted) {
						h_deltaK_shifted_tagh_cut->Fill(tag_counter, 
							deltaK_shifted, fill_weight);
					}
					if(k_cut_smeared) {
						h_deltaK_smeared_tagh_cut->Fill(tag_counter, 
							deltaK_smeared, fill_weight);
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
					h_deltaK_shifted_tagm->Fill(tag_counter, deltaK_shifted, 
						fill_weight);
					h_deltaK_smeared_tagm->Fill(tag_counter, deltaK_smeared, 
						fill_weight);
					
					if(k_cut) {
						h_deltaK_tagm_cut->Fill(tag_counter, 
							deltaK, fill_weight);
						if(bunch_val) { 
							h_deltaK_tagm_cut_main->Fill(tag_counter,
								deltaK);
						}
						h_beam_rf_dt_cut->Fill(brfdt);
					} 
					if(k_cut_shifted) {
						h_deltaK_shifted_tagm_cut->Fill(tag_counter, 
							deltaK_shifted, fill_weight);
					}
					if(k_cut_smeared) {
						h_deltaK_smeared_tagm_cut->Fill(tag_counter, 
							deltaK_smeared, fill_weight);
					}
				}
			}
		}
		
		if(e_cut && p_cut && k_cut) {
			h_fcal_xy->Fill(x1, y1);
			h_ccal_xy->Fill(x2, y2);
		}
	}
	
	return;
}



void JEventProcessor_compton_simulation_detailed::set_cuts(int32_t runnumber)
{
	if(runnumber > 60000 && runnumber < 61355) {
		
		// Phase I, Be Target
		
		deltaE_mu_p0_mc    = -3.64447e-02;
		deltaE_mu_p1_mc    = -1.13010e-02;
		deltaE_mu_p2_mc    =  9.76324e-04;
		deltaE_mu_p3_mc    =  0.;
		
		deltaE_sig_p0_mc   =  9.55909e-03;
		deltaE_sig_p1_mc   =  3.89604e-02;
		deltaE_sig_p2_mc   =  0.;
		//------------------------------//
		deltaE_mu_p0_data  = -3.92986e-03;
		deltaE_mu_p1_data  =  5.16391e-03;
		deltaE_mu_p2_data  = -2.87647e-04;
		deltaE_mu_p3_data  =  0.;
		
		deltaE_sig_p0_data =  8.14994e-03;
		deltaE_sig_p1_data =  4.71470e-02;
		deltaE_sig_p2_data =  0.;
		
		
		deltaPhi_mu_p0_mc    =  1.80081e+02;
		deltaPhi_mu_p1_mc    = -1.14527e-01;
		deltaPhi_mu_p2_mc    =  1.14052e-02;
		deltaPhi_mu_p3_mc    = -4.68967e-04;
		
		deltaPhi_sig_p0_mc   =  1.16289e+01;
		deltaPhi_sig_p1_mc   = -1.67600e+00;
		deltaPhi_sig_p2_mc   =  1.62980e-01;
		deltaPhi_sig_p3_mc   = -5.43148e-03;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.78360e+02;
		deltaPhi_mu_p1_data  =  5.30878e-01;
		deltaPhi_mu_p2_data  = -6.44321e-02;
		deltaPhi_mu_p3_data  =  2.49178e-03;
		
		deltaPhi_sig_p0_data =  1.22866e+01;
		deltaPhi_sig_p1_data = -1.83819e+00;
		deltaPhi_sig_p2_data =  1.75185e-01;
		deltaPhi_sig_p3_data = -5.65782e-03;
		
		
		deltaK_mu_p0_mc    =  1.02951e-01;
		deltaK_mu_p1_mc    = -3.08838e-02;
		deltaK_mu_p2_mc    =  2.19226e-03;
		deltaK_mu_p3_mc    = -8.52794e-05;
		
		deltaK_sig_p0_mc   =  2.31343e-02;
		deltaK_sig_p1_mc   =  1.87014e-02;
		deltaK_sig_p2_mc   =  5.84751e-04;
		deltaK_sig_p3_mc   = -4.64438e-05;
		//--------------------------------//
		deltaK_mu_p0_data  =  1.39527e-02;
		deltaK_mu_p1_data  =  1.04239e-02;
		deltaK_mu_p2_data  = -3.45410e-03;
		deltaK_mu_p3_data  =  1.47912e-04;
		
		deltaK_sig_p0_data = -1.41305e-02;
		deltaK_sig_p1_data =  3.20553e-02;
		deltaK_sig_p2_data = -8.37718e-04;
		deltaK_sig_p3_data =  1.15860e-05;
		
	} else if(runnumber>60000 && runnumber<61911) {
		
		// Phase I, He Target, 200nA
		
		deltaE_mu_p0_mc    = -5.22300e-02;
		deltaE_mu_p1_mc    = -7.73078e-03;
		deltaE_mu_p2_mc    =  7.99221e-04;
		deltaE_mu_p3_mc    =  0.;
		
		deltaE_sig_p0_mc   =  8.56900e-03;
		deltaE_sig_p1_mc   =  3.88883e-02;
		deltaE_sig_p2_mc   =  2.35416e-09;
		//------------------------------//
		deltaE_mu_p0_data  =  4.10435e-01;
		deltaE_mu_p1_data  = -4.20598e-02;
		deltaE_mu_p2_data  = -6.19015e-04;
		deltaE_mu_p3_data  =  0.;
		
		deltaE_sig_p0_data =  1.26554e-02;
		deltaE_sig_p1_data =  4.28404e-02;
		deltaE_sig_p2_data =  0.;
		
		
		deltaPhi_mu_p0_mc    =  1.79526e+02;
		deltaPhi_mu_p1_mc    =  9.69292e-02;
		deltaPhi_mu_p2_mc    = -1.29449e-02;
		deltaPhi_mu_p3_mc    =  4.38689e-04;
		
		deltaPhi_sig_p0_mc   =  1.18746e+01;
		deltaPhi_sig_p1_mc   = -1.87957e+00;
		deltaPhi_sig_p2_mc   =  1.90664e-01;
		deltaPhi_sig_p3_mc   = -6.56656e-03;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.79777e+02;
		deltaPhi_mu_p1_data  = -7.73534e-03;
		deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_data =  1.11023e+01;
		deltaPhi_sig_p1_data = -1.58498e+00;
		deltaPhi_sig_p2_data =  1.53958e-01;
		deltaPhi_sig_p3_data = -5.09573e-03;
		
		
		deltaK_mu_p0_mc    =  1.33908e-02;
		deltaK_mu_p1_mc    =  3.93957e-03;
		deltaK_mu_p2_mc    = -1.90372e-03;
		deltaK_mu_p3_mc    =  7.16575e-05;
		
		deltaK_sig_p0_mc   = -1.46091e-03;
		deltaK_sig_p1_mc   =  2.67449e-02;
		deltaK_sig_p2_mc   = -3.00888e-04;
		deltaK_sig_p3_mc   = -9.70243e-06;
		//--------------------------------//
		deltaK_mu_p0_data    = -4.23755e-02;
		deltaK_mu_p1_data    =  5.52975e-02;
		deltaK_mu_p2_data    = -9.11057e-03;
		deltaK_mu_p3_data    =  3.35119e-04;
		
		deltaK_sig_p0_data   = -3.52834e-02;
		deltaK_sig_p1_data   =  4.58032e-02;
		deltaK_sig_p2_data   = -2.87025e-03;
		deltaK_sig_p3_data   =  1.05750e-04;
		
	} else if(runnumber>60000 && runnumber<61940) {
		
		// Phase I, He Target, 50nA
		
		deltaE_mu_p0_mc    = -4.52241e-02;
		deltaE_mu_p1_mc    = -9.08425e-03;
		deltaE_mu_p2_mc    =  8.94119e-04;
		deltaE_mu_p3_mc    =  0.;
		
		deltaE_sig_p0_mc   =  8.69176e-03;
		deltaE_sig_p1_mc   =  3.88600e-02;
		deltaE_sig_p2_mc   =  2.15792e-09;
		//------------------------------//
		deltaE_mu_p0_data  =  3.03822e-02;
		deltaE_mu_p1_data  = -1.48003e-02;
		deltaE_mu_p2_data  =  3.31982e-04;
		deltaE_mu_p3_data  =  0.;
		
		deltaE_sig_p0_data =  1.43507e-02;
		deltaE_sig_p1_data =  3.61306e-02;
		deltaE_sig_p2_data =  6.05400e-07;
		
		
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
		
	} else if(runnumber>60000 && runnumber<61956) {
		
		// Phase I, He Target, 100nA
		
		deltaE_mu_p0_mc    = -4.52241e-02;
		deltaE_mu_p1_mc    = -9.08425e-03;
		deltaE_mu_p2_mc    =  8.94119e-04;
		deltaE_mu_p3_mc    =  0.;
		
		deltaE_sig_p0_mc   =  8.69176e-03;
		deltaE_sig_p1_mc   =  3.88600e-02;
		deltaE_sig_p2_mc   =  2.15792e-09;
		//------------------------------//
		deltaE_mu_p0_data  = -1.34178e-03;
		deltaE_mu_p1_data  = -5.79806e-03;
		deltaE_mu_p2_data  = -1.12057e-05;
		deltaE_mu_p3_data  =  0.;
		
		deltaE_sig_p0_data =  1.48973e-02;
		deltaE_sig_p1_data =  3.56447e-02;
		deltaE_sig_p2_data =  1.81881e-02;
		
		
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
		
	} else if(runnumber>80000 && runnumber<81400) {
		
		// Be Target, Phase-II
		
		deltaE_mu_p0_mc      = -3.31056e-02;
		deltaE_mu_p1_mc      = -1.29897e-02;
		deltaE_mu_p2_mc      =  1.12944e-03;
		deltaE_mu_p3_mc      =  0.;
		
		deltaE_sig_p0_mc     =  8.77889e-03;
		deltaE_sig_p1_mc     =  3.97123e-02;
		deltaE_sig_p2_mc     =  0.;
		//--------------------------------//
		deltaE_mu_p0_data    =  1.34256e+00;
		deltaE_mu_p1_data    = -5.06498e-01;
		deltaE_mu_p2_data    =  6.21370e-02;
		deltaE_mu_p3_data    = -2.54423e-03;
		
		deltaE_sig_p0_data   =  1.31333e-02;
		deltaE_sig_p1_data   =  3.91007e-02;
		deltaE_sig_p2_data   =  0.;
		
		
		deltaPhi_mu_p0_mc    =  1.79592e+02;
		deltaPhi_mu_p1_mc    =  1.74329e-01;
		deltaPhi_mu_p2_mc    = -3.36433e-02;
		deltaPhi_mu_p3_mc    =  1.58799e-03;
		
		deltaPhi_sig_p0_mc   =  1.41982e+01;
		deltaPhi_sig_p1_mc   = -2.59293e+00;
		deltaPhi_sig_p2_mc   =  2.72994e-01;
		deltaPhi_sig_p3_mc   = -9.92059e-03;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.79847e+02;
		deltaPhi_mu_p1_data  = -1.08261e-02;
		deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_data =  1.26581e+01;
		deltaPhi_sig_p1_data = -2.00145e+00;
		deltaPhi_sig_p2_data =  1.97636e-01;
		deltaPhi_sig_p3_data = -6.64585e-03;
		
		
		deltaK_mu_p0_mc      = -3.89233e-02;
		deltaK_mu_p1_mc      =  2.74883e-02;
		deltaK_mu_p2_mc      = -5.80026e-03;
		deltaK_mu_p3_mc      =  2.73530e-04;
		
		deltaK_sig_p0_mc     = -4.08550e-02;
		deltaK_sig_p1_mc     =  3.74660e-02;
		deltaK_sig_p2_mc     = -1.16770e-03;
		deltaK_sig_p3_mc     =  7.05775e-06;
		//--------------------------------//
		deltaK_mu_p0_data    = -2.55073e-02;
		deltaK_mu_p1_data    =  2.13147e-02;
		deltaK_mu_p2_data    = -4.24219e-03;
		deltaK_mu_p3_data    =  1.58222e-04;
		
		deltaK_sig_p0_data   = -7.45541e-02;
		deltaK_sig_p1_data   =  5.69673e-02;
		deltaK_sig_p2_data   = -3.75135e-03;
		deltaK_sig_p3_data   =  1.23297e-04;
		
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
		
	} else {
		
		// He Target, Field ON, Phase-II
		
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


double JEventProcessor_compton_simulation_detailed::get_event_weight(double vertex_z) {
	
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
