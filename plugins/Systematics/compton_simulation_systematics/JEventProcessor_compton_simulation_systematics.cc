// $Id$
//
//    File: JEventProcessor_compton_simulation_systematics.cc
// Created: Mon Nov 23 16:24:15 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_compton_simulation_systematics.h"

extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_compton_simulation_systematics());
}
} // "C"

JEventProcessor_compton_simulation_systematics::JEventProcessor_compton_simulation_systematics() {
	
	m_USE_REACTION_WEIGHT = 0;
	gPARMS->SetDefaultParameter("compton_simulation_systematics:USE_REACTION_WEIGHT", 
		m_USE_REACTION_WEIGHT);
	
	m_REACTION_CUT_WEIGHT = 1.e4;
	gPARMS->SetDefaultParameter("compton_simulation_systematics:REACTION_CUT_WEIGHT", 
		m_REACTION_CUT_WEIGHT);
}

//------------------
// init
//------------------
jerror_t JEventProcessor_compton_simulation_systematics::init(void)
{
	rndm_gen = new TRandom3(0);
	
	TDirectory *dir_systematics = new TDirectoryFile("compton_simulation_systematics", 
		"compton_simulation_systematics");
	dir_systematics->cd();
	
	h_beam           = new TH1F("beam", "Is there a beam photon?", 2, -0.5, 1.5);
	h_tagh_flux      = new TH1F("tagh_flux", "TAGH Flux", 274, 0.5, 274.5);
	h_tagm_flux      = new TH1F("tagm_flux", "TAGM Flux", 102, 0.5, 102.5);
	h_double_compton = new TH1F("double_compton", "Is e2 > 0?", 2, -0.5, 1.5);
	
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
	
	TDirectory *dir_fcalE = new TDirectoryFile("fcalE", "fcalE");
	dir_fcalE->cd();
	for(int ihist=0; ihist<20; ihist++) {
		
		int loc_cut_int = 5*ihist;
		double loc_cut  = 0.05 * (double)(ihist);
		
		h_deltaK_tagh_fcalE[ihist] = new TH2F(Form("deltaK_tagh_%02dfcalE", loc_cut_int), 
			Form("#DeltaK (E_{FCAL} > %.2f GeV)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_fcalE[ihist] = new TH2F(Form("deltaK_tagm_%02dfcalE", loc_cut_int), 
			Form("#DeltaK (E_{FCAL} > %.2f GeV)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_fcalE[ihist] = new TH2F(Form("deltaK_tagh_cut_%02dfcalE", loc_cut_int), 
			Form("#DeltaK (E_{FCAL} > %.2f GeV)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_fcalE[ihist] = new TH2F(Form("deltaK_tagm_cut_%02dfcalE", loc_cut_int), 
			Form("#DeltaK (E_{FCAL} > %.2f GeV)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	dir_fcalE->cd("../");
	
	//--------------------------------------------------------------------//
	
	TDirectory *dir_ccalE = new TDirectoryFile("ccalE", "ccalE");
	dir_ccalE->cd();
	for(int ihist=0; ihist<13; ihist++) {
		
		int loc_cut_int = 5*ihist;
		double loc_cut  = 0.5 * (double)(ihist);
		
		h_deltaK_tagh_ccalE[ihist] = new TH2F(Form("deltaK_tagh_%02dccalE", loc_cut_int), 
			Form("#DeltaK (E_{CCAL} > %.2f GeV)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_ccalE[ihist] = new TH2F(Form("deltaK_tagm_%02dccalE", loc_cut_int), 
			Form("#DeltaK (E_{CCAL} > %.2f GeV)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_ccalE[ihist] = new TH2F(Form("deltaK_tagh_cut_%02dccalE", loc_cut_int), 
			Form("#DeltaK (E_{CCAL} > %.2f GeV)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_ccalE[ihist] = new TH2F(Form("deltaK_tagm_cut_%02dccalE", loc_cut_int), 
			Form("#DeltaK (E_{CCAL} > %.2f GeV)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	dir_ccalE->cd("../");
	
	//--------------------------------------------------------------------//
	
	TDirectory *dir_deltaE = new TDirectoryFile("DeltaE", "DeltaE");
	dir_deltaE->cd();
	for(int ihist=0; ihist<16; ihist++) {
		
		double loc_cut  = cut_sigmas[ihist];
		int loc_cut_int = (int)(loc_cut*10.);
		
		h_deltaK_tagh_sigE[ihist] = new TH2F(Form("deltaK_tagh_%03dsigE", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #DeltaE Cut)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_sigE[ihist] = new TH2F(Form("deltaK_tagm_%03dsigE", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #DeltaE Cut)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_sigE[ihist] = new TH2F(Form("deltaK_tagh_cut_%03dsigE", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #DeltaE Cut)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_sigE[ihist] = new TH2F(Form("deltaK_tagm_cut_%03dsigE", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #DeltaE Cut)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	dir_deltaE->cd("../");
	
	//--------------------------------------------------------------------//
	
	TDirectory *dir_deltaPhi = new TDirectoryFile("DeltaPhi", "DeltaPhi");
	dir_deltaPhi->cd();
	for(int ihist=0; ihist<16; ihist++) {
		
		double loc_cut  = cut_sigmas[ihist];
		int loc_cut_int = (int)(loc_cut*10.);
		
		h_deltaK_tagh_sigPhi[ihist] = new TH2F(Form("deltaK_tagh_%03dsigPhi", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #Delta#phi Cut)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_sigPhi[ihist] = new TH2F(Form("deltaK_tagm_%03dsigPhi", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #Delta#phi Cut)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_sigPhi[ihist] = new TH2F(Form("deltaK_tagh_cut_%03dsigPhi", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #Delta#phi Cut)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_sigPhi[ihist] = new TH2F(Form("deltaK_tagm_cut_%03dsigPhi", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #Delta#phi Cut)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	dir_deltaPhi->cd("../");
	
	//--------------------------------------------------------------------//
	
	TDirectory *dir_deltaK = new TDirectoryFile("DeltaK", "DeltaK");
	dir_deltaK->cd();
	for(int ihist=0; ihist<16; ihist++) {
		
		double loc_cut  = cut_sigmas[ihist];
		int loc_cut_int = (int)(loc_cut*10.);
		
		h_deltaK_tagh_sigK[ihist] = new TH2F(Form("deltaK_tagh_%03dsigK", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #DeltaK Cut)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_sigK[ihist] = new TH2F(Form("deltaK_tagm_%03dsigK", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #DeltaK Cut)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	dir_deltaK->cd("../");
	
	//--------------------------------------------------------------------//
	
	TDirectory *dir_fcalT = new TDirectoryFile("fcalT", "fcalT");
	dir_fcalT->cd();
	for(int ihist=0; ihist<10; ihist++) {
		
		int loc_cut_int = 5*(ihist+1);
		double loc_cut  = 0.5 * (double)(ihist+1);
		
		h_deltaK_tagh_fcalT[ihist] = new TH2F(Form("deltaK_tagh_%02dfcalT", loc_cut_int), 
			Form("#DeltaK (|t_{FCAL}-t_{RF}| < %.1f ns)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_fcalT[ihist] = new TH2F(Form("deltaK_tagm_%02dfcalT", loc_cut_int), 
			Form("#DeltaK (|t_{FCAL}-t_{RF}| < %.1f ns)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_fcalT[ihist] = new TH2F(Form("deltaK_tagh_cut_%02dfcalT", loc_cut_int), 
			Form("#DeltaK (|t_{FCAL}-t_{RF}| < %.1f ns)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_fcalT[ihist] = new TH2F(Form("deltaK_tagm_cut_%02dfcalT", loc_cut_int), 
			Form("#DeltaK (|t_{FCAL}-t_{RF}| < %.1f ns)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	dir_fcalT->cd("../");
	
	//--------------------------------------------------------------------//
	
	TDirectory *dir_ccalT = new TDirectoryFile("ccalT", "ccalT");
	dir_ccalT->cd();
	for(int ihist=0; ihist<10; ihist++) {
		
		int loc_cut_int = 5*ihist;
		double loc_cut  = 0.5 * (double)(ihist);
		
		h_deltaK_tagh_ccalT[ihist] = new TH2F(Form("deltaK_tagh_%02dccalT", loc_cut_int), 
			Form("#DeltaK (|t_{CCAL}-t_{RF}| < %.1f ns)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_ccalT[ihist] = new TH2F(Form("deltaK_tagm_%02dccalT", loc_cut_int), 
			Form("#DeltaK (|t_{CCAL}-t_{RF}| < %.1f ns)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_ccalT[ihist] = new TH2F(Form("deltaK_tagh_cut_%02dccalT", loc_cut_int), 
			Form("#DeltaK (|t_{CCAL}-t_{RF}| < %.1f ns)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_ccalT[ihist] = new TH2F(Form("deltaK_tagm_cut_%02dccalT", loc_cut_int), 
			Form("#DeltaK (|t_{CCAL}-t_{RF}| < %.1f ns)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	dir_ccalT->cd("../");
	
	//--------------------------------------------------------------------//
	
	TDirectory *dir_fcal_fid = new TDirectoryFile("fcal_fid", "fcal_fid");
	dir_fcal_fid->cd();
	for(int icut = 0; icut < N_FID_CUTS; icut++) {
		
		h_deltaK_tagh_fcalfid[icut] = new TH2F(Form("deltaK_tagh_fcalfid_%02d",icut), 
			"#DeltaK", 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_fcalfid[icut] = new TH2F(Form("deltaK_tagm_fcalfid_%02d",icut), 
			"#DeltaK", 102, 0.5, 102.5, 2000, -8., 8.);
	}
	dir_fcal_fid->cd("../");
	
	//--------------------------------------------------------------------//
	
	TDirectory *dir_ccal_fid = new TDirectoryFile("ccal_fid", "ccal_fid");
	dir_ccal_fid->cd();
	for(int icut = 0; icut < N_FID_CUTS; icut++) {
		
		h_deltaK_tagh_ccalfid[icut] = new TH2F(Form("deltaK_tagh_ccalfid_%02d",icut), 
			"#DeltaK", 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_ccalfid[icut] = new TH2F(Form("deltaK_tagm_ccalfid_%02d",icut), 
			"#DeltaK", 102, 0.5, 102.5, 2000, -8., 8.);
	}
	dir_ccal_fid->cd( "../" );
	
	//--------------------------------------------------------------------//
	
	TDirectory *dir_fcal_phi = new TDirectoryFile( "fcal_phi", "fcal_phi" );
	dir_fcal_phi->cd();
	for(int ip = 0; ip < 8; ip++) {
		
		h_deltaK_tagh_fcal_phi[ip]     = new TH2F(Form("deltaK_tagh_fcal_phi_%d",ip), 
			Form("#DeltaK (FCAL Phi Sect. %d)",ip), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_fcal_phi[ip]     = new TH2F(Form("deltaK_tagm_fcal_phi_%d",ip), 
			Form("#DeltaK (FCAL Phi Sect. %d)",ip), 102, 0.5, 102.5, 2000, -8., 8.);
		h_deltaK_tagh_cut_fcal_phi[ip] = new TH2F(Form("deltaK_tagh_cut_fcal_phi_%d",ip), 
			Form("#DeltaK (FCAL Phi Sect. %d)",ip), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_fcal_phi[ip] = new TH2F(Form("deltaK_tagm_cut_fcal_phi_%d",ip), 
			Form("#DeltaK (FCAL Phi Sect. %d)",ip), 102, 0.5, 102.5, 2000, -8., 8.);
		h_xy_fcal_phi[ip] = new TH2F(Form("xy_fcal_phi_%d",ip), 
			Form("FCAL Y vs. X (Phi Sect. %d)",ip), 500, -100., 100., 500, -100., 100.);
	}
	dir_fcal_phi->cd("../");
	
	//--------------------------------------------------------------------//
	
	TDirectory *dir_ccal_phi = new TDirectoryFile("ccal_phi", "ccal_phi");
	dir_ccal_phi->cd();
	for(int ip = 0; ip < 8; ip++) {
		
		h_deltaK_tagh_ccal_phi[ip]     = new TH2F(Form("deltaK_tagh_ccal_phi_%d",ip), 
			Form("#DeltaK (CCAL Phi Sect. %d)",ip), 274, 0.5, 274.5, 2000, -8.,  8.);
		h_deltaK_tagm_ccal_phi[ip]     = new TH2F(Form("deltaK_tagm_ccal_phi_%d",ip), 
			Form("#DeltaK (CCAL Phi Sect. %d)",ip), 102, 0.5, 102.5, 2000, -8.,  8.);
		h_deltaK_tagh_cut_ccal_phi[ip] = new TH2F(Form("deltaK_tagh_cut_ccal_phi_%d",ip), 
			Form("#DeltaK (CCAL Phi Sect. %d)",ip), 274, 0.5, 274.5, 2000, -8.,  8.);
		h_deltaK_tagm_cut_ccal_phi[ip] = new TH2F(Form("deltaK_tagm_cut_ccal_phi_%d",ip), 
			Form("#DeltaK (CCAL Phi Sect. %d)",ip), 102, 0.5, 102.5, 2000, -8.,  8.);
		h_xy_ccal_phi[ip] = new TH2F(Form("xy_ccal_phi_%d",ip), 
			Form("CCAL Y vs. X (Phi Sect. %d)",ip), 500, -13., 13., 500, -13., 13.);
	}
	dir_ccal_phi->cd("../");
	
	//--------------------------------------------------------------------//
	
	TDirectory *dir_fcal_layer = new TDirectoryFile("fcal_layer", "fcal_layer");
	dir_fcal_layer->cd();
	
	for(int ip = 0; ip < 8; ip++) {
		
		h_deltaK_tagh_fcal_layer[ip]     = new TH2F(Form("deltaK_tagh_fcal_layer_%d",ip+1), 
			Form("#DeltaK (FCAL Layer %d)",ip+1), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_fcal_layer[ip]     = new TH2F(Form("deltaK_tagm_fcal_layer_%d",ip+1), 
			Form("#DeltaK (FCAL Layer %d)",ip+1), 102, 0.5, 102.5, 2000, -8., 8.);
		h_deltaK_tagh_cut_fcal_layer[ip] = new TH2F(Form("deltaK_tagh_cut_fcal_layer_%d",ip+1), 
			Form("#DeltaK (FCAL Layer %d)",ip+1), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_fcal_layer[ip] = new TH2F(Form("deltaK_tagm_cut_fcal_layer_%d",ip+1), 
			Form("#DeltaK (FCAL Layer %d)",ip+1), 102, 0.5, 102.5, 2000, -8., 8.);
		h_xy_fcal_layer[ip] = new TH2F(Form("xy_fcal_layer_%d",ip+1), 
			Form("FCAL Y vs. X (Layer %d)",ip+1), 500, -100., 100., 500, -100., 100.);
	}
	dir_fcal_layer->cd("../");
	
	//--------------------------------------------------------------------//
	
	TDirectory *dir_ccal_layer = new TDirectoryFile("ccal_layer", "ccal_layer");
	dir_ccal_layer->cd();
	for(int ip = 0; ip < 5; ip++) {
		
		h_deltaK_tagh_ccal_layer[ip]     = new TH2F(Form("deltaK_tagh_ccal_layer_%d",ip+1), 
			Form("#DeltaK (CCAL Layer %d)",ip+1), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_ccal_layer[ip]     = new TH2F(Form("deltaK_tagm_ccal_layer_%d",ip+1), 
			Form("#DeltaK (CCAL Layer %d)",ip+1), 102, 0.5, 102.5, 2000, -8., 8.);
		h_deltaK_tagh_cut_ccal_layer[ip] = new TH2F(Form("deltaK_tagh_cut_ccal_layer_%d",ip+1), 
			Form("#DeltaK (CCAL Layer %d)",ip+1), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_ccal_layer[ip] = new TH2F(Form("deltaK_tagm_cut_ccal_layer_%d",ip+1), 
			Form("#DeltaK (CCAL Layer %d)",ip+1), 102, 0.5, 102.5, 2000, -8., 8.);
		h_xy_ccal_layer[ip] = new TH2F(Form("xy_ccal_layer_%d",ip+1), 
			Form("CCAL Y vs. X (Layer %d)",ip+1), 500, -13., 13., 500, -13., 13.);
	}
	dir_ccal_layer->cd("../");
	
	//--------------------------------------------------------------------//
	
	h_deltaK_tagh_single = new TH2F("deltaK_tagh_single", 
		"#DeltaK (only one CCAL Shower, FCAL Shower, and DBeamPhoton); TAGH Counter; [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm_single = new TH2F("deltaK_tagm_single", 
		"#DeltaK (only one CCAL Shower, FCAL Shower, and DBeamPhoton); TAGM Counter; [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	
	dir_systematics->cd("../");
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_compton_simulation_systematics::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	
	DGeometry*   dgeom = NULL;
	DApplication* dapp = dynamic_cast<DApplication*>(eventLoop->GetJApplication());
	if(dapp)     dgeom = dapp->GetDGeometry(runnumber);
   	
	if(dgeom){
		dgeom->GetTargetZ(m_beamZ);
		dgeom->GetFCALPosition(m_fcalX, m_fcalY, m_fcalZ);
		dgeom->GetCCALPosition(m_ccalX, m_ccalY, m_ccalZ);
	} else{
		cerr << "No geometry accessbile to compton_simulation_reweight plugin." << endl;
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
jerror_t JEventProcessor_compton_simulation_systematics::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
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
	
	double e_extra_gam = mc_thrown[0]->energy();
	if(e_extra_gam > 0.) 
		h_double_compton->Fill(1);
	else {
		h_double_compton->Fill(0);
	}
	
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
	}
	
	if(locRFBunch->dNumParticleVotes < 2) {
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	
	//----------     Check FCAL-CCAL Pairs     ----------//
	
	vector<ComptonCandidate_t> candidates;
	
	for(vector<const DFCALShower*>::const_iterator show1 = fcal_showers.begin(); 
		show1 != fcal_showers.end(); show1++) {
		
		double e1     = (*show1)->getEnergy();
		DVector3 pos1 = (*show1)->getPosition_log() - vertex + fcal_correction;
		
		double t1     = (*show1)->getTime() - (pos1.Mag()/c);
		double phi1   = pos1.Phi() * (180. / TMath::Pi());
		double theta1 = pos1.Theta();
		
		int fcal_phi_sect;
		if(     -180. < phi1 && phi1 <= -135.) fcal_phi_sect = 0;
		else if(-135. < phi1 && phi1 <=  -90.) fcal_phi_sect = 1;
		else if( -90. < phi1 && phi1 <=  -45.) fcal_phi_sect = 2;
		else if( -45. < phi1 && phi1 <=    0.) fcal_phi_sect = 3;
		else if(   0. < phi1 && phi1 <=   45.) fcal_phi_sect = 4;
		else if(  45. < phi1 && phi1 <=   90.) fcal_phi_sect = 5;
		else if(  90. < phi1 && phi1 <=  135.) fcal_phi_sect = 6;
		else if( 135. < phi1 && phi1 <=  180.) fcal_phi_sect = 7;
		else fcal_phi_sect = 8;
		
		double fcal_face_x = vertex.X() + (pos1.X() * (m_fcalZ - vertex.Z())/pos1.Z());
		double fcal_face_y = vertex.Y() + (pos1.Y() * (m_fcalZ - vertex.Z())/pos1.Z());
		
		fcal_face_x -= m_fcalX_new;
		fcal_face_y -= m_fcalY_new;
		
		if(phase_val==1) {
			if((-32. < fcal_face_y && fcal_face_y < -20.) 
				&& (-8. < fcal_face_x && fcal_face_x < 4.)) continue;
		}
		
		int fcal_layer;
		if(     (-4.0157*(2.5) < fcal_face_x && fcal_face_x < 4.0157*(2.5)) 
			&& (-4.0157*(2.5) < fcal_face_y && fcal_face_y < 4.0157*(2.5))) 
			fcal_layer = 1;
		else if((-4.0157*(3.5) < fcal_face_x && fcal_face_x < 4.0157*(3.5)) 
			&& (-4.0157*(3.5) < fcal_face_y && fcal_face_y < 4.0157*(3.5))) 
			fcal_layer = 2;
		else if((-4.0157*(4.5) < fcal_face_x && fcal_face_x < 4.0157*(4.5)) 
			&& (-4.0157*(4.5) < fcal_face_y && fcal_face_y < 4.0157*(4.5))) 
			fcal_layer = 3;
		else if((-4.0157*(5.5) < fcal_face_x && fcal_face_x < 4.0157*(5.5)) 
			&& (-4.0157*(5.5) < fcal_face_y && fcal_face_y < 4.0157*(5.5))) 
			fcal_layer = 4;
		else if((-4.0157*(6.5) < fcal_face_x && fcal_face_x < 4.0157*(6.5)) 
			&& (-4.0157*(6.5) < fcal_face_y && fcal_face_y < 4.0157*(6.5))) 
			fcal_layer = 5;
		else if((-4.0157*(7.5) < fcal_face_x && fcal_face_x < 4.0157*(7.5)) 
			&& (-4.0157*(7.5) < fcal_face_y && fcal_face_y < 4.0157*(7.5))) 
			fcal_layer = 6;
		else if((-4.0157*(8.5) < fcal_face_x && fcal_face_x < 4.0157*(8.5)) 
			&& (-4.0157*(8.5) < fcal_face_y && fcal_face_y < 4.0157*(8.5))) 
			fcal_layer = 7;
		else fcal_layer = 8;
		
		for(vector<const DCCALShower*>::const_iterator show2 = ccal_showers.begin(); 
			show2 != ccal_showers.end(); show2++) {
			
			double e2 = (*show2)->E;
			DVector3 pos2((*show2)->x1, (*show2)->y1, (*show2)->z);
			pos2      = pos2 - vertex + ccal_correction;
			
			double t2     = (*show2)->time - (pos2.Mag()/c);
			double phi2   = pos2.Phi() * (180. / TMath::Pi());
			double theta2 = pos2.Theta();
			
			// calculate deltaPhi and deltaT:
			
			double deltaPhi = fabs(phi2 - phi1);
			double deltaT   = t2 - t1;
			
			// Get CCAL Layer and Phi Section:
			
			int ccal_phi_sect;
			if(     -180. < phi2 && phi2 <=-135.) ccal_phi_sect = 0;
			else if(-135. < phi2 && phi2 <= -90.) ccal_phi_sect = 1;
			else if( -90. < phi2 && phi2 <= -45.) ccal_phi_sect = 2;
			else if( -45. < phi2 && phi2 <=   0.) ccal_phi_sect = 3;
			else if(   0. < phi2 && phi2 <=  45.) ccal_phi_sect = 4;
			else if(  45. < phi2 && phi2 <=  90.) ccal_phi_sect = 5;
			else if(  90. < phi2 && phi2 <= 135.) ccal_phi_sect = 6;
			else if( 135. < phi2 && phi2 <= 180.) ccal_phi_sect = 7;
			else ccal_phi_sect = 8;
			
			double ccal_face_x = vertex.X() + (pos2.X() * (m_ccalZ - vertex.Z())/pos2.Z());
			double ccal_face_y = vertex.Y() + (pos2.Y() * (m_ccalZ - vertex.Z())/pos2.Z());
			
			ccal_face_x -= m_ccalX_new;
			ccal_face_y -= m_ccalY_new;
			
			int ccal_layer;
			if( (-2.09*(1.+1.) < ccal_face_x && ccal_face_x < 2.09*(1.+1.)) 
				&& (-2.09*(1.+1.) < ccal_face_y && ccal_face_y < 2.09*(1.+1.)))
				ccal_layer = 1;
			else if((-2.09*(1.+2.) < ccal_face_x && ccal_face_x < 2.09*(1.+2.)) 
				&& (-2.09*(1.+2.) < ccal_face_y && ccal_face_y < 2.09*(1.+2.))) 
				ccal_layer = 2;
			else if((-2.09*(1.+3.) < ccal_face_x && ccal_face_x < 2.09*(1.+3.)) 
				&& (-2.09*(1.+3.) < ccal_face_y && ccal_face_y < 2.09*(1.+3.))) 
				ccal_layer = 3;
			else if((-2.09*(1.+4.) < ccal_face_x && ccal_face_x < 2.09*(1.+4.)) 
				&& (-2.09*(1.+4.) < ccal_face_y && ccal_face_y < 2.09*(1.+4.))) 
				ccal_layer = 4;
			else ccal_layer = 5;
			
			// loop over beam photons:
			
			for(vector<const DBeamPhoton*>::const_iterator gam = beam_photons.begin();
				gam != beam_photons.end(); gam++) {
				
				double eb = (*gam)->lorentzMomentum().E();
				double tb = (*gam)->time();
				
				double brfdt = tb - rfTime;
				
				int bunch_val;
				
				if(fabs(brfdt) < 2.004)
					bunch_val = 1;
				else if((-(2.004 + 5.*4.008)<=brfdt && brfdt<=-(2.004 + 3.*4.008))
					||((2.004 + 3.*4.008)<=brfdt && brfdt<=(2.004 + 5.*4.008)))
					bunch_val = 0;
				else 
					continue;
				
				if(eb < BEAM_min_energy_cut) continue;
				
				double deltaE = (e1 + e2) - (eb + m_e);
				/*
				double ecomp1 =  1. / ((1./eb) + (1./m_e)*(1.-cos(theta1)));
				double ecomp2 =  1. / ((1./eb) + (1./m_e)*(1.-cos(theta2)));
				double deltaK = (ecomp1 + ecomp2) - (eb + m_e);
				*/
				double deltaK = m_e * sin(theta1+theta2) / 
					(sin(theta1) + sin(theta2) - sin(theta1+theta2));
				deltaK       -= eb;
				
				ComptonCandidate_t loc_Cand;
				
				loc_Cand.bunch_val     = bunch_val;
				
				loc_Cand.e1            = e1;
				loc_Cand.t1            = t1;
				loc_Cand.x1            = pos1.X();
				loc_Cand.y1            = pos1.Y();
				loc_Cand.z1            = pos1.Z();
				loc_Cand.e2            = e2;
				loc_Cand.t2            = t2;
				loc_Cand.x2            = pos2.X();
				loc_Cand.y2            = pos2.Y();
				loc_Cand.z2            = pos2.Z();
				
				loc_Cand.rfTime        = rfTime;
				
				loc_Cand.deltaPhi      = deltaPhi;
				loc_Cand.deltaT        = deltaT;
				loc_Cand.deltaE        = deltaE;
				loc_Cand.deltaK        = deltaK;
				
				loc_Cand.fcal_phi_sect = fcal_phi_sect;
				loc_Cand.ccal_phi_sect = ccal_phi_sect;
				
				loc_Cand.fcal_layer    = fcal_layer;
				loc_Cand.ccal_layer    = ccal_layer;
				
				loc_Cand.eb            = eb;
				loc_Cand.tag_counter   = (*gam)->dCounter;
				
				DetectorSystem_t sys   = (*gam)->dSystem;
				if(sys==SYS_TAGH)      loc_Cand.tag_sys = 0;
				else if(sys==SYS_TAGM) loc_Cand.tag_sys = 1;
				
				loc_Cand.event_weight  = reaction_weight;
				
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
jerror_t JEventProcessor_compton_simulation_systematics::erun(void)
{
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_compton_simulation_systematics::fini(void)
{
	
	return NOERROR;
}


int JEventProcessor_compton_simulation_systematics::fcal_fiducial_cut(DVector3 pos, DVector3 vertex)
{
	int fid_cut = 0;
	
	double fcal_inner_layer_cut = 2.5 * 4.0157;
	
	double fcal_face_x = vertex.X() + (pos.X() * (m_fcalZ - vertex.Z())/pos.Z());
	double fcal_face_y = vertex.Y() + (pos.Y() * (m_fcalZ - vertex.Z())/pos.Z());
	
	fcal_face_x -= m_fcalX;
	fcal_face_y -= m_fcalY;
	
	if((-1.*fcal_inner_layer_cut < fcal_face_x && fcal_face_x < fcal_inner_layer_cut)
		&& (-1.*fcal_inner_layer_cut < fcal_face_y 
		&& fcal_face_y < fcal_inner_layer_cut)) fid_cut = 1;
	
	/*
	if((-32.<fcal_face_y && fcal_face_y<-20.) && (-8.<fcal_face_x && fcal_face_x<4.))
			fid_cut = 1;
	*/
	
	return fid_cut;
}


int JEventProcessor_compton_simulation_systematics::ccal_fiducial_cut(DVector3 pos, DVector3 vertex)
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



void JEventProcessor_compton_simulation_systematics::fill_histograms(
	vector<ComptonCandidate_t> Comp_Candidates, DVector3 vertex) {
	
	int n_candidates = static_cast<int>(Comp_Candidates.size());
	
	int n_good_cands = 0;
	
	int loc_tag_sys, loc_tag_counter;
	double loc_deltaK;
	
	for(int ic = 0; ic < n_candidates; ic++) {
				
		ComptonCandidate_t loc_Cand = Comp_Candidates[ic];
				
		//-------------------------------------------//
		
		int bunch_val     = loc_Cand.bunch_val;
		
		double e1         = loc_Cand.e1;
		double x1         = loc_Cand.x1;
		double y1         = loc_Cand.y1;
		double z1         = loc_Cand.z1;
		
		double e2         = loc_Cand.e2;
		double x2         = loc_Cand.x2;
		double y2         = loc_Cand.y2;
		double z2         = loc_Cand.z2;
		
		double eb         = loc_Cand.eb;
		int tag_sys       = loc_Cand.tag_sys;
		int tag_counter   = loc_Cand.tag_counter;
			
		double deltaPhi   = loc_Cand.deltaPhi;
		double deltaE     = loc_Cand.deltaE;
		double deltaK     = loc_Cand.deltaK;
		
		int fcal_phi_sect = loc_Cand.fcal_phi_sect;
		int ccal_phi_sect = loc_Cand.ccal_phi_sect;
		
		int fcal_layer    = loc_Cand.fcal_layer;
		int ccal_layer    = loc_Cand.ccal_layer;
		
		double event_weight = loc_Cand.event_weight;
		
		//--------------     Cuts      --------------//
		
		double fcal_face_x = vertex.X() + (x1 * (m_fcalZ - vertex.Z())/z1);
		double fcal_face_y = vertex.Y() + (y1 * (m_fcalZ - vertex.Z())/z1);
		
		fcal_face_x -= m_fcalX_new;
		fcal_face_y -= m_fcalY_new;
		
		double loc_onelayer_cut = 2.5*4.0157;
		
		int fcal_fid_cut = 0;
		if((-1.*loc_onelayer_cut < fcal_face_x 
			&& fcal_face_x < loc_onelayer_cut) 
			&& (-1.*loc_onelayer_cut < fcal_face_y 
			&& fcal_face_y < loc_onelayer_cut)) 
			fcal_fid_cut = 1;
		
		double ccal_face_x = vertex.X() + (x2 * (m_ccalZ - vertex.Z())/z2);
		double ccal_face_y = vertex.Y() + (y2 * (m_ccalZ - vertex.Z())/z2);
		
		ccal_face_x -= m_ccalX_new;
		ccal_face_y -= m_ccalY_new;
		
		int ccal_fid_cut = 0;
		if((-4.18 < ccal_face_x && ccal_face_x < 4.18) 
			&& (-4.18 < ccal_face_y && ccal_face_y < 4.18)) 
			ccal_fid_cut = 1;
		if(ccal_face_x < -8.36 || ccal_face_x > 10.45 
			|| ccal_face_y < -10.45 || ccal_face_y > 10.45) 
			ccal_fid_cut = 1;
		
		
		double deltaE_smeared, deltaPhi_smeared, deltaK_smeared;
		
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
		
		double deltaK_mu_mc    = f_deltaK_mu_mc->Eval(eb);
		double deltaK_mu_data  = f_deltaK_mu_data->Eval(eb);
		double deltaK_sig_mc   = f_deltaK_sig_mc->Eval(eb);
		double deltaK_sig_data = f_deltaK_sig_data->Eval(eb);
		
		deltaK_smeared = deltaK + (deltaK_mu_data - deltaK_mu_mc);
		if(deltaK_sig_data > deltaK_sig_mc) {
			double loc_deltaK_smear = sqrt(pow(deltaK_sig_data,2.0)
				- pow(deltaK_sig_mc,2.0));
			deltaK_smeared += rndm_gen->Gaus(0.,loc_deltaK_smear);
		}
		
		int phi_cut[16], e_cut[16], k_cut[16];
		for(int icut=0; icut<16; icut++) {
			double loc_cut = cut_sigmas[icut];
			
			if(fabs(deltaPhi_smeared-deltaPhi_mu_data) < loc_cut*deltaPhi_sig_data) {
				phi_cut[icut] = 1;
			} else {
				phi_cut[icut] = 0;
			}
			
			if(fabs(deltaE_smeared-deltaE_mu_data) < loc_cut*deltaE_sig_data) {
				e_cut[icut] = 1;
			} else {
				e_cut[icut] = 0;
			}
			
			if(fabs(deltaK_smeared-deltaK_mu_data) < loc_cut*deltaK_sig_data) {
				k_cut[icut] = 1;
			} else {
				k_cut[icut] = 0;
			}
		}
		
		int e5_cut =   e_cut[6];
		int k5_cut =   k_cut[6];
		int p5_cut = phi_cut[6];
		
		//---------------------------------------------------------------------------//
		// RF Timing Cuts:
		
		double t1     = loc_Cand.t1;
		double t2     = loc_Cand.t2;
		double rfTime = loc_Cand.rfTime;
		
		int fcal_t_cut = 0, ccal_t_cut = 0;
		if(fabs(t1 - rfTime) < FCAL_RF_time_cut) fcal_t_cut = 1;
		if(fabs(t2 - rfTime) < CCAL_RF_time_cut) ccal_t_cut = 1;
		
		int fcal_t_cut_vec[10], ccal_t_cut_vec[10];
		for(int icut=0; icut<10; icut++) {
			double loc_cut = 0.5 * (double)(icut+1);
			
			if(fabs(t1-rfTime) < loc_cut) {
				fcal_t_cut_vec[icut] = 1;
			} else {
				fcal_t_cut_vec[icut] = 0;
			}
			
			if(fabs(t2-rfTime) < loc_cut) {
				ccal_t_cut_vec[icut] = 1;
			} else {
				ccal_t_cut_vec[icut] = 0;
			}
		}
		
		//---------------------------------------------------------------------------//
		
		
		double fill_weight;
		if(bunch_val) fill_weight =  1.0;
		else          fill_weight = -0.25;
		
		if(m_USE_REACTION_WEIGHT) {
			fill_weight *= event_weight;
		}
		
		if(e5_cut && p5_cut && !fcal_fid_cut && !ccal_fid_cut 
			&& e1 > FCAL_min_energy_cut && e2 > CCAL_min_energy_cut
			&& fcal_t_cut && ccal_t_cut) {
			
			n_good_cands++;
			loc_deltaK      = deltaK_smeared;
			loc_tag_sys     = tag_sys;
			loc_tag_counter = tag_counter;
			
			if(tag_sys==0) {
				h_deltaK_tagh_fcal_phi[fcal_phi_sect]->Fill(
					tag_counter, deltaK_smeared, fill_weight);
				h_deltaK_tagh_ccal_phi[ccal_phi_sect]->Fill(
					tag_counter, deltaK_smeared, fill_weight);
				h_deltaK_tagh_fcal_layer[fcal_layer-1]->Fill(
					tag_counter, deltaK_smeared, fill_weight);
				h_deltaK_tagh_ccal_layer[ccal_layer-1]->Fill(
					tag_counter, deltaK_smeared, fill_weight);
			} else { 
				h_deltaK_tagm_fcal_phi[fcal_phi_sect]->Fill(
					tag_counter, deltaK_smeared, fill_weight);
				h_deltaK_tagm_ccal_phi[ccal_phi_sect]->Fill(
					tag_counter, deltaK_smeared, fill_weight);
				h_deltaK_tagm_fcal_layer[fcal_layer-1]->Fill(
					tag_counter, deltaK_smeared, fill_weight);
				h_deltaK_tagm_ccal_layer[ccal_layer-1]->Fill(
					tag_counter, deltaK_smeared, fill_weight);
			}
			
			if(k5_cut) {
				
				if( tag_sys==0 ) {
					h_deltaK_tagh_cut_fcal_phi[fcal_phi_sect]->Fill(
						tag_counter, deltaK_smeared, fill_weight);
					h_deltaK_tagh_cut_ccal_phi[ccal_phi_sect]->Fill(
						tag_counter, deltaK_smeared, fill_weight);
					h_deltaK_tagh_cut_fcal_layer[fcal_layer-1]->Fill(
						tag_counter, deltaK_smeared, fill_weight);
					h_deltaK_tagh_cut_ccal_layer[ccal_layer-1]->Fill(
						tag_counter, deltaK_smeared, fill_weight);
				} else { 
					h_deltaK_tagm_cut_fcal_phi[fcal_phi_sect]->Fill(
						tag_counter, deltaK_smeared, fill_weight);
					h_deltaK_tagm_cut_ccal_phi[ccal_phi_sect]->Fill(
						tag_counter, deltaK_smeared, fill_weight);
					h_deltaK_tagm_cut_fcal_layer[fcal_layer-1]->Fill(
						tag_counter, deltaK_smeared, fill_weight);
					h_deltaK_tagm_cut_ccal_layer[ccal_layer-1]->Fill(
						tag_counter, deltaK_smeared, fill_weight);
				}
				
				h_xy_fcal_phi[fcal_phi_sect]->Fill(x1, y1, fill_weight);
				h_xy_ccal_phi[ccal_phi_sect]->Fill(x2, y2, fill_weight);
				h_xy_fcal_layer[fcal_layer-1]->Fill(x1, y1, fill_weight);
				h_xy_ccal_layer[ccal_layer-1]->Fill(x2, y2, fill_weight);
			}
		}
		
		//-------------------------------------------//
		
		if(tag_sys==0) { // TAGH
			
			if(!fcal_fid_cut && !ccal_fid_cut) {
				
				if(e5_cut && p5_cut && e2 > CCAL_min_energy_cut 
					&& fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<20; icut++) {
						double loc_cut = 0.05 * (double)(icut);
						if(e1 > loc_cut) {
							h_deltaK_tagh_fcalE[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							if(k5_cut) {
								h_deltaK_tagh_cut_fcalE[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && p5_cut && e1 > FCAL_min_energy_cut 
					&& fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<13; icut++) {
						double loc_cut = 0.5 * (double)(icut);
						if(e2 > loc_cut) {
							h_deltaK_tagh_ccalE[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							if(k5_cut) {
								h_deltaK_tagh_cut_ccalE[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && e1 > FCAL_min_energy_cut && e2 > CCAL_min_energy_cut 
					 && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<16; icut++) {
						if(phi_cut[icut]) {
							h_deltaK_tagh_sigPhi[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							if(k5_cut) {
								h_deltaK_tagh_cut_sigPhi[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							}
						}
					}
				}
				
				if(p5_cut && e1 > FCAL_min_energy_cut && e2 > CCAL_min_energy_cut 
					 && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<16; icut++) {
						if(e_cut[icut]) {
							h_deltaK_tagh_sigE[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							if(k5_cut) {
								h_deltaK_tagh_cut_sigE[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && p5_cut && e1 > FCAL_min_energy_cut 
					&& e2 > CCAL_min_energy_cut && ccal_t_cut) {
					
					for(int icut=0; icut<10; icut++) {
						if(fcal_t_cut_vec[icut]) {
							h_deltaK_tagh_fcalT[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							if(k5_cut) {
								h_deltaK_tagh_cut_fcalT[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && p5_cut && e1 > FCAL_min_energy_cut 
					&& e2 > CCAL_min_energy_cut && fcal_t_cut) {
					
					for(int icut=0; icut<10; icut++) {
						if(ccal_t_cut_vec[icut]) {
							h_deltaK_tagh_ccalT[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							if(k5_cut) {
								h_deltaK_tagh_cut_ccalT[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							}
						}
					}
				}
				
				if(p5_cut && e5_cut && e1 > FCAL_min_energy_cut 
					&& e2 > CCAL_min_energy_cut && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<16; icut++) {
						if(k_cut[icut]) {
							h_deltaK_tagh_sigK[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
						}
					}
				}
				
			} // if(!fcal_fid_cut && !ccal_fid_cut)
			
		} else { // TAGM
			
			if(!fcal_fid_cut && !ccal_fid_cut) {
				
				if(e5_cut && p5_cut && e2 > CCAL_min_energy_cut 
					&& fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<20; icut++) {
						double loc_cut = 0.05 * (double)(icut);
						if(e1 > loc_cut) {
							h_deltaK_tagm_fcalE[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							if(k5_cut) {
								h_deltaK_tagm_cut_fcalE[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && p5_cut && e1 > FCAL_min_energy_cut 
					&& fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<13; icut++) {
						double loc_cut = 0.5 * (double)(icut);
						if(e2 > loc_cut) {
							h_deltaK_tagm_ccalE[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							if(k5_cut) {
								h_deltaK_tagm_cut_ccalE[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && e1 > FCAL_min_energy_cut && e2 > CCAL_min_energy_cut 
					 && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<16; icut++) {
						if(phi_cut[icut]) {
							h_deltaK_tagm_sigPhi[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							if(k5_cut) {
								h_deltaK_tagm_cut_sigPhi[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							}
						}
					}
				}
				
				if(p5_cut && e1 > FCAL_min_energy_cut && e2 > CCAL_min_energy_cut 
					 && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<16; icut++) {
						if(e_cut[icut]) {
							h_deltaK_tagm_sigE[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							if(k5_cut) {
								h_deltaK_tagm_cut_sigE[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && p5_cut && e1 > FCAL_min_energy_cut 
					&& e2 > CCAL_min_energy_cut && ccal_t_cut) {
					
					for(int icut=0; icut<10; icut++) {
						if(fcal_t_cut_vec[icut]) {
							h_deltaK_tagm_fcalT[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							if(k5_cut) {
								h_deltaK_tagm_cut_fcalT[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && p5_cut && e1 > FCAL_min_energy_cut 
					&& e2 > CCAL_min_energy_cut && fcal_t_cut) {
					
					for(int icut=0; icut<10; icut++) {
						if(ccal_t_cut_vec[icut]) {
							h_deltaK_tagm_ccalT[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							if(k5_cut) {
								h_deltaK_tagm_cut_ccalT[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
							}
						}
					}
				}
				
				if(p5_cut && e5_cut && e1 > FCAL_min_energy_cut 
					&& e2 > CCAL_min_energy_cut && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<16; icut++) {
						if(k_cut[icut]) {
							h_deltaK_tagm_sigK[icut]->Fill(tag_counter, deltaK_smeared, fill_weight);
						}
					}
				}
				
			} // if(!fcal_fid_cut && !ccal_fid_cut)
		}
		
		
		if(e5_cut && p5_cut && k5_cut && e1 > FCAL_min_energy_cut 
			&& e2 > CCAL_min_energy_cut && fcal_t_cut && ccal_t_cut) 
		{
			for(int icut = 0; icut < N_FID_CUTS; icut++) {
				
				int loc_fcal_fid_cut = 0;
				double loc_fcal_inner_layer_cut = 2.0*4.0157 + 
					(4.0157/10.)*(double)icut;
				
				if( (-1.*loc_fcal_inner_layer_cut < fcal_face_x
					&& fcal_face_x < loc_fcal_inner_layer_cut)
					&& (-1.*loc_fcal_inner_layer_cut < fcal_face_y 
					&& fcal_face_y < loc_fcal_inner_layer_cut) ) 
						loc_fcal_fid_cut = 1;
				
				if(!loc_fcal_fid_cut && !ccal_fid_cut) {
					
					if(tag_sys==0) 
						h_deltaK_tagh_fcalfid[icut]->Fill(tag_counter, 
							deltaK_smeared, fill_weight);
					else 
						h_deltaK_tagm_fcalfid[icut]->Fill(tag_counter, 
							deltaK_smeared, fill_weight);
				}
			}
			
			for(int icut = 0; icut < N_FID_CUTS; icut++) {
				
				int loc_ccal_fid_cut = 0;
				double loc_ccal_inner_layer_cut = 1.0*2.09 + 
					(2.09/10.)*(double)icut;
				
				if( (-1.*loc_ccal_inner_layer_cut < ccal_face_x
					&& ccal_face_x < loc_ccal_inner_layer_cut)
					&& (-1.*loc_ccal_inner_layer_cut < ccal_face_y 
					&& ccal_face_y < loc_ccal_inner_layer_cut) ) 
						loc_ccal_fid_cut = 1;
				
				if(!loc_ccal_fid_cut && !fcal_fid_cut) {
					
					if(tag_sys==0) 
						h_deltaK_tagh_ccalfid[icut]->Fill(tag_counter, 
							deltaK_smeared, fill_weight);
					else 
						h_deltaK_tagm_ccalfid[icut]->Fill(tag_counter, 
							deltaK_smeared, fill_weight);
				}
			}
		}
	}
	
	if(n_good_cands==1) {
		if(loc_tag_sys==0) 
			h_deltaK_tagh_single->Fill(loc_tag_counter, loc_deltaK, 1.0);
		else
			h_deltaK_tagm_single->Fill(loc_tag_counter, loc_deltaK, 1.0);
	}
	
	return;
}



void JEventProcessor_compton_simulation_systematics::set_cuts(int32_t runnumber)
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
		
		// Phase I, He Target (stand-in values, adjust later)
		
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
		
		deltaE_mu_p0_mc      = -4.85724e-02;  deltaE_mu_p0_data    = -2.85756e-01;
		deltaE_mu_p1_mc      = -6.69318e-03;  deltaE_mu_p1_data    =  8.91005e-02;
		deltaE_mu_p2_mc      =  8.18727e-04;  deltaE_mu_p2_data    = -8.73318e-03;
		deltaE_mu_p3_mc      =  0.; 			    deltaE_mu_p3_data    =  2.45143e-04;
		
		deltaE_sig_p0_mc     =  9.25478e-03;  deltaE_sig_p0_data   =  1.81080e-02;
		deltaE_sig_p1_mc     =  3.95553e-02;  deltaE_sig_p1_data   =  1.38460e-02;
		deltaE_sig_p2_mc     =  7.63579e-08;  deltaE_sig_p2_data   =  5.76009e-02;
		//----------------------------------------------------------------------//
		deltaPhi_mu_p0_mc    =  1.79431e+02;  deltaPhi_mu_p0_data  =  1.80011e+02;
		deltaPhi_mu_p1_mc    =  1.10701e-01;  deltaPhi_mu_p1_data  = -8.18407e-03;
		deltaPhi_mu_p2_mc    = -1.31671e-02;  deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_mc    =  4.13213e-04;  deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_mc   =  1.24358e+01;  deltaPhi_sig_p0_data =  1.22089e+01;
		deltaPhi_sig_p1_mc   = -1.92471e+00;  deltaPhi_sig_p1_data = -1.81467e+00;
		deltaPhi_sig_p2_mc   =  1.86377e-01;  deltaPhi_sig_p2_data =  1.72870e-01;
		deltaPhi_sig_p3_mc   = -6.17908e-03;  deltaPhi_sig_p3_data = -5.55077e-03;
		//----------------------------------------------------------------------//
		deltaK_mu_p0_mc      =  1.87111e-01;  deltaK_mu_p0_data    =  3.22415e-01;
		deltaK_mu_p1_mc      = -5.01049e-02;  deltaK_mu_p1_data    = -8.24914e-02;
		deltaK_mu_p2_mc      =  6.19433e-04;  deltaK_mu_p2_data    =  4.25586e-03;
		deltaK_mu_p3_mc      = -2.58358e-05;  deltaK_mu_p3_data    = -1.88916e-04;
		
		deltaK_sig_p0_mc     =  7.66563e-01;  deltaK_sig_p0_data   =  5.70095e-01;
		deltaK_sig_p1_mc     = -1.20835e-01;  deltaK_sig_p1_data   = -4.93059e-02;
		deltaK_sig_p2_mc     =  1.93384e-02;  deltaK_sig_p2_data   =  1.19542e-02;
		deltaK_sig_p3_mc     = -6.91971e-04;  deltaK_sig_p3_data   = -4.17058e-04;
		
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


double JEventProcessor_compton_simulation_systematics::get_vertex_weight(double vertex_z) {
	
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
