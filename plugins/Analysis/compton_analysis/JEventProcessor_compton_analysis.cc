// $Id$
//
//    File: JEventProcessor_compton_analysis.cc
// Created: Thu Feb 10 21:33:24 EST 2022
// Creator: andrsmit (on Linux ifarm1802.jlab.org 3.10.0-1160.11.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_compton_analysis.h"

// Routine used to create our JEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_compton_analysis());
}
} // "C"

//------------------
// init
//------------------
jerror_t JEventProcessor_compton_analysis::init(void)
{
	TDirectory *dir_compton = new TDirectoryFile("compton_analysis", 
		"compton_analysis");
	dir_compton->cd();
	
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
	
	h_ccalE     = new TH2F("ccalE", "Energy of CCAL Shower; E_{Beam} [GeV]; E_{CCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	h_fcalE     = new TH2F("fcalE", "Energy of FCAL Shower; E_{Beam} [GeV]; E_{FCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	h_ccalE_cut = new TH2F("ccalE_cut", "Energy of CCAL Shower; E_{Beam} [GeV]; E_{CCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	h_fcalE_cut = new TH2F("fcalE_cut", "Energy of FCAL Shower; E_{Beam} [GeV]; E_{FCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	
	//--------------------------------------------------------------------//
	
	h_deltaE_vs_deltaK = new TH2F("deltaE_vs_deltaK", 
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
	
	h_deltaPhi_tagh = new TH2F("deltaPhi_tagh", 
		"#Delta#phi; TAGH Counter; |#phi_{1}-#phi_{2}| [deg.]", 
		274, 0.5, 274.5, 3600, 0.0, 360.0);
	h_deltaPhi_tagm = new TH2F("deltaPhi_tagm", 
		"#Delta#phi; TAGM Counter; |#phi_{1}-#phi_{2}| [deg.]", 
		102, 0.5, 102.5, 3600, 0.0, 360.0);
	
	h_deltaK_tagh = new TH2F("deltaK_tagh", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm = new TH2F("deltaK_tagm", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	
	h_deltaK_tagh_cut = new TH2F("deltaK_tagh_cut", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm_cut = new TH2F("deltaK_tagm_cut", 
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
jerror_t JEventProcessor_compton_analysis::brun(JEventLoop *eventLoop, int32_t runnumber)
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
		
		if(runnumber<110622) {
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
		m_ccalX_new = 0.135;
		m_ccalY_new = 0.135;
		
		m_beamX     = 0.146;
		m_beamY     = 0.017;
	}
	
	fcal_correction.SetXYZ(m_fcalX_new-m_fcalX, m_fcalY_new-m_fcalY, 0.);
	ccal_correction.SetXYZ(m_ccalX_new-m_ccalX, m_ccalY_new-m_ccalY, 0.);
	
	set_cuts(runnumber);
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_compton_analysis::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
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
	
	vector<const DFCALShower*> good_fcal_showers;
	vector<const DCCALShower*> good_ccal_showers;
	
	//-----   RF Bunch   -----//
	
	const DEventRFBunch *locRFBunch = NULL;
	try { 
		eventLoop->GetSingle(locRFBunch, "CalorimeterOnly");
	} catch (...) { 
		return NOERROR;
	}
	double rfTime = locRFBunch->dTime;
	
	japp->RootFillLock(this);  // Acquire root lock
	
	for(vector< const DBeamPhoton* >::const_iterator gam = beam_photons.begin();
		gam != beam_photons.end(); gam++) {
		
		int counter = (*gam)->dCounter;
		
		h_beam_rf_dt->Fill((*gam)->time() - rfTime);
		
		DetectorSystem_t sys = (*gam)->dSystem;
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
			good_fcal_showers.push_back((*show));
			if((*show)->getEnergy() > FCAL_min_energy_cut) {
				n_fcal_showers++;
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
			good_ccal_showers.push_back((*show));
			if((*show)->E > CCAL_min_energy_cut) {
				n_ccal_showers++;
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
				else if((-(2.004 + 5.*4.008)<=brfdt && brfdt<=-(2.004 + 3.*4.008))
					||((2.004 + 3.*4.008)<=brfdt && brfdt<=(2.004 + 5.*4.008)))
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
				
				ComptonCandidate_t loc_Cand;
				
				loc_Cand.e1 = e1;
				loc_Cand.x1 = pos1.X();
				loc_Cand.y1 = pos1.Y();
				loc_Cand.z1 = pos1.Z();
				loc_Cand.e2 = e2;
				loc_Cand.x2 = pos2.X();
				loc_Cand.y2 = pos2.Y();
				loc_Cand.z2 = pos2.Z();
				
				loc_Cand.deltaPhi    = deltaPhi;
				loc_Cand.deltaT      = deltaT;
				loc_Cand.deltaE      = deltaE;
				loc_Cand.deltaK      = deltaK;
				
				loc_Cand.bunch_val   = bunch_val;
				loc_Cand.eb          = eb;
				loc_Cand.brfdt       = brfdt;
				loc_Cand.tag_counter = (*gam)->dCounter;
				
				DetectorSystem_t sys = (*gam)->dSystem;
				if(sys==SYS_TAGH)      loc_Cand.tag_sys = 0;
				else if(sys==SYS_TAGM) loc_Cand.tag_sys = 1;
				
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
jerror_t JEventProcessor_compton_analysis::erun(void)
{
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_compton_analysis::fini(void)
{
	
	return NOERROR;
}


int JEventProcessor_compton_analysis::fcal_fiducial_cut(DVector3 pos, DVector3 vertex)
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


int JEventProcessor_compton_analysis::ccal_fiducial_cut(DVector3 pos, DVector3 vertex)
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



void JEventProcessor_compton_analysis::fill_histograms(
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
		
		double e1 = loc_Cand.e1;
		double x1 = loc_Cand.x1;
		double y1 = loc_Cand.y1;
		double e2 = loc_Cand.e2;
		double x2 = loc_Cand.x2;
		double y2 = loc_Cand.y2;
		
		//--------------     Cuts      --------------//
		
		double deltaE_mu_data  = f_deltaE_mu_data->Eval(eb);
		double deltaE_sig_data = eb * f_deltaE_sig_data->Eval(eb);
		
		int e_cut = 0;
		if(fabs(deltaE - deltaE_mu_data) < 5.0*deltaE_sig_data) e_cut = 1;
		
		double deltaPhi_mu_data  = f_deltaPhi_mu_data->Eval(eb);
		double deltaPhi_sig_data = f_deltaPhi_sig_data->Eval(eb);
		
		int p_cut = 0;
		if(bfield_val) {
			if(fabs(deltaPhi - 180.) < 50.) {
				p_cut = 1;
			}
		} else {
			if(fabs(deltaPhi - deltaPhi_mu_data) < 5.0*deltaPhi_sig_data) p_cut = 1;
		}
		
		double deltaK_mu_data  = f_deltaK_mu_data->Eval(eb);
		double deltaK_sig_data = f_deltaK_sig_data->Eval(eb);
		
		int k_cut = 0;
		if(fabs(deltaK - deltaK_mu_data) < 5.0*deltaK_sig_data) k_cut = 1;
		
		int fcal_e_cut = 0;
		if(e1 > FCAL_min_energy_cut) fcal_e_cut = 1;
		
		int ccal_e_cut = 0;
		if(e2 > CCAL_min_energy_cut) ccal_e_cut = 1;
		
		//-------------------------------------------//
		
		double fill_weight;
		if(bunch_val) fill_weight =  1.0;
		else          fill_weight = -0.25;
		
		h_fcalE->Fill(eb, e1);
		h_ccalE->Fill(eb, e2);
		if(e_cut && p_cut && k_cut && fcal_e_cut && ccal_e_cut) {
			h_fcalE_cut->Fill(eb, e1);
			h_ccalE_cut->Fill(eb, e2);
			h_fcal_xy->Fill(x1, y1, fill_weight);
			h_ccal_xy->Fill(x2, y2, fill_weight);
		}
		if(!fcal_e_cut || !ccal_e_cut) continue;
		
		h_deltaE_vs_deltaK->Fill(deltaK, deltaE, fill_weight);
		
		if(tag_sys==0) {
			h_deltaE_tagh->Fill(tag_counter, deltaE, fill_weight);
			if(e_cut) {
				h_deltaPhi_tagh->Fill(tag_counter, deltaPhi, fill_weight);
				if(p_cut) {
					h_deltaK_tagh->Fill(tag_counter, deltaK, fill_weight);
					if(k_cut) {
						h_deltaK_tagh_cut->Fill(tag_counter, deltaK, fill_weight);
						if(bunch_val) { 
							h_deltaK_tagh_cut_main->Fill(tag_counter, deltaK);
						} else {
							h_deltaK_tagh_cut_acc->Fill(tag_counter, deltaK, 0.5);
						}
						h_beam_rf_dt_cut->Fill(brfdt);
					}
				}
			}
		} else {
			h_deltaE_tagm->Fill(tag_counter, deltaE, fill_weight);
			if(e_cut) {
				h_deltaPhi_tagm->Fill(tag_counter, deltaPhi, fill_weight);
				if(p_cut) {
					h_deltaK_tagm->Fill(tag_counter, deltaK, fill_weight);
					if(k_cut) {
						h_deltaK_tagm_cut->Fill(tag_counter, deltaK, fill_weight);
						if(bunch_val) { 
							h_deltaK_tagm_cut_main->Fill(tag_counter, deltaK);
						} else {
							h_deltaK_tagm_cut_acc->Fill(tag_counter, deltaK, 0.5);
						}
						h_beam_rf_dt_cut->Fill(brfdt);
					}
				}
			}
		}
	}
	
	return;
}


void JEventProcessor_compton_analysis::set_cuts(int32_t runnumber)
{
	if(runnumber>60000 && runnumber<61355) {
	  
		// Phase I, Be Target
		
		deltaE_mu_p0_data    =  8.33517e-03;
		deltaE_mu_p1_data    =  2.09025e-03;
		deltaE_mu_p2_data    = -1.09342e-04;
		deltaE_mu_p3_data    =  0.;
		
		deltaE_sig_p0_data   =  8.37004e-03;
		deltaE_sig_p1_data   =  4.56259e-02;
		deltaE_sig_p2_data   =  0.;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.79943e+02;
		deltaPhi_mu_p1_data  = -2.11766e-02;
		deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_data =  1.20139e+01;
		deltaPhi_sig_p1_data = -1.75486e+00;
		deltaPhi_sig_p2_data =  1.67515e-01;
		deltaPhi_sig_p3_data = -5.48316e-03;
		//--------------------------------//
		deltaK_mu_p0_data    = -9.36095e-02;
		deltaK_mu_p1_data    =  5.48923e-02;
		deltaK_mu_p2_data    = -1.19844e-02;
		deltaK_mu_p3_data    =  4.38188e-04;
		
		deltaK_sig_p0_data   =  6.68283e-01;
		deltaK_sig_p1_data   = -8.45642e-02;
		deltaK_sig_p2_data   =  1.61255e-02;
		deltaK_sig_p3_data   = -5.93363e-04;
		
	} else if(runnumber>60000 && runnumber<69999) {
		
		// Phase I, He Target (stand-in values, adjust later)
		
		deltaE_mu_p0_data    =  8.33517e-03;
		deltaE_mu_p1_data    =  2.09025e-03;
		deltaE_mu_p2_data    = -1.09342e-04;
		deltaE_mu_p3_data    =  0.;
		
		deltaE_sig_p0_data   =  8.37004e-03;
		deltaE_sig_p1_data   =  4.56259e-02;
		deltaE_sig_p2_data   =  0.;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.79943e+02;
		deltaPhi_mu_p1_data  = -2.11766e-02;
		deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_data =  1.20139e+01;
		deltaPhi_sig_p1_data = -1.75486e+00;
		deltaPhi_sig_p2_data =  1.67515e-01;
		deltaPhi_sig_p3_data = -5.48316e-03;
		//--------------------------------//
		deltaK_mu_p0_data    = -9.36095e-02;
		deltaK_mu_p1_data    =  5.48923e-02;
		deltaK_mu_p2_data    = -1.19844e-02;
		deltaK_mu_p3_data    =  4.38188e-04;
		
		deltaK_sig_p0_data   =  6.68283e-01;
		deltaK_sig_p1_data   = -8.45642e-02;
		deltaK_sig_p2_data   =  1.61255e-02;
		deltaK_sig_p3_data   = -5.93363e-04;
		
	} else if(runnumber>80000 && runnumber<81400) {
		
		// Phase II, Be Target
		
		deltaE_mu_p0_data    = -0.05;
		deltaE_mu_p1_data    =  0.0;
		deltaE_mu_p2_data    =  0.0;
		deltaE_mu_p3_data    =  0.0;
		
		deltaE_sig_p0_data   =  1.40025e-02;
		deltaE_sig_p1_data   =  3.55640e-02;
		deltaE_sig_p2_data   =  5.15193e-04;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.79816e+02;
		deltaPhi_mu_p1_data  = -7.06279e-03;
		deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_data =  1.21214e+01;
		deltaPhi_sig_p1_data = -1.84405e+00;
		deltaPhi_sig_p2_data =  1.81923e-01;
		deltaPhi_sig_p3_data = -6.13418e-03;
		//--------------------------------//
		deltaK_mu_p0_data    = -1.74619e-01;
		deltaK_mu_p1_data    =  1.01328e-01;
		deltaK_mu_p2_data    = -1.84593e-02;
		deltaK_mu_p3_data    =  7.13390e-04;
		
		deltaK_sig_p0_data   =  3.88777e-01;
		deltaK_sig_p1_data   =  2.50588e-02;
		deltaK_sig_p2_data   =  2.20573e-03;
		deltaK_sig_p3_data   = -2.76677e-05;
		
	} else if(runnumber>80000 && runnumber<81473) {
		
		// Phase II, He Target, Field OFF (stand-in values, adjust later)
		
		deltaE_mu_p0_data    = -0.05;
		deltaE_mu_p1_data    =  0.0;
		deltaE_mu_p2_data    =  0.0;
		deltaE_mu_p3_data    =  0.0;
		
		deltaE_sig_p0_data   =  1.40025e-02;
		deltaE_sig_p1_data   =  3.55640e-02;
		deltaE_sig_p2_data   =  5.15193e-04;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.79816e+02;
		deltaPhi_mu_p1_data  = -7.06279e-03;
		deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_data =  1.21214e+01;
		deltaPhi_sig_p1_data = -1.84405e+00;
		deltaPhi_sig_p2_data =  1.81923e-01;
		deltaPhi_sig_p3_data = -6.13418e-03;
		//--------------------------------//
		deltaK_mu_p0_data    = -1.74619e-01;
		deltaK_mu_p1_data    =  1.01328e-01;
		deltaK_mu_p2_data    = -1.84593e-02;
		deltaK_mu_p3_data    =  7.13390e-04;
		
		deltaK_sig_p0_data   =  3.88777e-01;
		deltaK_sig_p1_data   =  2.50588e-02;
		deltaK_sig_p2_data   =  2.20573e-03;
		deltaK_sig_p3_data   = -2.76677e-05;
		
	} else if(runnumber>80000 && runnumber<89999) {
		
		// Phase II, He Target, Field ON (stand-in values, adjust later)
		
		deltaE_mu_p0_data    = -0.05;
		deltaE_mu_p1_data    =  0.0;
		deltaE_mu_p2_data    =  0.0;
		deltaE_mu_p3_data    =  0.0;
		
		deltaE_sig_p0_data   =  1.40025e-02;
		deltaE_sig_p1_data   =  3.55640e-02;
		deltaE_sig_p2_data   =  5.15193e-04;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.79816e+02;
		deltaPhi_mu_p1_data  = -7.06279e-03;
		deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_data =  1.21214e+01;
		deltaPhi_sig_p1_data = -1.84405e+00;
		deltaPhi_sig_p2_data =  1.81923e-01;
		deltaPhi_sig_p3_data = -6.13418e-03;
		//--------------------------------//
		deltaK_mu_p0_data    = -1.74619e-01;
		deltaK_mu_p1_data    =  1.01328e-01;
		deltaK_mu_p2_data    = -1.84593e-02;
		deltaK_mu_p3_data    =  7.13390e-04;
		
		deltaK_sig_p0_data   =  3.88777e-01;
		deltaK_sig_p1_data   =  2.50588e-02;
		deltaK_sig_p2_data   =  2.20573e-03;
		deltaK_sig_p3_data   = -2.76677e-05;
		
	} else if(runnumber>110000 && runnumber<110584) {
		
		// Phase III, Be Target, Field ON (stand-in values, adjust later)
		
		deltaE_mu_p0_data    = -2.85756e-01;
		deltaE_mu_p1_data    =  8.91005e-02;
		deltaE_mu_p2_data    = -8.73318e-03;
		deltaE_mu_p3_data    =  2.45143e-04;
		
		deltaE_sig_p0_data   =  1.81080e-02;
		deltaE_sig_p1_data   =  1.38460e-02;
		deltaE_sig_p2_data   =  5.76009e-02;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.80011e+02;
		deltaPhi_mu_p1_data  = -8.18407e-03;
		deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_data =  1.22089e+01;
		deltaPhi_sig_p1_data = -1.81467e+00;
		deltaPhi_sig_p2_data =  1.72870e-01;
		deltaPhi_sig_p3_data = -5.55077e-03;
		//--------------------------------//
		deltaK_mu_p0_data    =  3.22415e-01;
		deltaK_mu_p1_data    = -8.24914e-02;
		deltaK_mu_p2_data    =  4.25586e-03;
		deltaK_mu_p3_data    = -1.88916e-04;
		
		deltaK_sig_p0_data   =  5.70095e-01;
		deltaK_sig_p1_data   = -4.93059e-02;
		deltaK_sig_p2_data   =  1.19542e-02;
		deltaK_sig_p3_data   = -4.17058e-04;
		
	} else if(runnumber>110000 && runnumber<110622) {
		
		// Phase III, Be Target, Field OFF
		
		deltaE_mu_p0_data    = -2.85756e-01;
		deltaE_mu_p1_data    =  8.91005e-02;
		deltaE_mu_p2_data    = -8.73318e-03;
		deltaE_mu_p3_data    =  2.45143e-04;
		
		deltaE_sig_p0_data   =  1.81080e-02;
		deltaE_sig_p1_data   =  1.38460e-02;
		deltaE_sig_p2_data   =  5.76009e-02;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.80011e+02;
		deltaPhi_mu_p1_data  = -8.18407e-03;
		deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_data =  1.22089e+01;
		deltaPhi_sig_p1_data = -1.81467e+00;
		deltaPhi_sig_p2_data =  1.72870e-01;
		deltaPhi_sig_p3_data = -5.55077e-03;
		//--------------------------------//
		deltaK_mu_p0_data    =  3.22415e-01;
		deltaK_mu_p1_data    = -8.24914e-02;
		deltaK_mu_p2_data    =  4.25586e-03;
		deltaK_mu_p3_data    = -1.88916e-04;
		
		deltaK_sig_p0_data   =  5.70095e-01;
		deltaK_sig_p1_data   = -4.93059e-02;
		deltaK_sig_p2_data   =  1.19542e-02;
		deltaK_sig_p3_data   = -4.17058e-04;
		
	} else if(runnumber>110000 && runnumber<111969) {
		
		// Phase III, He Target, Field ON (stand-in values, adjust later)
		
		deltaE_mu_p0_data    = -2.85756e-01;
		deltaE_mu_p1_data    =  8.91005e-02;
		deltaE_mu_p2_data    = -8.73318e-03;
		deltaE_mu_p3_data    =  2.45143e-04;
		
		deltaE_sig_p0_data   =  1.81080e-02;
		deltaE_sig_p1_data   =  1.38460e-02;
		deltaE_sig_p2_data   =  5.76009e-02;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.80011e+02;
		deltaPhi_mu_p1_data  = -8.18407e-03;
		deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_data =  1.22089e+01;
		deltaPhi_sig_p1_data = -1.81467e+00;
		deltaPhi_sig_p2_data =  1.72870e-01;
		deltaPhi_sig_p3_data = -5.55077e-03;
		//--------------------------------//
		deltaK_mu_p0_data    =  3.22415e-01;
		deltaK_mu_p1_data    = -8.24914e-02;
		deltaK_mu_p2_data    =  4.25586e-03;
		deltaK_mu_p3_data    = -1.88916e-04;
		
		deltaK_sig_p0_data   =  5.70095e-01;
		deltaK_sig_p1_data   = -4.93059e-02;
		deltaK_sig_p2_data   =  1.19542e-02;
		deltaK_sig_p3_data   = -4.17058e-04;
		
	} else if(runnumber>110000 && runnumber<119999) {
		
		// Phase III, He Target, Field OFF (stand-in values, adjust later)
		
		deltaE_mu_p0_data    = -2.85756e-01;
		deltaE_mu_p1_data    =  8.91005e-02;
		deltaE_mu_p2_data    = -8.73318e-03;
		deltaE_mu_p3_data    =  2.45143e-04;
		
		deltaE_sig_p0_data   =  1.81080e-02;
		deltaE_sig_p1_data   =  1.38460e-02;
		deltaE_sig_p2_data   =  5.76009e-02;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.80011e+02;
		deltaPhi_mu_p1_data  = -8.18407e-03;
		deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_data =  1.22089e+01;
		deltaPhi_sig_p1_data = -1.81467e+00;
		deltaPhi_sig_p2_data =  1.72870e-01;
		deltaPhi_sig_p3_data = -5.55077e-03;
		//--------------------------------//
		deltaK_mu_p0_data    =  3.22415e-01;
		deltaK_mu_p1_data    = -8.24914e-02;
		deltaK_mu_p2_data    =  4.25586e-03;
		deltaK_mu_p3_data    = -1.88916e-04;
		
		deltaK_sig_p0_data   =  5.70095e-01;
		deltaK_sig_p1_data   = -4.93059e-02;
		deltaK_sig_p2_data   =  1.19542e-02;
		deltaK_sig_p3_data   = -4.17058e-04;
		
	} else {
		
		// Placeholder for non-PrimEx runs:
		
		deltaE_mu_p0_data    = -2.85756e-01;
		deltaE_mu_p1_data    =  8.91005e-02;
		deltaE_mu_p2_data    = -8.73318e-03;
		deltaE_mu_p3_data    =  2.45143e-04;
		
		deltaE_sig_p0_data   =  1.81080e-02;
		deltaE_sig_p1_data   =  1.38460e-02;
		deltaE_sig_p2_data   =  5.76009e-02;
		//--------------------------------//
		deltaPhi_mu_p0_data  =  1.80011e+02;
		deltaPhi_mu_p1_data  = -8.18407e-03;
		deltaPhi_mu_p2_data  =  0.;
		deltaPhi_mu_p3_data  =  0.;
		
		deltaPhi_sig_p0_data =  1.22089e+01;
		deltaPhi_sig_p1_data = -1.81467e+00;
		deltaPhi_sig_p2_data =  1.72870e-01;
		deltaPhi_sig_p3_data = -5.55077e-03;
		//--------------------------------//
		deltaK_mu_p0_data    =  3.22415e-01;
		deltaK_mu_p1_data    = -8.24914e-02;
		deltaK_mu_p2_data    =  4.25586e-03;
		deltaK_mu_p3_data    = -1.88916e-04;
		
		deltaK_sig_p0_data   =  5.70095e-01;
		deltaK_sig_p1_data   = -4.93059e-02;
		deltaK_sig_p2_data   =  1.19542e-02;
		deltaK_sig_p3_data   = -4.17058e-04;
	}
	
	
	f_deltaE_mu_data = new TF1("f_deltaE_mu_data", "pol3", 3.0, 12.0);
	f_deltaE_mu_data->SetParameters(deltaE_mu_p0_data, deltaE_mu_p1_data, 
		deltaE_mu_p2_data, deltaE_mu_p3_data);
	
	f_deltaE_sig_data = new TF1("f_deltaE_sig_data", 
		"sqrt([0]*[0] + ([1]/sqrt(x))*([1]/sqrt(x)) + ([2]/x)*([2]/x))", 3.0, 12.0);
	f_deltaE_sig_data->SetParameters(deltaE_sig_p0_data, deltaE_sig_p1_data, 
		deltaE_sig_p2_data);
	
	//--------------------------------------------------------------------------------------//
	
	f_deltaPhi_mu_data = new TF1("f_deltaPhi_mu_data", "pol3", 3.0, 12.0);
	f_deltaPhi_mu_data->SetParameters(deltaPhi_mu_p0_data, deltaPhi_mu_p1_data, 
		deltaPhi_mu_p2_data, deltaPhi_mu_p3_data);
	
	f_deltaPhi_sig_data = new TF1("f_deltaPhi_sig_data", "pol3", 3.0, 12.0);
	f_deltaPhi_sig_data->SetParameters(deltaPhi_sig_p0_data, deltaPhi_sig_p1_data, 
		deltaPhi_sig_p2_data, deltaPhi_sig_p3_data);
	
	//--------------------------------------------------------------------------------------//
	
	f_deltaK_mu_data = new TF1("f_deltaK_mu_data", "pol3", 3.0, 12.0);
	f_deltaK_mu_data->SetParameters(deltaK_mu_p0_data, deltaK_mu_p1_data, 
		deltaK_mu_p2_data, deltaK_mu_p3_data);
	
	f_deltaK_sig_data = new TF1("f_deltaK_sig_data", "pol3", 3.0, 12.0);
	f_deltaK_sig_data->SetParameters(deltaK_sig_p0_data, deltaK_sig_p1_data, 
		deltaK_sig_p2_data, deltaK_sig_p3_data);
	
	
	return;
}
