// $Id$
//
//    File: JEventProcessor_compton_tree.cc
// Created: Wed Dec  9 16:09:48 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_compton_tree.h"


// Routine used to create our JEventProcessor
extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_compton_tree());
	}
} // "C"

thread_local DTreeFillData JEventProcessor_compton_tree::dTreeFillData;

JEventProcessor_compton_tree::JEventProcessor_compton_tree() {
	
	m_RFTimeCut = 6.0;
	gPARMS->SetDefaultParameter("compton_tree:RFTimeCut", m_RFTimeCut);
	
	m_BeamEnergyCut = 6.0;
	gPARMS->SetDefaultParameter("compton_tree:BeamEnergyCut", m_BeamEnergyCut);
	
	m_DeltaECut = 3.0;
	gPARMS->SetDefaultParameter("compton_tree:BeamEnergyCut", m_DeltaECut);
}

//------------------
// init
//------------------
jerror_t JEventProcessor_compton_tree::init(void)
{
	dTreeInterface = DTreeInterface::Create_DTreeInterface("primex_compton","primex_compton.root");
	
	DTreeBranchRegister locTreeBranchRegister;
	
	locTreeBranchRegister.Register_Single<Int_t>("eventNum");
	locTreeBranchRegister.Register_Single<Double_t>("rfTime");
	
	// Beam Photon:
	locTreeBranchRegister.Register_Single<Int_t>("nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("tag_counter","nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("tag_sys","nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("e_beam","nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("t_beam","nbeam");
	
	// FCAL Showers:
	locTreeBranchRegister.Register_Single<Int_t>("nfcal");
	// reconstructed quantities:
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("e_fcal","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("x_fcal","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("y_fcal","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("z_fcal","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("t_fcal","nfcal");
	// info about shower:
	locTreeBranchRegister.Register_FundamentalArray<Int_t>( "nblocks_fcal","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>( "e1e9_fcal","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("e9e25_fcal","nfcal");
	
	// CCAL Showers:
	locTreeBranchRegister.Register_Single<Int_t>("nccal");
	// reconstructed quantities:
	locTreeBranchRegister.Register_FundamentalArray<Double_t>( "e_ccal","nccal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>( "x_ccal","nccal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>( "y_ccal","nccal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("x1_ccal","nccal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("y1_ccal","nccal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>( "z_ccal","nccal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>( "t_ccal","nccal");
	// info about shower:
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("nblocks_ccal","nccal");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>(  "idmax_ccal","nccal");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>(     "id_ccal","nccal");
	
	// TOF Points:
	locTreeBranchRegister.Register_Single<Int_t>("ntof");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("xtof","ntof");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("ytof","ntof");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("ztof","ntof");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("ttof","ntof");
	
	dTreeInterface->Create_Branches(locTreeBranchRegister);
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_compton_tree::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	DGeometry*   dgeom = NULL;
	DApplication* dapp = dynamic_cast<DApplication*>(eventLoop->GetJApplication());
	if(dapp)     dgeom = dapp->GetDGeometry(runnumber);
	
	if(dgeom){
		dgeom->GetTargetZ(m_beamZ);
		dgeom->GetFCALPosition(m_fcalX, m_fcalY, m_fcalZ);
		dgeom->GetCCALPosition(m_ccalX, m_ccalY, m_ccalZ);
	} else{
		cerr << "No geometry accessbile to compton_analysis plugin." << endl;
		return RESOURCE_UNAVAILABLE;
	}
	
	double m_fcalX_new = m_fcalX;
	double m_fcalY_new = m_fcalY;
	double m_ccalX_new = m_ccalX;
	double m_ccalY_new = m_ccalY;
	
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
		
	} else if(runnumber>80000 && runnumber<89999) {
		
		// PhaseII:
		m_fcalX_new =  0.408;
		m_fcalY_new =  0.027;
		m_ccalX_new =  0.108;
		m_ccalY_new =  0.130;
		
	} else if(runnumber>110000) {
		
		// Phase III:
		/*
		m_fcalX_new =  0.408;
		m_fcalY_new =  0.027;
		m_ccalX_new =  0.135;
		m_ccalY_new =  0.135;
		
		// beam parameters should ideally come from Simon's beam spot calibration using tracks.
		// however, in the meantime, I assume the FCAL position is the same from periods II to III.
		// Then, the offset between the beam and the FCAL from the compton alignment procedure 
		// gives the center of the beam spot. For simulations, I will also use this value as a 
		// placeholder, and will just use the same size of the beam spot from phase II.
		
		m_beamX     =  0.146;
		m_beamY     =  0.017;
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
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_compton_tree::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
	//--------------------------------------------------------------------------------------//
	// Check Trigger:
	
	const DL1Trigger *locTrig = NULL;
	try {
		eventLoop->GetSingle(locTrig);
	} catch (...) {}
	if(locTrig == NULL) { return NOERROR; }
	
	uint32_t locTrigMask   = locTrig->trig_mask;
	uint32_t locTrigFPMask = locTrig->fp_trig_mask;
	
	// skip FP trigger events:
	if(locTrigFPMask) return NOERROR;
	
	// check if main physics trigger is set:
	if(!(locTrigMask & 1)) return NOERROR;
	
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
	
	//--------------------------------------------------------------------------------------//
	// RF Time
	
	const DEventRFBunch *locRFBunch = NULL;
	try { 
		eventLoop->GetSingle(locRFBunch, "CalorimeterOnly");
	} catch (...) { return NOERROR; }
	double locRFTime = locRFBunch->dTime;
	if(locRFBunch->dNumParticleVotes < 2) return NOERROR;
	
	//--------------------------------------------------------------------------------------//
	// Criteria for writing events to Trees:
	//  0. Trigger bit 0 is set
	//  1. 1 DCCALShower & 1 DFCALShower within +/-6ns of RF time
	//  2. Corresponding >6GeV beam photon within +/-6ns of RF time
	//  3. Energy conservation within 3 GeV
	
	bool locEventSelector = false;
	
	for(vector<const DFCALShower*>::const_iterator show1 = locFCALShowers.begin(); 
		show1 != locFCALShowers.end(); show1++) {
		
		DVector3 pos1 = (*show1)->getPosition_log() - locVertex + m_fcal_correction;
		double e1     = (*show1)->getEnergy();
		double t1     = (*show1)->getTime() - (pos1.Mag()/c);
		
		if(fabs(t1-locRFTime) > m_RFTimeCut) continue;
		
		for(vector<const DCCALShower*>::const_iterator show2 = locCCALShowers.begin(); 
			show2 != locCCALShowers.end(); show2++) {
			
			DVector3 pos2((*show2)->x1, (*show2)->y1, (*show2)->z);
			pos2 += (m_ccal_correction - locVertex);
			double t2 = (*show2)->time - (pos2.Mag()/c);
			double e2 = (*show2)->E;
			
			if(fabs(t2-locRFTime) > m_RFTimeCut) continue;
			
			for(vector<const DBeamPhoton*>::const_iterator gam = locBeamPhotons.begin();
				gam != locBeamPhotons.end(); gam++) {
				
				double eb = (*gam)->lorentzMomentum().E();
				double tb = (*gam)->time();
				
				if(eb < m_BeamEnergyCut) continue;
				if(fabs(tb-locRFTime) < m_RFTimeCut) continue;
				
				double locDeltaE = e1+e2 - eb;
				if(fabs(locDeltaE) < m_DeltaECut) {
					locEventSelector = true;
				}
				
			} // end DBeamPhoton loop
			
		} // end DCCALShower loop
		
	} // end DFCALShower loop
	
	if(locEventSelector) {
		write_events(eventnumber, locRFTime, locBeamPhotons, locFCALShowers, locCCALShowers, locTOFPoints);
	}
	
	dTreeInterface->Fill(dTreeFillData);
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_compton_tree::fini(void)
{
	delete dTreeInterface;
	return NOERROR;
}

void JEventProcessor_compton_tree::write_events(uint64_t eventnumber, double rfTime,
	vector<const DBeamPhoton*> beam_photons, 
	vector<const DFCALShower*> fcal_showers,
	vector<const DCCALShower*> ccal_showers,
	vector<const DTOFPoint*> tof_points) {
	
	dTreeFillData.Fill_Single<Int_t>("eventNum", eventnumber);
	dTreeFillData.Fill_Single<Double_t>("rfTime", rfTime);
	
	// Beam Photons:
	size_t n_beam_photon = 0;
	for(vector<const DBeamPhoton*>::const_iterator gam = beam_photons.begin();
		gam != beam_photons.end(); gam++) {
		
		int loc_counter = (*gam)->dCounter;
		int loc_sys = -1;
		if((*gam)->dSystem == SYS_TAGH) loc_sys = 0;
		else if((*gam)->dSystem == SYS_TAGM) loc_sys = 1;
		
		dTreeFillData.Fill_Array<Int_t>("tag_counter", loc_counter, n_beam_photon);
		dTreeFillData.Fill_Array<Int_t>("tag_sys",     loc_sys,     n_beam_photon);
		dTreeFillData.Fill_Array<Double_t>("e_beam", (*gam)->lorentzMomentum().E(), n_beam_photon);
		dTreeFillData.Fill_Array<Double_t>("t_beam", (*gam)->time(),                n_beam_photon);
		
		n_beam_photon++;
	}
	dTreeFillData.Fill_Single<Int_t>("nbeam", n_beam_photon);
	
	// FCAL Showers:
	size_t n_fcal_shower = 0;
	for(vector<const DFCALShower*>::const_iterator show = fcal_showers.begin(); 
		show != fcal_showers.end(); show++) {
		
		dTreeFillData.Fill_Array<Double_t>("e_fcal", (*show)->getEnergy(),           n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("x_fcal", (*show)->getPosition_log().X(), n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("y_fcal", (*show)->getPosition_log().Y(), n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("z_fcal", (*show)->getPosition_log().Z(), n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("t_fcal", (*show)->getTime(),             n_fcal_shower);
		
		dTreeFillData.Fill_Array<Int_t>( "nblocks_fcal", (*show)->getNumBlocks(), n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>( "e1e9_fcal", (*show)->getE1E9(),      n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("e9e25_fcal", (*show)->getE9E25(),     n_fcal_shower);
		
		n_fcal_shower++;
	}
	dTreeFillData.Fill_Single<Int_t>("nfcal", n_fcal_shower);
	
	
	// CCAL Showers:
	size_t n_ccal_shower = 0;
	for(vector<const DCCALShower*>::const_iterator show = ccal_showers.begin(); 
		show != ccal_showers.end(); show++) {
		
		dTreeFillData.Fill_Array<Double_t>( "e_ccal", (*show)->E,    n_ccal_shower);
		dTreeFillData.Fill_Array<Double_t>( "x_ccal", (*show)->x,    n_ccal_shower);
		dTreeFillData.Fill_Array<Double_t>( "y_ccal", (*show)->y,    n_ccal_shower);
		dTreeFillData.Fill_Array<Double_t>("x1_ccal", (*show)->x1,   n_ccal_shower);
		dTreeFillData.Fill_Array<Double_t>("y1_ccal", (*show)->y1,   n_ccal_shower);
		dTreeFillData.Fill_Array<Double_t>( "z_ccal", (*show)->z,    n_ccal_shower);
		dTreeFillData.Fill_Array<Double_t>( "t_ccal", (*show)->time, n_ccal_shower);
		
		dTreeFillData.Fill_Array<Int_t>("nblocks_ccal", (*show)->dime,  n_ccal_shower);
		dTreeFillData.Fill_Array<Int_t>(  "idmax_ccal", (*show)->idmax, n_ccal_shower);
		dTreeFillData.Fill_Array<Int_t>(     "id_ccal", (*show)->id,    n_ccal_shower);
		
		n_ccal_shower++;
	}
	dTreeFillData.Fill_Single<Int_t>("nccal", n_ccal_shower);
	
	// TOF Points:
	size_t n_tof_points = 0;
	for(vector<const DTOFPoint*>::const_iterator tof = tof_points.begin(); 
		tof != tof_points.end(); tof++) {
		
		dTreeFillData.Fill_Array<Double_t>("x_tof", (*tof)->pos.X(), n_tof_points);
		dTreeFillData.Fill_Array<Double_t>("y_tof", (*tof)->pos.Y(), n_tof_points);
		dTreeFillData.Fill_Array<Double_t>("z_tof", (*tof)->pos.Z(), n_tof_points);
		dTreeFillData.Fill_Array<Double_t>("t_tof", (*tof)->t,       n_tof_points);
		
		n_tof_points++;
	}
	dTreeFillData.Fill_Single<Int_t>("ntof", n_tof_points);
	
	return;
}
