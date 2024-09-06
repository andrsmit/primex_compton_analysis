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
	
	/*
	There are 2 options for writing out events to TTrees:
	1. Apply very minimal cuts and write out as much information as needed for performing systematics studies.
		This should be applied for select-run used for systematics studies, but doesn't need to be done
		for all runs.
	2. Apply more selective cuts to reduce the output file size.
	
	These two options are selected by the parameter 'm_LOOSE_CUTS' which is either 0 (default) or 1, and 
	can be changed at runtime in a jana config file.
	*/
	
	m_LOOSE_CUTS = 0;
	gPARMS->SetDefaultParameter("compton_tree:LOOSE_CUTS", m_LOOSE_CUTS);
	
	// set up event selection based:
	if(m_LOOSE_CUTS) {
		
		// minimum energy cuts [GeV]:
		m_cut_fcalE = 0.0;
		m_cut_ccalE = 0.0;
		m_cut_beamE = 6.0;
		
		// timing cuts [ns]:
		m_cut_fcalrfdt = 1.e3;
		m_cut_ccalrfdt = 1.e3;
		m_cut_beamrfdt = 1.e3;
		
		// deltaE cut width [GeV]:
		m_cut_deltaE = 8.0;
		
	} else {
		
		// minimum energy cuts [GeV]:
		m_cut_fcalE = 0.35;
		m_cut_ccalE = 3.0;
		m_cut_beamE = 6.0;
		
		// timing cuts [ns]:
		m_cut_fcalrfdt = 2.004;
		m_cut_ccalrfdt = 2.004;
		m_cut_beamrfdt = 2.004;
		
		// deltaE cut width [GeV]:
		m_cut_deltaE = 3.0;
	}
	
	m_SAVE_MC_NOHITS = 0;
	gPARMS->SetDefaultParameter("compton_tree:m_SAVE_MC_NOHITS", m_SAVE_MC_NOHITS);
}

//------------------
// init
//------------------
jerror_t JEventProcessor_compton_tree::init(void)
{
	// Histogram to keep track of the number of events in each file:
	
	h_gen_flux   = new TH1F("gen_flux",   "Generated Beam Energy; [GeV]", 12000, 0., 12.);
	h_gen_weight = new TH2F("gen_weight", "Generated Reaction Weight vs. Beam Energy; [GeV]; [mb]", 
		12000, 0., 12., 1000, 0., 500.);
	
	dTreeInterface = DTreeInterface::Create_DTreeInterface("primex_compton","primex_compton.root");
	
	DTreeBranchRegister locTreeBranchRegister;
	
	locTreeBranchRegister.Register_Single<Int_t>("eventNum");
	locTreeBranchRegister.Register_Single<Double_t>("rfTime");
	
	// Beam Photon:
	locTreeBranchRegister.Register_Single<Int_t>("nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("tag_counter","nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("tag_sys","nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("beam_e","nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("beam_t","nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("acc_scale_factor","nbeam");
	
	// FCAL Showers:
	locTreeBranchRegister.Register_Single<Int_t>("nfcal");
	// reconstructed quantities:
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_e","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_x","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_y","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_z","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_t","nfcal");
	// info about shower:
	locTreeBranchRegister.Register_FundamentalArray<Int_t>(   "fcal_nblocks","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_e1e9",   "nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_e9e25",  "nfcal");
	
	// CCAL Showers:
	locTreeBranchRegister.Register_Single<Int_t>("nccal");
	// reconstructed quantities:
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("ccal_e", "nccal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("ccal_x", "nccal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("ccal_y", "nccal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("ccal_x1","nccal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("ccal_y1","nccal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("ccal_z", "nccal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("ccal_t", "nccal");
	// info about shower:
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("ccal_nblocks","nccal");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("ccal_idmax",  "nccal");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("ccal_id",     "nccal");
	
	// TOF Points:
	locTreeBranchRegister.Register_Single<Int_t>("ntof");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("tof_x","ntof");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("tof_y","ntof");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("tof_z","ntof");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("tof_t","ntof");
	
	// MC Thrown:
	locTreeBranchRegister.Register_Single<Int_t>("nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_pdgtype","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_x","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_y","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_z","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_t","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_e","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_p","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_theta","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_phi","nmc");
	
	locTreeBranchRegister.Register_Single<Int_t>("mc_reaction_type");
	locTreeBranchRegister.Register_Single<Double_t>("mc_reaction_weight");
	locTreeBranchRegister.Register_Single<Double_t>("mc_reaction_energy");
	
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
	
	// Code to obtain the scaling factors for accidental beam bunches (copied from DAnalysisUtilities.cc in gluex_root_analysis)
	
	ostringstream locCommandStream;
	locCommandStream << "ccdb dump ANALYSIS/accidental_scaling_factor -r " << runnumber;
	FILE* locInputFile = gSystem->OpenPipe(locCommandStream.str().c_str(), "r");
	if(locInputFile == NULL) {
		
		m_HodoscopeHiFactor    = 1.00;
		m_HodoscopeHiFactorErr = 0.01;
		m_HodoscopeLoFactor    = 1.00;
		m_HodoscopeLoFactorErr = 0.01;
		m_MicroscopeFactor     = 1.00;
		m_MicroscopeFactorErr  = 0.01;
		m_TAGMEnergyBoundHi    = 9.00;
		m_TAGMEnergyBoundLo    = 8.00;
		
		return NOERROR;
		/*
		cerr << "Could not load ANALYSIS/accidental_scaling_factor from CCDB !" << endl;
		gSystem->Exit(1);        // make sure we don't fail silently
		return RESOURCE_UNAVAILABLE;    // sanity check, this shouldn't be executed!
		*/
	}
	
	//get the first line
	char buff[1024]; // I HATE char buffers
	if(fgets(buff, sizeof(buff), locInputFile) == NULL)
	{
		m_HodoscopeHiFactor    = 1.00;
		m_HodoscopeHiFactorErr = 0.01;
		m_HodoscopeLoFactor    = 1.00;
		m_HodoscopeLoFactorErr = 0.01;
		m_MicroscopeFactor     = 1.00;
		m_MicroscopeFactorErr  = 0.01;
		m_TAGMEnergyBoundHi    = 9.00;
		m_TAGMEnergyBoundLo    = 8.00;
		
		gSystem->ClosePipe(locInputFile);
		return NOERROR;
		/*
		//vector<double> locCachedValues = { -1., -1., -1., -1., -1., -1., -1., -1. };
		//dAccidentalScalingFactor_Cache[runnumber] = locCachedValues;   // give up for this run
		gSystem->ClosePipe(locInputFile);
		cerr << "Could not parse ANALYSIS/accidental_scaling_factor from CCDB !" << endl;
		gSystem->Exit(1);        // make sure we don't fail silently
		return RESOURCE_UNAVAILABLE;    // sanity check, this shouldn't be executed!
		*/
	}
	
	//get the second line (where the # is)
	if(fgets(buff, sizeof(buff), locInputFile) == NULL)
	{
		m_HodoscopeHiFactor    = 1.00;
		m_HodoscopeHiFactorErr = 0.01;
		m_HodoscopeLoFactor    = 1.00;
		m_HodoscopeLoFactorErr = 0.01;
		m_MicroscopeFactor     = 1.00;
		m_MicroscopeFactorErr  = 0.01;
		m_TAGMEnergyBoundHi    = 9.00;
		m_TAGMEnergyBoundLo    = 8.00;
		
		gSystem->ClosePipe(locInputFile);
		return NOERROR;
		/*
		//vector<double> locCachedValues = { -1., -1., -1., -1., -1., -1., -1., -1. };
		//dAccidentalScalingFactor_Cache[runnumber] = locCachedValues;   // give up for this run
		gSystem->ClosePipe(locInputFile);
		cerr << "Could not parse ANALYSIS/accidental_scaling_factor from CCDB !" << endl;
		gSystem->Exit(1);        // make sure we don't fail silently
		return RESOURCE_UNAVAILABLE;    // sanity check, this shouldn't be executed!
		*/
	}
	
	// catch some CCDB error conditions
	if(strncmp(buff, "Cannot", 6) == 0) 
	{
		m_HodoscopeHiFactor    = 1.00;
		m_HodoscopeHiFactorErr = 0.01;
		m_HodoscopeLoFactor    = 1.00;
		m_HodoscopeLoFactorErr = 0.01;
		m_MicroscopeFactor     = 1.00;
		m_MicroscopeFactorErr  = 0.01;
		m_TAGMEnergyBoundHi    = 9.00;
		m_TAGMEnergyBoundLo    = 8.00;
		
		gSystem->ClosePipe(locInputFile);
		return NOERROR;
		/*
		// no assignment for this run
		//vector<double> locCachedValues = { -1., -1., -1., -1., -1., -1., -1., -1. };
		//dAccidentalScalingFactor_Cache[runnumber] = locCachedValues;   // give up for this run
		gSystem->ClosePipe(locInputFile);
		cerr << "No data available for ANALYSIS/accidental_scaling_factor, run " << runnumber << " from CCDB !" << endl;
		gSystem->Exit(1);        // make sure we don't fail silently
		return RESOURCE_UNAVAILABLE;    // sanity check, this shouldn't be executed!
		*/
	}
	
	istringstream locStringStream(buff);
	
	double locHodoscopeHiFactor    = -1.0;
	double locHodoscopeHiFactorErr = -1.0;
	double locHodoscopeLoFactor    = -1.0;
	double locHodoscopeLoFactorErr = -1.0;
	double locMicroscopeFactor     = -1.0;
	double locMicroscopeFactorErr  = -1.0;
	double locTAGMEnergyBoundHi    = -1.0;
	double locTAGMEnergyBoundLo    = -1.0;
	
	//extract it
	locStringStream >> locHodoscopeHiFactor >> locHodoscopeHiFactorErr >> locHodoscopeLoFactor
		>> locHodoscopeLoFactorErr >> locMicroscopeFactor >> locMicroscopeFactorErr
		>> locTAGMEnergyBoundHi >> locTAGMEnergyBoundLo;
	
	//Close the pipe
	gSystem->ClosePipe(locInputFile);
	
	m_HodoscopeHiFactor    = locHodoscopeHiFactor;
	m_HodoscopeHiFactorErr = locHodoscopeHiFactorErr;
	m_HodoscopeLoFactor    = locHodoscopeLoFactor;
	m_HodoscopeLoFactorErr = locHodoscopeLoFactorErr;
	m_MicroscopeFactor     = locMicroscopeFactor;
	m_MicroscopeFactorErr  = locMicroscopeFactorErr;
	m_TAGMEnergyBoundHi    = locTAGMEnergyBoundHi;
	m_TAGMEnergyBoundLo    = locTAGMEnergyBoundLo;
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_compton_tree::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
	//--------------------------------------------------------------------------------------//
	// Check for thrown MC information:
	
	vector<const DMCThrown*> locMCThrown;
	eventLoop->Get(locMCThrown);
	
	vector<const DMCReaction*> locMCReaction;
	eventLoop->Get(locMCReaction);
	
	int locIsMC = 0;
	if(locMCThrown.size()) locIsMC = 1;
	
	//--------------------------------------------------------------------------------------//
	// Check Trigger:
	
	const DL1Trigger *locTrig = NULL;
	try {
		eventLoop->GetSingle(locTrig);
	} catch (...) {}
	if(locTrig == NULL && locIsMC==0) { return NOERROR; }
	
	if(locTrig != NULL) {
		uint32_t locTrigMask   = locTrig->trig_mask;
		uint32_t locTrigFPMask = locTrig->fp_trig_mask;
		
		// skip FP trigger events:
		if(locTrigFPMask) return NOERROR;
		
		// check if main physics trigger is set:
		if(!(locTrigMask & 1)) return NOERROR;
	}
	
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
	// Save information on weights and number of simulated events for MC:
	
	if(locIsMC) {
		japp->RootFillLock(this);
		h_gen_flux->Fill(locMCReaction[0]->beam.energy());
		h_gen_weight->Fill(locMCReaction[0]->beam.energy(), locMCReaction[0]->weight*1.e-3);
		japp->RootFillUnLock(this);
	}
	
	//--------------------------------------------------------------------------------------//
	// RF Time
	
	const DEventRFBunch *locRFBunch = NULL;
	try { 
		eventLoop->GetSingle(locRFBunch, "CalorimeterOnly");
	} catch (...) { 
		if(locIsMC && m_SAVE_MC_NOHITS) {
			write_events(eventnumber, 0.0, locMCThrown, locMCReaction);
			dTreeInterface->Fill(dTreeFillData);
		}
		return NOERROR;
	}
	double locRFTime = locRFBunch->dTime;
	if(locRFBunch->dNumParticleVotes < 2) {
		if(locIsMC && m_SAVE_MC_NOHITS) {
			write_events(eventnumber, locRFTime, locMCThrown, locMCReaction);
			dTreeInterface->Fill(dTreeFillData);
		}
		return NOERROR;
	}
	
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
		
		if(fabs(t1-locRFTime) >= m_cut_fcalrfdt) continue;
		if(e1 <= m_cut_fcalE) continue;
		
		for(vector<const DCCALShower*>::const_iterator show2 = locCCALShowers.begin(); 
			show2 != locCCALShowers.end(); show2++) {
			
			DVector3 pos2((*show2)->x1, (*show2)->y1, (*show2)->z);
			pos2 += (m_ccal_correction - locVertex);
			double t2 = (*show2)->time - (pos2.Mag()/c);
			double e2 = (*show2)->E;
			
			if(fabs(t2-locRFTime) >= m_cut_ccalrfdt) continue;
			if(e2 <= m_cut_ccalE) continue;
			
			for(vector<const DBeamPhoton*>::const_iterator gam = locBeamPhotons.begin();
				gam != locBeamPhotons.end(); gam++) {
				
				double eb = (*gam)->lorentzMomentum().E();
				if(eb < m_cut_beamE) continue;
				
				double brfdt = fabs((*gam)->time() - locRFTime);
				
				// only save events where there is a beam photon passing the cuts in the main RF bunch 
				// OR in one of the accidental sidebands that will be used for subtraction. 
				// NOTE: we need to use the exact same sidebands when analyzing the trees and doing the 
				// accidental subtraction that are used in these cuts.
				
				int bunch_val = 0;
				if(brfdt < (m_beam_bunches_main*m_cut_beamrfdt)) bunch_val = 1;
				else if(
					((m_beam_bunches_main+5.5)*4.008 < brfdt) &&
					(brfdt < (m_beam_bunches_main+5.5+m_beam_bunches_acc)*4.008)
				) bunch_val = 2;
				
				if(bunch_val==0) continue;
				
				double locDeltaE = e1+e2 - eb;
				if(fabs(locDeltaE) < m_cut_deltaE) {
					locEventSelector = true;
				}
				
			} // end DBeamPhoton loop
			
		} // end DCCALShower loop
		
	} // end DFCALShower loop
	
	if(locEventSelector) {
		write_events(eventnumber, locRFTime, locBeamPhotons, locFCALShowers, locCCALShowers, locTOFPoints, 
			locMCThrown, locMCReaction);
		dTreeInterface->Fill(dTreeFillData);
	}
	else if(locIsMC && m_SAVE_MC_NOHITS) {
		write_events(eventnumber, locRFTime, locMCThrown, locMCReaction);
		dTreeInterface->Fill(dTreeFillData);
	}
	
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
	vector<const DMCThrown*> mc_thrown, vector<const DMCReaction*> mc_reaction) {
	
	dTreeFillData.Fill_Single<Int_t>("eventNum", eventnumber);
	dTreeFillData.Fill_Single<Double_t>("rfTime", rfTime);
	
	dTreeFillData.Fill_Single<Int_t>("nbeam", 0);
	dTreeFillData.Fill_Single<Int_t>("nfcal", 0);
	dTreeFillData.Fill_Single<Int_t>("nccal", 0);
	dTreeFillData.Fill_Single<Int_t>("ntof",  0);
	
	// MC Thrown:
	size_t n_mc_thrown = 0;
	for(vector<const DMCThrown*>::const_iterator mc = mc_thrown.begin(); mc != mc_thrown.end(); mc++) {
		
		dTreeFillData.Fill_Array<Double_t>("mc_pdgtype", (*mc)->pdgtype, n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_x",       (*mc)->position().X(),                n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_y",       (*mc)->position().Y(),                n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_z",       (*mc)->position().Z(),                n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_t",       (*mc)->time(),                        n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_e",       (*mc)->energy(),                      n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_p",       (*mc)->momentum().Mag(),              n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_theta",   (*mc)->momentum().Theta()*180.0/M_PI, n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_phi",     (*mc)->momentum().Phi()*180.0/M_PI,   n_mc_thrown);
		
		n_mc_thrown++;
	}
	dTreeFillData.Fill_Single<Int_t>("nmc", n_mc_thrown);
	
	if(mc_reaction.size() == 1) {
		dTreeFillData.Fill_Single<Int_t>("mc_reaction_type", mc_reaction[0]->type);
		dTreeFillData.Fill_Single<Double_t>("mc_reaction_weight", mc_reaction[0]->weight);
		dTreeFillData.Fill_Single<Double_t>("mc_reaction_energy", mc_reaction[0]->beam.energy());
	}
	
	return;
}

void JEventProcessor_compton_tree::write_events(uint64_t eventnumber, double rfTime,
	vector<const DBeamPhoton*> beam_photons, 
	vector<const DFCALShower*> fcal_showers,
	vector<const DCCALShower*> ccal_showers,
	vector<const DTOFPoint*> tof_points,
	vector<const DMCThrown*> mc_thrown,
	vector<const DMCReaction*> mc_reaction) {
	
	dTreeFillData.Fill_Single<Int_t>("eventNum", eventnumber);
	dTreeFillData.Fill_Single<Double_t>("rfTime", rfTime);
	
	// Beam Photons:
	size_t n_beam_photon = 0;
	for(vector<const DBeamPhoton*>::const_iterator gam = beam_photons.begin();
		gam != beam_photons.end(); gam++) {
		
		// 7.23.24: 
		// Only write out beam photons in the prompt RF peak, and the sidebands used for accidental subtraction:
		//
		// NOTE: we need to use the exact same sidebands when analyzing the trees and doing the 
		// accidental subtraction that are used in these cuts.
		
		double brfdt = fabs((*gam)->time() - rfTime);
		
		int bunch_val = 0;
		//if(brfdt < (m_beam_bunches_main*m_cut_beamrfdt)) bunch_val = 1;
		if(brfdt < (m_beam_bunches_main*2.004)) bunch_val = 1;
		else if(
			((m_beam_bunches_main+5.5)*4.008 < brfdt) &&
			(brfdt < (m_beam_bunches_main+5.5+m_beam_bunches_acc)*4.008)
		) bunch_val = 2;
		
		if(bunch_val==0) continue;
		
		int loc_counter = (*gam)->dCounter;
		int loc_sys = -1;
		if((*gam)->dSystem == SYS_TAGH) loc_sys = 0;
		else if((*gam)->dSystem == SYS_TAGM) loc_sys = 1;
		
		double loc_eb = (*gam)->lorentzMomentum().E();
		dTreeFillData.Fill_Array<Int_t>("tag_counter", loc_counter,    n_beam_photon);
		dTreeFillData.Fill_Array<Int_t>("tag_sys",     loc_sys,        n_beam_photon);
		dTreeFillData.Fill_Array<Double_t>("beam_e",   loc_eb,         n_beam_photon);
		dTreeFillData.Fill_Array<Double_t>("beam_t",  (*gam)->time(),  n_beam_photon);
		dTreeFillData.Fill_Array<Double_t>("acc_scale_factor", get_acc_scaling_factor(loc_eb), n_beam_photon);
		
		n_beam_photon++;
	}
	dTreeFillData.Fill_Single<Int_t>("nbeam", n_beam_photon);
	
	// FCAL Showers:
	size_t n_fcal_shower = 0;
	for(vector<const DFCALShower*>::const_iterator show = fcal_showers.begin(); 
		show != fcal_showers.end(); show++) {
		
		dTreeFillData.Fill_Array<Double_t>("fcal_e", (*show)->getEnergy(),           n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("fcal_x", (*show)->getPosition_log().X(), n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("fcal_y", (*show)->getPosition_log().Y(), n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("fcal_z", (*show)->getPosition_log().Z(), n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("fcal_t", (*show)->getTime(),             n_fcal_shower);
		
		dTreeFillData.Fill_Array<Int_t>(   "fcal_nblocks", (*show)->getNumBlocks(), n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("fcal_e1e9",    (*show)->getE1E9(),      n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("fcal_e9e25",   (*show)->getE9E25(),     n_fcal_shower);
		
		n_fcal_shower++;
	}
	dTreeFillData.Fill_Single<Int_t>("nfcal", n_fcal_shower);
	
	
	// CCAL Showers:
	size_t n_ccal_shower = 0;
	for(vector<const DCCALShower*>::const_iterator show = ccal_showers.begin(); 
		show != ccal_showers.end(); show++) {
		
		dTreeFillData.Fill_Array<Double_t>("ccal_e",  (*show)->E,    n_ccal_shower);
		dTreeFillData.Fill_Array<Double_t>("ccal_x",  (*show)->x,    n_ccal_shower);
		dTreeFillData.Fill_Array<Double_t>("ccal_y",  (*show)->y,    n_ccal_shower);
		dTreeFillData.Fill_Array<Double_t>("ccal_x1", (*show)->x1,   n_ccal_shower);
		dTreeFillData.Fill_Array<Double_t>("ccal_y1", (*show)->y1,   n_ccal_shower);
		dTreeFillData.Fill_Array<Double_t>("ccal_z",  (*show)->z,    n_ccal_shower);
		dTreeFillData.Fill_Array<Double_t>("ccal_t",  (*show)->time, n_ccal_shower);
		
		dTreeFillData.Fill_Array<Int_t>("ccal_nblocks", (*show)->dime,  n_ccal_shower);
		dTreeFillData.Fill_Array<Int_t>("ccal_idmax",   (*show)->idmax, n_ccal_shower);
		dTreeFillData.Fill_Array<Int_t>("ccal_id",      (*show)->id,    n_ccal_shower);
		
		n_ccal_shower++;
	}
	dTreeFillData.Fill_Single<Int_t>("nccal", n_ccal_shower);
	
	// TOF Points:
	size_t n_tof_points = 0;
	for(vector<const DTOFPoint*>::const_iterator tof = tof_points.begin(); 
		tof != tof_points.end(); tof++) {
		
		dTreeFillData.Fill_Array<Double_t>("tof_x", (*tof)->pos.X(), n_tof_points);
		dTreeFillData.Fill_Array<Double_t>("tof_y", (*tof)->pos.Y(), n_tof_points);
		dTreeFillData.Fill_Array<Double_t>("tof_z", (*tof)->pos.Z(), n_tof_points);
		dTreeFillData.Fill_Array<Double_t>("tof_t", (*tof)->t,       n_tof_points);
		
		n_tof_points++;
	}
	dTreeFillData.Fill_Single<Int_t>("ntof", n_tof_points);
	
	// MC Thrown:
	size_t n_mc_thrown = 0;
	for(vector<const DMCThrown*>::const_iterator mc = mc_thrown.begin(); mc != mc_thrown.end(); mc++) {
		
		dTreeFillData.Fill_Array<Double_t>("mc_pdgtype", (*mc)->pdgtype, n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_x",       (*mc)->position().X(),                n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_y",       (*mc)->position().Y(),                n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_z",       (*mc)->position().Z(),                n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_t",       (*mc)->time(),                        n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_e",       (*mc)->energy(),                      n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_p",       (*mc)->momentum().Mag(),              n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_theta",   (*mc)->momentum().Theta()*180.0/M_PI, n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_phi",     (*mc)->momentum().Phi()*180.0/M_PI,   n_mc_thrown);
		
		n_mc_thrown++;
	}
	dTreeFillData.Fill_Single<Int_t>("nmc", n_mc_thrown);
	
	if(mc_reaction.size() == 1) {
		dTreeFillData.Fill_Single<Int_t>("mc_reaction_type", mc_reaction[0]->type);
		dTreeFillData.Fill_Single<Double_t>("mc_reaction_weight", mc_reaction[0]->weight);
		dTreeFillData.Fill_Single<Double_t>("mc_reaction_energy", mc_reaction[0]->beam.energy());
	}
	
	return;
}

double JEventProcessor_compton_tree::get_acc_scaling_factor(double eb)
{
	if(eb > m_TAGMEnergyBoundHi)
		return m_HodoscopeHiFactor;
	else if(eb > m_TAGMEnergyBoundLo)
		return m_MicroscopeFactor;
	else
		return m_HodoscopeLoFactor;
}
