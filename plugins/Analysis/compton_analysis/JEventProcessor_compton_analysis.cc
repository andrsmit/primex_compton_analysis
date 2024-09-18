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
// JEventProcessor_compton_analysis (Constructor)
//------------------
JEventProcessor_compton_analysis::JEventProcessor_compton_analysis()
{
	// default values for the RF timing cuts for each sub-detector:
	m_cut_fcalrfdt = 2.004;
	m_cut_ccalrfdt = 2.004;
	m_cut_beamrfdt = 2.004;
	
	// default values for the minimum energy cuts:
	m_cut_fcalE = 0.35;
	m_cut_ccalE = 3.00;
	m_cut_beamE = 6.00;
	
	// default values for the widths of the Compton cuts:
	m_cut_deltaE   = 5.0;
	m_cut_deltaPhi = 5.0;
	m_cut_deltaK   = 5.0;
	
	// how many inner layers of the FCAL to cut:
	m_cut_fcal_layers = 1.0;
	
	// how many beam bunches to use for accidental subtraction:
	m_beam_bunches_acc = 5;
	
	//-------------------------------------------------------------------------------------//
	// allow for command-line overriding of the default values:
	
	gPARMS->SetDefaultParameter("compton_analysis:cut_fcalrfdt",     m_cut_fcalrfdt);
	gPARMS->SetDefaultParameter("compton_analysis:cut_ccalrfdt",     m_cut_ccalrfdt);
	gPARMS->SetDefaultParameter("compton_analysis:cut_beamrfdt",     m_cut_beamrfdt);
	gPARMS->SetDefaultParameter("compton_analysis:cut_fcalE",        m_cut_fcalE);
	gPARMS->SetDefaultParameter("compton_analysis:cut_ccalE",        m_cut_ccalE);
	gPARMS->SetDefaultParameter("compton_analysis:cut_beamE",        m_cut_beamE);
	gPARMS->SetDefaultParameter("compton_analysis:cut_fcal_layers",  m_cut_fcal_layers);
	gPARMS->SetDefaultParameter("compton_analysis:cut_deltaE",       m_cut_deltaE);
	gPARMS->SetDefaultParameter("compton_analysis:cut_deltaPhi",     m_cut_deltaPhi);
	gPARMS->SetDefaultParameter("compton_analysis:cut_deltaK",       m_cut_deltaK);
	gPARMS->SetDefaultParameter("compton_analysis:beam_bunches_acc", m_beam_bunches_acc);
}

//------------------
// init
//------------------
jerror_t JEventProcessor_compton_analysis::init(void)
{
	TDirectory *dir_compton = new TDirectoryFile("compton_analysis", "compton_analysis");
	dir_compton->cd();
	
	h_trig   = new TH1F("trig",   "GTP Trigger Bits", 34, -0.5, 33.5);
	h_fptrig = new TH1F("fp_trig", "FP Trigger Bits", 34, -0.5, 33.5);
	
	TDirectory *dir_timing = new TDirectoryFile("timing", "timing");
	dir_timing->cd();
	
	h_fcal_rf_dt     = new TH1F("fcal_rf_dt", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	h_ccal_rf_dt     = new TH1F("ccal_rf_dt", "t_{CCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	h_beam_rf_dt     = new TH1F("beam_rf_dt", "t_{Beam} - t_{RF}; [ns]", 2000, -100., 100.);
	h_beam_rf_dt_cut = new TH1F("beam_rf_dt_cut", "t_{Beam} - t_{RF} (Compton Cuts); [ns]", 
		2000, -100., 100.);
	h_beam_rf_dt_tagh = new TH2F("beam_rf_dt_tagh", "t_{Beam} - t_{RF}; TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -100., 100.);
	h_beam_rf_dt_tagm = new TH2F("beam_rf_dt_tagm", "t_{Beam} - t_{RF}; TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -100., 100.);
	
	dir_timing->cd("../");
	
	//--------------------------------------------------------------------//
	
	TDirectory *dir_multiplicity = new TDirectoryFile("multiplicity", "multiplicity");
	dir_multiplicity->cd();
	
	h_nccal = new TH1F("nccal", "Number of CCAL Showers", 10, -0.5, 9.5);
	h_nfcal = new TH1F("nfcal", "Number of FCAL Showers", 10, -0.5, 9.5);
	
	h_n_ccal_fcal     = new TH1F("n_ccal_fcal",     "Is there a shower in both CCAL+FCAL", 2, -0.5, 1.5);
	h_n_ccal_fcal_cut = new TH1F("n_ccal_fcal_cut", "Is there a shower in both CCAL+FCAL", 2, -0.5, 1.5);
	
	h_ccalE     = new TH2F("ccalE", "Energy of CCAL Shower; E_{Beam} [GeV]; E_{CCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	h_fcalE     = new TH2F("fcalE", "Energy of FCAL Shower; E_{Beam} [GeV]; E_{FCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	h_ccalE_cut = new TH2F("ccalE_cut", "Energy of CCAL Shower; E_{Beam} [GeV]; E_{CCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	h_fcalE_cut = new TH2F("fcalE_cut", "Energy of FCAL Shower; E_{Beam} [GeV]; E_{FCAL} [GeV]", 
		120, 0., 12., 1200, 0., 12.);
	
	dir_multiplicity->cd("../");
	
	//--------------------------------------------------------------------//
	
	h_deltaK_vs_deltaE = new TH2F("deltaK_vs_deltaE", 
		"#DeltaK vs. #DeltaE", 500, -8.0, 8.0, 500, -8.0, 8.0);
	h_deltaK_vs_deltaE->GetXaxis()->SetTitle("#left(E_{1}+E_{2}#right) - #left(E_{#gamma}+m_{e}#right) [GeV]");
	h_deltaK_vs_deltaE->GetYaxis()->SetTitle("E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]");
	
	h_deltaE_tagh = new TH2F("deltaE_tagh", 
		"#DeltaE; TAGH Counter; #left(E_{1}+E_{2}#right) - #left(E_{#gamma}+m_{e}#right) [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0);
	h_deltaE_tagm = new TH2F("deltaE_tagm", 
		"#DeltaE; TAGM Counter; #left(E_{1}+E_{2}#right) - #left(E_{#gamma}+m_{e}#right) [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0);
	h_deltaE_tagh->Sumw2();
	h_deltaE_tagm->Sumw2();
	
	h_deltaPhi_tagh = new TH2F("deltaPhi_tagh", 
		"#Delta#phi; TAGH Counter; #left|#phi_{1}-#phi_{2}#right| [deg.]", 
		274, 0.5, 274.5, 3600, 0.0, 360.0);
	h_deltaPhi_tagm = new TH2F("deltaPhi_tagm", 
		"#Delta#phi; TAGM Counter; #left|#phi_{1}-#phi_{2}#right| [deg.]", 
		102, 0.5, 102.5, 3600, 0.0, 360.0);
	h_deltaPhi_tagh->Sumw2();
	h_deltaPhi_tagm->Sumw2();
	
	h_deltaK_tagh = new TH2F("deltaK_tagh", 
		"#DeltaK; TAGH Counter; E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm = new TH2F("deltaK_tagm", 
		"#DeltaK; TAGM Counter; E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	h_deltaK_tagh->Sumw2();
	h_deltaK_tagm->Sumw2();
	
	h_deltaK_tagh_cut = new TH2F("deltaK_tagh_cut", 
		"#DeltaK; TAGH Counter; E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm_cut = new TH2F("deltaK_tagm_cut", 
		"#DeltaK; TAGM Counter; E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	h_deltaK_tagh_cut->Sumw2();
	h_deltaK_tagm_cut->Sumw2();
	
	h_deltaK_tagh_main = new TH2F("deltaK_tagh_main", 
		"#DeltaK; TAGH Counter; E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm_main = new TH2F("deltaK_tagm_main", 
		"#DeltaK; TAGM Counter; E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	h_deltaK_tagh_main->Sumw2();
	h_deltaK_tagm_main->Sumw2();
	
	h_deltaK_tagh_acc = new TH2F("deltaK_tagh_acc", 
		"#DeltaK; TAGH Counter; E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm_acc = new TH2F("deltaK_tagm_acc", 
		"#DeltaK; TAGM Counter; E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	h_deltaK_tagh_acc->Sumw2();
	h_deltaK_tagm_acc->Sumw2();
	
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
	
	double loc_targetZ = 65.0;
	
	if(dgeom){
		
		dgeom->GetTargetZ(loc_targetZ);
		
		double loc_x, loc_y, loc_z;
		
		dgeom->GetFCALPosition(loc_x, loc_y, loc_z);
		m_fcalFace.SetXYZ(loc_x, loc_y, loc_z);
		
		dgeom->GetCCALPosition(loc_x, loc_y, loc_z);
		m_ccalFace.SetXYZ(loc_x, loc_y, loc_z);
		
	} else{
		cerr << "No geometry accessbile to PrimExComptonAnalysis plugin." << endl;
		return RESOURCE_UNAVAILABLE;
	}
	
	jana::JCalibration *jcalib = japp->GetJCalibration(runnumber);
	std::map<string, float> beam_spot;
	jcalib->Get("PHOTON_BEAM/beam_spot", beam_spot);
	m_beamSpot.SetXYZ(beam_spot.at("x"), beam_spot.at("y"), loc_targetZ);
	
	if(runnumber>60000 && runnumber<69999) {
		
		m_phase_val  = 1;
		m_bfield_val = 0;
		
		if(runnumber<61355) {
			m_Target         = Be9;
			m_target_length  = 1.7755;
			m_target_density = 1.85;
			m_atten          = 0.01172;
		} else {
			m_Target         = Helium;
			m_target_length  = 29.535;
			m_target_density = 0.1217;
			m_atten          = 0.00821;
		}
		
		//--------------------------------------------------------------------//
		// For phase 1, update the geometry from Compton alignment studies:
		
		// (2/4/2024): Correction to alignment after Simon updated beam spot with new CDC alignment:
		
		m_fcal_correction.SetXYZ(0.455 - m_fcalFace.X(), -0.032 - m_fcalFace.Y(), 0.0);
		m_fcalFace += m_fcal_correction;
		
		m_ccal_correction.SetXYZ(-0.082 - m_ccalFace.X(), 0.061 - m_ccalFace.Y(), 0.0);
		if(runnumber>=61483) m_ccal_correction.SetY(0.051 - m_ccalFace.Y());
		m_ccalFace += m_ccal_correction;
		
		if(runnumber<61483) 
		{
			m_beamSpot.SetX( 0.027);
			m_beamSpot.SetY(-0.128);
		} else if(runnumber<61774) 
		{
			m_beamSpot.SetX( 0.001);
			m_beamSpot.SetY(-0.077);
		} else 
		{
			m_beamSpot.SetX( 0.038);
			m_beamSpot.SetY(-0.095);
		}
		//--------------------------------------------------------------------//
		
	} else if(runnumber>80000 && runnumber<89999) {
		
		m_phase_val = 2;
		
		if(runnumber<81384) {
			m_Target         = Be9;
			m_target_length  = 1.7755;
			m_target_density = 1.85;
			m_atten          = 0.01172;
		} else {
			m_Target         = Helium;
			m_target_length  = 29.535;
			m_target_density = 0.1217;
			m_atten          = 0.00821;
		}
		
		if(runnumber<81471) m_bfield_val = 0;
		else m_bfield_val = 1;
		
		m_fcal_correction.SetXYZ(0.0, 0.0, 0.0);
		m_ccal_correction.SetXYZ(0.0, 0.0, 0.0);
		
	} else if(runnumber>110000 && runnumber<119999) {
		
		m_phase_val = 3;
		
		if(runnumber<110622) {
			m_Target         = Be9;
			m_target_length  = 1.7755;
			m_target_density = 1.85;
			m_atten          = 0.01172;
		} else {
			m_Target         = Helium;
			m_target_length  = 29.535;
			m_target_density = 0.1217;
			m_atten          = 0.00821;
		}
		
		if(
			((110476<=runnumber) && (runnumber<=110482)) || // Be empty
			((110600<=runnumber) && (runnumber<=110621)) || // Be target
			((111969<=runnumber) && (runnumber<=112001))    // He target
		)    m_bfield_val = 0;
		else m_bfield_val = 1;
		
		m_fcal_correction.SetXYZ(0.0, 0.0, 0.0);
		m_ccal_correction.SetXYZ(0.0, 0.0, 0.0);
		
	} else {
		
		// non-PrimEx runs - use defaults:
		
		m_phase_val = 0;
		
		// these are for Helium, we should really change the default to Proton, but for now it doesn't matter:
		m_Target         = Helium;
		m_target_length  = 29.535;
		m_target_density = 0.1217;
		m_atten          = 0.00821;
		
		// assume solenoid is on:
		m_bfield_val = 1;
		
		m_fcal_correction.SetXYZ(0.0, 0.0, 0.0);
		m_ccal_correction.SetXYZ(0.0, 0.0, 0.0);
	}
	
	/*------------------------------------------------------------------------------------------------------*/
	// Code to obtain the scaling factors for accidental beam bunches 
	//     (copied from DAnalysisUtilities.cc in gluex_root_analysis)
	
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
	
	/*------------------------------------------------------------------------------------------------------*/
	
	set_cuts(runnumber);
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_compton_analysis::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
	//-----------------------------------------------------//
	// Get RF Time
	
	const DEventRFBunch *locRFBunch = NULL;
	try {
		eventLoop->GetSingle(locRFBunch, "CalorimeterOnly");
	} catch(...) {return NOERROR;}
	double locRFTime = locRFBunch->dTime;
	
	//-----------------------------------------------------//
	// Data objects
	
	vector<const DBeamPhoton*> locBeamPhotons;
	eventLoop->Get(locBeamPhotons);
	
	vector<const DCCALShower*> locCCALShowers;
	eventLoop->Get(locCCALShowers);
	
	vector<const DFCALShower*> locFCALShowers;
	eventLoop->Get(locFCALShowers);
	
	// filtered vectors of objects used for analysis:
	
	vector<pair<double,const DBeamPhoton*>> locGoodBeamPhotons;
	vector<const DCCALShower*> locGoodCCALShowers;
	vector<const DFCALShower*> locGoodFCALShowers;
	
	//-----------------------------------------------------//
	// Trigger information
	
	uint32_t trigmask, fp_trigmask;
	const DL1Trigger *trig = NULL;
	try {
		eventLoop->GetSingle(trig);
	} catch (...) {}
	if (trig == NULL) { 
		trigmask    = 0;
		fp_trigmask = 0;
	} else {
		trigmask    = trig->trig_mask;
		fp_trigmask = trig->fp_trig_mask;
	}
	
	//-----------------------------------------------------//
	// Apply fill lock for multi-threaded running:
	
	japp->RootFillLock(this);
	
	//-----------------------------------------------------//
	// Fill trigger histograms:
	
	for(int ibit = 0; ibit < 34; ibit++) {
		if(trigmask & (1 << ibit))    h_trig->Fill(ibit);
		if(fp_trigmask & (1 << ibit)) h_fptrig->Fill(ibit);
	}
	
	//-----------------------------------------------------//
	// Fill timing histograms:
	
	int locNBeamPhotons = 0;
	for(vector<const DBeamPhoton*>::const_iterator gam = locBeamPhotons.begin();
		gam != locBeamPhotons.end(); gam++) {
		
		double eb = (*gam)->lorentzMomentum().E();
		if(eb < m_cut_beamE) continue;
		
		int counter = (*gam)->dCounter;
		double loc_t = (*gam)->time() - locRFTime;
		
		h_beam_rf_dt->Fill(loc_t);
		
		DetectorSystem_t sys = (*gam)->dSystem;
		if(sys==SYS_TAGH) h_beam_rf_dt_tagh->Fill(counter, loc_t);
		if(sys==SYS_TAGM) h_beam_rf_dt_tagm->Fill(counter, loc_t);
		
		double loc_weight = 0.0;
		if(fabs(loc_t) < m_cut_beamrfdt) loc_weight = 1.0;
		else if(
			(6.5*4.008 < fabs(loc_t)) &&
			(fabs(loc_t) < (6.5+(double)m_beam_bunches_acc)*4.008)
		) loc_weight = -1.0*get_acc_scaling_factor(eb)/(2.0*(double)m_beam_bunches_acc);
		else continue;
		
		locNBeamPhotons++;
		locGoodBeamPhotons.push_back({loc_weight, (*gam)});
	}
	
	int locNFCALShowers = 0, locNGoodFCALShowers = 0;	
	for(vector<const DFCALShower*>::const_iterator show = locFCALShowers.begin(); 
		show != locFCALShowers.end(); show++) {
		
		DVector3 loc_pos = (*show)->getPosition() - m_beamSpot + m_fcal_correction;
		double loc_t     = (*show)->getTime() - (loc_pos.Mag()/m_c) - locRFTime;
		
		int loc_fid_cut  = fcal_fiducial_cut(loc_pos);
		
		if(fabs(loc_t) < m_cut_fcalrfdt) {
			locNFCALShowers++;
			if((*show)->getEnergy() > m_cut_fcalE) {
				locNGoodFCALShowers++;
				if(!loc_fid_cut) {
					locGoodFCALShowers.push_back((*show));
				}
			}
		}
		h_fcal_rf_dt->Fill(loc_t);
	}
	
	int locNCCALShowers = 0, locNGoodCCALShowers = 0;
	for(vector<const DCCALShower*>::const_iterator show = locCCALShowers.begin();
		show != locCCALShowers.end(); show++) {
		
		DVector3 loc_pos((*show)->x1, (*show)->y1, (*show)->z);
		loc_pos = loc_pos - m_beamSpot + m_ccal_correction;
		double loc_t = (*show)->time - (loc_pos.Mag()/m_c) - locRFTime;
		
		int loc_fid_cut = ccal_fiducial_cut(loc_pos);
		
		if(fabs(loc_t) < m_cut_ccalrfdt) {
			locNCCALShowers++;
			if((*show)->E > m_cut_ccalE) {
				locNGoodCCALShowers++;
				if(!loc_fid_cut) {
					locGoodCCALShowers.push_back((*show));
				}
			}
		}
		h_ccal_rf_dt->Fill(loc_t);
	}
	
	//-----------------------------------------------------//
	// Fill multiplicity histograms:
	
	h_nccal->Fill(locNCCALShowers);
	h_nfcal->Fill(locNFCALShowers);
	
	if((locNCCALShowers>0) && (locNFCALShowers>0))
		h_n_ccal_fcal->Fill(1.);
	else 
		h_n_ccal_fcal->Fill(0.);
	
	if((locNGoodCCALShowers>0) && (locNGoodFCALShowers>0))
		h_n_ccal_fcal_cut->Fill(1.);
	else 
		h_n_ccal_fcal_cut->Fill(0.);
	
	if(locRFBunch->dNumParticleVotes < 2) {
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	
	//----------     Check FCAL-CCAL Pairs     ----------//
	
	vector<ComptonCandidate_t> locCandidates;
	
	for(vector<const DFCALShower*>::const_iterator show1 = locGoodFCALShowers.begin(); 
		show1 != locGoodFCALShowers.end(); show1++) {
		
		double e1     = (*show1)->getEnergy();
		DVector3 pos1 = (*show1)->getPosition() - m_beamSpot + m_fcal_correction;
		
		double phi1   = pos1.Phi() * (180./TMath::Pi());
		double theta1 = pos1.Theta();
		
		for(vector<const DCCALShower*>::const_iterator show2 = locGoodCCALShowers.begin(); 
			show2 != locGoodCCALShowers.end(); show2++) {
			
			double e2  = (*show2)->E;
			DVector3 pos2( (*show2)->x1, (*show2)->y1, (*show2)->z );
			pos2       = pos2 - m_beamSpot + m_ccal_correction;
			
			double phi2   = pos2.Phi() * (180./TMath::Pi());
			double theta2 = pos2.Theta();
			
			// calculate deltaPhi and deltaT:
			
			double deltaPhi = fabs(phi2 - phi1);
			
			// loop over beam photons:
			
			for(auto gam_it = locGoodBeamPhotons.begin(); gam_it != locGoodBeamPhotons.end(); gam_it++) {
				
				double fill_weight     = gam_it->first;
				const DBeamPhoton *gam = gam_it->second;
				
				double eb = gam->lorentzMomentum().E();
				double tb = gam->time();
				
				double brfdt = tb - locRFTime;
				
				double deltaE = (e1 + e2) - (eb + m_e);
				double deltaK = m_e*(sin(theta1) / (2.*pow(sin(theta1/2.),2.)*tan(theta2))) - m_e - eb;
				
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
				loc_Cand.deltaE      = deltaE;
				loc_Cand.deltaK      = deltaK;
				
				loc_Cand.weight      = fill_weight;
				loc_Cand.eb          = eb;
				loc_Cand.brfdt       = brfdt;
				loc_Cand.tag_counter = gam->dCounter;
				
				DetectorSystem_t sys = gam->dSystem;
				if(sys==SYS_TAGH)      loc_Cand.tag_sys = 0;
				else if(sys==SYS_TAGM) loc_Cand.tag_sys = 1;
				
				locCandidates.push_back(loc_Cand);
				
			} // end DBeamPhoton loop
		} // end DCCALShower loop
	} // end DFCALShower loop
	
	fill_histograms(locCandidates);
	
	japp->RootFillUnLock(this);  // Release root lock
	
	return NOERROR;
}


int JEventProcessor_compton_analysis::fcal_fiducial_cut(DVector3 pos)
{
	int fid_cut = 0;
	
	double fcal_inner_layer_cut = (1.5 + m_cut_fcal_layers) * 4.0157;
	
	double fcal_face_x = m_beamSpot.X() + (pos.X() * (m_fcalFace.Z() - m_beamSpot.Z())/pos.Z());
	double fcal_face_y = m_beamSpot.Y() + (pos.Y() * (m_fcalFace.Z() - m_beamSpot.Z())/pos.Z());
	
	fcal_face_x -= m_fcalFace.X();
	fcal_face_y -= m_fcalFace.Y();
	
	if((-1.*fcal_inner_layer_cut < fcal_face_x && fcal_face_x < fcal_inner_layer_cut)
		&& (-1.*fcal_inner_layer_cut < fcal_face_y 
		&& fcal_face_y < fcal_inner_layer_cut)) fid_cut = 1;
	
	// only apply the next fiducial cut for runs from phase-I:
	
	if(m_phase_val==1) {
		if((-32.<fcal_face_y && fcal_face_y<-20.) && (-8.<fcal_face_x && fcal_face_x<4.)) fid_cut = 1;
	}
	
	return fid_cut;
}


int JEventProcessor_compton_analysis::ccal_fiducial_cut(DVector3 pos)
{
	int fid_cut = 0;
	
	double ccal_inner_layer_cut = 2.0 * 2.09;
	
	double ccal_face_x = m_beamSpot.X() + (pos.X() * (m_ccalFace.Z() - m_beamSpot.Z())/pos.Z());
	double ccal_face_y = m_beamSpot.Y() + (pos.Y() * (m_ccalFace.Z() - m_beamSpot.Z())/pos.Z());
	
	ccal_face_x -= m_ccalFace.X();
	ccal_face_y -= m_ccalFace.Y();
	
	if((-1.*ccal_inner_layer_cut < ccal_face_x && ccal_face_x < ccal_inner_layer_cut)
		&& (-1.*ccal_inner_layer_cut < ccal_face_y 
		&& ccal_face_y < ccal_inner_layer_cut)) fid_cut = 1;
	
	if(ccal_face_x<-8.36 || ccal_face_x>10.45 
		|| ccal_face_y<-10.45 || ccal_face_y>10.45) fid_cut = 1;
	
	return fid_cut;
}


double JEventProcessor_compton_analysis::get_acc_scaling_factor(double eb)
{
	if(eb > m_TAGMEnergyBoundHi)
		return m_HodoscopeHiFactor;
	else if(eb > m_TAGMEnergyBoundLo)
		return m_MicroscopeFactor;
	else
		return m_HodoscopeLoFactor;
}


void JEventProcessor_compton_analysis::fill_histograms(vector<ComptonCandidate_t> Comp_Candidates) 
{
	int n_candidates = static_cast<int>(Comp_Candidates.size());
	
	for(int ic = 0; ic < n_candidates; ic++) {
		
		ComptonCandidate_t loc_Cand = Comp_Candidates[ic];
		
		//-------------------------------------------//
		
		double weight   = loc_Cand.weight;
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
		
		int   e_cut = cut_deltaE(  deltaE,   eb, m_cut_deltaE,   1.e2);
		int phi_cut = cut_deltaPhi(deltaPhi, eb, m_cut_deltaPhi, m_cut_deltaPhi);
		int   k_cut = cut_deltaK(  deltaK,   eb, m_cut_deltaK,   m_cut_deltaK);
		
		//-------------------------------------------//
		
		h_fcalE->Fill(eb, e1, weight);
		h_ccalE->Fill(eb, e2, weight);
		if(e_cut && phi_cut && k_cut) {
			h_fcalE_cut->Fill(eb, e1, weight);
			h_ccalE_cut->Fill(eb, e2, weight);
			h_fcal_xy->Fill(x1, y1, weight);
			h_ccal_xy->Fill(x2, y2, weight);
			h_beam_rf_dt_cut->Fill(brfdt);
		}
		
		h_deltaK_vs_deltaE->Fill(deltaE, deltaK, weight);
		
		if(tag_sys==0) {
			h_deltaE_tagh->Fill(tag_counter, deltaE, weight);
			if(e_cut) {
				h_deltaPhi_tagh->Fill(tag_counter, deltaPhi, weight);
				if(phi_cut) {
					h_deltaK_tagh->Fill(tag_counter, deltaK, weight);
					if(weight > 0.0) h_deltaK_tagh_main->Fill(tag_counter, deltaK);
					else h_deltaK_tagh_acc->Fill(tag_counter, deltaK, -1.0*weight);
					if(k_cut) {
						h_deltaK_tagh_cut->Fill(tag_counter, deltaK, weight);
					}
				}
			}
		} else {
			h_deltaE_tagm->Fill(tag_counter, deltaE, weight);
			if(e_cut) {
				h_deltaPhi_tagm->Fill(tag_counter, deltaPhi, weight);
				if(phi_cut) {
					h_deltaK_tagm->Fill(tag_counter, deltaK, weight);
					if(weight > 0.0) h_deltaK_tagh_main->Fill(tag_counter, deltaK);
					else h_deltaK_tagm_acc->Fill(tag_counter, deltaK, -1.0*weight);
					if(k_cut) {
						h_deltaK_tagm_cut->Fill(tag_counter, deltaK, weight);
					}
				}
			}
		}
	}
	
	return;
}


int JEventProcessor_compton_analysis::cut_deltaE(double deltaE, double eb, double n_sigma_left, double n_sigma_right) 
{
	double loc_mu = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_mu += (m_deltaE_mu_pars[ipar]*pow(eb,(double)ipar));
	
	double loc_sigma = sqrt(pow(m_deltaE_sigma_pars[0],2.0) + pow(m_deltaE_sigma_pars[1]/sqrt(eb),2.0) 
		+ pow(m_deltaE_sigma_pars[2]/eb,2.0));
	loc_sigma *= eb;
	
	double loc_diff = deltaE - loc_mu;
	if((-1.*n_sigma_left*loc_sigma < loc_diff) && (loc_diff < n_sigma_right*loc_sigma)) return 1;
	else return 0;
}


int JEventProcessor_compton_analysis::cut_deltaPhi(double deltaPhi, double eb, double n_sigma_left, double n_sigma_right) 
{
	if(m_bfield_val==1) {
		if(fabs(deltaPhi-180.) < 50.) return 1;
		else return 0;
	}
	
	double loc_mu = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_mu += (m_deltaPhi_mu_pars[ipar]*pow(eb,(double)ipar));
	
	double loc_sigma = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_sigma += (m_deltaPhi_sigma_pars[ipar]*pow(eb,(double)ipar));
	
	double loc_diff = deltaPhi - loc_mu;
	if((-1.*n_sigma_left*loc_sigma < loc_diff) && (loc_diff < n_sigma_right*loc_sigma)) return 1;
	else return 0;
}


int JEventProcessor_compton_analysis::cut_deltaK(double deltaK, double eb, double n_sigma_left, double n_sigma_right) 
{
	double loc_mu = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_mu += (m_deltaE_mu_pars[ipar]*pow(eb,(double)ipar));
	
	double loc_sigma = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_sigma += (m_deltaK_sigma_pars[ipar]*pow(eb,(double)ipar));
	
	double loc_diff = deltaK - loc_mu;
	if((-1.*n_sigma_left*loc_sigma < loc_diff) && (loc_diff < n_sigma_right*loc_sigma)) return 1;
	else return 0;
}


void JEventProcessor_compton_analysis::set_cuts(int32_t runnumber)
{
	if(m_Target==Be9) {
		
		m_deltaE_mu_pars[0] =  2.80447e-02;
		m_deltaE_mu_pars[1] = -9.73792e-03;
		m_deltaE_mu_pars[2] =  6.27229e-04;
		m_deltaE_mu_pars[3] = -1.11857e-05;
		
		m_deltaE_sigma_pars[0] = 9.19528e-03;
		m_deltaE_sigma_pars[1] = 4.44202e-02;
		m_deltaE_sigma_pars[2] = 1.70059e-06;
		//----------------------------------//
		m_deltaPhi_mu_pars[0] = 180.0;
		m_deltaPhi_mu_pars[1] = 0.0;
		m_deltaPhi_mu_pars[2] = 0.0;
		m_deltaPhi_mu_pars[3] = 0.0;
		
		m_deltaPhi_sigma_pars[0] = 6.0;
		m_deltaPhi_sigma_pars[1] = 0.0;
		m_deltaPhi_sigma_pars[2] = 0.0;
		m_deltaPhi_sigma_pars[3] = 0.0;
		//----------------------------------//
		m_deltaK_mu_pars[0] =  4.11025e-01;
		m_deltaK_mu_pars[1] = -8.98612e-02;
		m_deltaK_mu_pars[2] =  1.90051e-03;
		m_deltaK_mu_pars[3] =  0.0;
		
		m_deltaK_sigma_pars[0] =  5.05212e-01;
		m_deltaK_sigma_pars[1] = -4.77652e-02;
		m_deltaK_sigma_pars[2] =  1.25422e-02;
		m_deltaK_sigma_pars[3] = -4.75279e-04;
		
	} else if(m_Target==Helium) {
		
		m_deltaE_mu_pars[0] =  3.74024e-02;
		m_deltaE_mu_pars[1] = -1.12178e-02;
		m_deltaE_mu_pars[2] =  5.70029e-04;
		m_deltaE_mu_pars[3] =  5.02307e-07;
		
		m_deltaE_sigma_pars[0] = 1.10937e-02;
		m_deltaE_sigma_pars[1] = 4.09888e-02;
		m_deltaE_sigma_pars[2] = 2.79658e-02;
		//----------------------------------//
		m_deltaPhi_mu_pars[0] = 180.0;
		m_deltaPhi_mu_pars[1] = 0.0;
		m_deltaPhi_mu_pars[2] = 0.0;
		m_deltaPhi_mu_pars[3] = 0.0;
		
		m_deltaPhi_sigma_pars[0] = 6.0;
		m_deltaPhi_sigma_pars[1] = 0.0;
		m_deltaPhi_sigma_pars[2] = 0.0;
		m_deltaPhi_sigma_pars[3] = 0.0;
		//----------------------------------//
		m_deltaK_mu_pars[0] =  4.17999e-01;
		m_deltaK_mu_pars[1] = -8.43377e-02;
		m_deltaK_mu_pars[2] =  1.59661e-03;
		m_deltaK_mu_pars[3] =  0.0;
		
		m_deltaK_sigma_pars[0] =  3.72050e-01;
		m_deltaK_sigma_pars[1] =  6.44563e-03;
		m_deltaK_sigma_pars[2] =  5.72365e-03;
		m_deltaK_sigma_pars[3] = -1.73985e-04;
		
	} else {
		
		// copy Helium target runs from above:
		
		m_deltaE_mu_pars[0] =  3.74024e-02;
		m_deltaE_mu_pars[1] = -1.12178e-02;
		m_deltaE_mu_pars[2] =  5.70029e-04;
		m_deltaE_mu_pars[3] =  5.02307e-07;
		
		m_deltaE_sigma_pars[0] = 1.10937e-02;
		m_deltaE_sigma_pars[1] = 4.09888e-02;
		m_deltaE_sigma_pars[2] = 2.79658e-02;
		//----------------------------------//
		m_deltaPhi_mu_pars[0] = 180.0;
		m_deltaPhi_mu_pars[1] = 0.0;
		m_deltaPhi_mu_pars[2] = 0.0;
		m_deltaPhi_mu_pars[3] = 0.0;
		
		m_deltaPhi_sigma_pars[0] = 6.0;
		m_deltaPhi_sigma_pars[1] = 0.0;
		m_deltaPhi_sigma_pars[2] = 0.0;
		m_deltaPhi_sigma_pars[3] = 0.0;
		//----------------------------------//
		m_deltaK_mu_pars[0] =  4.17999e-01;
		m_deltaK_mu_pars[1] = -8.43377e-02;
		m_deltaK_mu_pars[2] =  1.59661e-03;
		m_deltaK_mu_pars[3] =  0.0;
		
		m_deltaK_sigma_pars[0] =  3.72050e-01;
		m_deltaK_sigma_pars[1] =  6.44563e-03;
		m_deltaK_sigma_pars[2] =  5.72365e-03;
		m_deltaK_sigma_pars[3] = -1.73985e-04;
		
	}
	
	return;
}
