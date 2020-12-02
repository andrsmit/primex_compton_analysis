// $Id$
//
//    File: JEventProcessor_compton_analysis2.cc
// Created: Mon Nov 23 16:24:15 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_compton_analysis2.h"


extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_compton_analysis2());
}
} // "C"



//------------------
// Constructor
//------------------
JEventProcessor_compton_analysis2::JEventProcessor_compton_analysis2()
{
	// Set defaults:
	RUN_GROUP = 0;
	gPARMS->SetDefaultParameter( "COMPTON_ANALYSIS2:RUN_GROUP", RUN_GROUP );
	
	sprintf( cut_pathName, "/work/halld/home/andrsmit/primex_compton_analysis/cuts" );
	
	
	// read in cut values:
	
	read_cuts();
	
	
	f_deltaT_mu  = new TF1( "f_deltaT_mu",  "pol2", BEAM_min_energy_cut, 12.0 );
	f_deltaT_mu->SetParameters(  deltaT_mu_p0,  deltaT_mu_p1,  deltaT_mu_p2 );
	
	f_deltaT_sig = new TF1( "f_deltaT_sig", "pol2", BEAM_min_energy_cut, 12.0 );
	f_deltaT_sig->SetParameters( deltaT_sig_p0, deltaT_sig_p1, deltaT_sig_p2 );
}



//------------------
// init
//------------------
jerror_t JEventProcessor_compton_analysis2::init(void)
{
	
	
	TDirectory *dir_compton_analysis2 = new TDirectoryFile( "compton_analysis2", "compton_analysis2" );
	dir_compton_analysis2->cd();
	
	
	hTrig   = new TH1F( "hTrig",   "L1 Trigger Bits",          33, -0.5, 32.5 );
	hfpTrig = new TH1F( "hfpTrig", "Front Panel Trigger Bits", 33, -0.5, 32.5 );
	
	
	h_fcal_rf_dt = new TH1F( "fcal_rf_dt", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100. );
	h_ccal_rf_dt = new TH1F( "ccal_rf_dt", "t_{CCAL} - t_{RF}; [ns]", 2000, -100., 100. );
	h_beam_rf_dt = new TH1F( "beam_rf_dt", "t_{Beam} - t_{RF}; [ns]", 2000, -100., 100. );
	
	
	
	TDirectory *dir_deltaPhi = new TDirectoryFile( "DeltaPhi", "DeltaPhi" );
	dir_deltaPhi->cd();
	
	h_deltaPhi[0]      = new TH2F( "deltaPhi_tagh",      "#Delta#phi; TAGH Counter; [deg]",                                                274, 0.5, 274.5, 3600, 0., 360. );
	h_deltaPhi_t[0]    = new TH2F( "deltaPhi_tagh_t",    "#Delta#phi (#DeltaT Cut); TAGH Counter; [deg]",                                  274, 0.5, 274.5, 3600, 0., 360. );
	h_deltaPhi_te[0]   = new TH2F( "deltaPhi_tagh_te",   "#Delta#phi (#DeltaT + #DeltaE Cut); TAGH Counter; [deg]",                        274, 0.5, 274.5, 3600, 0., 360. );
	h_deltaPhi_tek[0]  = new TH2F( "deltaPhi_tagh_tek",  "#Delta#phi (#DeltaT + #DeltaE + #DeltaK Cut); TAGH Counter; [deg]",              274, 0.5, 274.5, 3600, 0., 360. );
	h_deltaPhi_tekp[0] = new TH2F( "deltaPhi_tagh_tekp", "#Delta#phi (#DeltaT + #DeltaE + #DeltaK + #Delta#phi Cut); TAGH Counter; [deg]", 274, 0.5, 274.5, 3600, 0., 360. );
	
	h_deltaPhi[1]      = new TH2F( "deltaPhi_tagm",      "#Delta#phi; TAGM Counter; [deg]",                                                102, 0.5, 102.5, 3600, 0., 360. );
	h_deltaPhi_t[1]    = new TH2F( "deltaPhi_tagm_t",    "#Delta#phi (#DeltaT Cut); TAGM Counter; [deg]",                                  102, 0.5, 102.5, 3600, 0., 360. );
	h_deltaPhi_te[1]   = new TH2F( "deltaPhi_tagm_te",   "#Delta#phi (#DeltaT + #DeltaE Cut); TAGM Counter; [deg]",                        102, 0.5, 102.5, 3600, 0., 360. );
	h_deltaPhi_tek[1]  = new TH2F( "deltaPhi_tagm_tek",  "#Delta#phi (#DeltaT + #DeltaE + #DeltaK Cut); TAGM Counter; [deg]",              102, 0.5, 102.5, 3600, 0., 360. );
	h_deltaPhi_tekp[1] = new TH2F( "deltaPhi_tagm_tekp", "#Delta#phi (#DeltaT + #DeltaE + #DeltaK + #Delta#phi Cut); TAGM Counter; [deg]", 102, 0.5, 102.5, 3600, 0., 360. );
	
	dir_deltaPhi->cd("../");
	
	
	TDirectory *dir_deltaT = new TDirectoryFile( "DeltaT", "DeltaT" );
	dir_deltaT->cd();
	
	h_deltaT[0]      = new TH2F( "deltaT_tagh",      "#DeltaT; TAGH Counter; [ns]",                                                274, 0.5, 274.5, 2000, -100., 100. );
	h_deltaT_e[0]    = new TH2F( "deltaT_tagh_e",    "#DeltaT (#DeltaE Cut); TAGH Counter; [ns]",                                  274, 0.5, 274.5, 2000, -100., 100. );
	h_deltaT_ep[0]   = new TH2F( "deltaT_tagh_ep",   "#DeltaT (#DeltaE + #Delta#phi Cut); TAGH Counter; [ns]",                     274, 0.5, 274.5, 2000, -100., 100. );
	h_deltaT_epk[0]  = new TH2F( "deltaT_tagh_epk",  "#DeltaT (#DeltaE + #Delta#phi + #DeltaK Cut); TAGH Counter; [ns]",           274, 0.5, 274.5, 2000, -100., 100. );
	h_deltaT_epkt[0] = new TH2F( "deltaT_tagh_epkt", "#DeltaT (#DeltaE + #Delta#phi + #DeltaK + #DeltaT Cut); TAGH Counter; [ns]", 274, 0.5, 274.5, 2000, -100., 100. );
	
	h_deltaT[1]      = new TH2F( "deltaT_tagm",      "#DeltaT; TAGM Counter; [ns]",                                                102, 0.5, 102.5, 2000, -100., 100. );
	h_deltaT_e[1]    = new TH2F( "deltaT_tagm_e",    "#DeltaT (#DeltaE Cut); TAGM Counter; [ns]",                                  102, 0.5, 102.5, 2000, -100., 100. );
	h_deltaT_ep[1]   = new TH2F( "deltaT_tagm_ep",   "#DeltaT (#DeltaE + #Delta#phi Cut); TAGM Counter; [ns]",                     102, 0.5, 102.5, 2000, -100., 100. );
	h_deltaT_epk[1]  = new TH2F( "deltaT_tagm_epk",  "#DeltaT (#DeltaE + #Delta#phi + #DeltaK Cut); TAGM Counter; [ns]",           102, 0.5, 102.5, 2000, -100., 100. );
	h_deltaT_epkt[1] = new TH2F( "deltaT_tagm_epkt", "#DeltaT (#DeltaE + #Delta#phi + #DeltaK + #DeltaT Cut); TAGM Counter; [ns]", 102, 0.5, 102.5, 2000, -100., 100. );
	
	dir_deltaT->cd("../");
	
	
	TDirectory *dir_deltaR = new TDirectoryFile( "DeltaR", "DeltaR" );
	dir_deltaR->cd();
	
	h_deltaR[0]      = new TH2F( "deltaR_tagh",      "#DeltaR; TAGH Counter; [ns]",                                                274, 0.5, 274.5, 1000, 0., 100. );
	h_deltaR_t[0]    = new TH2F( "deltaR_tagh_t",    "#DeltaR (#DeltaT Cut); TAGH Counter; [ns]",                                  274, 0.5, 274.5, 1000, 0., 100. );
	h_deltaR_te[0]   = new TH2F( "deltaR_tagh_te",   "#DeltaR (#DeltaT + #DeltaE Cut); TAGH Counter; [ns]",                        274, 0.5, 274.5, 1000, 0., 100. );
	h_deltaR_tep[0]  = new TH2F( "deltaR_tagh_tep",  "#DeltaR (#DeltaT + #DeltaE + #Delta#phi Cut); TAGH Counter; [ns]",           274, 0.5, 274.5, 1000, 0., 100. );
	h_deltaR_tepk[0] = new TH2F( "deltaR_tagh_tepk", "#DeltaR (#DeltaT + #DeltaE + #Delta#phi + #DeltaK Cut); TAGH Counter; [ns]", 274, 0.5, 274.5, 1000, 0., 100. );
	
	h_deltaR[1]      = new TH2F( "deltaR_tagm",      "#DeltaR; TAGM Counter; [ns]",                                                102, 0.5, 102.5, 1000, 0., 100. );
	h_deltaR_t[1]    = new TH2F( "deltaR_tagm_t",    "#DeltaR (#DeltaT Cut); TAGM Counter; [ns]",                                  102, 0.5, 102.5, 1000, 0., 100. );
	h_deltaR_te[1]   = new TH2F( "deltaR_tagm_te",   "#DeltaR (#DeltaT + #DeltaE Cut); TAGM Counter; [ns]",                        102, 0.5, 102.5, 1000, 0., 100. );
	h_deltaR_tep[1]  = new TH2F( "deltaR_tagm_tep",  "#DeltaR (#DeltaT + #DeltaE + #Delta#phi Cut); TAGM Counter; [ns]",           102, 0.5, 102.5, 1000, 0., 100. );
	h_deltaR_tepk[1] = new TH2F( "deltaR_tagm_tepk", "#DeltaR (#DeltaT + #DeltaE + #Delta#phi + #DeltaK Cut); TAGM Counter; [ns]", 102, 0.5, 102.5, 1000, 0., 100. );
	
	dir_deltaR->cd("../");
	
	
	TDirectory *dir_deltaE = new TDirectoryFile( "DeltaE", "DeltaE" );
	dir_deltaE->cd();
	
	h_deltaE[0]      = new TH2F( "deltaE_tagh",      "#DeltaE; TAGH Counter; [ns]",                                                274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaE_t[0]    = new TH2F( "deltaE_tagh_t",    "#DeltaE (#DeltaT Cut); TAGH Counter; [ns]",                                  274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaE_tp[0]   = new TH2F( "deltaE_tagh_tp",   "#DeltaE (#DeltaT + #Delta#phi Cut); TAGH Counter; [ns]",                     274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaE_tpk[0]  = new TH2F( "deltaE_tagh_tpk",  "#DeltaE (#DeltaT + #Delta#phi + #DeltaK Cut); TAGH Counter; [ns]",           274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaE_tpke[0] = new TH2F( "deltaE_tagh_tpke", "#DeltaE (#DeltaT + #Delta#phi + #DeltaK + #DeltaE Cut); TAGH Counter; [ns]", 274, 0.5, 274.5, 2000, -4.0, 4.0 );
	
	h_deltaE[1]      = new TH2F( "deltaE_tagm",      "#DeltaE; TAGM Counter; [ns]",                                                102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaE_t[1]    = new TH2F( "deltaE_tagm_t",    "#DeltaE (#DeltaT Cut); TAGM Counter; [ns]",                                  102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaE_tp[1]   = new TH2F( "deltaE_tagm_tp",   "#DeltaE (#DeltaT + #Delta#phi Cut); TAGM Counter; [ns]",                     102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaE_tpk[1]  = new TH2F( "deltaE_tagm_tpk",  "#DeltaE (#DeltaT + #Delta#phi + #DeltaK Cut); TAGM Counter; [ns]",           102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaE_tpke[1] = new TH2F( "deltaE_tagm_tpke", "#DeltaE (#DeltaT + #Delta#phi + #DeltaK + #DeltaE Cut); TAGM Counter; [ns]", 102, 0.5, 102.5, 2000, -4.0, 4.0 );
	
	dir_deltaE->cd("../");
	
	
	TDirectory *dir_deltaK = new TDirectoryFile( "DeltaK", "DeltaK" );
	dir_deltaK->cd();
	
	h_deltaK[0]      = new TH2F( "deltaK_tagh",      "#DeltaK; TAGH Counter; [ns]",                                                274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaK_t[0]    = new TH2F( "deltaK_tagh_t",    "#DeltaK (#DeltaT Cut); TAGH Counter; [ns]",                                  274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaK_te[0]   = new TH2F( "deltaK_tagh_te",   "#DeltaK (#DeltaT + #DeltaE Cut); TAGH Counter; [ns]",                        274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaK_tep[0]  = new TH2F( "deltaK_tagh_tep",  "#DeltaK (#DeltaT + #DeltaE + #Delta#phi Cut); TAGH Counter; [ns]",           274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaK_tepk[0] = new TH2F( "deltaK_tagh_tepk", "#DeltaK (#DeltaT + #DeltaE + #Delta#phi + #DeltaK Cut); TAGH Counter; [ns]", 274, 0.5, 274.5, 2000, -4.0, 4.0 );
	
	h_deltaK[1]      = new TH2F( "deltaK_tagm",      "#DeltaK; TAGM Counter; [ns]",                                                102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaK_t[1]    = new TH2F( "deltaK_tagm_t",    "#DeltaK (#DeltaT Cut); TAGM Counter; [ns]",                                  102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaK_te[1]   = new TH2F( "deltaK_tagm_te",   "#DeltaK (#DeltaT + #DeltaE Cut); TAGM Counter; [ns]",                        102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaK_tep[1]  = new TH2F( "deltaK_tagm_tep",  "#DeltaK (#DeltaT + #DeltaE + #Delta#phi Cut); TAGM Counter; [ns]",           102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaK_tepk[1] = new TH2F( "deltaK_tagm_tepk", "#DeltaK (#DeltaT + #DeltaE + #Delta#phi + #DeltaK Cut); TAGM Counter; [ns]", 102, 0.5, 102.5, 2000, -4.0, 4.0 );
	
	dir_deltaK->cd("../");
	
	
	TDirectory *dir_deltaK2 = new TDirectoryFile( "DeltaK2", "DeltaK2" );
	dir_deltaK2->cd();
	
	h_deltaK2[0]      = new TH2F( "deltaK2_tagh",      "#DeltaK2; TAGH Counter; [ns]",                                                274, 0.5, 274.5, 2000, -8.0, 8.0 );
	h_deltaK2_t[0]    = new TH2F( "deltaK2_tagh_t",    "#DeltaK2 (#DeltaT Cut); TAGH Counter; [ns]",                                  274, 0.5, 274.5, 2000, -8.0, 8.0 );
	h_deltaK2_te[0]   = new TH2F( "deltaK2_tagh_te",   "#DeltaK2 (#DeltaT + #DeltaE Cut); TAGH Counter; [ns]",                        274, 0.5, 274.5, 2000, -8.0, 8.0 );
	h_deltaK2_tep[0]  = new TH2F( "deltaK2_tagh_tep",  "#DeltaK2 (#DeltaT + #DeltaE + #Delta#phi Cut); TAGH Counter; [ns]",           274, 0.5, 274.5, 2000, -8.0, 8.0 );
	h_deltaK2_tepk[0] = new TH2F( "deltaK2_tagh_tepk", "#DeltaK2 (#DeltaT + #DeltaE + #Delta#phi + #DeltaK Cut); TAGH Counter; [ns]", 274, 0.5, 274.5, 2000, -8.0, 8.0 );
	
	h_deltaK2[1]      = new TH2F( "deltaK2_tagm",      "#DeltaK2; TAGM Counter; [ns]",                                                102, 0.5, 102.5, 2000, -8.0, 8.0 );
	h_deltaK2_t[1]    = new TH2F( "deltaK2_tagm_t",    "#DeltaK2 (#DeltaT Cut); TAGM Counter; [ns]",                                  102, 0.5, 102.5, 2000, -8.0, 8.0 );
	h_deltaK2_te[1]   = new TH2F( "deltaK2_tagm_te",   "#DeltaK2 (#DeltaT + #DeltaE Cut); TAGM Counter; [ns]",                        102, 0.5, 102.5, 2000, -8.0, 8.0 );
	h_deltaK2_tep[1]  = new TH2F( "deltaK2_tagm_tep",  "#DeltaK2 (#DeltaT + #DeltaE + #Delta#phi Cut); TAGM Counter; [ns]",           102, 0.5, 102.5, 2000, -8.0, 8.0 );
	h_deltaK2_tepk[1] = new TH2F( "deltaK2_tagm_tepk", "#DeltaK2 (#DeltaT + #DeltaE + #Delta#phi + #DeltaK Cut); TAGM Counter; [ns]", 102, 0.5, 102.5, 2000, -8.0, 8.0 );
	
	dir_deltaK2->cd("../");
	
	
	
	h_fcal_xy = new TH2F( "fcal_xy", "FCAL Shower Position; x_{FCAL} [cm]; y_{FCAL} [cm]", 1000, -60., 60., 1000, -60., 60. );
	h_ccal_xy = new TH2F( "ccal_xy", "CCAL Shower Position; x_{CCAL} [cm]; y_{CCAL} [cm]", 1000, -13., 13., 1000, -13., 13. );
	
	dir_compton_analysis2->cd("../");
	
	
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_compton_analysis2::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	
	DGeometry*   dgeom = NULL;
  	DApplication* dapp = dynamic_cast< DApplication* >( eventLoop->GetJApplication() );
  	if( dapp )   dgeom = dapp->GetDGeometry( runnumber );
   	
	if( dgeom ){
    	  	dgeom->GetTargetZ( m_beamZ );
		dgeom->GetFCALPosition( m_fcalX, m_fcalY, m_fcalZ );
		dgeom->GetCCALPosition( m_ccalX, m_ccalY, m_ccalZ );
  	} else{
    	  	cerr << "No geometry accessbile to compton_analysis2 plugin." << endl;
    	  	return RESOURCE_UNAVAILABLE;
  	}
	
	jana::JCalibration *jcalib = japp->GetJCalibration(runnumber);
  	std::map<string, float> beam_spot;
  	jcalib->Get("PHOTON_BEAM/beam_spot", beam_spot);
  	m_beamX  =  beam_spot.at("x");
  	m_beamY  =  beam_spot.at("y");
	
	
	// Update beam and ccal position offsets:
	
	m_ccalX_new = -0.029;
	m_ccalY_new =  0.094;
	
	m_beamX_new =  0.0784;
	m_beamY_new = -0.0532;
	
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_compton_analysis2::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
	
	//-----   Check Trigger   -----//
	
	const DL1Trigger *trig = NULL;
  	try {
      	  	eventLoop->GetSingle(trig);
  	} catch (...) {}
	if (trig == NULL) { return NOERROR; }
	
	
	uint32_t trigmask    = trig->trig_mask;	
	uint32_t fp_trigmask = trig->fp_trig_mask;
	
	for( int ibit = 0; ibit < 33; ibit++ ) {
	  	if(trigmask & (1 << ibit))    hTrig->Fill(ibit);
	  	if(fp_trigmask & (1 << ibit)) hfpTrig->Fill(ibit);
	}
	
	//if( !trig->Get_IsPhysicsEvent() ) return NOERROR;
	if( trigmask==8 ) return NOERROR;
	if( fp_trigmask ) return NOERROR;
	
	
	//-----   RF Bunch   -----//
	
	const DEventRFBunch *locRFBunch = NULL;
	try { 
	  	eventLoop->GetSingle( locRFBunch, "CalorimeterOnly" );
	} catch (...) { return NOERROR; }
	double rfTime = locRFBunch->dTime;
	if( locRFBunch->dNumParticleVotes < 2 ) return NOERROR;
	
	
	
	//-----   Geometry Objects   -----//
	
	vector< const DFCALGeometry* > fcalGeomVec;
  	vector< const DCCALGeometry* > ccalGeomVec;
  	eventLoop->Get( fcalGeomVec );
	eventLoop->Get( ccalGeomVec );
	
  	if( fcalGeomVec.size() != 1 ) {
    		cerr << "No FCAL geometry accessbile." << endl;
    		return RESOURCE_UNAVAILABLE;
  	}
	
	if( ccalGeomVec.size() != 1 ) {
    		cerr << "No CCAL geometry accessbile." << endl;
    		return RESOURCE_UNAVAILABLE;
  	}
	
	const DFCALGeometry *fcalGeom = fcalGeomVec[0];
	const DCCALGeometry *ccalGeom = ccalGeomVec[0];
	
	
	
	DVector3 vertex;
	vertex.SetXYZ( m_beamX_new, m_beamY_new, m_beamZ );
	
	
	
	//-----   Data Objects   -----//
	
	
	vector< const DBeamPhoton* > beam_photons;
	vector< const DCCALShower* > ccal_showers;
	vector< const DFCALShower* > fcal_showers;
	
	eventLoop->Get( beam_photons );
	eventLoop->Get( ccal_showers );
	eventLoop->Get( fcal_showers );
	
	vector< const DFCALShower* > fcal_candidates;
	vector< const DCCALShower* > ccal_candidates;
	
	
	
	
	
	japp->RootFillLock(this);  // Acquire root lock
	
	
	int n_fcal_showers = 0;	
	for( vector< const DFCALShower* >::const_iterator show = fcal_showers.begin(); 
		show != fcal_showers.end(); show++ ) {
		
		DVector3 loc_pos  =  (*show)->getPosition_log()  -  vertex;
		double loc_theta  =  loc_pos.Theta() * (180. / TMath::Pi());
		double loc_r      =  loc_pos.Mag();
		
		double loc_t      =  (*show)->getTime() - (loc_r/c) - rfTime;
		
		double face_x     =  vertex.X()  +  (loc_pos.X() * (m_fcalZ - vertex.Z())/loc_pos.Z());
		double face_y     =  vertex.Y()  +  (loc_pos.Y() * (m_fcalZ - vertex.Z())/loc_pos.Z());
		
		int row   = fcalGeom->row(    static_cast<float>(face_y) );
		int col   = fcalGeom->column( static_cast<float>(face_x) );
		int layer = fcalLayer( row, col );
		
		if( ((*show)->getEnergy() > FCAL_min_energy_cut) 
			&& (layer > 1) && (loc_theta < 4.0) ) {
		  	n_fcal_showers++;
		  	fcal_candidates.push_back( (*show) );
		}
		
		h_fcal_rf_dt->Fill( loc_t );
	}
	
	
	int n_ccal_showers = 0;	
	for( vector< const DCCALShower* >::const_iterator show = ccal_showers.begin();
		show != ccal_showers.end(); show++ ) {
		
		double loc_x = (*show)->x1  -  vertex.X() - (m_ccalX + m_ccalX_new);
		double loc_y = (*show)->y1  -  vertex.Y() - (m_ccalY + m_ccalY_new);
		double loc_z = (*show)->z   -  vertex.Z();
		
		double loc_r = sqrt( loc_x*loc_x  +  loc_y*loc_y  +  loc_z*loc_z );
		double loc_t = (*show)->time - (loc_r/c) - rfTime;
		
		double xface = vertex.X()  +  (loc_x * (m_ccalZ - vertex.Z())/loc_z);
		double yface = vertex.Y()  +  (loc_y * (m_ccalZ - vertex.Z())/loc_z);
		
		int row   = ccalGeom->row(    static_cast<float>(yface-m_ccalY_new) );
		int col   = ccalGeom->column( static_cast<float>(xface-m_ccalX_new) );
		int layer = ccalLayer( row, col );
		
		if( (fabs(loc_t) < CCAL_RF_time_cut) && ((*show)->E > CCAL_min_energy_cut) 
			&& (layer > 1) && (layer < 5) && (col != 1) ) {
			n_ccal_showers++;
			ccal_candidates.push_back( (*show) );
		}
		
		h_ccal_rf_dt->Fill( loc_t );
	}
	
	
	
	
	//----------     Check FCAL-CCAL Pairs     ----------//
	
	vector< ComptonCandidate_t > candidates;
	
	for( vector< const DFCALShower* >::const_iterator show1 = fcal_candidates.begin(); 
		show1 != fcal_candidates.end(); show1++ ) {
		
		
		double e1         =  (*show1)->getEnergy();
		DVector3 pos1     =  (*show1)->getPosition_log()  -  vertex;
		
		double r1         =  pos1.Mag();
		double t1         =  (*show1)->getTime()  -  (r1/c);
		double phi1       =  pos1.Phi() * (180. / TMath::Pi());
		double theta1     =  pos1.Theta();
		
		
		for( vector< const DCCALShower* >::const_iterator show2 = ccal_candidates.begin(); 
			show2 != ccal_candidates.end(); show2++ ) {
			
			
			double e2     =  (*show2)->E;
			DVector3 pos2( (*show2)->x1, (*show2)->y1, (*show2)->z );
			pos2          =  pos2  -  vertex;
			
			DVector3 ccal_offset( m_ccalX - m_ccalX_new, m_ccalY - m_ccalY_new, 0.0 );
			pos2          =  pos2  -  ccal_offset;
			
			double r2     =  pos2.Mag();
			double t2     =  (*show2)->time  -  (r2 / c);
			double phi2   =  pos2.Phi() * (180. / TMath::Pi());
			double theta2 =  pos2.Theta();
			
			int loc_phi_slice;
			if(      -180. < phi2  &&  phi2 <= -135. ) loc_phi_slice = 0;
			else if( -135. < phi2  &&  phi2 <=  -90. ) loc_phi_slice = 1;
			else if(  -90. < phi2  &&  phi2 <=  -45. ) loc_phi_slice = 2;
			else if(  -45. < phi2  &&  phi2 <=    0. ) loc_phi_slice = 3;
			else if(    0. < phi2  &&  phi2 <=   45. ) loc_phi_slice = 4;
			else if(   45. < phi2  &&  phi2 <=   90. ) loc_phi_slice = 5;
			else if(   90. < phi2  &&  phi2 <=  135. ) loc_phi_slice = 6;
			else if(  135. < phi2  &&  phi2 <=  180. ) loc_phi_slice = 7;
			else loc_phi_slice = 8;
			
			// calculate deltaPhi and deltaT:
			
			double deltaPhi  =  fabs(phi2 - phi1);
			double deltaT    =  t2 - t1;
			
			// calculate separation distance of particles at FCAL z position:
			
			double x1_face   =  ((m_fcalZ - vertex.Z()) / pos1.Z()) * pos1.X();
			double y1_face   =  ((m_fcalZ - vertex.Z()) / pos1.Z()) * pos1.Y();
			double x2_face   =  ((m_fcalZ - vertex.Z()) / pos2.Z()) * pos2.X();
			double y2_face   =  ((m_fcalZ - vertex.Z()) / pos2.Z()) * pos2.Y();
			double deltaR    =  sqrt( pow(x1_face-x2_face,2.0) 
				+ pow(y1_face-y2_face,2.0) );
			
			// loop over beam photons:
			
			for( vector< const DBeamPhoton* >::const_iterator gam = beam_photons.begin();
				gam != beam_photons.end(); gam++ ) {
				
				
				double eb    = (*gam)->lorentzMomentum().E();
				double tb    = (*gam)->time();
				
				double brfdt = tb - rfTime;
				
				int bunch_val;
				
				if( fabs(brfdt) < 2.004 )
					bunch_val = 1;
				else if( (  -(2.004 + 3.*4.008) <= brfdt && brfdt <= -(2.004 + 2.*4.008) )
					|| ( (2.004 + 2.*4.008) <= brfdt && brfdt <=  (2.004 + 3.*4.008) ) )
					bunch_val = 0;
				else 
					continue;
				
				if( eb < BEAM_min_energy_cut ) continue;
				
				double ecomp1   =   1. / ( (1./eb)  +  (1./m_e)*(1. - cos(theta1)) );
				double ecomp2   =   1. / ( (1./eb)  +  (1./m_e)*(1. - cos(theta2)) );
				double deltaE   =  (e1     + e2    ) - (eb + m_e);
				double deltaK   =  (ecomp1 + ecomp2) - (eb + m_e);
				
				double deltaK2  =  m_e * sin(theta1+theta2) / 
					(sin(theta1) + sin(theta2) - sin(theta1+theta2));
				deltaK2        -=  eb;
				
				ComptonCandidate_t loc_Cand;
				
				loc_Cand.bunch_val = bunch_val;
				
				loc_Cand.e1        = e1;
				loc_Cand.x1        = pos1.X();
				loc_Cand.y1        = pos1.Y();
				loc_Cand.z1        = pos1.Z();
				loc_Cand.e2        = e2;
				loc_Cand.x2        = pos2.X();
				loc_Cand.y2        = pos2.Y();
				loc_Cand.z2        = pos2.Z();
				
				loc_Cand.phi_slice = loc_phi_slice;
				
				loc_Cand.deltaPhi  = deltaPhi;
				loc_Cand.deltaT    = deltaT;
				loc_Cand.deltaR    = deltaR;
				loc_Cand.deltaE    = deltaE;
				loc_Cand.deltaK    = deltaK;
				loc_Cand.deltaK2   = deltaK2;
				
				loc_Cand.eb          = eb;
				loc_Cand.tag_counter = (*gam)->dCounter;
				
				DetectorSystem_t sys = (*gam)->dSystem;
				if( sys==SYS_TAGH )      loc_Cand.tag_sys = 0;
				else if( sys==SYS_TAGM ) loc_Cand.tag_sys = 1;
				
				candidates.push_back( loc_Cand );
				
				
			} // end DBeamPhoton loop
			
		} // end DCCALShower loop
		
	} // end DFCALShower loop
	
	
	fill_histograms( candidates );
	
	
	japp->RootFillUnLock(this);  // Release root lock
	
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_compton_analysis2::erun(void)
{
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_compton_analysis2::fini(void)
{
	
	return NOERROR;
}


//------------------
// read in cut values
//------------------
void JEventProcessor_compton_analysis2::read_cuts()
{
	
	for( int i=0; i<274; i++ ) {
		  deltaE_mu_tagh[i] = 0.;   deltaE_sig_tagh[i] = 0.;
		deltaPhi_mu_tagh[i] = 0.; deltaPhi_sig_tagh[i] = 0.;
		  deltaK_mu_tagh[i] = 0.;   deltaK_sig_tagh[i] = 0.;
	}
	for( int i=0; i<102; i++ ) {
		  deltaE_mu_tagm[i] = 0.;   deltaE_sig_tagm[i] = 0.;
		deltaPhi_mu_tagm[i] = 0.; deltaPhi_sig_tagm[i] = 0.;
		  deltaK_mu_tagm[i] = 0.;   deltaK_sig_tagm[i] = 0.;
	}
	
	char run_group_name[256];
	
	if( RUN_GROUP==0 ) {
		sprintf( run_group_name, "%s/Be",    cut_pathName );
	} else if( RUN_GROUP==1 ) {
		sprintf( run_group_name, "%s/He50",  cut_pathName );
	} else if( RUN_GROUP==2 ) {
		sprintf( run_group_name, "%s/He100", cut_pathName );
	} else {
		sprintf( run_group_name, "%s/Be",    cut_pathName );
	}
	
	char fname[256];
	
	
	//-----   DeltaE   -----//
	
	sprintf( fname, "%s/deltaE_tagh.dat", run_group_name );
	ifstream infile_deltaE_tagh(fname);
	if( infile_deltaE_tagh.good() ) {
		int a; double b, c;
		for( int i=0; i<274; i++ ) {
			infile_deltaE_tagh >> a >> b >> c;
			deltaE_mu_tagh[i]  = b;
			deltaE_sig_tagh[i] = c;
		}
	}
	infile_deltaE_tagh.close();
	
	sprintf( fname, "%s/deltaE_tagm.dat", run_group_name );
	ifstream infile_deltaE_tagm(fname);
	if( infile_deltaE_tagm.good() ) {
		int a; double b, c;
		for( int i=0; i<102; i++ ) {
			infile_deltaE_tagm >> a >> b >> c;
			deltaE_mu_tagm[i]  = b;
			deltaE_sig_tagm[i] = c;
		}
	}
	infile_deltaE_tagm.close();
	
	
	
	//-----   DeltaPhi   -----//
	
	sprintf( fname, "%s/deltaPhi_tagh.dat", run_group_name );
	ifstream infile_deltaPhi_tagh(fname);
	if( infile_deltaPhi_tagh.good() ) {
		int a; double b, c;
		for( int i=0; i<274; i++ ) {
			infile_deltaPhi_tagh >> a >> b >> c;
			deltaPhi_mu_tagh[i]  = b;
			deltaPhi_sig_tagh[i] = c;
		}
	}
	infile_deltaPhi_tagh.close();
	
	sprintf( fname, "%s/deltaPhi_tagm.dat", run_group_name );
	ifstream infile_deltaPhi_tagm(fname);
	if( infile_deltaPhi_tagm.good() ) {
		int a; double b, c;
		for( int i=0; i<102; i++ ) {
			infile_deltaPhi_tagm >> a >> b >> c;
			deltaPhi_mu_tagm[i]  = b;
			deltaPhi_sig_tagm[i] = c;
		}
	}
	infile_deltaPhi_tagm.close();
	
	
	
	//-----   DeltaK   -----//
	
	sprintf( fname, "%s/deltaK_tagh.dat", run_group_name );
	ifstream infile_deltaK_tagh(fname);
	if( infile_deltaK_tagh.good() ) {
		int a; double b, c;
		for( int i=0; i<274; i++ ) {
			infile_deltaK_tagh >> a >> b >> c;
			deltaK_mu_tagh[i]  = b;
			deltaK_sig_tagh[i] = c;
		}
	}
	infile_deltaK_tagh.close();
	
	sprintf( fname, "%s/deltaK_tagm.dat", run_group_name );
	ifstream infile_deltaK_tagm(fname);
	if( infile_deltaK_tagm.good() ) {
		int a; double b, c;
		for( int i=0; i<102; i++ ) {
			infile_deltaK_tagm >> a >> b >> c;
			deltaK_mu_tagm[i]  = b;
			deltaK_sig_tagm[i] = c;
		}
	}
	infile_deltaK_tagm.close();
	
	
	return;
}


//--------------------------------------------
// Get CCAL layer from row and column
//--------------------------------------------
int JEventProcessor_compton_analysis2::ccalLayer(int row, int col) {
	
	
	int layer;
	
	if( (row > 3 && row < 8) && (col > 3 && col < 8) ) layer = 1;
	else if( (row > 2 && row <  9) && (col > 2 && col <  9) ) layer = 2;
	else if( (row > 1 && row < 10) && (col > 1 && col < 10) ) layer = 3;
	else if( (row > 0 && row < 11) && (col > 0 && col < 11) ) layer = 4;
	else layer = 5;	
	
	return layer;
}


//--------------------------------------------
// Get FCAL layer from row and column
//--------------------------------------------
int JEventProcessor_compton_analysis2::fcalLayer(int row, int col) {
	
	
	int layer;
	
	if( (row > 26 && row < 32) && (col > 26 && col < 32) ) layer = 1;
	else if( (row > 25 && row < 33) && (col > 25 && col < 33) ) layer =  2;
	else if( (row > 24 && row < 34) && (col > 24 && col < 34) ) layer =  3;
	else if( (row > 23 && row < 35) && (col > 23 && col < 35) ) layer =  4;
	else if( (row > 22 && row < 36) && (col > 22 && col < 36) ) layer =  5;
	else if( (row > 21 && row < 37) && (col > 21 && col < 37) ) layer =  6;
	else if( (row > 20 && row < 38) && (col > 20 && col < 38) ) layer =  7;
	else if( (row > 19 && row < 39) && (col > 19 && col < 39) ) layer =  8;
	else if( (row > 18 && row < 40) && (col > 18 && col < 40) ) layer =  9;
	else if( (row > 17 && row < 41) && (col > 17 && col < 41) ) layer = 10;
	else layer = 11;
	
	
	return layer;
}



//--------------------------------------------
// Fill Histograms
//--------------------------------------------
void JEventProcessor_compton_analysis2::fill_histograms( vector< ComptonCandidate_t > Comp_Candidates ) {
	
	
	
	int n_candidates  =  static_cast<int>( Comp_Candidates.size() );
	
	for( int ic = 0; ic < n_candidates; ic++ ) {
				
		ComptonCandidate_t loc_Cand = Comp_Candidates[ic];
				
		//-------------------------------------------//
		
		int bunch_val        =  loc_Cand.bunch_val;
		
		double eb            =  loc_Cand.eb;
		int tag_sys          =  loc_Cand.tag_sys;
		int tag_counter      =  loc_Cand.tag_counter;
		
		double deltaPhi      =  loc_Cand.deltaPhi;
		double deltaT        =  loc_Cand.deltaT;
		double deltaR        =  loc_Cand.deltaR;
		double deltaE        =  loc_Cand.deltaE;
		double deltaK        =  loc_Cand.deltaK;
		double deltaK2       =  loc_Cand.deltaK2;
		
		double x1            =  loc_Cand.x1;
		double y1            =  loc_Cand.y1;
		double x2            =  loc_Cand.x2;
		double y2            =  loc_Cand.y2;
		
		//--------------     Cuts      --------------//
		
		double   deltaE_mu,   deltaE_sig;
		double deltaPhi_mu, deltaPhi_sig;
		double   deltaK_mu,   deltaK_sig;
		
		if( tag_sys==0 ) {
			
			deltaE_mu    = deltaE_mu_tagh[tag_counter-1];
			deltaE_sig   = deltaE_sig_tagh[tag_counter-1];
			
			deltaPhi_mu  = deltaPhi_mu_tagh[tag_counter-1];
			deltaPhi_sig = deltaPhi_sig_tagh[tag_counter-1];
			
			deltaK_mu    = deltaK_mu_tagh[tag_counter-1];
			deltaK_sig   = deltaK_sig_tagh[tag_counter-1];
			
		} else {
			
			deltaE_mu    = deltaE_mu_tagm[tag_counter-1];
			deltaE_sig   = deltaE_sig_tagm[tag_counter-1];
			
			deltaPhi_mu  = deltaPhi_mu_tagm[tag_counter-1];
			deltaPhi_sig = deltaPhi_sig_tagm[tag_counter-1];
			
			deltaK_mu    = deltaK_mu_tagm[tag_counter-1];
			deltaK_sig   = deltaK_sig_tagm[tag_counter-1];
			
		}
		
		double deltaT_mu    = f_deltaT_mu->Eval(eb);
		double deltaT_sig   = f_deltaT_sig->Eval(eb);
		
		int e5_cut = 0;
		if( fabs(deltaE  -  deltaE_mu) < 5.0*  deltaE_sig ) e5_cut = 1;
		
		int p5_cut = 0;
		if( fabs(deltaPhi-deltaPhi_mu) < 5.0*deltaPhi_sig ) p5_cut = 1;
		
		int k5_cut = 0;
		if( fabs(deltaK  -  deltaK_mu) < 5.0*  deltaK_sig ) k5_cut = 1;
		
		int t5_cut = 0;
		if( fabs(deltaT  -  deltaT_mu) < 5.0*  deltaT_sig ) t5_cut = 1;
		
		//-------------------------------------------//
		
		
		double fill_weight;
		if( bunch_val ) fill_weight =  1.0;
		else            fill_weight = -0.5;
		
		
		
		h_deltaPhi[tag_sys]->Fill( tag_counter, deltaPhi, fill_weight );
		if( t5_cut ) {
			h_deltaPhi_t[tag_sys]->Fill( tag_counter, deltaPhi, fill_weight );
			if( e5_cut ) {
				h_deltaPhi_te[tag_sys]->Fill( tag_counter, deltaPhi, fill_weight );
				if( k5_cut ) {
					h_deltaPhi_tek[tag_sys]->Fill( tag_counter, deltaPhi, fill_weight );
					if( p5_cut ) {
						h_deltaPhi_tekp[tag_sys]->Fill( tag_counter, deltaPhi, fill_weight );
					}
				}
			}
		}
		
		h_deltaT[tag_sys]->Fill( tag_counter, deltaT, fill_weight );
		if( e5_cut ) {
			h_deltaT_e[tag_sys]->Fill( tag_counter, deltaT, fill_weight );
			if( p5_cut ) {
				h_deltaT_ep[tag_sys]->Fill( tag_counter, deltaT, fill_weight );
				if( k5_cut ) {
					h_deltaT_epk[tag_sys]->Fill( tag_counter, deltaT, fill_weight );
					if( t5_cut ) {
						h_deltaT_epkt[tag_sys]->Fill( tag_counter, deltaT, fill_weight );
					}
				}
			}
		}
		
		h_deltaR[tag_sys]->Fill( tag_counter, deltaR, fill_weight );
		if( t5_cut ) {
			h_deltaR_t[tag_sys]->Fill( tag_counter, deltaR, fill_weight );
			if( e5_cut ) {
				h_deltaR_te[tag_sys]->Fill( tag_counter, deltaR, fill_weight );
				if( p5_cut ) {
					h_deltaR_tep[tag_sys]->Fill( tag_counter, deltaR, fill_weight );
					if( k5_cut ) {
						h_deltaR_tepk[tag_sys]->Fill( tag_counter, deltaR, fill_weight );
					}
				}
			}
		}
		
		h_deltaE[tag_sys]->Fill( tag_counter, deltaE, fill_weight );
		if( t5_cut ) {
			h_deltaE_t[tag_sys]->Fill( tag_counter, deltaE, fill_weight );
			if( p5_cut ) {
				h_deltaE_tp[tag_sys]->Fill( tag_counter, deltaE, fill_weight );
				if( k5_cut ) {
					h_deltaE_tpk[tag_sys]->Fill( tag_counter, deltaE, fill_weight );
					if( e5_cut ) {
						h_deltaE_tpke[tag_sys]->Fill( tag_counter, deltaE, fill_weight );
					}
				}
			}
		}
		
		h_deltaK[tag_sys]->Fill( tag_counter, deltaK, fill_weight );
		if( t5_cut ) {
			h_deltaK_t[tag_sys]->Fill( tag_counter, deltaK, fill_weight );
			if( e5_cut ) {
				h_deltaK_te[tag_sys]->Fill( tag_counter, deltaK, fill_weight );
				if( p5_cut ) {
					h_deltaK_tep[tag_sys]->Fill( tag_counter, deltaK, fill_weight );
					if( k5_cut ) {
						h_deltaK_tepk[tag_sys]->Fill( tag_counter, deltaK, fill_weight );
					}
				}
			}
		}
		
		h_deltaK2[tag_sys]->Fill( tag_counter, deltaK2, fill_weight );
		if( t5_cut ) {
			h_deltaK2_t[tag_sys]->Fill( tag_counter, deltaK2, fill_weight );
			if( e5_cut ) {
				h_deltaK2_te[tag_sys]->Fill( tag_counter, deltaK2, fill_weight );
				if( p5_cut ) {
					h_deltaK2_tep[tag_sys]->Fill( tag_counter, deltaK2, fill_weight );
					if( k5_cut ) {
						h_deltaK2_tepk[tag_sys]->Fill( tag_counter, deltaK2, fill_weight );
					}
				}
			}
		}
		
		
		
		if( e5_cut && p5_cut && k5_cut ) {
			
			h_fcal_xy->Fill( x1,y1 );
			h_ccal_xy->Fill( x2,y2 );
			
		}
		
		
	}
	
	
	
	
	
	
	return;
}

