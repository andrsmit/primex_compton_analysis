// $Id$
//
//    File: JEventProcessor_compton_analysis.cc
// Created: Mon Nov 23 16:24:15 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_compton_analysis.h"


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
	
	
	hTrig   = new TH1F( "hTrig",   "L1 Trigger Bits",          33, -0.5, 32.5 );
	hfpTrig = new TH1F( "hfpTrig", "Front Panel Trigger Bits", 33, -0.5, 32.5 );
	
	
	h_fcal_rf_dt = new TH1F( "fcal_rf_dt", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100. );
	h_ccal_rf_dt = new TH1F( "ccal_rf_dt", "t_{CCAL} - t_{RF}; [ns]", 2000, -100., 100. );
	h_beam_rf_dt = new TH1F( "beam_rf_dt", "t_{Beam} - t_{RF}; [ns]", 2000, -100., 100. );
	
	
	
	TDirectory *dir_deltaPhi = new TDirectoryFile( "DeltaPhi", "DeltaPhi" );
	dir_deltaPhi->cd();
	
	h_deltaPhi[0]     = new TH2F( "deltaPhi_tagh",     
		"#Delta#phi; TAGH Counter; [deg]", 
		274, 0.5, 274.5, 3600, 0., 360. );
	h_deltaPhi_e[0]   = new TH2F( "deltaPhi_tagh_e", 
		"#Delta#phi (#DeltaE Cut); TAGH Counter; [deg]", 
		274, 0.5, 274.5, 3600, 0., 360. );
	h_deltaPhi_ek[0]  = new TH2F( "deltaPhi_tagh_ek", 
		"#Delta#phi (#DeltaE + #DeltaK Cut); TAGH Counter; [deg]", 
		274, 0.5, 274.5, 3600, 0., 360. );
	h_deltaPhi_ekp[0] = new TH2F( "deltaPhi_tagh_ekp", 
		"#Delta#phi (#DeltaE + #DeltaK + #Delta#phi Cut); TAGH Counter; [deg]", 
		274, 0.5, 274.5, 3600, 0., 360. );
	
	h_deltaPhi[1]     = new TH2F( "deltaPhi_tagm", 
		"#Delta#phi; TAGM Counter; [deg]", 
		102, 0.5, 102.5, 3600, 0., 360. );
	h_deltaPhi_e[1]   = new TH2F( "deltaPhi_tagm_e", 
		"#Delta#phi (#DeltaE Cut); TAGM Counter; [deg]", 
		102, 0.5, 102.5, 3600, 0., 360. );
	h_deltaPhi_ek[1]  = new TH2F( "deltaPhi_tagm_ek", 
		"#Delta#phi (#DeltaE + #DeltaK Cut); TAGM Counter; [deg]", 
		102, 0.5, 102.5, 3600, 0., 360. );
	h_deltaPhi_ekp[1] = new TH2F( "deltaPhi_tagm_ekp", 
		"#Delta#phi (#DeltaE + #DeltaK + #Delta#phi Cut); TAGM Counter; [deg]", 
		102, 0.5, 102.5, 3600, 0., 360. );
	
	dir_deltaPhi->cd("../");
	
	
	TDirectory *dir_deltaT = new TDirectoryFile( "DeltaT", "DeltaT" );
	dir_deltaT->cd();
	
	h_deltaT[0]     = new TH2F( "deltaT_tagh", 
		"#DeltaT; TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -100., 100. );
	h_deltaT_e[0]   = new TH2F( "deltaT_tagh_e", 
		"#DeltaT (#DeltaE Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -100., 100. );
	h_deltaT_ep[0]  = new TH2F( "deltaT_tagh_ep", 
		"#DeltaT (#DeltaE + #Delta#phi Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -100., 100. );
	h_deltaT_epk[0] = new TH2F( "deltaT_tagh_epk", 
		"#DeltaT (#DeltaE + #Delta#phi + #DeltaK Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -100., 100. );
	
	h_deltaT[1]     = new TH2F( "deltaT_tagm", 
		"#DeltaT; TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -100., 100. );
	h_deltaT_e[1]   = new TH2F( "deltaT_tagm_e", 
		"#DeltaT (#DeltaE Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -100., 100. );
	h_deltaT_ep[1]  = new TH2F( "deltaT_tagm_ep", 
		"#DeltaT (#DeltaE + #Delta#phi Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -100., 100. );
	h_deltaT_epk[1] = new TH2F( "deltaT_tagm_epk", 
		"#DeltaT (#DeltaE + #Delta#phi + #DeltaK Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -100., 100. );
	
	dir_deltaT->cd("../");
	
	
	TDirectory *dir_deltaR = new TDirectoryFile( "DeltaR", "DeltaR" );
	dir_deltaR->cd();
	
	h_deltaR[0]     = new TH2F( "deltaR_tagh", 
		"#DeltaR; TAGH Counter; [ns]", 
		274, 0.5, 274.5, 1000, 0., 100. );
	h_deltaR_e[0]   = new TH2F( "deltaR_tagh_e", 
		"#DeltaR (#DeltaE Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 1000, 0., 100. );
	h_deltaR_ep[0]  = new TH2F( "deltaR_tagh_ep", 
		"#DeltaR (#DeltaE + #Delta#phi Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 1000, 0., 100. );
	h_deltaR_epk[0] = new TH2F( "deltaR_tagh_epk", 
		"#DeltaR (#DeltaE + #Delta#phi + #DeltaK Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 1000, 0., 100. );
	
	h_deltaR[1]     = new TH2F( "deltaR_tagm", 
		"#DeltaR; TAGM Counter; [ns]", 
		102, 0.5, 102.5, 1000, 0., 100. );
	h_deltaR_e[1]   = new TH2F( "deltaR_tagm_e", 
		"#DeltaR (#DeltaE Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 1000, 0., 100. );
	h_deltaR_ep[1]  = new TH2F( "deltaR_tagm_ep", 
		"#DeltaR (#DeltaE + #Delta#phi Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 1000, 0., 100. );
	h_deltaR_epk[1] = new TH2F( "deltaR_tagm_epk", 
		"#DeltaR (#DeltaE + #Delta#phi + #DeltaK Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 1000, 0., 100. );
	
	dir_deltaR->cd("../");
	
	
	TDirectory *dir_deltaE = new TDirectoryFile( "DeltaE", "DeltaE" );
	dir_deltaE->cd();
	
	h_deltaE[0]     = new TH2F( "deltaE_tagh", 
		"#DeltaE; TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaE_p[0]   = new TH2F( "deltaE_tagh_p", 
		"#DeltaE (#Delta#phi Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaE_pk[0]  = new TH2F( "deltaE_tagh_pk", 
		"#DeltaE (#Delta#phi + #DeltaK Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaE_pke[0] = new TH2F( "deltaE_tagh_pke", 
		"#DeltaE (#Delta#phi + #DeltaK + #DeltaE Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	
	h_deltaE[1]     = new TH2F( "deltaE_tagm", 
		"#DeltaE; TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaE_p[1]   = new TH2F( "deltaE_tagm_p", 
		"#DeltaE (#Delta#phi Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaE_pk[1]  = new TH2F( "deltaE_tagm_pk", 
		"#DeltaE (#Delta#phi + #DeltaK Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaE_pke[1] = new TH2F( "deltaE_tagm_pke", 
		"#DeltaE (#Delta#phi + #DeltaK + #DeltaE Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	
	dir_deltaE->cd("../");
	
	
	TDirectory *dir_deltaK = new TDirectoryFile( "DeltaK", "DeltaK" );
	dir_deltaK->cd();
	
	h_deltaK[0]     = new TH2F( "deltaK_tagh", 
		"#DeltaK; TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaK_e[0]   = new TH2F( "deltaK_tagh_e", 
		"#DeltaK (#DeltaE Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaK_ep[0]  = new TH2F( "deltaK_tagh_ep", 
		"#DeltaK (#DeltaE + #Delta#phi Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaK_epk[0] = new TH2F( "deltaK_tagh_epk", 
		"#DeltaK (#DeltaE + #Delta#phi + #DeltaK Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	
	h_deltaK[1]     = new TH2F( "deltaK_tagm", 
		"#DeltaK; TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaK_e[1]   = new TH2F( "deltaK_tagm_e", 
		"#DeltaK (#DeltaE Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaK_ep[1]  = new TH2F( "deltaK_tagm_ep", 
		"#DeltaK (#DeltaE + #Delta#phi Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaK_epk[1] = new TH2F( "deltaK_tagm_epk", 
		"#DeltaK (#DeltaE + #Delta#phi + #DeltaK Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	
	dir_deltaK->cd("../");
	
	
	TDirectory *dir_deltaK2 = new TDirectoryFile( "DeltaK2", "DeltaK2" );
	dir_deltaK2->cd();
	
	h_deltaK2[0]     = new TH2F( "deltaK2_tagh", 
		"#DeltaK2; TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0 );
	h_deltaK2_e[0]   = new TH2F( "deltaK2_tagh_e", 
		"#DeltaK2 (#DeltaE Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0 );
	h_deltaK2_ep[0]  = new TH2F( "deltaK2_tagh_ep", 
		"#DeltaK2 (#DeltaE + #Delta#phi Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0 );
	h_deltaK2_epk[0] = new TH2F( "deltaK2_tagh_epk", 
		"#DeltaK2 (#DeltaE + #Delta#phi + #DeltaK Cut); TAGH Counter; [ns]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0 );
	
	h_deltaK2[1]     = new TH2F( "deltaK2_tagm", 
		"#DeltaK2; TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0 );
	h_deltaK2_e[1]   = new TH2F( "deltaK2_tagm_e", 
		"#DeltaK2 (#DeltaE Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0 );
	h_deltaK2_ep[1]  = new TH2F( "deltaK2_tagm_ep", 
		"#DeltaK2 (#DeltaE + #Delta#phi Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0 );
	h_deltaK2_epk[1] = new TH2F( "deltaK2_tagm_epk", 
		"#DeltaK2 (#DeltaE + #Delta#phi + #DeltaK Cut); TAGM Counter; [ns]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0 );
	
	dir_deltaK2->cd("../");
	
	
	
	h_fcal_xy = new TH2F( "fcal_xy", 
		"FCAL Shower Position; x_{FCAL} [cm]; y_{FCAL} [cm]", 
		1000, -60., 60., 1000, -60., 60. );
	h_ccal_xy = new TH2F( "ccal_xy", 
		"CCAL Shower Position; x_{CCAL} [cm]; y_{CCAL} [cm]", 
		1000, -13., 13., 1000, -13., 13. );
	
	
	
	
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_compton_analysis::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	
	DGeometry*   dgeom = NULL;
  	DApplication* dapp = dynamic_cast< DApplication* >( eventLoop->GetJApplication() );
  	if( dapp )   dgeom = dapp->GetDGeometry( runnumber );
   	
	if( dgeom ){
    	  	dgeom->GetTargetZ( m_beamZ );
		dgeom->GetFCALPosition( m_fcalX, m_fcalY, m_fcalZ );
		dgeom->GetCCALPosition( m_ccalX, m_ccalY, m_ccalZ );
  	} else{
    	  	cerr << "No geometry accessbile to compton_analysis plugin." << endl;
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
	
	
	set_cuts( runnumber );
	
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_compton_analysis::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
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
	/*
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
	*/
	
	
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
		
		face_x -= m_fcalX;
		face_y -= m_fcalY;
		
		double inner_layer_cut = 2.5*4.0157;
		
		int fid_cut = 0;
		if( (-1.*inner_layer_cut < face_x && face_x < inner_layer_cut) 
			&& (-1.*inner_layer_cut < face_y && face_y < inner_layer_cut) ) fid_cut = 1;
		
		if( (fabs(loc_t) < FCAL_RF_time_cut) && ((*show)->getEnergy() > FCAL_min_energy_cut) 
			&& !fid_cut && (loc_theta < 4.0) ) {
		  	n_fcal_showers++;
		  	fcal_candidates.push_back( (*show) );
		}
		
		h_fcal_rf_dt->Fill( loc_t );
	}
	
	
	int n_ccal_showers = 0;	
	for( vector< const DCCALShower* >::const_iterator show = ccal_showers.begin();
		show != ccal_showers.end(); show++ ) {
		
		/* 
		In the simulation, the CCAL position was shifted to (x,y) = (m_ccalX_new,m_ccalY_new)
		because this was the position found using the offline alignment procedure.
		
		In the reconstruction software (DCCALShower_factory), however, the showers are 
		applied the (now incorrect) offset of x1 -> x1 + m_ccalX && y1 -> y1 + m_ccalY.
		
		Our goal with this bit of code is to apply a fiducial cut on the CCAL using the 
		x,y position that the electrons/photons entered the CCAL.
		
		In my own build of halld_recon, I changed the DCCALShower_factory to give x,y,z 
		inside the calorimeter instead of projecting it to the surface. 
		So what we need to do hear is apply the correct offset, then project to 
		the face of the CCAL.
		*/
		
		double loc_x = (*show)->x1  -  vertex.X() -  m_ccalX  +  m_ccalX_new;
		double loc_y = (*show)->y1  -  vertex.Y() -  m_ccalY  +  m_ccalY_new;
		double loc_z = (*show)->z   -  vertex.Z();
		
		/* 
		Now, loc_x, loc_y, and loc_z give x,y,z correctly in a coordinate system
		where the target center (assumed interaction vertex) is the origin.
		Now, we project the coordinates to the face of the CCAL.
		*/
		
		double xface = vertex.X()  +  (loc_x * (m_ccalZ - vertex.Z())/loc_z);
		double yface = vertex.Y()  +  (loc_y * (m_ccalZ - vertex.Z())/loc_z);
		
		/*
		Right now, xface and yface are the x,y coordinates of the CCAL shower
		projected to the z-position of the CCAL face. However, they are still
		in the coordinate system where the origin is the target center.
		For more straightforward application of the fiducial cut, let's convert
		them to the coordinate system where the center of CCAL is (0,0):
		*/
		
		xface -= m_ccalX_new;
		yface -= m_ccalY_new;
		
		/*
		The fiducial cut removes the first inner layer of the CCAL, the most outer
		layer of the CCAL, and one additional column in the negative x direction
		that is shadowed by the FCAL:
		*/
		
		int fid_cut = 0;
		if( (-4.18 < xface && xface < 4.18) && (-4.18 < yface && yface < 4.18) ) fid_cut = 1;
		if(  xface < -8.36 || xface > 10.45 || yface < -10.45 || yface > 10.45 ) fid_cut = 1;
		
		
		
		double loc_r = sqrt( loc_x*loc_x  +  loc_y*loc_y  +  loc_z*loc_z );
		double loc_t = (*show)->time - (loc_r/c) - rfTime;
		
		if( (fabs(loc_t) < CCAL_RF_time_cut) && ((*show)->E > CCAL_min_energy_cut) 
			&& !fid_cut ) {
			n_ccal_showers++;
			ccal_candidates.push_back( (*show) );
		}
		
		h_ccal_rf_dt->Fill( loc_t );
	}
	
	
	for( vector< const DBeamPhoton* >::const_iterator gam = beam_photons.begin();
		gam != beam_photons.end(); gam++ ) {
		
		double loc_t = (*gam)->time() - rfTime;
		h_beam_rf_dt->Fill( loc_t );
	}
	
	
	//----------     Check FCAL-CCAL Pairs     ----------//
	
	vector< ComptonCandidate_t > candidates;
	
	for( vector< const DFCALShower* >::const_iterator show1 = fcal_candidates.begin(); 
		show1 != fcal_candidates.end(); show1++ ) {
		
		
		double e1         =  (*show1)->getEnergy();
		DVector3 pos1     =  (*show1)->getPosition_log()  -  vertex;
		
		if( (-32. < pos1.Y() && pos1.Y() < -20.) && (-8. < pos1.X() && pos1.X() < 4.) )
			continue;
		
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
			
			/*
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
			*/
			
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
				
				//loc_Cand.phi_slice = loc_phi_slice;
				
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


//------------------
// read in cut values
//------------------
void JEventProcessor_compton_analysis::set_cuts( int32_t runnumber )
{
	
	if( runnumber < 61355 ) { // Be Target (200 nA)
		
		deltaE_mu_p0    =  6.35736e-02;
		deltaE_mu_p1    = -1.51078e-02;
		deltaE_mu_p2    =  3.15066e-04;
		deltaE_mu_p3    =  8.28735e-06;
		
		deltaE_sig_p0   =  9.48484e-03;
		deltaE_sig_p1   =  2.49121e-02;
		deltaE_sig_p2   =  3.84999e-03;
		
		
		deltaPhi_mu_p0  =  1.78907e+02;
		deltaPhi_mu_p1  =  3.59216e-01;
		deltaPhi_mu_p2  = -4.88098e-02;
		deltaPhi_mu_p3  =  2.03322e-03;
		
		deltaPhi_sig_p0 =  1.20540e+01;
		deltaPhi_sig_p1 = -1.41635e+00;
		deltaPhi_sig_p2 =  1.15506e-01;
		deltaPhi_sig_p3 = -3.23715e-03;
		
				
		deltaK_mu_p0    =  5.14497e-02;
		deltaK_mu_p1    =  8.53598e-03;
		deltaK_mu_p2    = -4.43561e-03;
		deltaK_mu_p3    =  2.15182e-04;
		
		deltaK_sig_p0   =  9.33422e-02;
		deltaK_sig_p1   = -5.89212e-04;
		deltaK_sig_p2   =  3.18354e-03;
		deltaK_sig_p3   = -1.60047e-04;
		
		
	} else if( runnumber < 61911 ) { // He Target (200 nA)
		
		deltaE_mu_p0    =  1.25368e-01;
		deltaE_mu_p1    = -4.08206e-02;
		deltaE_mu_p2    =  3.62775e-03;
		deltaE_mu_p3    = -1.36018e-04;
		
		deltaE_sig_p0   =  1.59421e-02;
		deltaE_sig_p1   =  4.19103e-12;
		deltaE_sig_p2   =  3.32391e-02;
		
		
		deltaPhi_mu_p0  =  1.80422e+02;
		deltaPhi_mu_p1  = -2.70566e-01;
		deltaPhi_mu_p2  =  2.96910e-02;
		deltaPhi_mu_p3  = -1.14780e-03;
		
		deltaPhi_sig_p0 =  1.05340e+01;
		deltaPhi_sig_p1 = -1.13232e+00;
		deltaPhi_sig_p2 =  9.35420e-02;
		deltaPhi_sig_p3 = -2.69885e-03;
		
		
		deltaK_mu_p0    =  6.88062e-02;
		deltaK_mu_p1    =  6.88754e-04;
		deltaK_mu_p2    = -3.16666e-03;
		deltaK_mu_p3    =  1.50777e-04;
		
		deltaK_sig_p0   =  9.95112e-02;
		deltaK_sig_p1   = -3.22801e-05;
		deltaK_sig_p2   =  2.82362e-03;
		deltaK_sig_p3   = -1.35114e-04;
		
	} else if( runnumber < 61940 ) { // He Target (50 nA)
		
		deltaE_mu_p0    =  6.61443e-02;
		deltaE_mu_p1    = -2.01740e-02;
		deltaE_mu_p2    =  3.06065e-04;
		deltaE_mu_p3    =  1.92091e-05;
		
		deltaE_sig_p0   =  8.20626e-03;
		deltaE_sig_p1   =  3.30092e-02;
		deltaE_sig_p2   =  7.09972e-10;
		
		
		deltaPhi_mu_p0  =  1.79917e+02;
		deltaPhi_mu_p1  = -4.91278e-02;
		deltaPhi_mu_p2  = -1.95718e-04;
		deltaPhi_mu_p3  =  1.75895e-04;
		
		deltaPhi_sig_p0 =  1.12987e+01;
		deltaPhi_sig_p1 = -1.41973e+00;
		deltaPhi_sig_p2 =  1.27952e-01;
		deltaPhi_sig_p3 = -4.04766e-03;
		
		
		deltaK_mu_p0    =  5.52678e-02;
		deltaK_mu_p1    =  6.75644e-03;
		deltaK_mu_p2    = -4.00935e-03;
		deltaK_mu_p3    =  1.88147e-04;
		
		deltaK_sig_p0   =  2.29331e-01;
		deltaK_sig_p1   = -4.96555e-02;
		deltaK_sig_p2   =  9.01915e-03;
		deltaK_sig_p3   = -3.91743e-04;
		
	} else { // He Target (100 nA)
		
		deltaE_mu_p0    =  1.96300e-01;
		deltaE_mu_p1    = -6.91539e-02;
		deltaE_mu_p2    =  7.11088e-03;
		deltaE_mu_p3    = -2.69863e-04;
		
		deltaE_sig_p0   =  1.61072e-02;
		deltaE_sig_p1   =  3.19930e-12;
		deltaE_sig_p2   =  3.20052e-02;
		
		
		deltaPhi_mu_p0  =  1.80227e+02;
		deltaPhi_mu_p1  = -2.00610e-01;
		deltaPhi_mu_p2  =  2.07285e-02;
		deltaPhi_mu_p3  = -7.50372e-04;
		
		deltaPhi_sig_p0 =  1.38342e+01;
		deltaPhi_sig_p1 = -2.34462e+00;
		deltaPhi_sig_p2 =  2.38851e-01;
		deltaPhi_sig_p3 = -8.45802e-03;
		
		
		deltaK_mu_p0    =  2.35848e-02;
		deltaK_mu_p1    =  1.66183e-02;
		deltaK_mu_p2    = -5.16862e-03;
		deltaK_mu_p3    =  2.37393e-04;
		
		deltaK_sig_p0   =  2.60072e-01;
		deltaK_sig_p1   = -5.69117e-02;
		deltaK_sig_p2   =  9.42118e-03;
		deltaK_sig_p3   = -3.89880e-04;
		
	}
	
	
	f_deltaE_mu = new TF1(   "f_deltaE_mu", "pol3", 3.0, 12.0 );
	f_deltaE_mu->SetParameters( deltaE_mu_p0, deltaE_mu_p1, deltaE_mu_p2, deltaE_mu_p3 );
	
	f_deltaE_sig = new TF1( "f_deltaE_sig", "[0] + [1]/sqrt(x) + [2]/x", 3.0, 12.0 );
	f_deltaE_sig->SetParameters( deltaE_sig_p0, deltaE_sig_p1, deltaE_sig_p2 );
	
	
	
	f_deltaPhi_mu = new TF1(   "f_deltaPhi_mu", "pol3", 3.0, 12.0 );
	f_deltaPhi_mu->SetParameters(   deltaPhi_mu_p0,  deltaPhi_mu_p1,  deltaPhi_mu_p2, 
		deltaPhi_mu_p3 );
	
	f_deltaPhi_sig = new TF1( "f_deltaPhi_sig", "pol3", 3.0, 12.0 );
	f_deltaPhi_sig->SetParameters( deltaPhi_sig_p0, deltaPhi_sig_p1, deltaPhi_sig_p2, 
		deltaPhi_sig_p3 );
	
	
	
	f_deltaK_mu = new TF1(   "f_deltaK_mu", "pol3", 3.0, 12.0 );
	f_deltaK_mu->SetParameters(   deltaK_mu_p0,  deltaK_mu_p1,  deltaK_mu_p2,  deltaK_mu_p3 );
	
	f_deltaK_sig = new TF1( "f_deltaK_sig", "pol3", 3.0, 12.0 );
	f_deltaK_sig->SetParameters( deltaK_sig_p0, deltaK_sig_p1, deltaK_sig_p2, deltaK_sig_p3 );
	
	
	return;
}


//--------------------------------------------
// Get CCAL layer from row and column
//--------------------------------------------
int JEventProcessor_compton_analysis::ccalLayer(int row, int col) {
	
	
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
int JEventProcessor_compton_analysis::fcalLayer(int row, int col) {
	
	
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
void JEventProcessor_compton_analysis::fill_histograms( vector< ComptonCandidate_t > Comp_Candidates ) {
	
	
	
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
		
		deltaE_mu  = f_deltaE_mu->Eval(eb);
		deltaE_sig = eb * f_deltaE_sig->Eval(eb);
		
		deltaPhi_mu   = f_deltaPhi_mu->Eval(eb);
		deltaPhi_sig  = f_deltaPhi_sig->Eval(eb);
		
		deltaK_mu     = f_deltaK_mu->Eval(eb);
		deltaK_sig    = f_deltaK_sig->Eval(eb);
		
		/*
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
		*/
		
		int e5_cut = 0;
		if( fabs(deltaE  -  deltaE_mu) < 5.0*  deltaE_sig ) e5_cut = 1;
		
		int p5_cut = 0;
		if( fabs(deltaPhi-deltaPhi_mu) < 5.0*deltaPhi_sig ) p5_cut = 1;
		
		int k5_cut = 0;
		if( fabs(deltaK  -  deltaK_mu) < 5.0*  deltaK_sig ) k5_cut = 1;
		
		//-------------------------------------------//
		
		
		double fill_weight;
		if( bunch_val ) fill_weight =  1.0;
		else            fill_weight = -0.5;
		
		
		
		h_deltaPhi[tag_sys]->Fill( tag_counter, deltaPhi, fill_weight );
		if( e5_cut ) {
			h_deltaPhi_e[tag_sys]->Fill( tag_counter, deltaPhi, fill_weight );
			if( k5_cut ) {
				h_deltaPhi_ek[tag_sys]->Fill( tag_counter, deltaPhi, fill_weight );
				if( p5_cut ) {
					h_deltaPhi_ekp[tag_sys]->Fill( tag_counter, deltaPhi, fill_weight );
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
				}
			}
		}
		
		h_deltaR[tag_sys]->Fill( tag_counter, deltaR, fill_weight );
		if( e5_cut ) {
			h_deltaR_e[tag_sys]->Fill( tag_counter, deltaR, fill_weight );
			if( p5_cut ) {
				h_deltaR_ep[tag_sys]->Fill( tag_counter, deltaR, fill_weight );
				if( k5_cut ) {
					h_deltaR_epk[tag_sys]->Fill( tag_counter, deltaR, fill_weight );
				}
			}
		}
		
		h_deltaE[tag_sys]->Fill( tag_counter, deltaE, fill_weight );
		if( p5_cut ) {
			h_deltaE_p[tag_sys]->Fill( tag_counter, deltaE, fill_weight );
			if( k5_cut ) {
				h_deltaE_pk[tag_sys]->Fill( tag_counter, deltaE, fill_weight );
				if( e5_cut ) {
					h_deltaE_pke[tag_sys]->Fill( tag_counter, deltaE, fill_weight );
				}
			}
		}
		
		h_deltaK[tag_sys]->Fill( tag_counter, deltaK, fill_weight );
		if( e5_cut ) {
			h_deltaK_e[tag_sys]->Fill( tag_counter, deltaK, fill_weight );
			if( p5_cut ) {
				h_deltaK_ep[tag_sys]->Fill( tag_counter, deltaK, fill_weight );
				if( k5_cut ) {
					h_deltaK_epk[tag_sys]->Fill( tag_counter, deltaK, fill_weight );
				}
			}
		}
		
		h_deltaK2[tag_sys]->Fill( tag_counter, deltaK2, fill_weight );
		if( e5_cut ) {
			h_deltaK2_e[tag_sys]->Fill( tag_counter, deltaK2, fill_weight );
			if( p5_cut ) {
				h_deltaK2_ep[tag_sys]->Fill( tag_counter, deltaK2, fill_weight );
				if( k5_cut ) {
					h_deltaK2_epk[tag_sys]->Fill( tag_counter, deltaK2, fill_weight );
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

