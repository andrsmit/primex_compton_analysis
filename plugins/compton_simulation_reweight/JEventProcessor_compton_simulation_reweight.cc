// $Id$
//
//    File: JEventProcessor_compton_simulation_reweight.cc
// Created: Mon Nov 23 16:24:15 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_compton_simulation_reweight.h"


extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_compton_simulation_reweight());
}
} // "C"



//------------------
// init
//------------------
jerror_t JEventProcessor_compton_simulation_reweight::init(void)
{
	
	
	h_fcal_rf_dt = new TH1F( "fcal_rf_dt", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100. );
	h_ccal_rf_dt = new TH1F( "ccal_rf_dt", "t_{CCAL} - t_{RF}; [ns]", 2000, -100., 100. );
	h_beam_rf_dt = new TH1F( "beam_rf_dt", "t_{Beam} - t_{RF}; [ns]", 2000, -100., 100. );
	
	
	h_beam           = new TH1F( "beam", "Is there a beam photon?", 2, -0.5, 1.5 );
	h_tagh_flux      = new TH1F( "tagh_flux", "TAGH Flux", 274, 0.5, 274.5 );
	h_tagm_flux      = new TH1F( "tagm_flux", "TAGM Flux", 102, 0.5, 102.5 );
	h_double_compton = new TH1F( "double_compton", "Is e2 > 0?", 2, -0.5, 1.5 );
	
	
	h_vertex          = new TH1F( "vertex",          
		"Vertex Z Position (unweighted)", 5000, 0., 100. );
	h_vertex_weighted = new TH1F( "vertex_weighted", 
		"Vertex Z Position (weighted)",   5000, 0., 100. );
	
	h_vertex_cut          = new TH1F( "vertex_cut", 
		"Vertex Z Position (unweighted)", 5000, 0., 100. );
	h_vertex_cut_weighted = new TH1F( "vertex_cut_weighted", 
		"Vertex Z Position (weighted)",   5000, 0., 100. );
	
	
	h_deltaE_tagh = new TH2F( "deltaE_tagh", 
		"#DeltaE; TAGH Counter; E_{1}+E_{2} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaE_tagm = new TH2F( "deltaE_tagm", 
		"#DeltaE; TAGM Counter; E_{1}+E_{2} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaE_tagh_weighted = new TH2F( "deltaE_tagh_weighted", 
		"#DeltaE; TAGH Counter; E_{1}+E_{2} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaE_tagm_weighted = new TH2F( "deltaE_tagm_weighted", 
		"#DeltaE; TAGM Counter; E_{1}+E_{2} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	
	h_deltaPhi_tagh = new TH2F( "deltaPhi_tagh", 
		"#Delta#phi; TAGH Counter; |#phi_{1}-#phi_{2}| [deg.]", 
		274, 0.5, 274.5, 3600, 0.0, 360.0 );
	h_deltaPhi_tagm = new TH2F( "deltaPhi_tagm", 
		"#Delta#phi; TAGM Counter; |#phi_{1}-#phi_{2}| [deg.]", 
		102, 0.5, 102.5, 3600, 0.0, 360.0 );
	h_deltaPhi_tagh_weighted = new TH2F( "deltaPhi_tagh_weighted", 
		"#Delta#phi; TAGH Counter; |#phi_{1}-#phi_{2}| [deg.]", 
		274, 0.5, 274.5, 3600, 0.0, 360.0 );
	h_deltaPhi_tagm_weighted = new TH2F( "deltaPhi_tagm_weighted", 
		"#Delta#phi; TAGM Counter; |#phi_{1}-#phi_{2}| [deg.]", 
		102, 0.5, 102.5, 3600, 0.0, 360.0 );
	
	h_deltaK_tagh = new TH2F( "deltaK_tagh", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaK_tagm = new TH2F( "deltaK_tagm", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaK_tagh_weighted = new TH2F( "deltaK_tagh_weighted", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaK_tagm_weighted = new TH2F( "deltaK_tagm_weighted", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	
	h_deltaK_tagh_cut = new TH2F( "deltaK_tagh_cut", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaK_tagm_cut = new TH2F( "deltaK_tagm_cut", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	h_deltaK_tagh_cut_weighted = new TH2F( "deltaK_tagh_cut_weighted", 
		"#DeltaK; TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0 );
	h_deltaK_tagm_cut_weighted = new TH2F( "deltaK_tagm_cut_weighted", 
		"#DeltaK; TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0 );
	
	h_fcal_xy          = new TH2F( "fcal_xy", 
		"FCAL Y vs. X (unweighed); x_{FCAL} [cm]; y_{FCAL} [cm]", 
		1000, -100., 100., 1000, -100., 100. );
	h_fcal_xy_weighted = new TH2F( "fcal_xy_weighted", 
		"FCAL Y vs. X (weighted); x_{FCAL} [cm]; y_{FCAL} [cm]", 
		1000, -100., 100., 1000, -100., 100. );
	
	h_ccal_xy          = new TH2F( "ccal_xy", 
		"CCAL Y vs. X (unweighed); x_{CCAL} [cm]; y_{CCAL} [cm]", 
		1000, -15., 15., 1000, -15., 15. );
	h_ccal_xy_weighted = new TH2F( "ccal_xy_weighted", 
		"CCAL Y vs. X (weighted); x_{CCAL} [cm]; y_{CCAL} [cm]", 
		1000, -15., 15., 1000, -15., 15. );
	
	
	
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_compton_simulation_reweight::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	
	DGeometry*   dgeom = NULL;
  	DApplication* dapp = dynamic_cast< DApplication* >( eventLoop->GetJApplication() );
  	if( dapp )   dgeom = dapp->GetDGeometry( runnumber );
   	
	if( dgeom ){
    	  	dgeom->GetTargetZ( m_beamZ );
		dgeom->GetFCALPosition( m_fcalX, m_fcalY, m_fcalZ );
		dgeom->GetCCALPosition( m_ccalX, m_ccalY, m_ccalZ );
  	} else{
    	  	cerr << "No geometry accessbile to compton_simulation_reweight plugin." << endl;
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
jerror_t JEventProcessor_compton_simulation_reweight::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
	
	//-----   RF Bunch   -----//
	
	const DEventRFBunch *locRFBunch = NULL;
	try { 
	  	eventLoop->GetSingle( locRFBunch, "CalorimeterOnly" );
	} catch (...) { return NOERROR; }
	double rfTime = locRFBunch->dTime;
	
	
	
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
	vector< const DMCThrown*   >   mc_thrown;
	
	eventLoop->Get( beam_photons );
	eventLoop->Get( ccal_showers );
	eventLoop->Get( fcal_showers );
	eventLoop->Get(   mc_thrown  );
	
	vector< const DFCALShower* > fcal_candidates;
	vector< const DCCALShower* > ccal_candidates;
	
	
	
	
	
	japp->RootFillLock(this);  // Acquire root lock
	
	
	double e_extra_gam = mc_thrown[0]->energy();
	if( e_extra_gam > 0. ) 
		h_double_compton->Fill( 1 );
	else {
		h_double_compton->Fill( 0 );
		//japp->RootFillUnLock(this);
		//return NOERROR;
	}
	
	
	int n_beam_photons = (int)beam_photons.size();
	if( n_beam_photons ) 
		h_beam->Fill( 1 );
	else 
		h_beam->Fill( 0 );
	
	
	
	double vertex_z     = mc_thrown[0]->position().Z();
	double event_weight = 1.;
	
	h_vertex->Fill( vertex_z );
	h_vertex_weighted->Fill( vertex_z, event_weight );
	
	
	if( locRFBunch->dNumParticleVotes < 2 ) {
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	
	
	
	for( vector< const DBeamPhoton* >::const_iterator gam = beam_photons.begin();
		gam != beam_photons.end(); gam++ ) {
		
		int counter = (*gam)->dCounter;
		
		DetectorSystem_t sys = (*gam)->dSystem;
		if( sys==SYS_TAGH )      h_tagh_flux->Fill( counter );
		else if( sys==SYS_TAGM ) h_tagm_flux->Fill( counter );
	}
	
	
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
		
		double inner_layer_cut = 1.5*4.0157;
		
		int fid_cut = 0;
		if( (-1.*inner_layer_cut < face_x && face_x < inner_layer_cut) 
			&& (-1.*inner_layer_cut < face_y && face_y < inner_layer_cut) ) fid_cut = 1;
		
		if( (fabs(loc_t) < FCAL_RF_time_cut) && ((*show)->getEnergy() > FCAL_min_energy_cut) 
			&& !fid_cut && (loc_theta < 4.0) ) {
		  	n_fcal_showers++;
		  	fcal_candidates.push_back( (*show) );
		}
		
	}
	
	
	int n_ccal_showers = 0;	
	for( vector< const DCCALShower* >::const_iterator show = ccal_showers.begin();
		show != ccal_showers.end(); show++ ) {
		
		/* 
		In the simulation, the CCAL position was shifted to (x,y) = (m_ccalX_new,m_ccalY_new)
		because this was the position found using the offline alignment procedure.
		
		In the reconstruction software (DCCALShower_factory), however, the showers are applied the 
		(now incorrect) offset of x1 -> x1 + m_ccalX && y1 -> y1 + m_ccalY.
		
		Our goal with this bit of code is to apply a fiducial cut on the CCAL using the x,y position
		that the electrons/photons entered the CCAL.
		
		In my own build of halld_recon, I changed the DCCALShower_factory to give x,y,z inside the calorimeter
		instead of projecting it to the surface. 
		So what we need to do hear is apply the correct offset, then project to the face of the CCAL.
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
				
				loc_Cand.bunch_val    = bunch_val;
				
				loc_Cand.e1           = e1;
				loc_Cand.x1           = pos1.X();
				loc_Cand.y1           = pos1.Y();
				loc_Cand.z1           = pos1.Z();
				loc_Cand.e2           = e2;
				loc_Cand.x2           = pos2.X();
				loc_Cand.y2           = pos2.Y();
				loc_Cand.z2           = pos2.Z();
				
				loc_Cand.phi_slice    = loc_phi_slice;
				
				loc_Cand.deltaPhi     = deltaPhi;
				loc_Cand.deltaT       = deltaT;
				loc_Cand.deltaR       = deltaR;
				loc_Cand.deltaE       = deltaE;
				loc_Cand.deltaK       = deltaK;
				loc_Cand.deltaK2      = deltaK2;
				
				loc_Cand.vz           = vertex_z;
				loc_Cand.event_weight = event_weight;
				
				loc_Cand.eb           = eb;
				loc_Cand.tag_counter  = (*gam)->dCounter;
				
				DetectorSystem_t sys  = (*gam)->dSystem;
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
jerror_t JEventProcessor_compton_simulation_reweight::erun(void)
{
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_compton_simulation_reweight::fini(void)
{
	
	return NOERROR;
}


//--------------------------------------------
// Get CCAL layer from row and column
//--------------------------------------------
int JEventProcessor_compton_simulation_reweight::ccalLayer(int row, int col) {
	
	
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
int JEventProcessor_compton_simulation_reweight::fcalLayer(int row, int col) {
	
	
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
void JEventProcessor_compton_simulation_reweight::fill_histograms( vector< ComptonCandidate_t > Comp_Candidates ) {
	
	
	
	int n_candidates  =  static_cast<int>( Comp_Candidates.size() );
	
	for( int ic = 0; ic < n_candidates; ic++ ) {
				
		ComptonCandidate_t loc_Cand = Comp_Candidates[ic];
				
		//-------------------------------------------//
		
		int bunch_val        =  loc_Cand.bunch_val;
		
		double eb            =  loc_Cand.eb;
		int tag_sys          =  loc_Cand.tag_sys;
		int tag_counter      =  loc_Cand.tag_counter;
		
		double deltaPhi      =  loc_Cand.deltaPhi;
		double deltaE        =  loc_Cand.deltaE;
		double deltaK        =  loc_Cand.deltaK;
		
		double x1            =  loc_Cand.x1;
		double y1            =  loc_Cand.y1;
		double x2            =  loc_Cand.x2;
		double y2            =  loc_Cand.y2;
		
		double vertex_z      =  loc_Cand.vz;
		double event_weight  =  loc_Cand.event_weight;
		
		//--------------     Cuts      --------------//
		
		double   deltaE_mu,   deltaE_sig;
		double deltaPhi_mu, deltaPhi_sig;
		double   deltaK_mu,   deltaK_sig;
		
		deltaE_mu     = f_deltaE_mu->Eval(eb);
		deltaE_sig    = eb * f_deltaE_sig->Eval(eb);
		
		deltaPhi_mu   = f_deltaPhi_mu->Eval(eb);
		deltaPhi_sig  = f_deltaPhi_sig->Eval(eb);
		
		deltaK_mu     = f_deltaK_mu->Eval(eb);
		deltaK_sig    = f_deltaK_sig->Eval(eb);
		
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
		
		
		
		if( tag_sys==0 ) {
			
			h_deltaE_tagh->Fill( tag_counter, deltaE, fill_weight );
			h_deltaE_tagh_weighted->Fill( tag_counter, deltaE, event_weight*fill_weight );
			
			if( e5_cut ) {
				h_deltaPhi_tagh->Fill( tag_counter, deltaPhi, fill_weight );
				h_deltaPhi_tagh_weighted->Fill( tag_counter, deltaPhi, 
					event_weight*fill_weight );
				if( p5_cut ) {
					h_deltaK_tagh->Fill( tag_counter, deltaK, fill_weight );
					h_deltaK_tagh_weighted->Fill( tag_counter, deltaK, 
						fill_weight*event_weight );
					if( k5_cut ) {
						h_deltaK_tagh_cut->Fill( tag_counter, deltaK, 
							fill_weight );
						h_deltaK_tagh_cut_weighted->Fill( tag_counter, deltaK, 
							fill_weight*event_weight );
					}
				}
			}
			
		} else {
			
			h_deltaE_tagm->Fill( tag_counter, deltaE, fill_weight );
			h_deltaE_tagm_weighted->Fill( tag_counter, deltaE, event_weight*fill_weight );
			
			if( e5_cut ) {
				h_deltaPhi_tagm->Fill( tag_counter, deltaPhi, fill_weight );
				h_deltaPhi_tagm_weighted->Fill( tag_counter, deltaPhi, 
					event_weight*fill_weight );
				if( p5_cut ) {
					h_deltaK_tagm->Fill( tag_counter, deltaK, fill_weight );
					h_deltaK_tagm_weighted->Fill( tag_counter, deltaK, 
						fill_weight*event_weight );
					if( k5_cut ) {
						h_deltaK_tagm_cut->Fill( tag_counter, deltaK, 
							fill_weight );
						h_deltaK_tagm_cut_weighted->Fill( tag_counter, deltaK, 
							fill_weight*event_weight );
					}
				}
			}
			
		}
		
		
		if( e5_cut && p5_cut && k5_cut ) {
			
			h_fcal_xy->Fill( x1, y1 );
			h_fcal_xy_weighted->Fill( x1, y1, event_weight );
			
			h_ccal_xy->Fill( x2, y2 );
			h_ccal_xy_weighted->Fill( x2, y2, event_weight );
			
			h_vertex_cut->Fill( vertex_z );
			h_vertex_cut_weighted->Fill( vertex_z, event_weight );
			
		}
		
	}
	
	
	
	
	
	
	return;
}



void JEventProcessor_compton_simulation_reweight::set_cuts( int32_t runnumber )
{
	
	if( runnumber > 61320 && runnumber < 61355 ) {
		
		deltaE_mu_p0    =  1.14092e-01;
		deltaE_mu_p1    = -6.54721e-02;
		deltaE_mu_p2    =  7.10876e-03;
		deltaE_mu_p3    = -2.28944e-04;
		
		deltaE_sig_p0   =  4.95969e-03;
		deltaE_sig_p1   =  3.33127e-02;
		deltaE_sig_p2   =  3.06011e-13;
		
		
		deltaPhi_mu_p0  =  1.81299e+02;
		deltaPhi_mu_p1  = -5.17968e-01;
		deltaPhi_mu_p2  =  5.35586e-02;
		deltaPhi_mu_p3  = -1.94653e-03;
		
		deltaPhi_sig_p0 =  1.05813e+01;
		deltaPhi_sig_p1 = -1.01581e+00;
		deltaPhi_sig_p2 =  8.43296e-02;
		deltaPhi_sig_p3 = -2.24102e-03;
		
		
		deltaK_mu_p0    =  2.44163e-01;
		deltaK_mu_p1    = -7.22041e-02;
		deltaK_mu_p2    =  6.14512e-03;
		deltaK_mu_p3    = -2.26866e-04;
		
		deltaK_sig_p0   =  1.46801e-01;
		deltaK_sig_p1   = -1.40138e-02;
		deltaK_sig_p2   =  3.70946e-03;
		deltaK_sig_p3   = -1.21212e-04;
		
		
	} else {
		
		deltaE_mu_p0    = -5.78515e-04;
		deltaE_mu_p1    = -2.37929e-02;
		deltaE_mu_p2    =  2.81620e-03;
		deltaE_mu_p3    = -3.99247e-05;
		
		deltaE_sig_p0   =  4.16177e-03;
		deltaE_sig_p1   =  3.59818e-02;
		deltaE_sig_p2   =  1.58440e-04;
		
		
		deltaPhi_mu_p0  =  1.78868e+02;
		deltaPhi_mu_p1  =  3.22387e-01;
		deltaPhi_mu_p2  = -4.06228e-02;
		deltaPhi_mu_p3  =  1.52118e-03;
		
		deltaPhi_sig_p0 =  9.68719e+00;
		deltaPhi_sig_p1 = -8.32833e-01;
		deltaPhi_sig_p2 =  6.66080e-02;
		deltaPhi_sig_p3 = -1.60904e-03;
		
		
		deltaK_mu_p0    =  2.23855e-01;
		deltaK_mu_p1    = -6.59566e-02;
		deltaK_mu_p2    =  5.43829e-03;
		deltaK_mu_p3    = -1.92947e-04;
		
		deltaK_sig_p0   =  1.52679e-01;
		deltaK_sig_p1   = -1.57996e-02;
		deltaK_sig_p2   =  3.80448e-03;
		deltaK_sig_p3   = -1.27476e-04;
		
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

