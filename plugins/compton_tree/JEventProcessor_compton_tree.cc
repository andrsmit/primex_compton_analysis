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


//------------------
// Constructor
//------------------
JEventProcessor_compton_tree::JEventProcessor_compton_tree()
{
	// Set defaults:
	
    	RUN_GROUP = 0;
	gPARMS->SetDefaultParameter( "COMPTON_TREE:RUN_GROUP", RUN_GROUP );
	
	sprintf( cut_pathName, "/work/halld/home/andrsmit/primex_compton_analysis/cuts" );
	
	read_cuts();
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
	
	locTreeBranchRegister.Register_Single<Int_t>("nComp");
	
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("tag_counter","nComp");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("tag_sys",    "nComp");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("bunch_val",  "nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("eb","nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("tb","nComp");
	
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("e1","nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("x1","nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("y1","nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("z1","nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("t1","nComp");
	
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("e2","nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("x2","nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("y2","nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("z2","nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("t2","nComp");
	
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("DeltaE",  "nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("DeltaPhi","nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("DeltaK",  "nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("DeltaK2", "nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("DeltaT",  "nComp");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("DeltaR",  "nComp");
	
	dTreeInterface->Create_Branches(locTreeBranchRegister);
	

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_compton_tree::brun(JEventLoop *eventLoop, int32_t runnumber)
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
	
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_compton_tree::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
	
	//-----   Check Trigger   -----//
	
	const DL1Trigger *trig = NULL;
  	try {
      	  	eventLoop->GetSingle(trig);
  	} catch (...) {}
	if (trig == NULL) { return NOERROR; }
	
	
	uint32_t trigmask    = trig->trig_mask;	
	uint32_t fp_trigmask = trig->fp_trig_mask;
	
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
	
	
	
	
	
	
	int n_fcal_showers = 0;	
	for( vector< const DFCALShower* >::const_iterator show = fcal_showers.begin(); 
		show != fcal_showers.end(); show++ ) {
		
		DVector3 loc_pos  =  (*show)->getPosition_log()  -  vertex;
		double loc_theta  =  loc_pos.Theta() * (180. / TMath::Pi());
		
		double face_x     =  vertex.X()  +  (loc_pos.X() * (m_fcalZ - vertex.Z())/loc_pos.Z());
		double face_y     =  vertex.Y()  +  (loc_pos.Y() * (m_fcalZ - vertex.Z())/loc_pos.Z());
		
		int row   = fcalGeom->row(    static_cast<float>(face_y) );
		int col   = fcalGeom->column( static_cast<float>(face_x) );
		int layer = fcalLayer( row, col );
		
		if( ((*show)->getEnergy() > FCAL_min_energy_cut) && (layer > 1) && (loc_theta < 4.0) ) {
		  	n_fcal_showers++;
		  	fcal_candidates.push_back( (*show) );
		}
		
	}
	
	
	int n_ccal_showers = 0;	
	for( vector< const DCCALShower* >::const_iterator show = ccal_showers.begin();
		show != ccal_showers.end(); show++ ) {
		
		double loc_x = (*show)->x1  -  vertex.X() - (m_ccalX + m_ccalX_new);
		double loc_y = (*show)->y1  -  vertex.Y() - (m_ccalY + m_ccalY_new);
		double loc_z = (*show)->z   -  vertex.Z();
		
		double xface = vertex.X()  +  (loc_x * (m_ccalZ - vertex.Z())/loc_z);
		double yface = vertex.Y()  +  (loc_y * (m_ccalZ - vertex.Z())/loc_z);
		
		int row   = ccalGeom->row(    static_cast<float>(yface-m_ccalY_new) );
		int col   = ccalGeom->column( static_cast<float>(xface-m_ccalX_new) );
		int layer = ccalLayer( row, col );
		
		if( ((*show)->E > CCAL_min_energy_cut) && (layer > 1) && (layer < 5) && (col != 1) ) {
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
				
				double eb       = (*gam)->lorentzMomentum().E();
				double tb       = (*gam)->time();
				int tag_counter = (*gam)->dCounter;
				
				double brfdt = tb - rfTime;
				
				int bunch_val;
				
				if( fabs(brfdt) < 6.012 )
					bunch_val = 1;
				else if( (  -(6.012 + 3.*4.008) <= brfdt && brfdt <= -(6.012 + 2.*4.008) )
					|| ( (6.012 + 2.*4.008) <= brfdt && brfdt <=  (6.012 + 3.*4.008) ) )
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
				
				
				DetectorSystem_t sys = (*gam)->dSystem;
				
				double   deltaE_mu,   deltaE_sig;
				double deltaPhi_mu, deltaPhi_sig;
				double   deltaK_mu,   deltaK_sig;
				
				if( sys==SYS_TAGH ) {
					
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
				
				int e_cut = 0;
				if( fabs(deltaE  -  deltaE_mu) < 10.0*  deltaE_sig ) e_cut = 1;
				
				int p_cut = 0;
				if( fabs(deltaPhi-deltaPhi_mu) < 10.0*deltaPhi_sig ) p_cut = 1;
				
				int k_cut = 0;
				if( fabs(deltaK  -  deltaK_mu) < 10.0*  deltaK_sig ) k_cut = 1;
				
				if( e_cut && p_cut && k_cut ) {
					
					ComptonCandidate_t loc_Cand;
					
					loc_Cand.bunch_val = bunch_val;
					
					loc_Cand.e1        = e1;
					loc_Cand.x1        = pos1.X();
					loc_Cand.y1        = pos1.Y();
					loc_Cand.z1        = pos1.Z();
					loc_Cand.t1        = t1;
					loc_Cand.e2        = e2;
					loc_Cand.x2        = pos2.X();
					loc_Cand.y2        = pos2.Y();
					loc_Cand.z2        = pos2.Z();
					loc_Cand.t2        = t2;
					
					loc_Cand.deltaPhi  = deltaPhi;
					loc_Cand.deltaT    = deltaT;
					loc_Cand.deltaR    = deltaR;
					loc_Cand.deltaE    = deltaE;
					loc_Cand.deltaK    = deltaK;
					loc_Cand.deltaK2   = deltaK2;
					
					loc_Cand.eb          = eb;
					loc_Cand.tb          = tb;
					loc_Cand.tag_counter = tag_counter;
					
					if( sys==SYS_TAGH )      loc_Cand.tag_sys = 0;
					else if( sys==SYS_TAGM ) loc_Cand.tag_sys = 1;
					
					candidates.push_back( loc_Cand );
				
				}
				
			} // end DBeamPhoton loop
			
		} // end DCCALShower loop
		
	} // end DFCALShower loop
	
	
	int n_cands = (int)candidates.size();
	
	dTreeFillData.Fill_Single<Int_t>("eventNum", eventnumber);
	dTreeFillData.Fill_Single<Double_t>("rfTime", rfTime);
	
	size_t candIndex = 0;
	for( int ic = 0; ic < n_cands; ic++ ) {
		
		dTreeFillData.Fill_Array<Double_t>("eb", candidates[ic].eb, candIndex);
		dTreeFillData.Fill_Array<Double_t>("tb", candidates[ic].tb, candIndex);
		dTreeFillData.Fill_Array<Int_t>("tag_counter", candidates[ic].tag_counter, candIndex);
		dTreeFillData.Fill_Array<Int_t>("tag_sys",     candidates[ic].tag_sys,     candIndex);
		dTreeFillData.Fill_Array<Int_t>("bunch_val",   candidates[ic].bunch_val,   candIndex);
		
		
		dTreeFillData.Fill_Array<Double_t>("e1", candidates[ic].e1, candIndex);
		dTreeFillData.Fill_Array<Double_t>("x1", candidates[ic].x1, candIndex);
		dTreeFillData.Fill_Array<Double_t>("y1", candidates[ic].y1, candIndex);
		dTreeFillData.Fill_Array<Double_t>("z1", candidates[ic].z1, candIndex);
		dTreeFillData.Fill_Array<Double_t>("t1", candidates[ic].t1, candIndex);
		
		dTreeFillData.Fill_Array<Double_t>("e2", candidates[ic].e2, candIndex);
		dTreeFillData.Fill_Array<Double_t>("x2", candidates[ic].x2, candIndex);
		dTreeFillData.Fill_Array<Double_t>("y2", candidates[ic].y2, candIndex);
		dTreeFillData.Fill_Array<Double_t>("z2", candidates[ic].z2, candIndex);
		dTreeFillData.Fill_Array<Double_t>("t2", candidates[ic].t2, candIndex);
		
		dTreeFillData.Fill_Array<Double_t>("DeltaE",   candidates[ic].deltaE,   candIndex);
		dTreeFillData.Fill_Array<Double_t>("DeltaPhi", candidates[ic].deltaPhi, candIndex);
		dTreeFillData.Fill_Array<Double_t>("DeltaK",   candidates[ic].deltaK,   candIndex);
		dTreeFillData.Fill_Array<Double_t>("DeltaK2",  candidates[ic].deltaK2,  candIndex);
		dTreeFillData.Fill_Array<Double_t>("DeltaT",   candidates[ic].deltaT,   candIndex);
		dTreeFillData.Fill_Array<Double_t>("DeltaR",   candidates[ic].deltaR,   candIndex);
		
		candIndex++;
	}
	dTreeFillData.Fill_Single<Int_t>("nComp", candIndex);
	
	dTreeInterface->Fill(dTreeFillData);
	


	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_compton_tree::erun(void)
{
	
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


//------------------
// read in cut values
//------------------
void JEventProcessor_compton_tree::read_cuts()
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
int JEventProcessor_compton_tree::ccalLayer(int row, int col) {
	
	
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
int JEventProcessor_compton_tree::fcalLayer(int row, int col) {
	
	
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

