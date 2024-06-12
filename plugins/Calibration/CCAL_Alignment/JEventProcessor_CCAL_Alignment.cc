// $Id$
//
//    File: JEventProcessor_CCAL_Alignment.cc
// Created: Fri Jun  7 10:39:01 AM EDT 2024
// Creator: andrsmit (on Linux ifarm180302.jlab.org 5.14.0-362.24.2.el9_3.x86_64 x86_64)
//

#include "JEventProcessor_CCAL_Alignment.h"

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_CCAL_Alignment());
}
} // "C"

thread_local DTreeFillData JEventProcessor_CCAL_Alignment::dTreeFillData;

//------------------
// init
//------------------
jerror_t JEventProcessor_CCAL_Alignment::init(void)
{
	
	TDirectory *dir_CCAL_Alignment = new TDirectoryFile("CCAL_Alignment", 
		"CCAL_Alignment");
	dir_CCAL_Alignment->cd();
	
	for(int ihist=0; ihist<n_cuts; ihist++) {
		
		h_deltaE[ihist] = new TH1F("deltaE_cut_%d", 
			"#DeltaE; E_{1} + E_{2} - E_{#gamma} [GeV]", 2000, -8., 8.);
		
		h_deltaK[ihist] = new TH1F("deltaK_cut_%d", 
			"#DeltaK; E_{Compton} - E_{#gamma} [GeV]", 2000, -8., 8.);
		
		h_deltaPhi[ihist] = new TH1F("deltaPhi_cut_%d", 
			"#Delta#phi; |#phi_{1} - #phi_{2}| [#circ]", 3600, 0., 360.);
		
		h_e1[ihist] = new TH1F("e1_cut_%d",
			"Energy of Shower 1; E_{1} [GeV]", 6000, 0., 6.);
		h_e2[ihist] = new TH1F("e2_cut_%d",
			"Energy of Shower 2; E_{2} [GeV]", 6000, 0., 6.);
		h_eb[ihist] = new TH1F("eb_cut_%d",
			"Energy of Beam Photon; E_{#gamma} [GeV]", 1200, 0., 12.);
		
		h_xy1[ihist] = new TH2F("xy1_cut_%d",
			"Position of Shower 1; x_{1} [cm]; y_{1} [cm]", 500, -50., 50., 500, -50., 50.);
		h_xy2[ihist] = new TH2F("xy2_cut_%d",
			"Position of Shower 1; x_{1} [cm]; y_{1} [cm]", 500, -50., 50., 500, -50., 50.);
	}
	
	dir_CCAL_Alignment->cd("../");
	
	dTreeInterface = DTreeInterface::Create_DTreeInterface(
		"ccal_compton", "ccal_compton.root");
	DTreeBranchRegister locTreeBranchRegister;
	
	locTreeBranchRegister.Register_Single<Double_t>("eb");
	locTreeBranchRegister.Register_Single<Double_t>("tb");
	locTreeBranchRegister.Register_Single<Double_t>("rfTime");
	
	locTreeBranchRegister.Register_Single<Double_t>("e1");
	locTreeBranchRegister.Register_Single<Double_t>("x1");
	locTreeBranchRegister.Register_Single<Double_t>("y1");
	locTreeBranchRegister.Register_Single<Double_t>("z1");
	locTreeBranchRegister.Register_Single<Double_t>("t1");
	locTreeBranchRegister.Register_Single<Int_t>("tof_match");
	
	locTreeBranchRegister.Register_Single<Double_t>("e2");
	locTreeBranchRegister.Register_Single<Double_t>("x2");
	locTreeBranchRegister.Register_Single<Double_t>("y2");
	locTreeBranchRegister.Register_Single<Double_t>("z2");
	locTreeBranchRegister.Register_Single<Double_t>("t2");
	
	dTreeInterface->Create_Branches(locTreeBranchRegister);
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_CCAL_Alignment::brun(JEventLoop *eventLoop, int32_t runnumber)
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
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_CCAL_Alignment::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	//--------------------------------------------------------------------//
	// Check trigger:
	
	const DL1Trigger *locTrig = NULL;
	try {
		loop->GetSingle(locTrig);
	} catch (...) {}
	if(locTrig == NULL) return NOERROR;
	if(locTrig->fp_trig_mask) return NOERROR;
	
	//--------------------------------------------------------------------//
	// Check RF time:
	
	const DEventRFBunch *locRFBunch = NULL;
	try { 
		loop->GetSingle(locRFBunch, "CalorimeterOnly");
	} catch (...) { return NOERROR; }
	double locRFTime = locRFBunch->dTime;
	if(locRFBunch->dNumParticleVotes < 2) return NOERROR;
	
	//--------------------------------------------------------------------//
	// Get releveant data objects:
	
	DVector3 locVertex;
	locVertex.SetXYZ(m_beamX, m_beamY, m_beamZ);
	
	vector<const DBeamPhoton*> locBeamPhotons;
	loop->Get(locBeamPhotons);
	
	vector<const DFCALShower*> locFCALShowers;
	loop->Get(locFCALShowers);
	
	vector<const DCCALShower*> locCCALShowers;
	loop->Get(locCCALShowers);
	
	vector<const DBCALShower*> locBCALShowers;
	loop->Get(locBCALShowers);
	
	vector<const DTOFPoint*> locTOFPoints;
	loop->Get(locTOFPoints);
	
	//--------------------------------------------------------------------//
	// Fill lock for multi-threaded running:
	
	japp->RootFillLock(this);
	
	//--------------------------------------------------------------------//
	// Select events with correct multiplicities:
	
	int locNFCALShowers = 0, locNCCALShowers = 0, locNBCALShowers = 0;
	
	vector<const DFCALShower*> locGoodFCALShowers;
	locGoodFCALShowers.clear();
	vector<const DCCALShower*> locGoodCCALShowers;
	locGoodCCALShowers.clear();
	
	for(vector<const DFCALShower*>::const_iterator show = locFCALShowers.begin(); 
		show != locFCALShowers.end(); show++) {
		
		DVector3 pos = (*show)->getPosition_log() - locVertex;
		double rf_dt = (*show)->getTime() - (pos.Mag()/c) - locRFTime;
		
		if(fabs(rf_dt)<2.0) {
			locNFCALShowers++;
			locGoodFCALShowers.push_back((*show));
		}
	}
	
	// Loop over CCAL Showers:
	
	for(vector<const DCCALShower*>::const_iterator show = locCCALShowers.begin();
		show != locCCALShowers.end(); show++) {
		
		DVector3 pos((*show)->x1, (*show)->y1, (*show)->z);
		pos = pos - locVertex;
		double rf_dt = (*show)->time - (pos.Mag()/c) - locRFTime;
		
		if(fabs(rf_dt)<2.0) { 
			locNCCALShowers++;
			locGoodCCALShowers.push_back((*show));
		}
	}
	
	// Loop over BCAL Shower:
	
	for(vector<const DBCALShower*>::const_iterator show = locBCALShowers.begin();
		show != locBCALShowers.end(); show++) {
		
		DVector3 pos((*show)->x, (*show)->y, (*show)->z);
		pos = pos - locVertex;
		
		double rf_dt = (*show)->t - (pos.Mag()/c) - locRFTime;
		
		if(fabs(rf_dt)<6.0) { locNBCALShowers++; }
	}
	
	if(locNBCALShowers>0 || locNCCALShowers!=1 || locNFCALShowers!=1) {
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	
	//--------------------------------------------------------------------//
	// Evaluate 2-cluster pair in the FCAL:
	
	const DFCALShower *show1 = locGoodFCALShowers[0];
	const DCCALShower *show2 = locGoodCCALShowers[0];
	
	DVector3 pos1 = show1->getPosition_log() - locVertex;
	double e1     = show1->getEnergy();
	double phi1   = pos1.Phi() * (180. / TMath::Pi());
	double theta1 = pos1.Theta();
	
	DVector3 pos2(show2->x1, show2->y1, show2->z);
	pos2 -= locVertex;
	double e2     = show2->E;
	double phi2   = pos2.Phi() * (180. / TMath::Pi());
	double theta2 = pos2.Theta();
	
	double deltaPhi = fabs(phi2 - phi1);
	
	int   locNGoodCands = 0;
	double      eb_best = 0.;
	double      tb_best = 0.;
	
	for(vector< const DBeamPhoton* >::const_iterator gam = locBeamPhotons.begin();
		gam != locBeamPhotons.end(); gam++) {
		
		double eb    = (*gam)->lorentzMomentum().E();
		double tb    = (*gam)->time();
		double brfdt = tb - locRFTime;
		
		int bunch_val = -1;
		
		if(fabs(brfdt) < 2.004)
			bunch_val = 1;
		else if(( -(6.012 + 5.*4.008) <= brfdt && brfdt <= -(6.012 + 3.*4.008))
			|| ((6.012 + 3.*4.008) <= brfdt && brfdt <=  (6.012 + 5.*4.008)))
			bunch_val = 0;
		else 
			continue;
		
		double ecomp1 =  1. / ((1./eb) + (1./m_e)*(1.-cos(theta1)));
		double ecomp2 =  1. / ((1./eb) + (1./m_e)*(1.-cos(theta2)));
		double deltaE = (e1     + e2    ) - (eb + m_e);
		double deltaK = (ecomp1 + ecomp2) - (eb + m_e);
		
		double fill_weight = 0.;
		if(bunch_val==1)      fill_weight =  1.0;
		else if(bunch_val==0) fill_weight = -0.25;
		
		int cut_vals[n_cuts];
		for(int icut = 0; icut < n_cuts; icut++) cut_vals[icut] = 0;
		
		// first cut is trivial:
		cut_vals[0] = 1;
		
		// second cut is a fiducial cut on FCAL:
		if(!fcal_fiducial_cut(pos1,locVertex,1.0)) {
			
			cut_vals[1] = 1;
			
			// next cut is coplanarity:
			if(fabs(deltaPhi-180) < 50.) {
				
				cut_vals[2] = 1;
				
				// next cut is shower energy and elasticity:
				if(e1 > 0.5 && e2 > 0.5 && fabs(deltaE) < 1.0) {
					
					cut_vals[3] = 1;
				}
			}
		}
		
		// the final cut is everything above, but with a stricter fiducial cut:
		if(!fcal_fiducial_cut(pos1,locVertex,2.0)) {
			if(fabs(deltaPhi-180) < 50.) {
				if(e1 > 0.5 && e2 > 0.5 && fabs(deltaE) < 1.0) {
					cut_vals[4] = 1;
				}
			}
		}
		
		//-----------------------------------------------------------------//
		// Fill Histograms:
		
		for(int icut=0; icut<n_cuts; icut++) {
			if(cut_vals[icut]) {
				h_deltaE[icut]->Fill(deltaE, fill_weight);
				h_deltaK[icut]->Fill(deltaK, fill_weight);
				h_deltaPhi[icut]->Fill(deltaPhi, fill_weight);
				h_e1[icut]->Fill(e1, fill_weight);
				h_e2[icut]->Fill(e2, fill_weight);
				h_eb[icut]->Fill(eb, fill_weight);
				h_xy1[icut]->Fill(pos1.X(), pos1.Y());
				h_xy2[icut]->Fill(pos2.X(), pos2.Y());
			}
		}
		
		if(bunch_val==1) {
			if(fabs(deltaPhi-180)<50. && fabs(deltaE)<1. && e1>0.5 && e2>0.5) {
				locNGoodCands++;
				eb_best = eb;
				tb_best = tb;
			}
		}
		
	} // end DBeamPhoton loop
	
	//-----------------------------------------------------------------//
	// Fill TTree Data for the event:
	
	if(locNGoodCands==1) {
		dTreeFillData.Fill_Single<Double_t>("eb", eb_best);
		dTreeFillData.Fill_Single<Double_t>("tb", tb_best);
		dTreeFillData.Fill_Single<Double_t>("rfTime",locRFTime);
		dTreeFillData.Fill_Single<Double_t>("e1", e1);
		dTreeFillData.Fill_Single<Double_t>("x1", show1->getPosition_log().X());
		dTreeFillData.Fill_Single<Double_t>("y1", show1->getPosition_log().Y());
		dTreeFillData.Fill_Single<Double_t>("z1", show1->getPosition_log().Z());
		dTreeFillData.Fill_Single<Double_t>("t1", show1->getTime());
		dTreeFillData.Fill_Single<Double_t>("e2", e2);
		dTreeFillData.Fill_Single<Double_t>("x2", show2->x1);
		dTreeFillData.Fill_Single<Double_t>("y2", show2->y1);
		dTreeFillData.Fill_Single<Double_t>("z2", show2->z);
		dTreeFillData.Fill_Single<Double_t>("t2", show2->time);
		
		//-----------------------------------------------------------------//
		// Check for matches between FCAL shower and TOF:
		
		int tof_match = 0;
		double tof_dx, tof_dy, tof_dt;
		check_TOF_match(pos1, locRFTime, locVertex, locTOFPoints, tof_dx, tof_dy, tof_dt, 1.0);
		if(sqrt(pow(tof_dx,2.0) + pow(tof_dy,2.0)) < 8.0) {
			tof_match = 1;
		}
		dTreeFillData.Fill_Single<Int_t>("tof_match", tof_match);
		
		dTreeInterface->Fill(dTreeFillData);
	}
	
	japp->RootFillUnLock(this);
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_CCAL_Alignment::fini(void)
{
	delete dTreeInterface;
	return NOERROR;
}

//------------------
// fcal_fiducial_cut
//------------------

int JEventProcessor_CCAL_Alignment::fcal_fiducial_cut(DVector3 pos, DVector3 vertex, 
	double cut_layers) {
	
  int fid_cut = 0;
  
  double fcal_inner_layer_cut = (1.5 + cut_layers) * 4.0157;
  
  double fcal_face_x = vertex.X() + (pos.X() * (m_fcalZ - vertex.Z())/pos.Z());
  double fcal_face_y = vertex.Y() + (pos.Y() * (m_fcalZ - vertex.Z())/pos.Z());
  
  fcal_face_x -= m_fcalX;
  fcal_face_y -= m_fcalY;
  
  if((-1.*fcal_inner_layer_cut < fcal_face_x && fcal_face_x < fcal_inner_layer_cut)
     && (-1.*fcal_inner_layer_cut < fcal_face_y 
	 && fcal_face_y < fcal_inner_layer_cut)) fid_cut = 1;
  
  return fid_cut;
}

//------------------
// ccal_fiducial_cut
//------------------

int JEventProcessor_CCAL_Alignment::ccal_fiducial_cut(DVector3 pos, DVector3 vertex)
{
	int fid_cut = 0;
	
	double ccal_inner_layer_cut = 2.0 * 2.09;
	
	double ccal_face_x = vertex.X() + (pos.X() * (m_ccalZ - vertex.Z())/pos.Z());
	double ccal_face_y = vertex.Y() + (pos.Y() * (m_ccalZ - vertex.Z())/pos.Z());
	
	ccal_face_x -= m_ccalX;
	ccal_face_y -= m_ccalY;
	
	if((-1.*ccal_inner_layer_cut < ccal_face_x && ccal_face_x < ccal_inner_layer_cut)
		&& (-1.*ccal_inner_layer_cut < ccal_face_y 
		&& ccal_face_y < ccal_inner_layer_cut)) fid_cut = 1;
	
	return fid_cut;
}

//------------------
// check_TOF_match
//------------------

void JEventProcessor_CCAL_Alignment::check_TOF_match(DVector3 pos1, double rfTime, 
	DVector3 vertex, vector<const DTOFPoint*> tof_points, double &dx_min, double &dy_min, 
	double &dt_min, double rf_time_cut) {
  
  dx_min = 1000.;
  dy_min = 1000.;
  dt_min = 1000.;
  
  for(vector<const DTOFPoint*>::const_iterator tof = tof_points.begin(); 
      tof != tof_points.end(); tof++) {
    
    double xt = (*tof)->pos.X() - vertex.X();
    double yt = (*tof)->pos.Y() - vertex.Y();
    double zt = (*tof)->pos.Z() - vertex.Z();
    double rt = sqrt(xt*xt + yt*yt + zt*zt);
    double tt = (*tof)->t - (rt/c);
    double dt = tt - rfTime;
    xt *= pos1.Z() / zt;
    yt *= pos1.Z() / zt;
    double dx = pos1.X() - xt;
    double dy = pos1.Y() - yt;
    
    if(fabs(dt) < rf_time_cut) {
      if((dx*dx + dy*dy) < (dx_min*dx_min + dy_min*dy_min)) {
				dx_min = dx;
				dy_min = dy;
				dt_min = dt;
      }
    }
  }
  
  return;
}
