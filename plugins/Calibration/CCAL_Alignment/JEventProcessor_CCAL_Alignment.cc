// $Id$
//
//    File: JEventProcessor_CCAL_Alignment.cc
// Created: Sat Jan 29 00:56:12 EST 2022
// Creator: andrsmit (on Linux ifarm1801.jlab.org 3.10.0-1160.11.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_CCAL_Alignment.h"
using namespace jana;

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
	
	h_beam_rf_dt = new TH1F("beam_rf_dt", 
		"Beam Photon Time (propogated to target) - Event RF Time; [ns]", 
		2000, -100., 100.);
	h_fcal_rf_dt = new TH1F("fcal_rf_dt", 
		"FCAL Shower Time (propogated to target) - Event RF Time; [ns]", 
		2000, -100., 100.);
	h_ccal_rf_dt = new TH1F("ccal_rf_dt", 
		"CCAL Shower Time (propogated to target) - Event RF Time; [ns]", 
		2000, -100., 100.);
	
	h_deltaE      = new TH1F("deltaE", 
		"E_{#gamma} - (E_{1}+E_{2}); #DeltaE [GeV]", 2000, -8., 8.);
	h_deltaE_cut  = new TH1F("deltaE_cut", 
		"E_{#gamma} - (E_{1}+E_{2}); #DeltaE [GeV]", 2000, -8., 8.);
	h_deltaE_cut2 = new TH1F("deltaE_cut2", 
		"E_{#gamma} - (E_{1}+E_{2}); #DeltaE [GeV]", 2000, -8., 8.);
	h_deltaE_cut3 = new TH1F("deltaE_cut3", 
		"E_{#gamma} - (E_{1}+E_{2}); #DeltaE [GeV]", 2000, -8., 8.);
	h_deltaE_cut4 = new TH1F("deltaE_cut4", 
		"E_{#gamma} - (E_{1}+E_{2}); #DeltaE [GeV]", 2000, -8., 8.);
	
	h_deltaK      = new TH1F("deltaK", 
		"E_{#gamma} - (E_{Compton}); #DeltaK [GeV]", 2000, -8., 8.);
	h_deltaK_cut  = new TH1F("deltaK_cut", 
		"E_{#gamma} - (E_{Compton}); #DeltaK [GeV]", 2000, -8., 8.);
	h_deltaK_cut2 = new TH1F("deltaK_cut2", 
		"E_{#gamma} - (E_{Compton}); #DeltaK [GeV]", 2000, -8., 8.);
	h_deltaK_cut3 = new TH1F("deltaK_cut3", 
		"E_{#gamma} - (E_{Compton}); #DeltaK [GeV]", 2000, -8., 8.);
	h_deltaK_cut4 = new TH1F("deltaK_cut4", 
		"E_{#gamma} - (E_{Compton}); #DeltaK [GeV]", 2000, -8., 8.);
	
	h_deltaPhi      = new TH1F("deltaPhi", 
		"|#phi_{1} - #phi_{2}|; #Delta#phi [#circ]", 3600, 0., 360.);
	h_deltaPhi_cut  = new TH1F("deltaPhi_cut", 
		"|#phi_{1} - #phi_{2}|; #Delta#phi [#circ]", 3600, 0., 360.);
	h_deltaPhi_cut2 = new TH1F("deltaPhi_cut2", 
		"|#phi_{1} - #phi_{2}|; #Delta#phi [#circ]", 3600, 0., 360.);
	h_deltaPhi_cut3 = new TH1F("deltaPhi_cut3", 
		"|#phi_{1} - #phi_{2}|; #Delta#phi [#circ]", 3600, 0., 360.);
	h_deltaPhi_cut4 = new TH1F("deltaPhi_cut4", 
		"|#phi_{1} - #phi_{2}|; #Delta#phi [#circ]", 3600, 0., 360.);
	
	h_e1      = new TH1F("e1",      
		"Energy of Shower 1; E_{1} [GeV]", 12000, 0., 12.);
	h_e1_cut  = new TH1F("e1_cut",  
		"Energy of Shower 1; E_{1} [GeV]", 12000, 0., 12.);
	h_e1_cut2 = new TH1F("e1_cut2", 
		"Energy of Shower 1; E_{1} [GeV]", 12000, 0., 12.);
	h_e1_cut3 = new TH1F("e1_cut3", 
		"Energy of Shower 1; E_{1} [GeV]", 12000, 0., 12.);
	h_e1_cut4 = new TH1F("e1_cut4", 
		"Energy of Shower 1; E_{1} [GeV]", 12000, 0., 12.);
	
	h_e2      = new TH1F("e2",      
		"Energy of Shower 2; E_{2} [GeV]", 12000, 0., 12.);
	h_e2_cut  = new TH1F("e2_cut",  
		"Energy of Shower 2; E_{2} [GeV]", 12000, 0., 12.);
	h_e2_cut2 = new TH1F("e2_cut2", 
		"Energy of Shower 2; E_{2} [GeV]", 12000, 0., 12.);
	h_e2_cut3 = new TH1F("e2_cut3", 
		"Energy of Shower 2; E_{2} [GeV]", 12000, 0., 12.);
	h_e2_cut4 = new TH1F("e2_cut4", 
		"Energy of Shower 1; E_{1} [GeV]", 12000, 0., 12.);
	
	h_eb      = new TH1F("eb",      
		"Energy of Beam Photon; E_{#gamma} [GeV]", 12000, 0., 12.);
	h_eb_cut  = new TH1F("eb_cut",  
		"Energy of Beam Photon; E_{#gamma} [GeV]", 12000, 0., 12.);
	h_eb_cut2 = new TH1F("eb_cut2", 
		"Energy of Beam Photon; E_{#gamma} [GeV]", 12000, 0., 12.);
	h_eb_cut3 = new TH1F("eb_cut3", 
		"Energy of Beam Photon; E_{#gamma} [GeV]", 12000, 0., 12.);
	h_eb_cut4 = new TH1F("eb_cut4", 
		"Energy of Beam Photon; E_{#gamma} [GeV]", 12000, 0., 12.);
	
	h_xy1      = new TH2F("xy1",      "Position of Shower 1; x_{1} [cm]; y_{1} [cm]", 
		500, -100., 100., 500, -100., 100.);
	h_xy1_cut  = new TH2F("xy1_cut",  "Position of Shower 1; x_{1} [cm]; y_{1} [cm]", 
		500, -100., 100., 500, -100., 100.);
	h_xy1_cut2 = new TH2F("xy1_cut2", "Position of Shower 1; x_{1} [cm]; y_{1} [cm]", 
		500, -100., 100., 500, -100., 100.);
	h_xy1_cut3 = new TH2F("xy1_cut3", "Position of Shower 1; x_{1} [cm]; y_{1} [cm]", 
		500, -100., 100., 500, -100., 100.);
	h_xy1_cut4 = new TH2F("xy1_cut4", "Position of Shower 1; x_{1} [cm]; y_{1} [cm]", 
		500, -100., 100., 500, -100., 100.);
	
	h_xy2      = new TH2F("xy2",      "Position of Shower 2; x_{2} [cm]; y_{2} [cm]", 
		500, -100., 100., 500, -100., 100.);
	h_xy2_cut  = new TH2F("xy2_cut",  "Position of Shower 2; x_{2} [cm]; y_{2} [cm]", 
		500, -100., 100., 500, -100., 100.);
	h_xy2_cut2 = new TH2F("xy2_cut2", "Position of Shower 2; x_{2} [cm]; y_{2} [cm]", 
		500, -100., 100., 500, -100., 100.);
	h_xy2_cut3 = new TH2F("xy2_cut3", "Position of Shower 2; x_{2} [cm]; y_{2} [cm]", 
		500, -100., 100., 500, -100., 100.);
	h_xy2_cut4 = new TH2F("xy2_cut4", "Position of Shower 2; x_{2} [cm]; y_{2} [cm]", 
		500, -100., 100., 500, -100., 100.);
	
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
	//-----   Check Trigger   -----//
	
	const DL1Trigger *trig = NULL;
  	try {
      	  	loop->GetSingle(trig);
  	} catch (...) {}
	if (trig == NULL) { return NOERROR; }
	if(trig->fp_trig_mask) return NOERROR;
	
	//-----   RF Bunch   -----//
	
	const DEventRFBunch *locRFBunch = NULL;
	try { 
	  	loop->GetSingle(locRFBunch, "CalorimeterOnly");
	} catch (...) { return NOERROR; }
	double rfTime = locRFBunch->dTime;
	if(locRFBunch->dNumParticleVotes < 2) return NOERROR;
	
	//-----   Data Objects   -----//
	
	DVector3 vertex;
	vertex.SetXYZ(m_beamX, m_beamY, m_beamZ);
	
	vector<const DBeamPhoton*> beam_photons;
	loop->Get(beam_photons);
	
	vector<const DFCALShower*> fcal_showers;
	loop->Get(fcal_showers);
	
	vector<const DCCALShower*> ccal_showers;
	loop->Get(ccal_showers);
	
	vector<const DTOFPoint*> tof_points;
	loop->Get(tof_points);
	
	japp->RootFillLock(this);
	
	int n_fcal_showers = 0;	
	vector<const DFCALShower*> good_fcal_showers;
	good_fcal_showers.clear();
	
	for(vector<const DFCALShower*>::const_iterator show = fcal_showers.begin(); 
		show != fcal_showers.end(); show++) {
		
		DVector3 pos = (*show)->getPosition_log() - vertex;
		
		double loc_rf_dt = (*show)->getTime() - (pos.Mag()/c) - rfTime;
		h_fcal_rf_dt->Fill(loc_rf_dt);
		
		if(fabs(loc_rf_dt)<6.0) {
		  	n_fcal_showers++;
		  	good_fcal_showers.push_back((*show));
		}
	}
	
	int n_ccal_showers = 0;	
	vector<const DCCALShower*> good_ccal_showers;
	good_ccal_showers.clear();
	
	for(vector<const DCCALShower*>::const_iterator show = ccal_showers.begin();
		show != ccal_showers.end(); show++) {
		
		DVector3 pos((*show)->x1, (*show)->y1, (*show)->z);
		pos = pos - vertex;
		
		double loc_rf_dt = (*show)->time - (pos.Mag()/c) - rfTime;
		h_ccal_rf_dt->Fill(loc_rf_dt);
		
		if(fabs(loc_rf_dt)<6.0) { 
			n_ccal_showers++;
			good_ccal_showers.push_back((*show));
		}
	}
	if(n_ccal_showers!=1 || n_fcal_showers!=1) {
		japp->RootFillUnLock(this);
		return NOERROR;
	}
	
	const DFCALShower *show1 = good_fcal_showers[0];
	const DCCALShower *show2 = good_ccal_showers[0];
	
	double e1     = show1->getEnergy();
	DVector3 pos1 = show1->getPosition_log() - vertex;
	
	double phi1   = pos1.Phi() * (180. / TMath::Pi());
	double theta1 = pos1.Theta();
	
	double e2     = show2->E;
	DVector3 pos2(show2->x1, show2->y1, show2->z);
	pos2 = pos2 - vertex;
	
	double phi2   = pos2.Phi() * (180. / TMath::Pi());
	double theta2 = pos2.Theta();
	
	double deltaPhi = fabs(phi2 - phi1);
	
	int   n_good_cands = 0;
	double     eb_best = 0.;
	double     tb_best = 0.;
	
	for(vector< const DBeamPhoton* >::const_iterator gam = beam_photons.begin();
		gam != beam_photons.end(); gam++) {
		
		double eb = (*gam)->lorentzMomentum().E();
		double tb = (*gam)->time();
		
		double brfdt = tb - rfTime;
		h_beam_rf_dt->Fill(brfdt);
		
		int bunch_val = -1;
		
		if(fabs(brfdt) < 2.004)
			bunch_val = 1;
		else if(( -(6.012 + 5.*4.008) <= brfdt && brfdt <= -(6.012 + 3.*4.008))
			|| ((6.012 + 3.*4.008) <= brfdt && brfdt <=  (6.012 + 5.*4.008)))
			bunch_val = 0;
		else 
			continue;
		
		//if(eb > 5.0) continue;
		
		double ecomp1 =  1. / ((1./eb) + (1./m_e)*(1.-cos(theta1)));
		double ecomp2 =  1. / ((1./eb) + (1./m_e)*(1.-cos(theta2)));
		double deltaE = (e1     + e2    ) - (eb + m_e);
		double deltaK = (ecomp1 + ecomp2) - (eb + m_e);
		
		double fill_weight = 0.;
		if(bunch_val==1)      fill_weight =  1.0;
		else if(bunch_val==0) fill_weight = -0.25;
		else                  fill_weight =  0.0;
		
		h_deltaK->Fill(deltaK, fill_weight);
		h_deltaPhi->Fill(deltaPhi, fill_weight);
		h_deltaE->Fill(deltaE, fill_weight);
		h_e1->Fill(e1, fill_weight);
		h_e2->Fill(e2, fill_weight);
		h_eb->Fill(eb, fill_weight);
		h_xy1->Fill(pos1.X(), pos1.Y());
		h_xy2->Fill(pos2.X(), pos2.Y());
		
		if(!fcal_fiducial_cut(pos1,vertex)) {
			
			h_deltaK_cut->Fill(deltaK, fill_weight);
			h_deltaPhi_cut->Fill(deltaPhi, fill_weight);
			h_deltaE_cut->Fill(deltaE, fill_weight);
			h_e1_cut->Fill(e1, fill_weight);
			h_e2_cut->Fill(e2, fill_weight);
			h_eb_cut->Fill(eb, fill_weight);
			h_xy1_cut->Fill(pos1.X(), pos1.Y());
			h_xy2_cut->Fill(pos2.X(), pos2.Y());
			
			if(fabs(deltaPhi-180)<50.) {
				
				h_deltaK_cut3->Fill(deltaK, fill_weight);
				h_deltaPhi_cut3->Fill(deltaPhi, fill_weight);
				h_deltaE_cut3->Fill(deltaE, fill_weight);
				h_e1_cut3->Fill(e1, fill_weight);
				h_e2_cut3->Fill(e2, fill_weight);
				h_eb_cut3->Fill(eb, fill_weight);
				h_xy1_cut3->Fill(pos1.X(), pos1.Y());
				h_xy2_cut3->Fill(pos2.X(), pos2.Y());
				
				if(e1>0.5 && e2>0.5 && fabs(deltaE)<1.0) {
					
					h_deltaK_cut4->Fill(deltaK, fill_weight);
					h_deltaPhi_cut4->Fill(deltaPhi, fill_weight);
					h_deltaE_cut4->Fill(deltaE, fill_weight);
					h_e1_cut4->Fill(e1, fill_weight);
					h_e2_cut4->Fill(e2, fill_weight);
					h_eb_cut4->Fill(eb, fill_weight);
					h_xy1_cut4->Fill(pos1.X(), pos1.Y());
					h_xy2_cut4->Fill(pos2.X(), pos2.Y());
				}
			}
		}
		
		if(!fcal_fiducial_cut2(pos1,vertex)) {
			h_deltaK_cut2->Fill(deltaK, fill_weight);
			h_deltaPhi_cut2->Fill(deltaPhi, fill_weight);
			h_deltaE_cut2->Fill(deltaE, fill_weight);
			h_e1_cut2->Fill(e1, fill_weight);
			h_e2_cut2->Fill(e2, fill_weight);
			h_eb_cut2->Fill(eb, fill_weight);
			h_xy1_cut2->Fill(pos1.X(), pos1.Y());
			h_xy2_cut2->Fill(pos2.X(), pos2.Y());
		}
		if(bunch_val==1) {
			if(fabs(deltaPhi-180)<50. && fabs(deltaE)<1. && e1>0.5 && e2>0.5) {
				n_good_cands++;
				eb_best = eb;
				tb_best = tb;
			}
		}
		
	} // end DBeamPhoton loop
	
	if(n_good_cands==1) {
		dTreeFillData.Fill_Single<Double_t>("eb", eb_best);
		dTreeFillData.Fill_Single<Double_t>("tb", tb_best);
		dTreeFillData.Fill_Single<Double_t>("rfTime",rfTime);
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
		
		int tof_match = check_TOF_match(pos1, rfTime, vertex, tof_points);
		dTreeFillData.Fill_Single<Int_t>("tof_match", tof_match);
		
		dTreeInterface->Fill(dTreeFillData);
	}
	
	japp->RootFillUnLock(this);
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_CCAL_Alignment::erun(void)
{
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


int JEventProcessor_CCAL_Alignment::fcal_fiducial_cut(DVector3 pos, DVector3 vertex)
{
	int fid_cut = 0;
	
	double fcal_inner_layer_cut = 2.5 * 4.0157;
	
	double fcal_face_x = vertex.X() + (pos.X() * (m_fcalZ - vertex.Z())/pos.Z());
	double fcal_face_y = vertex.Y() + (pos.Y() * (m_fcalZ - vertex.Z())/pos.Z());
	
	fcal_face_x       -= m_fcalX;
	fcal_face_y       -= m_fcalY;
	
	if((-1.*fcal_inner_layer_cut < fcal_face_x && fcal_face_x < fcal_inner_layer_cut)
		&& (-1.*fcal_inner_layer_cut < fcal_face_y 
		&& fcal_face_y < fcal_inner_layer_cut)) fid_cut = 1;
	
	return fid_cut;
}


int JEventProcessor_CCAL_Alignment::fcal_fiducial_cut2(DVector3 pos, DVector3 vertex)
{
	int fid_cut = 0;
	
	double fcal_inner_layer_cut = 3.5 * 4.0157;
	
	double fcal_face_x = vertex.X() + (pos.X() * (m_fcalZ - vertex.Z())/pos.Z());
	double fcal_face_y = vertex.Y() + (pos.Y() * (m_fcalZ - vertex.Z())/pos.Z());
	
	fcal_face_x       -= m_fcalX;
	fcal_face_y       -= m_fcalY;
	
	if((-1.*fcal_inner_layer_cut < fcal_face_x && fcal_face_x < fcal_inner_layer_cut)
		&& (-1.*fcal_inner_layer_cut < fcal_face_y 
		&& fcal_face_y < fcal_inner_layer_cut)) fid_cut = 1;
	
	return fid_cut;
}


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


int JEventProcessor_CCAL_Alignment::check_TOF_match(DVector3 pos1, double rfTime, 
	DVector3 vertex, vector<const DTOFPoint*> tof_points) {
	
	int global_tof_match = 0;
	
	for(vector<const DTOFPoint*>::const_iterator tof = tof_points.begin(); 
		tof != tof_points.end(); tof++) {
		
		double xt = (*tof)->pos.X() - vertex.X();
		double yt = (*tof)->pos.Y() - vertex.Y();
		double zt = (*tof)->pos.Z() - vertex.Z();
		double rt = sqrt(xt*xt + yt*yt + zt*zt);
		
		double tt = (*tof)->t;
		double dt = tt - (rt/c) - rfTime;
		
		//---------------------------------------------------------------------//
		
		xt *= pos1.Z() / zt;
		yt *= pos1.Z() / zt;
		
		double dx = pos1.X() - xt;
		double dy = pos1.Y() - yt;
		
		if(fabs(dx)<10.0 && fabs(dy)<10.0 && fabs(dt)<6.012) 
			if(fabs(dt)<6.012) global_tof_match++;
		
	}
	
	return global_tof_match;
}
