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
	DApplication* dapp = dynamic_cast< DApplication* >(eventLoop->GetJApplication());
	if(dapp)     dgeom = dapp->GetDGeometry(runnumber);
   	
	if(dgeom){
		dgeom->GetTargetZ(m_beamZ);
		dgeom->GetFCALPosition(m_fcalX, m_fcalY, m_fcalZ);
		dgeom->GetCCALPosition(m_ccalX, m_ccalY, m_ccalZ);
	} else{
		cerr << "No geometry accessbile to compton_analysis plugin." << endl;
		return RESOURCE_UNAVAILABLE;
	}
	
	jana::JCalibration *jcalib = japp->GetJCalibration(runnumber);
	std::map<string, float> beam_spot;
	jcalib->Get("PHOTON_BEAM/beam_spot", beam_spot);
	m_beamX  =  beam_spot.at("x");
	m_beamY  =  beam_spot.at("y");
	
	
	// apply correction to FCAL and CCAL positions from Compton-scattering alignment procedure:
	
	if(runnumber>60000 && runnumber<69999) {
		
		// Phase I:
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
		
		// PhaseII:
		m_fcalX_new =  0.408;
		m_fcalY_new =  0.027;
		m_ccalX_new =  0.108;
		m_ccalY_new =  0.130;
		
	} else if(runnumber>110000) {
		
		// Phase III:
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
	}
	
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
	
	//if(!trig->Get_IsPhysicsEvent()) return NOERROR;
	if(trigmask==8) return NOERROR;
	if(fp_trigmask) return NOERROR;
	
	
	//-----   RF Bunch   -----//
	
	const DEventRFBunch *locRFBunch = NULL;
	try { 
	  	eventLoop->GetSingle(locRFBunch, "CalorimeterOnly");
	} catch (...) { return NOERROR; }
	double rfTime = locRFBunch->dTime;
	if(locRFBunch->dNumParticleVotes < 2) return NOERROR;
	
	
	//-----   Geometry Objects   -----//
	
	DVector3 vertex;
	vertex.SetXYZ(m_beamX, m_beamY, m_beamZ);
	
	DVector3 fcal_correction(m_fcalX_new-m_fcalX, m_fcalY_new-m_fcalY, 0.);
	DVector3 ccal_correction(m_ccalX_new-m_ccalX, m_ccalY_new-m_ccalY, 0.);
	
	//-----   Data Objects   -----//
	
	vector<const DBeamPhoton*> beam_photons;
	vector<const DCCALShower*> ccal_showers;
	vector<const DFCALShower*> fcal_showers;
	
	eventLoop->Get(beam_photons);
	eventLoop->Get(ccal_showers);
	eventLoop->Get(fcal_showers);
	
	vector<const DFCALShower*> fcal_candidates;
	vector<const DCCALShower*> ccal_candidates;
	
	int n_fcal_showers = 0;	
	for(vector<const DFCALShower*>::const_iterator show = fcal_showers.begin(); 
		show != fcal_showers.end(); show++) {
		
		if(((*show)->getEnergy() > FCAL_min_energy_cut)) {
		  	n_fcal_showers++;
		  	fcal_candidates.push_back((*show));
		}
		
	}
	
	int n_ccal_showers = 0;	
	for(vector<const DCCALShower*>::const_iterator show = ccal_showers.begin();
		show != ccal_showers.end(); show++) {
		
		if(((*show)->E > CCAL_min_energy_cut)) {
			n_ccal_showers++;
			ccal_candidates.push_back((*show));
		}
		
	}
	
	//----------     Check FCAL-CCAL Pairs     ----------//
	
	vector<ComptonCandidate_t> candidates;
	
	for(vector<const DFCALShower*>::const_iterator show1 = fcal_candidates.begin(); 
		show1 != fcal_candidates.end(); show1++) {
		
		double e1     = (*show1)->getEnergy();
		DVector3 pos1 = (*show1)->getPosition_log() - vertex;
		pos1          += fcal_correction;
		
		double t1     =  (*show1)->getTime()  -  (pos1.Mag()/c);
		double phi1   =  pos1.Phi() * (180. / TMath::Pi());
		double theta1 =  pos1.Theta();
		
		for(vector<const DCCALShower*>::const_iterator show2 = ccal_candidates.begin(); 
			show2 != ccal_candidates.end(); show2++) {
			
			double e2     = (*show2)->E;
			DVector3 pos2((*show2)->x1, (*show2)->y1, (*show2)->z);
			pos2         -= vertex;
			pos2         += ccal_correction;
			
			double t2     = (*show2)->time  -  (pos2.Mag()/c);
			double phi2   = pos2.Phi() * (180. / TMath::Pi());
			double theta2 = pos2.Theta();
			
			// calculate deltaPhi and deltaT:
			
			double deltaPhi = fabs(phi2 - phi1);
			double deltaT   = t2 - t1;
			
			// calculate separation distance of particles at FCAL z position:
			
			double x1_face  = ((m_fcalZ - vertex.Z()) / pos1.Z()) * pos1.X();
			double y1_face  = ((m_fcalZ - vertex.Z()) / pos1.Z()) * pos1.Y();
			double x2_face  = ((m_fcalZ - vertex.Z()) / pos2.Z()) * pos2.X();
			double y2_face  = ((m_fcalZ - vertex.Z()) / pos2.Z()) * pos2.Y();
			double deltaR   = sqrt(pow(x1_face-x2_face,2.0) 
				+ pow(y1_face-y2_face,2.0));
			
			// loop over beam photons:
			
			for(vector<const DBeamPhoton*>::const_iterator gam = beam_photons.begin();
				gam != beam_photons.end(); gam++) {
				
				double eb       = (*gam)->lorentzMomentum().E();
				double tb       = (*gam)->time();
				int tag_counter = (*gam)->dCounter;
				
				double brfdt = tb - rfTime;
				
				int bunch_val;
				
				if(fabs(brfdt) < 6.012)
					bunch_val = 1;
				else if((-(6.012 + 5.*4.008) <= brfdt && brfdt <= -(6.012 + 3.*4.008))
					|| ((6.012 + 3.*4.008) <= brfdt && brfdt <=  (6.012 + 5.*4.008)))
					bunch_val = 0;
				else 
					continue;
				
				//if(eb < BEAM_min_energy_cut) continue;
				
				double ecomp1 =  1. / ((1./eb)  +  (1./m_e)*(1.-cos(theta1)));
				double ecomp2 =  1. / ((1./eb)  +  (1./m_e)*(1.-cos(theta2)));
				double deltaE = (e1     + e2    ) - (eb+m_e);
				double deltaK = (ecomp1 + ecomp2) - (eb+m_e);
				
				double deltaK2 = m_e * sin(theta1+theta2) / 
					(sin(theta1) + sin(theta2) - sin(theta1+theta2));
				deltaK2       -= eb;
				
				DetectorSystem_t sys = (*gam)->dSystem;
				
				int e_cut = 1;
				//int e_cut = 0;
				//if(fabs(deltaE) < 3.0) e_cut = 1;
				
				if(e_cut) {
					
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
					
					if(sys==SYS_TAGH)      loc_Cand.tag_sys = 0;
					else if(sys==SYS_TAGM) loc_Cand.tag_sys = 1;
					
					candidates.push_back(loc_Cand);
				}
				
			} // end DBeamPhoton loop
			
		} // end DCCALShower loop
		
	} // end DFCALShower loop
	
	
	int n_cands = (int)candidates.size();
	
	dTreeFillData.Fill_Single<Int_t>("eventNum", eventnumber);
	dTreeFillData.Fill_Single<Double_t>("rfTime", rfTime);
	
	size_t candIndex = 0;
	for(int ic = 0; ic < n_cands; ic++) {
		
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
