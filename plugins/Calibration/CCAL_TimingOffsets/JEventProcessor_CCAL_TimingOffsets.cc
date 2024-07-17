// $Id$
//
//    File: JEventProcessor_CCAL_TimingOffsets.cc
// Created: Mon Nov 30 14:09:58 EST 2020
// Creator: andrsmit (on Linux ifarm1801.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_CCAL_TimingOffsets.h"

extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_CCAL_TimingOffsets());
}
} // "C"

//------------------
// init
//------------------
jerror_t JEventProcessor_CCAL_TimingOffsets::init(void)
{
	TDirectory *dir_CCAL_TimingOffsets = new TDirectoryFile( "CCAL_TimingOffsets", 
		"CCAL_TimingOffsets");
	dir_CCAL_TimingOffsets->cd();
	
	h_nccal_hits         = new TH1F( "nccal_hits", 
		"Number of DCCALHits",           50, -0.5, 49.5 );
	h_nccal_showerHits   = new TH1F( "nccal_showerHits", 
		"Number of DCCALHits in Shower", 50, -0.5, 49.5 );
	
	h_ccal_rf_dt          = new TH1F( "ccal_rf_dt", 
		"DCCALHit Time - RF Time; ns", 2000, -100., 100. );
	h_ccal_rf_dt_vs_chan  = new TH2F( "ccal_rf_dt_vs_chan", 
		"DCCALHit Time - RF Time; CCAL id.; [ns]", 144, -0.5, 143.5, 2000, -100., 100. );
	h_ccal_rf_dt_vs_E     = new TH2F( "ccal_rf_dt_vs_E", 
		"DCCALHit Time - RF Time; E_{CCAL} [GeV]; [ns]", 200, 0.0, 10.0, 2000, -10., 10. );
	
	h_ccal_beam_dt          = new TH1F( "ccal_beam_dt", 
		"DCCALHit Time - Beam Photon Time; ns", 2000, -100., 100. );
	h_ccal_beam_dt_vs_chan  = new TH2F( "ccal_beam_dt_vs_chan", 
		"DCCALHit Time - Beam Photon Time; CCAL id.; [ns]", 
		144, -0.5, 143.5, 2000, -100., 100. );
	h_ccal_beam_dt_vs_E     = new TH2F( "ccal_beam_dt_vs_E", 
		"DCCALHit Time - Beam Photon Time; E_{CCAL} [GeV]; [ns]", 
		200, 0.0, 10.0, 2000, -10., 10. );
	
	h_ccal_rf_dt_show         = new TH1F( "ccal_rf_dt_show", 
		"DCCALShower Time - RF Time; ns", 2000, -100., 100. );
	h_ccal_rf_dt_show_vs_chan = new TH2F( "ccal_rf_dt_show_vs_chan", 
		"DCCALShower Time - RF Time; CCAL id.; [ns]", 144, -0.5, 143.5, 2000, -100., 100. );
	h_ccal_rf_dt_show_vs_E    = new TH2F( "ccal_rf_dt_show_vs_E", 
		"DCCALShower Time - RF Time; E_{CCAL} [GeV]; [ns]", 200, 0.0, 10.0, 2000, -10., 10. );
	
	h_ccal_beam_dt_show         = new TH1F( "ccal_beam_dt_show", 
		"DCCALShower Time - Beam Photon Time; ns", 2000, -100., 100. );
	h_ccal_beam_dt_show_vs_chan = new TH2F( "ccal_beam_dt_show_vs_chan", 
		"DCCALShower Time - Beam Photon Time; CCAL id.; [ns]", 
		144, -0.5, 143.5, 2000, -100., 100. );
	h_ccal_beam_dt_show_vs_E    = new TH2F( "ccal_beam_dt_show_vs_E", 
		"DCCALShower Time - Beam Photon Time; E_{CCAL} [GeV]; [ns]", 
		200, 0.0, 10.0, 2000, -10., 10. );
	
	dir_CCAL_TimingOffsets->cd("../");
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_CCAL_TimingOffsets::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	DGeometry*   dgeom = NULL;
	DApplication* dapp = dynamic_cast< DApplication* >( eventLoop->GetJApplication() );
	if(dapp)     dgeom = dapp->GetDGeometry( runnumber );
	
	if(dgeom){
		dgeom->GetTargetZ( m_beamZ );
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
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_CCAL_TimingOffsets::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
	//----------------------------------------------------------//
	// Reject FP triggers:
	
	const DL1Trigger *trig = NULL;
	try{
		eventLoop->GetSingle(trig);
	} catch (...) {}
	if(trig == NULL) { return NOERROR; }
	
	uint32_t fp_trigmask = trig->fp_trig_mask;
	if(fp_trigmask) return NOERROR;
	
	//----------------------------------------------------------//
	// Get the RF time using "CalorimeterOnly":
	
	const DEventRFBunch *locRFBunch = NULL;
	try { 
		eventLoop->GetSingle(locRFBunch, "CalorimeterOnly");
	} catch (...) { return NOERROR; }
	double rfTime = locRFBunch->dTime;
	if(locRFBunch->dNumParticleVotes < 2) return NOERROR;
	
	//----------------------------------------------------------//
	// Get DCCALGeometry object:
	
	vector<const DCCALGeometry*> ccalGeomVec;
	eventLoop->Get(ccalGeomVec);
	if(ccalGeomVec.size() != 1) {
		cerr << "No CCAL geometry accessbile." << endl;
		return RESOURCE_UNAVAILABLE;
	}
	const DCCALGeometry* ccalGeom = ccalGeomVec[0];
	
	//----------------------------------------------------------//
	// Get CCAL showers and hits:
	
	vector<const DCCALHit*> ccal_hits;
	eventLoop->Get(ccal_hits);
	
	vector<const DCCALShower*> ccal_showers;
	eventLoop->Get(ccal_showers);
	
	vector<const DBeamPhoton*> beam_photons;
	eventLoop->Get(beam_photons);
	
	int n_ccal_hits = (int)ccal_hits.size();
	
	DVector3 vertex;
	vertex.SetXYZ(m_beamX, m_beamY, m_beamZ);
	
	//----------------------------------------------------------//
	
	japp->RootFillLock(this);  // Acquire root lock
	
	int isGoodEvent = 0;
	
	for(vector<const DCCALShower*>::const_iterator show = ccal_showers.begin(); show != ccal_showers.end(); show++) {
		
		double energy = (*show)->E;
		if(energy < MIN_CCAL_ENERGY_SHOW) continue;
		
		DVector3 show_pos;
		show_pos.SetXYZ((*show)->x1 - vertex.X(), (*show)->y1 - vertex.Y(), (*show)->z - vertex.Z());
		
		double show_t = (*show)->time - (show_pos.Mag()/c);
		
		h_ccal_rf_dt_show->Fill(show_t - rfTime);
		h_ccal_rf_dt_show_vs_chan->Fill((*show)->idmax, show_t - rfTime);
		h_ccal_rf_dt_show_vs_E->Fill(energy, show_t - rfTime);
		
		//----------   Loop over hits belonging to this shower   ----------//
		
		vector<DCCALHit> shower_hits = (*show)->hitsInCluster;
		
		int loc_n_hits = (int)shower_hits.size();
		if( loc_n_hits < 1 ) continue;
		
		// only use the hit with maximum energy deposition:
		
		int max_hit_index = -1;
		double max_hit_energy = 0.0;
		
		for(int ihit = 0; ihit < loc_n_hits; ihit++) {
			double loc_hit_energy = (double)shower_hits[ihit].E;
			if( loc_hit_energy > max_hit_energy ) {
				max_hit_energy = loc_hit_energy;
				max_hit_index  = ihit;
			}
		}
		
		if(max_hit_index < 0) continue;
		
		DCCALHit ccalHit  = shower_hits[max_hit_index];
		
		double  hitE      =  1.e-3 * ccalHit.E;
		if(hitE < MIN_CCAL_ENERGY_HIT) continue;
		
		int ChannelNumber =  ccalGeom->channel(ccalHit.row, ccalHit.column);
		double chanx      =  ccalHit.x - vertex.X() + m_ccalX;
		double chany      =  ccalHit.y - vertex.Y() + m_ccalY;
		double chanz      =  show_pos.Z();
		double hitTime    =  ccalHit.t;
		
		double dR = sqrt(chanx*chanx + chany*chany + chanz*chanz);
		
		// propagate hit time to the interaction vertex:
		
		double tCorr = (m_ccalZ + DCCALGeometry::blockLength() - (*show)->z) / CCAL_C_EFFECTIVE;
		
		hitTime = hitTime - tCorr - (dR/c);
		
		double ccal_rf_dt = hitTime - rfTime;
		
		h_ccal_rf_dt->Fill(ccal_rf_dt);
		h_ccal_rf_dt_vs_chan->Fill(ChannelNumber, ccal_rf_dt);
		h_ccal_rf_dt_vs_E->Fill(hitE, ccal_rf_dt);
		
		h_nccal_showerHits->Fill(loc_n_hits);
		
		// Also look at time difference between beam photon:
		
		for(unsigned int ib = 0; ib < beam_photons.size(); ib++) {
			
			double beam_time    = beam_photons[ib]->time();
			double ccal_beam_dt = hitTime - beam_time;
			
			h_ccal_beam_dt->Fill(ccal_beam_dt);
			h_ccal_beam_dt_vs_chan->Fill(ChannelNumber, ccal_beam_dt);
			h_ccal_beam_dt_vs_E->Fill(hitE, ccal_beam_dt);
			
			h_ccal_beam_dt_show->Fill(show_t - beam_time);
			h_ccal_beam_dt_show_vs_chan->Fill(ChannelNumber, show_t - beam_time);
			h_ccal_beam_dt_show_vs_E->Fill(energy, show_t - beam_time);
		}
		
		isGoodEvent++;
	}
	
	if(isGoodEvent) h_nccal_hits->Fill(n_ccal_hits);
	
	japp->RootFillUnLock(this);  // Release root lock
	
	return NOERROR;
}
