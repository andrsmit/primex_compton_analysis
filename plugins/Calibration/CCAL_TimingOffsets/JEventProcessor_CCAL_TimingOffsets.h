// $Id$
//
//    File: JEventProcessor_CCAL_TimingOffsets.h
// Created: Mon Nov 30 14:09:58 EST 2020
// Creator: andrsmit (on Linux ifarm1801.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_CCAL_TimingOffsets_
#define _JEventProcessor_CCAL_TimingOffsets_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

using namespace jana;

#include <CCAL/DCCALGeometry.h>
#include <CCAL/DCCALShower.h>
#include <CCAL/DCCALHit.h>
#include <PID/DBeamPhoton.h>
#include <PID/DEventRFBunch.h>
#include <TRIGGER/DL1Trigger.h>
#include <HDGEOMETRY/DGeometry.h>

#include "units.h"
#include "DLorentzVector.h"
#include "DVector3.h"
#include "TRandom3.h"

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <thread>
#include <mutex>

class JEventProcessor_CCAL_TimingOffsets:public jana::JEventProcessor{
	public:
		JEventProcessor_CCAL_TimingOffsets() {};
		~JEventProcessor_CCAL_TimingOffsets(){};
		const char* className(void){return "JEventProcessor_CCAL_TimingOffsets";}

	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void) {return NOERROR;};
		jerror_t fini(void) {return NOERROR;};
		
		double m_beamX, m_beamY, m_beamZ;
		double m_ccalX, m_ccalY, m_ccalZ;
		
		const double c   = 29.9792458;
		const double m_e = 0.510998928e-3;
		
		const double CCAL_C_EFFECTIVE      = 13.6;   // [cm/ns]
		
		const double CCAL_RADIATION_LENGTH = 0.86;   // [cm]
		const double CCAL_CRITICAL_ENERGY  = 1.1e-3; // [GeV]
		
		const double MIN_CCAL_ENERGY_SHOW  = 1.0;    // [GeV]
		const double MIN_CCAL_ENERGY_HIT   = 0.015;   // [GeV]
		
		TH1F *h_nccal_hits,  *h_nccal_showerHits;
		TH1F *h_ccal_rf_dt;
		TH2F *h_ccal_rf_dt_vs_chan,   *h_ccal_rf_dt_vs_E;
		TH1F *h_ccal_beam_dt;
		TH2F *h_ccal_beam_dt_vs_chan, *h_ccal_beam_dt_vs_E;
		
		TH1F *h_ccal_rf_dt_show;
		TH2F *h_ccal_rf_dt_show_vs_chan,   *h_ccal_rf_dt_show_vs_E;
		TH1F *h_ccal_beam_dt_show;
		TH2F *h_ccal_beam_dt_show_vs_chan, *h_ccal_beam_dt_show_vs_E;
		
};

#endif // _JEventProcessor_CCAL_TimingOffsets_

