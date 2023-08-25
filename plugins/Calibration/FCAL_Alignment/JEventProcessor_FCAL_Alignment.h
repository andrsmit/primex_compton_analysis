// $Id$
//
//    File: JEventProcessor_FCAL_Alignment.h
// Created: Fri Jan 21 20:10:28 EST 2022
// Creator: andrsmit (on Linux ifarm1802.jlab.org 3.10.0-1160.11.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_FCAL_Alignment_
#define _JEventProcessor_FCAL_Alignment_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

#include "TH1.h"
#include "TH2.h"

using namespace jana;
using namespace std;

#include <CCAL/DCCALShower.h>
#include <FCAL/DFCALShower.h>
#include <PID/DBeamPhoton.h>
#include <PID/DEventRFBunch.h>
#include <TOF/DTOFPoint.h>
#include <TRIGGER/DL1Trigger.h>
#include <HDGEOMETRY/DGeometry.h>
#include <ANALYSIS/DTreeInterface.h>

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

class JEventProcessor_FCAL_Alignment:public jana::JEventProcessor{
	public:
		JEventProcessor_FCAL_Alignment(){};
		~JEventProcessor_FCAL_Alignment(){};
		const char* className(void){return "JEventProcessor_FCAL_Alignment";}

	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void);
		jerror_t fini(void);
		
		int fcal_fiducial_cut(DVector3 pos, DVector3 vertex);
		int fcal_fiducial_cut2(DVector3 pos, DVector3 vertex);
		int ccal_fiducial_cut(DVector3 pos, DVector3 vertex);
		
		int check_TOF_match(DVector3 pos1, double rfTime, 
			DVector3 vertex, vector<const DTOFPoint*> tof_points);
		
		double m_fcalX, m_fcalY, m_fcalZ;
		double m_ccalX, m_ccalY, m_ccalZ;
		double m_beamX, m_beamY, m_beamZ;
		
		const double c   = 29.9792458;
		const double m_e = 0.510998928e-3;
		
		const double FCAL_min_energy_cut = 0.3;
		
		TH1F *h_beam_rf_dt, *h_fcal_rf_dt, *h_ccal_rf_dt;
		
		TH1F *h_deltaE,      *h_deltaK,      *h_deltaPhi;
		TH1F *h_deltaE_cut,  *h_deltaK_cut,  *h_deltaPhi_cut;
		TH1F *h_deltaE_cut2, *h_deltaK_cut2, *h_deltaPhi_cut2;
		TH1F *h_deltaE_cut3, *h_deltaK_cut3, *h_deltaPhi_cut3;
		TH1F *h_deltaE_cut4, *h_deltaK_cut4, *h_deltaPhi_cut4;
		
		TH1F *h_e1,      *h_e2,      *h_eb;
		TH1F *h_e1_cut,  *h_e2_cut,  *h_eb_cut;
		TH1F *h_e1_cut2, *h_e2_cut2, *h_eb_cut2;
		TH1F *h_e1_cut3, *h_e2_cut3, *h_eb_cut3;
		TH1F *h_e1_cut4, *h_e2_cut4, *h_eb_cut4;
		
		TH2F *h_xy1,      *h_xy2;
		TH2F *h_xy1_cut,  *h_xy2_cut;
		TH2F *h_xy1_cut2, *h_xy2_cut2;
		TH2F *h_xy1_cut3, *h_xy2_cut3;
		TH2F *h_xy1_cut4, *h_xy2_cut4;
		
		DTreeInterface *dTreeInterface;
		static thread_local DTreeFillData dTreeFillData;
};

#endif // _JEventProcessor_FCAL_Alignment_

