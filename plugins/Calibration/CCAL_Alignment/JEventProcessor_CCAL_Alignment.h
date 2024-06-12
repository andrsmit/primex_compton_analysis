// $Id$
//
//    File: JEventProcessor_CCAL_Alignment.h
// Created: Fri Jun  7 10:39:01 AM EDT 2024
// Creator: andrsmit (on Linux ifarm180302.jlab.org 5.14.0-362.24.2.el9_3.x86_64 x86_64)
//

#ifndef _JEventProcessor_CCAL_Alignment_
#define _JEventProcessor_CCAL_Alignment_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

#include "TH1.h"
#include "TH2.h"

using namespace jana;
using namespace std;

#include <CCAL/DCCALShower.h>
#include <FCAL/DFCALShower.h>
#include <BCAL/DBCALShower.h>
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

class JEventProcessor_CCAL_Alignment:public jana::JEventProcessor{
	public:
		JEventProcessor_CCAL_Alignment(){};
		~JEventProcessor_CCAL_Alignment(){};
		const char* className(void){return "JEventProcessor_CCAL_Alignment";}

	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void){return NOERROR;};
		jerror_t fini(void);
		
		int fcal_fiducial_cut(DVector3 pos, DVector3 vertex, double cut_layer);
		int ccal_fiducial_cut(DVector3 pos, DVector3 vertex);
		
		void check_TOF_match(DVector3 pos, double rfTime, DVector3 vertex, 
			vector<const DTOFPoint*> tof_points, double &dx_min, double &dy_min, 
			double &dt_min, double rf_time_cut);
		
		double m_fcalX, m_fcalY, m_fcalZ;
		double m_ccalX, m_ccalY, m_ccalZ;
		double m_beamX, m_beamY, m_beamZ;
		
		const double c   = 29.9792458;
		const double m_e = 0.510998928e-3;
		
		static const int n_cuts = 5;
		TH1F *h_deltaE[n_cuts];
		TH1F *h_deltaK[n_cuts];
		TH1F *h_deltaPhi[n_cuts];
		TH1F *h_e1[n_cuts];
		TH1F *h_e2[n_cuts];
		TH1F *h_eb[n_cuts];
		TH2F *h_xy1[n_cuts];
		TH2F *h_xy2[n_cuts];
		
		DTreeInterface *dTreeInterface;
		static thread_local DTreeFillData dTreeFillData;
};

#endif // _JEventProcessor_CCAL_Alignment_

