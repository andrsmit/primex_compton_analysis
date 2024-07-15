// $Id$
//
//    File: JEventProcessor_compton_tree.h
// Created: Wed Dec  9 16:09:48 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_compton_tree_
#define _JEventProcessor_compton_tree_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

using namespace jana;
using namespace std;

#include <CCAL/DCCALShower.h>
#include <FCAL/DFCALShower.h>
#include <TOF/DTOFPoint.h>
#include <PID/DBeamPhoton.h>
#include <PID/DEventRFBunch.h>
#include <PID/DMCReaction.h>
#include <TRACKING/DMCThrown.h>

#include <CCAL/DCCALGeometry.h>
#include <FCAL/DFCALGeometry.h>

#include <TRIGGER/DL1Trigger.h>
#include <TRIGGER/DTrigger.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df250PulseData.h>

#include <HDGEOMETRY/DGeometry.h>
#include <ANALYSIS/DTreeInterface.h>

#include "units.h"
#include "DLorentzVector.h"
#include "DVector3.h"
#include "TRandom3.h"
#include "TSystem.h"

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <thread>
#include <mutex>

class JEventProcessor_compton_tree:public jana::JEventProcessor{
	public:
		JEventProcessor_compton_tree();
		~JEventProcessor_compton_tree(){};
		const char* className(void){return "JEventProcessor_compton_tree";}

	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void){return NOERROR;};
		jerror_t fini(void);
		
		void write_events(uint64_t eventnumber, double rfTime, 
			vector<const DMCThrown*> mc_thrown, vector<const DMCReaction*> mc_reaction);
		
		void write_events(uint64_t eventnumber, double rfTime,
			vector<const DBeamPhoton*> beam_photons, 
			vector<const DFCALShower*> fcal_showers,
			vector<const DCCALShower*> ccal_showers,
			vector<const DTOFPoint*> tof_points,
			vector<const DMCThrown*> mc_thrown,
			vector<const DMCReaction*> mc_reaction);
		
		double get_acc_scaling_factor(double eb);
		
		//----------   Constants   ----------//
		
		double m_beamX, m_beamY, m_beamZ;
		double m_fcalX, m_fcalY, m_fcalZ;
		double m_ccalX, m_ccalY, m_ccalZ;
		
		DVector3 m_fcal_correction, m_ccal_correction;
		
		const double c = 29.9792458; // cm/ns
		
		double m_RFTimeCut, m_BeamEnergyCut, m_DeltaECut;
		
		TH1F *h_gen_flux;
		TH2F *h_gen_weight;
		
		double m_HodoscopeHiFactor    = 1.0;
		double m_HodoscopeHiFactorErr = 1.0;
		double m_HodoscopeLoFactor    = 1.0;
		double m_HodoscopeLoFactorErr = 1.0;
		double m_MicroscopeFactor     = 1.0;
		double m_MicroscopeFactorErr  = 1.0;
		double m_TAGMEnergyBoundHi    = 1.0;
		double m_TAGMEnergyBoundLo    = 1.0;
		
		//----------   TTree   ----------//
		
		DTreeInterface *dTreeInterface;
		static thread_local DTreeFillData dTreeFillData;
		
		//-----------------------------------//
};

#endif // _JEventProcessor_compton_tree_

