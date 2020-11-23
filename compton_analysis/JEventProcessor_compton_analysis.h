// $Id$
//
//    File: JEventProcessor_compton_analysis.h
// Created: Mon Nov 23 16:24:15 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_compton_analysis_
#define _JEventProcessor_compton_analysis_

using namespace jana;

#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>

class JEventProcessor_compton_analysis:public jana::JEventProcessor{
	public:
		JEventProcessor_compton_analysis()  {};
		~JEventProcessor_compton_analysis() {};
		const char* className(void){return "JEventProcessor_compton_analysis";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
		
		
		
};

#endif // _JEventProcessor_compton_analysis_

