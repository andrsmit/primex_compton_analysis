// $Id$
//
//    File: JEventProcessor_compton_analysis.cc
// Created: Mon Nov 23 16:24:15 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_compton_analysis.h"


extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_compton_analysis());
}
} // "C"


//------------------
// init
//------------------
jerror_t JEventProcessor_compton_analysis::init(void)
{
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_compton_analysis::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_compton_analysis::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_compton_analysis::erun(void)
{
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_compton_analysis::fini(void)
{
	
	return NOERROR;
}

