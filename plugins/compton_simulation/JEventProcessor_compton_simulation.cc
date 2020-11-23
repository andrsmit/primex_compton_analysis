// $Id$
//
//    File: JEventProcessor_compton_simulation.cc
// Created: Mon Nov 23 16:28:33 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_compton_simulation.h"


extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_compton_simulation());
}
} // "C"


//------------------
// init
//------------------
jerror_t JEventProcessor_compton_simulation::init(void)
{
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_compton_simulation::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_compton_simulation::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_compton_simulation::erun(void)
{
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_compton_simulation::fini(void)
{
	
	return NOERROR;
}

