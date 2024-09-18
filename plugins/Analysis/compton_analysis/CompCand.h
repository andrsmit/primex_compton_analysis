#ifndef _CompCand_
#define _CompCand_

typedef struct {
	
	double weight;
	
	double deltaE;
	double deltaPhi;
	double deltaK;
	
	double e1;
	double x1;
	double y1;
	double z1;
	
	double e2;
	double x2;
	double y2;
	double z2;
	
	double eb, brfdt;
	int tag_counter;
	int tag_sys;
	
	/*
	double deltaT;
	double dTheta;
	double deltaR;
	
	int ring;
	int pair_cut;
	int phi_slice;
	
	double vz;
	double event_weight;
	*/
} ComptonCandidate_t;

#endif
