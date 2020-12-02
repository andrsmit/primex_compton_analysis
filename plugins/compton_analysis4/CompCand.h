#ifndef _CompCand_
#define _CompCand_

typedef struct {
	
	int bunch_val;
	
	double deltaT;	
	double deltaPhi;
	double dTheta;
	double deltaR;
	double deltaE;
	double deltaK;
	double deltaK2;
	
	double e1;
	double x1;
	double y1;
	double z1;
	
	double e2;
	double x2;
	double y2;
	double z2;
	
	int ring;
	int pair_cut;
	int phi_slice;
	
	double eb;
	int tag_counter;
	int tag_sys;
	
} ComptonCandidate_t;

#endif
