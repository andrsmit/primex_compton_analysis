#ifndef COMPTON_INCLUDE
#define COMPTON_INCLUDE

using namespace std;

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "TVector3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TSystem.h"
#include "TDirectory.h"

//---------------------------------------------//
// Structure for the configuration settings:

struct genSettings_t {
	int run_number;       //
	int tag_sys;          // 0 -> TAGH; 1-> TAGM
	int tag_counter_low;  //
	int tag_counter_high; //
	int beam_current;
	string  input_fname;
	string output_fname;
};

void printUsage(genSettings_t, int goYes);

//----------   Numerical Constants   ----------//

const double c   = 29.9792; // [cm/ns]
const double pi  = 3.1415926535;
const double m_e = 0.510998928e-3;
/*
double m_fcalX = 0.617, m_fcalY = -0.002;
const double m_fcalZ = 624.906;

double m_ccalX = 0.083, m_ccalY =  0.148;
const double m_ccalZ = 1279.376;

double m_beamX, m_beamY, m_beamZ;
*/
//----------   Function Declarations   ----------//

int main(int argc, char **argv);

//----------   Data Objects   ----------//

char rootTree_pathName[256], rootFile_pathName[256];

#endif
