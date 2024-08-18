
#include "compton.h"
#include "ComptonAna.h"

int main(int argc, char **argv) {
	
	genSettings_t genSettings;
	genSettings.run_number       = 61321;
	genSettings.tag_sys          =     0;
	genSettings.tag_counter_low  =   274;
	genSettings.tag_counter_high =   274;
	genSettings.beam_current     =   200;
	genSettings.input_fname      = "none";
	genSettings.output_fname     = "compton_ana.root";
	
	// parse command line:
	char *argptr;
	for(int iarg=1; iarg<argc; iarg++) {
		argptr = argv[iarg];
		if(*argptr == '-') {
			argptr++;
			switch(*argptr) {
			case 'r':
				genSettings.run_number = atoi(++argptr);
				break;
			case 's':
				genSettings.tag_sys = atoi(++argptr);
				break;
			case 'b':
				genSettings.tag_counter_low  = atoi(++argptr);
				break;
			case 'e':
				genSettings.tag_counter_high = atoi(++argptr);
				break;
			case 'c':
				genSettings.beam_current = atoi(++argptr);
				break;
			case 'i':
				genSettings.input_fname = ++argptr;
				break;
			case 'o':
				genSettings.output_fname = ++argptr;
				break;
			case 'h':
				printUsage(genSettings,0);
				break;
			default:
				fprintf(stderr,"Unrecognized argument: [-%s]\n\n",argptr);
				printUsage(genSettings,0);
				break;
			}
		}
	}
	printUsage(genSettings, 1);
	
	char loc_path[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	if(genSettings.beam_current==0) {
		
		// Directory where ROOT Trees are stored:
		
		sprintf(rootTree_pathName, 
			"%s/compton_mc/phase1/default_geometry/Run%06d/recRootTrees_noBkgd", 
			loc_path, genSettings.run_number);
		
		// Directory where output ROOT files will be stored:
		
		sprintf(rootFile_pathName, 
			"%s/analyze_trees/phase1/analyze_mc/rootFiles/Run%06d/compton/noBkgd", 
			loc_path, genSettings.run_number);
		
	} else {
		
		// Directory where ROOT Trees are stored:
		
		sprintf(rootTree_pathName, 
			"%s/compton_mc/phase1/default_geometry/Run%06d/recRootTrees_%03dnA", 
			loc_path, genSettings.run_number, genSettings.beam_current);
		
		// Directory where output ROOT files will be stored:
		
		sprintf(rootFile_pathName, 
			"%s/analyze_trees/phase1/analyze_mc/rootFiles/Run%06d/compton/%03dnA", 
			loc_path, genSettings.run_number, genSettings.beam_current);
	}
	
	// Construct analysis object:
	
	ComptonAna locAna;
	
	locAna.m_SHIFT_DISTRIBUTIONS = 1;
	locAna.m_SMEAR_DISTRIBUTIONS = 1;
	
	locAna.setRunNumber(genSettings.run_number);
	
	// Initialize histograms to be filled:
	
	locAna.initHistograms();
	
	//
	// Check if an input filename was specificed at runtime. 
	// If not, we'll do a loop over TAGGER counters from the rootTree directory above:
	//
	if(genSettings.input_fname!="none") {
		
		TString input_fname = Form("%s", genSettings.input_fname.c_str());
		if(gSystem->AccessPathName(input_fname.Data())) {
			fprintf(stderr,"Specified input filename is inaccessible.\n");
			exit(0);
		}
		locAna.setOutputFileName(Form("%s",genSettings.output_fname.c_str()));
		locAna.runAnalysis(input_fname.Data());
		
	} else {
		
		string loc_tag_sys = genSettings.tag_sys==0 ? "tagh" : "tagm";
		for(int loc_counter = genSettings.tag_counter_low; loc_counter <= genSettings.tag_counter_high; loc_counter++) {
			
			TString input_fname = Form("%s/%s_%03d.root", rootTree_pathName, loc_tag_sys.c_str(), loc_counter);
			
			// check if file exists:
			if(gSystem->AccessPathName(input_fname.Data())) continue;
			
			TString output_fname = Form("%s/%s_%03d.root", rootFile_pathName, loc_tag_sys.c_str(), loc_counter);
			
			// skip this counter if the output file already exists:
			if(!gSystem->AccessPathName(output_fname.Data())) continue;
			
			locAna.setOutputFileName(output_fname.Data());
			locAna.resetHistograms();
			locAna.runAnalysis(input_fname.Data());
			locAna.writeHistograms();
		}
		cout << "Done." << endl;
	}
	
	return 0;
}

void printUsage(genSettings_t genSettings, int goYes) {
	
	if(goYes==0) {
		fprintf(stderr,"\nSWITCHES:\n");
		fprintf(stderr,"-h\tPrintthis message\n");
		fprintf(stderr,"-r<arg>\tRun number for simulated data\n");
		fprintf(stderr,"-s<arg>\tTagger System (0 for TAGH, 1 for TAGM)\n");
		fprintf(stderr,"-b<arg>\tLower end of tagger counter number range to process\n");
		fprintf(stderr,"-e<arg>\tUpper end of tagger counter number range to process\n");
		fprintf(stderr,"-c<arg>\tBeam current of random trigger events mixed into simulation (default is 200)\n");
		fprintf(stderr,"-i<arg>\tInput file name (default is none)\n");
		fprintf(stderr,"-o<arg>\tOutput file name (default is compton_ana.root)\n\n\n");
	}
	
	if(goYes==1) {
		if(genSettings.input_fname!="none") {
			printf("\nAnalyzing simulations for run %d, input file: %s\n", genSettings.run_number, 
				genSettings.input_fname.c_str());
			cout << "" << endl;
		} else {
			string tsys = genSettings.tag_sys==0 ? "TAGH" : "TAGM";
			printf("\nAnalyzing simulations for run %d, %s counters %d-%d:\n", genSettings.run_number, tsys.c_str(), 
				genSettings.tag_counter_low, genSettings.tag_counter_high);
			cout << "" << endl;
		}
	}
	
	if(goYes==0) exit(0);
	
	return;
}
