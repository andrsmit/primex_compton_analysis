
#include "pairs.h"
#include "ComptonAna.h"

int main(int argc, char **argv) {
	
	genSettings_t genSettings;
	genSettings.run_number    = 61321;
	genSettings.file_ext_low  =     0;
	genSettings.file_ext_high =   600;
	genSettings.input_fname   = "none";
	genSettings.output_fname  = "pairs_systematics.root";
	
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
			case 'b':
				genSettings.file_ext_low = atoi(++argptr);
				break;
			case 'e':
				genSettings.file_ext_high = atoi(++argptr);
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
	
	// Directory where ROOT Trees are stored:
	sprintf(rootTree_pathName, 
		"/work/halld/home/andrsmit/primex_compton_analysis/bhgen_test/recRootTrees/Run%06d/trees", 
		genSettings.run_number);
	
	// Directory where output ROOT files will be stored:
	sprintf(rootFile_pathName, 
		"/work/halld/home/andrsmit/primex_compton_analysis/analyze_trees/phase1/analyze_mc/rootFiles/Run%06d/pair/systematics", 
		genSettings.run_number);
	
	// Construct analysis object:
	
	ComptonAna locAna;
	
	locAna.m_SHIFT_DISTRIBUTIONS = 1;
	locAna.m_SMEAR_DISTRIBUTIONS = 1;
	
	locAna.setRunNumber(genSettings.run_number);
	
	// Initialize histograms to be filled:
	
	locAna.initHistograms_systematics();
	//
	// Check if an input filename was specificed at runtime. 
	// If not, we'll do a loop over files from the rootTree directory above:
	//
	if(genSettings.input_fname!="none") {
		
		TString input_fname = Form("%s", genSettings.input_fname.c_str());
		if(gSystem->AccessPathName(input_fname.Data())) {
			fprintf(stderr,"Specified input filename is inaccessible.\n");
			exit(0);
		}
		locAna.setOutputFileName(Form("%s",genSettings.output_fname.c_str()));
		locAna.runAnalysis_systematics(input_fname.Data());
		locAna.writeHistograms_systematics();
		
	} else {
		
		int min_processed = genSettings.file_ext_low;
		int max_processed = genSettings.file_ext_high;
		
		int first = 1;
		
		for(int loc_ext = genSettings.file_ext_low; loc_ext <= genSettings.file_ext_high; loc_ext++) {
			
			TString input_fname = Form("%s/pair_rec_%04d.root", rootTree_pathName, loc_ext);
			
			// check if file exists:
			if(gSystem->AccessPathName(input_fname.Data())) continue;
			cout << "  processing ext " << loc_ext << endl;
			
			if(first) {
				min_processed = loc_ext;
				first = 0;
			}
			max_processed = loc_ext;
			locAna.runAnalysis_systematics(input_fname.Data());
		}
		
		// write results to file:
		locAna.setOutputFileName(Form("%s/pair_rec_%04d_%04d.root", rootFile_pathName, min_processed, max_processed));
		locAna.writeHistograms_systematics();
	}
	
	return 0;
}

void printUsage(genSettings_t genSettings, int goYes) {
	
	if(goYes==0) {
		fprintf(stderr,"\nSWITCHES:\n");
		fprintf(stderr,"-h\tPrintthis message\n");
		fprintf(stderr,"-r<arg>\tRun number for simulated data\n");
		fprintf(stderr,"-b<arg>\tLower end of file extension number range to process\n");
		fprintf(stderr,"-e<arg>\tUpper end of file extension number range to process\n");
		fprintf(stderr,"-i<arg>\tInput file name (default is none)\n");
		fprintf(stderr,"-o<arg>\tOutput file name (default is compton_ana.root)\n\n\n");
	}
	
	if(goYes==1) {
		if(genSettings.input_fname!="none") {
			printf("\nAnalyzing simulations for run %d, input file: %s\n", genSettings.run_number, 
				genSettings.input_fname.c_str());
			cout << "" << endl;
		} else {
			printf("\nAnalyzing simulations for run %d, files %d-%d:\n", genSettings.run_number, genSettings.file_ext_low, 
				genSettings.file_ext_high);
			cout << "" << endl;
		}
	}
	
	if(goYes==0) exit(0);
	
	return;
}
