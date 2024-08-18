
#include "compton.h"
#include "ComptonAna.h"

int main(int argc, char **argv) {
	
	genSettings_t genSettings;
	genSettings.min_run      = 61321;
	genSettings.max_run      = 61321;
	genSettings.run_number   =     0;
	genSettings.input_fname  = "none";
	genSettings.output_fname = "compton_systematics.root";
	
	// parse command line:
	char *argptr;
	for(int iarg=1; iarg<argc; iarg++) {
		argptr = argv[iarg];
		if(*argptr == '-') {
			argptr++;
			switch(*argptr) {
			case 'b':
				genSettings.min_run  = atoi(++argptr);
				break;
			case 'e':
				genSettings.max_run = atoi(++argptr);
				break;
			case 'i':
				genSettings.input_fname = ++argptr;
				break;
			case 'o':
				genSettings.output_fname = ++argptr;
				break;
			case 'r':
				genSettings.run_number = atoi(++argptr);
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
	
	// Directory where output ROOT files will be stored:
	sprintf(rootFile_pathName, 
		"/work/halld/home/andrsmit/primex_compton_analysis/analyze_trees/analyze_data/rootFiles/systematics", 
		genSettings.run_number);
	
	// Initialize histograms to be filled:
	
	ComptonAna locAna;
	locAna.initHistograms_systematics();
	
	//
	// Check if an input filename was specificed at runtime. 
	// If not, we'll do a loop over run numbers from the rootTree directory above:
	//
	if(genSettings.input_fname!="none") {
		
		if(genSettings.run_number==0) {
			fprintf(stderr,"Need to specify the runnumber for the file you are processing.\n");
			exit(0);
		}
		TString input_fname = Form("%s", genSettings.input_fname.c_str());
		if(gSystem->AccessPathName(input_fname.Data())) {
			fprintf(stderr,"Specified input filename is inaccessible.\n");
			exit(0);
		}
		locAna.setRunNumber(genSettings.run_number);
		locAna.setOutputFileName(Form("%s",genSettings.output_fname.c_str()));
		locAna.runAnalysis_systematics(input_fname.Data());
		
	} else {
		
		for(int loc_run = genSettings.min_run; loc_run <= genSettings.max_run; loc_run++) {
			
			// set rootTree path according to phase:
			sprintf(rootTree_pathName, "/work/halld/home/andrsmit/primex_compton_analysis/data/rootTrees/phase%d/%06d", 
				locAna.getPrimexPhase(loc_run), loc_run);
			
			// check if root tree directory exists for this runnumber. If not, skip it:
			if(gSystem->AccessPathName(rootTree_pathName)) continue;
			
			// check if output file already exists for this runnumber:
			TString output_fname = Form("%s/%d.root", rootFile_pathName, loc_run);
			if(!gSystem->AccessPathName(output_fname.Data())) continue;
			
			locAna.setRunNumber(loc_run);
			locAna.setOutputFileName(output_fname.Data());
			
			auto start = chrono::high_resolution_clock::now();
			
			// loop over all extensions in this directory:
			for(int loc_ext = 0; loc_ext < 250; loc_ext++) {
				
				// check if file exists:
				TString input_fname = Form("%s/%06d_%03d.root", rootTree_pathName, loc_run, loc_ext);
				if(gSystem->AccessPathName(input_fname.Data())) continue;
				
				cout << "  processing ext " << loc_ext << endl;
				locAna.runAnalysis_systematics(input_fname.Data());
			}
			
			auto stop = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
			//cout << duration.count() << endl;
			
			locAna.writeHistograms_systematics();
			locAna.resetHistograms_systematics();
		}
		cout << "Done." << endl;
	}
	
	return 0;
}

void printUsage(genSettings_t genSettings, int goYes) {
	
	if(goYes==0) {
		fprintf(stderr,"\nSWITCHES:\n");
		fprintf(stderr,"-h\tPrintthis message\n");
		fprintf(stderr,"-b<arg>\tMinimum run number to process\n");
		fprintf(stderr,"-e<arg>\tMaximum run number to process\n");
		fprintf(stderr,"-r<arg>\tRun number for specified input file\n");
		fprintf(stderr,"-i<arg>\tInput file name (default is none)\n");
		fprintf(stderr,"-o<arg>\tOutput file name (default is compton_ana.root)\n\n\n");
	}
	
	if(goYes==1) {
		if(genSettings.input_fname!="none") {
			printf("\nAnalyzing simulations for run %d, input file: %s\n", genSettings.run_number, 
				genSettings.input_fname.c_str());
			cout << "" << endl;
		} else {
			printf("\nAnalyzing simulations for runs %d-%d:\n", genSettings.min_run, genSettings.max_run);
			cout << "" << endl;
		}
	}
	
	if(goYes==0) exit(0);
	
	return;
}
