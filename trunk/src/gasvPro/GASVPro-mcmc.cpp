/**
 * Copyright 2012 Benjamin Raphael, Suzanne Sindi, Anthony Cannistra, Hsin-Ta Wu, Luke Peng, Selim Onal
 *
 *  This file is part of the GASVPro code distribution.
 * 
 *  GASVPro is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  GASVPro is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with gasv.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// % given a set of clusters  - samples over ambigously mapped fragments 
#define DEBUGGING
#undef DEBUGGING   // comment out to enable outputs to console

#include <vector>
#include <limits.h>
#include <algorithm>
#include <list>

#include <ctime>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>

#include "problemInstance.h"
	
using namespace std;

//Global Parameters.

int total_num_of_sv_clusters, total_num_of_fragments, dist_num_of_fragments, num_of_bases,threshold_for_unique_initial;
string fn; // filename holder 
int temp_int;
string temp_str;
double mut_rate;
double threshold_for_sampling = 1;  // the code ignores connected components that have less than this number of fragments

int Tokenize(const string& str,vector<string>& tokens, const string& delimiters = " "){
	int retVal = 0;
	string::size_type lastPos = str.find_first_not_of(delimiters, 0); // Skip delimiters at beginning.
    string::size_type pos     = str.find_first_of(delimiters, lastPos); // Find first "non-delimiter".
    while (string::npos != pos || string::npos != lastPos){
        tokens.push_back(str.substr(lastPos, pos - lastPos));// Found a token, add it to the vector.
		retVal++;
		lastPos = str.find_first_not_of(delimiters, pos); // Skip delimiters.  Note the "not_of"
        pos = str.find_first_of(delimiters, lastPos);// Find next "non-delimiter"
    }
	return retVal; //This gives you the number of tokens!
}

int findLocation(vector<pairedDependency>& locations, int numDependencies, pairedDependency key){	
	int first = 0;
	int last = numDependencies - 1;
	int mid = -1;
	int FOUND = 0;
	while (first <= last && FOUND == 0) {
		mid = (first + last) / 2;  // compute mid point.
		if (key > locations[mid]) 
			first = mid + 1;  // repeat search in top half.
		else if (key < locations[mid]) 
			last = mid - 1; // repeat search in bottom half.
		else
			FOUND = 1;     // found it. return position 
   	}
	// failed to find key
	if(FOUND == 0){ return -(first + 1);   }
	else{  return mid; }
}

int findLocation(vector<int>& locations, int numLocations, int key){
	int first = 0;
	int last = numLocations-1;
	int mid = -1;
	int FOUND = 0;
	while (first <= last && FOUND == 0) {
		mid = (first + last) / 2;  // compute mid point.
		if (key > locations[mid]) 
			first = mid + 1;  // repeat search in top half.
		else if (key < locations[mid]) 
			last = mid - 1; // repeat search in bottom half.
		else
			FOUND = 1;     // found it. return position /////
   	}
	// failed to find key
	if(FOUND == 0){ return -(first + 1);  }
	else{ return mid; }
}

int main(int argc, char* argv[] ){
	int START = -1;
	int END = -1;
	int MAX_SIZE_OF_CLUSTER = -1;
	string SUFFIX;
	
	int DO_UNIQUE = 0;
	
	if(argc < 3){
		cerr << "GASVPro-mcmc: MCMC sampling of assignments for PR with multiple mappigns\n";
		cerr << "Version:      1.2.1\n\n";
		
		cerr << "Usage: ./GASVPro-mcmc {ParametersFile} {Dir} (Optional:{START_CLUSTER} {END_CLUSTER})\n";
		
		cerr << "Required Parameters:\n";
		cerr << "\t\tParametersFile			   See GASVPro parameter file format\n";
		cerr << "\t\t<dir>                     GASVPro-graph Output Directory\n";
				
		cerr << "Optional Parameters:\n";
		cerr << "\t\tStart                     Cluster component to start with (-1 to do all)\n";
		cerr << "\t\tEnd                       Cluster component to end with (-1 to do all)\n";
		cerr << "\t\tNote: Start and End allow for parallel processing.\n\n";
	
		return -1;
	}
	else{
		if(argc>=5){
			START = atoi(argv[3]);
			END   = atoi(argv[4]);
		}
		
		/*
		if(argc>=5){ 
			DO_UNIQUE = atoi(argv[4]);
			if(! (DO_UNIQUE == 0 || DO_UNIQUE == 1) ){ 
				cerr << "Unique Flag can only be 0 or 1.\n"; 	
				cerr << "Proposed Unique Flag is --> " << DO_UNIQUE << endl;
				return -1;
			}
		}
		*/
	}
	
	/*
	cout << "To Do Notes:\n";
	cout << "(1) Current moves are add/remove and naive. Swap move is not presently supported.\n";
	*/
		
	/*
	 Notes to Tony:
	 
	 When parsing parameters file there are a few differences:
	 
	 These need to be new parameters in the parameters file:
	 BURN_IN
	 SAMPLE 
	 */
	
	
	
	//(1) Open Parameters File
	string parameters_file = argv[1];
	ifstream infos;
	infos.open(parameters_file.c_str());
	
	//(2) Determine input path; where all the svs are
	string input_path = argv[2];
		
    char lastChar = input_path.at( input_path.length() - 1 );
    if(lastChar!= '/'){
        input_path = input_path + "/";
    }

    
	//(3) Get Probability Model Parameters:
	double COVERAGE;
	double COVERAGE_SCALED;
	int LAVG;
	int LDIS;
	int READLEN;
	
	double Tolerance;
	int PRINT_FLAG;
	double LRTHRESHOLD = 0;
	bool PRINTALL = false;
		
	//(4) Get the MCMC Iteration Steps;
	int DEFAULT_BURN_IN = 100000;
	int DEFAULT_SAMPLE  = 900000;
	int BURN_IN = DEFAULT_BURN_IN;
	int SAMPLE = DEFAULT_SAMPLE;
			
	//(5) Get Error for the probability model
	double p_err = 0.01; //Default error parameter
	
	//Step 0: Process the Parameters File
	cout << "--Input Parameters--" << endl;
	string temp;
	ifstream p_file(argv[1], ios::in);
	if(p_file.is_open()){cout << "   Parameters file found." << endl;}
	else{cout << "WARNING: Parameters file \"" << argv[1] << "\" cannot be opened." << endl; exit(1);}
	while(getline(p_file, temp))
	{
		if(temp[0] == '#')
			continue;
		else
		{
			string term;
			string value;
			int spacePos = temp.find(' ');
			if(spacePos == -1)
			{
				spacePos = temp.find(':');
				value = temp.substr(spacePos+1);
				term = temp.substr(0,spacePos+1);
				
			}
			else{
				value = temp.substr(spacePos+1);
				term = temp.substr(0,spacePos);
			}
			
			if(term == "Lavg:")
			{
				LAVG = atoi(value.c_str());
				cout << "   Lavg: " << LAVG << endl;
				continue;
			}
			if(term == "ReadLen:"){
				READLEN = atoi(value.c_str());
				cout << "   Read Length: " << READLEN << endl;
				continue;
			}
			if(term == "Lambda:"){
				COVERAGE = atof(value.c_str());
				cout << "   Lambda: " << COVERAGE << endl;
				continue;
			}
			if(term == "Perr:"){
				p_err = atof(value.c_str());
				cout << "   Perr: " << p_err << endl;
				continue;
			}
			if(term == "Tolerance:"){
				Tolerance = atof(value.c_str());
				cout << "   Tolerance: " << Tolerance << endl;
				continue;
			}

			if(term == "Burnin:"){
				BURN_IN = atoi(value.c_str());
				cout << "   Burnin: " << BURN_IN << endl;
				continue;
			}
			
			if(term == "Sample:"){
				SAMPLE = atoi(value.c_str());
				cout << "   Sample: " << SAMPLE << endl;
				continue;
			}
			
			if(term == "Verbose:"){
				if(value == "Y" || value == "y" || value == "yes" || value == "YES" || value == "Yes"){
					PRINT_FLAG = 1;
				cout << "   ***Verbose Mode Enabled***" << endl;}
				else{
					PRINT_FLAG = 0;
				}
				continue;
			}	
			if(term == "LRThreshold:"){
				if(value == "all" || value == "All" || value == "ALL"){
					PRINTALL = true;
					cout << "   LR Threshold: PRINT_ALL" << endl;
				}
				else{
					LRTHRESHOLD = atof(value.c_str());
					PRINTALL = false;
					cout << "   LR Threshold: " << LRTHRESHOLD << endl;
				}
				continue;
			}
		}
	}
	
	//Error Checking:
	if(BURN_IN <=0){ 
		BURN_IN = DEFAULT_BURN_IN;
		cerr << "Error: Parameter Burnin must be positive. Resetting Burnin to default of " << DEFAULT_BURN_IN << endl;
	}
	
	if(SAMPLE <=0){ 
		SAMPLE = DEFAULT_SAMPLE;
		cerr << "Error: Parameter Sample must be positive. Resetting Sample to default of " << DEFAULT_SAMPLE << endl;
	}
	
	//Set coverage to the previous value:
	COVERAGE = COVERAGE/LAVG;
	
	//Set remaining parameters;
	LDIS = LAVG - 2*READLEN;
	COVERAGE_SCALED = COVERAGE * (LDIS*1.0/LAVG*1.0);
	
	// Step 1: Read the info from p_star.summary (MAY NOT BE NEEDED ANYMORE!)
	ifstream sum_info;
	fn = input_path + "/p_star.summary"; // file name for summary info
	sum_info.open(fn.c_str());
	if(sum_info.fail()){
		cerr<<"Missing file " << fn << ".\n";
        cerr<<"Indicates GASVPro-graph may not have run correctly.\n";
		return -1;
	}
	sum_info>> total_num_of_sv_clusters >> total_num_of_fragments >> dist_num_of_fragments;
	sum_info.close();
	cout << "Successfully Read in the Summary Info:\t" << total_num_of_sv_clusters << "\t" << total_num_of_fragments << "\t" << dist_num_of_fragments << endl;
	
	//Step 2: Read in/process each connected component:
		
	//Set up random number generator 
	srand((unsigned)time(0));
	
	ifstream subset_inp, subset_sv;
	
	int startComponent = 0;
	int endComponent =total_num_of_sv_clusters;
	if(START >= 0 && END >= 0){ startComponent = START; endComponent = END; }
	
	cout << "Processing SV Components " << startComponent << " to " << endComponent << ".\n";
		
	string output_path = input_path +"/" + parameters_file;
	
	for(int subset = startComponent; subset <= endComponent ; subset++){
		
		//Hack for making an integer a string;
		std::stringstream subsetString; 
		subsetString << subset;
		fn = input_path + "sv_" + subsetString.str() + ".sv";

		//Output level -- Default: only output the MCMCThreshold.clusters
		//Output level -- Full:    output the espFile and varFile.
		
		string espFile = output_path + "_sv_" + subsetString.str() + ".PR.results";
		string varFile = output_path + "_sv_" + subsetString.str() + ".Variant.results";
		//string likeFile = output_path + "_sv_" + subsetString.str() + ".Likelihood.results";
		//string mleFile = output_path + "_sv_" + subsetString.str() + ".MCMCMLE.clusters";
		//string averageFile = output_path + "_sv_" + subsetString.str() + ".MCMCAverage.clusters";
		//string avgVarOnly =  output_path + "_sv_" + subsetString.str() + ".MCMCAverageVar.clusters";
		string thresholdFile = output_path + "_sv_" + subsetString.str() + ".MCMCThreshold.clusters";
		
		subset_inp.open(fn.c_str());

		cout << "Trying to process " << fn.c_str() << endl;
		
		if(!subset_inp.fail()){
			cout << "Component " << subset << " does exist, processing.\n";
			
			//Make output files:
			
			//(1) ESP Results:
			ofstream outESP(espFile.c_str(),ios::out);
			
			//(2) Variant Results:
			ofstream outVAR(varFile.c_str(),ios::out);
			
			//(3) Likelihood Results:
			//ofstream outLIKE(likeFile.c_str(),ios::out);
			
			//(4) MCMC Results;
			//ofstream outMLE(mleFile.c_str(),ios::out);
			
			//(5) Thresholding Results;
			ofstream outTHRESHOLD(thresholdFile.c_str(),ios::out);
			
			//ofstream outAVERAGE(averageFile.c_str(),ios::out);
			
			//ofstream outAVERAGEVAR(avgVarOnly.c_str(),ios::out);
								   
			//(4) Assignment Matrices: 
			//Step 1: Read in the file once to count the number of ESPs and Variants needed
			
			cout << "\tFirst read, counting ESPs and Variants.\n";
			
			fn = input_path + "sv_" + subsetString.str() + ".sv";
			subset_sv.open(fn.c_str());
			int sv_count_subset; //This is not many SVs we need to process in this file.
			subset_sv>>sv_count_subset;
			string temp; getline(subset_sv,temp);
			
			string clusterLine;
			string coverageLine;
			vector <int> observedESPs;
			//Note: Did not want a vector with 0's in it! cout << "ObservedESPs size:\t " << observedESPs.size() << endl << flush;
			
			for(int i = 0; i<sv_count_subset; i++){
				getline(subset_sv,clusterLine);
				getline(subset_sv,coverageLine);
			
				vector <string> clusterTokens;
				vector <string> coverageTokens;				
				Tokenize(clusterLine,clusterTokens,"\t,"); //<--Needed to also split on fragments and coordinates.
				Tokenize(coverageLine,coverageTokens);
			
				//ClusterLine: (Note: The ESPs and coordinates are each on their own line!)
				//0     1    2  3       4 -> 4+NumESPS
				//c1	1	189	D	1_1_1_chr17_41883_42105_0:0:0_2:0:0_15789d3 	17	17	41642, 42623, 41821, 42623, 41532, 42334, 41532, 42513
				int numESPs = atoi(clusterTokens[1].c_str());
				for(int j = 4; j<4+numESPs; j++){
					vector <string> ESPTokens;
					Tokenize(clusterTokens[j],ESPTokens,"_");
					int a = atoi(ESPTokens[0].c_str());
					//Do we even need b and c?
					//int b = atoi(ESPTokens[1].c_str());
					//int c = atoi(ESPTokens[2].c_str());
					observedESPs.push_back(a);
				}
			}

			//Set up the unique ESPs; (Note: We just read in order to not have to worry about finding in succession)
			sort( observedESPs.begin(), observedESPs.end() );
			observedESPs.erase( unique( observedESPs.begin(), observedESPs.end() ), observedESPs.end() );
						
			std::vector<int>::iterator it;
			//int i = 0;
			//cout << "Here are the observed ESPs.\n";
			//for(it = observedESPs.begin(); it != observedESPs.end(); it++){    cout << i << " : " <<  *it << "\n"; i++;} 

			//Get the size of the vector;
			unsigned numObservedESPs = observedESPs.size();
			unsigned numObservedVariants = sv_count_subset;
			
			cout << "\tFinished, Component " << subset << " has " << numObservedESPs << " ESPs and " << numObservedVariants << " variants.\n";
			subset_sv.close();
			
			//Second Read Through; get the variants
			subset_sv.open(fn.c_str());
			
			cout << "\tSecond read, determining possible assignment matrices.\n";

			vector <int> numberOfMappings(observedESPs.size(),0); //When do we use this?
			vector <ESP> componentESPs(numObservedESPs);
			vector <variant> componentVariants(numObservedVariants);
			vector <pairedDependency> componentPairs(1);
			
			//cout << "Beginning second read through.\n";
			
			//Read the header line;
			subset_sv >>sv_count_subset;
			getline(subset_sv,temp);
			
			//New Comment //vector<int> ESPsToVary;
			//New Comment //int numESPsToVary = 0;
			
			for(int i = 0; i<sv_count_subset; i++){
				vector <int> seenESPs;
				int numSeen = 0;
				
				getline(subset_sv,clusterLine);
				getline(subset_sv,coverageLine);
				
				vector <string> clusterTokens;
				Tokenize(clusterLine,clusterTokens,"\t,"); //<--Needed to also split on fragments and coordinates.
				
				//ClusterLine: (Note: The ESPs and coordinates are each on their own line!)
				//0     1    2  3       4 -> 4+NumESPS								5	6	7
				//c1	1	189	D	1_1_1_chr17_41883_42105_0:0:0_2:0:0_15789d3 	17	17	41642, 42623, 41821, 42623, 41532, 42334, 41532, 42513
				componentVariants[i].setName(clusterTokens[0].c_str());
                componentVariants[i].setType(clusterTokens[3].c_str());
				
				//cout << "Processing variant " << i << ":\t" << componentVariants[i].getName() << endl << flush;
				
				//Set up the rest of a variant, this is what is to print out.
				vector <string> tmpTokens;
				int retVal = Tokenize(clusterLine,tmpTokens,"\t");
				
				//cout << "\tTokenized the row --> " << clusterLine << endl << flush;
				//cout << "\tWe have " << retVal << " number of tokens.\n";
				
				componentVariants[i].setTheRest(tmpTokens[2] + "\t" + tmpTokens[3] + "\t" + tmpTokens[4] + "\t" + tmpTokens[5] + "\t" + tmpTokens[6] + "\t" + tmpTokens[7]); 
				
				//cout << "\tSet the rest " << clusterLine << endl << flush;
				
				int numESPs = atoi(clusterTokens[1].c_str()); 
								
				for(int j = 4; j<4+numESPs; j++){
					vector <string> ESPTokens;
					Tokenize(clusterTokens[j],ESPTokens,"_");
					int a = atoi(ESPTokens[0].c_str());
					// Do we even need b and c?
					//int b = atoi(ESPTokens[1].c_str());
					//int c = atoi(ESPTokens[2].c_str());
					//a = Variant; b = currentNumberOfMapping; c = GlobalCounterOfMapping
				
					//Now tokenize to see if we are low-quality/repetitive;
					//Note: The low-quality and ambiguous mappings are the only ones we 
					//      want to consider for ambiguity. The other ESPs will be added to the permanent
					/*
					vector <string> ESPSplit;
					int numTokens = Tokenize(clusterTokens[j],ESPSplit,".");
					if(numTokens>1){
						vector <string> lowQualityValues;
						int numSplit = Tokenize(ESPSplit[numTokens-1],lowQualityValues,"_");
						if(numSplit == 3){
							//cout << "Low_Quality_Mapping:\t" << clusterTokens[j] << endl;
							int local = findLocation(ESPsToVary,numESPsToVary,a);
							if(local<0){ 
								int trueLocal = -1*local-1;
								ESPsToVary.insert(ESPsToVary.begin() + trueLocal,a);
								numESPsToVary++;
							}
						}
					}	
					*/
					int ESPIndex = findLocation(observedESPs,numObservedESPs, a);
					int alreadySeen = findLocation(seenESPs,numSeen,a);
					
					//cout << "Seeing esp " << a << " with index " << ESPIndex <<  " in variant " << clusterTokens[0].c_str() << endl;
				
					if(ESPIndex>=0 && alreadySeen < 0){ 
						seenESPs.push_back(a); //Mark as seen
						numSeen++;
						sort( seenESPs.begin(), seenESPs.end() ); //Resort the seen vector;
						numberOfMappings[ESPIndex]++; //Update Number of Views;
						componentVariants[i].addESP(ESPIndex); 	//Add ESP to the variant;
						componentESPs[ESPIndex].addVariant(i); //Add variant to the ESP:
					}
					else if(ESPIndex>=0){ } //We found it, but have already added it.
					else{ 
						cerr << clusterLine<< endl;
						cerr << "Could not find " << a << " in the ESP list.\n";
						exit(-1); 
					}					
				
				}
				
				//cout << "Finished variant (before sorting) " << clusterTokens[0].c_str() << componentVariants[i] << endl;

				
				//Sort the ESPs that have been recorded;
				componentVariants[i].sortPossibleESPs();
				
				//cout << "Finished variant " << clusterTokens[0].c_str() << componentVariants[i] << endl;
				
				//Note: The coverage line is NOW different!
				//0     1   2       3           4               5         6         7          8    9   10 11  
				//c1	D	1	-234.852	    -947.093	   41532	42334	    0	       0	0	0	0

				//c28	D	1	-1.72921e+08	-1.48652e+10	40688	39764630	37524588	0	-1	0	1

				
				vector <string> coverageTokens;	
				int numCoverageTokens = Tokenize(coverageLine,coverageTokens,"\t");
				componentVariants[i].setValues(atoi(coverageTokens[5].c_str()), atoi(coverageTokens[6].c_str()), atoi(coverageTokens[7].c_str()),
											   atoi(coverageTokens[8].c_str()), atof(coverageTokens[9].c_str()), atof(coverageTokens[10].c_str()),
											   atoi(coverageTokens[11].c_str()));
				cout << "\tFinished cluster --> " << coverageTokens[0].c_str() << endl;
				
			}
			subset_sv.close();
						
			//Need to set up the MCMC Moves?
			for(unsigned int i = 0; i<numObservedESPs; i++){
				//cout << "\tProcessing ESP " << i << "/" << numObservedESPs << " with " << componentESPs[i].getPossibleMappings() << " mappings\n" << flush;
				
				int numMappings = componentESPs[i].getPossibleMappings();
				componentESPs[i].setUpForMCMC();
				componentESPs[i].setName(observedESPs[i]);	
			}
			
			
			//Ready for the MCMC Moves; 
			vector<int> assignedVariants;
			vector<int> emptyVariants;

			int numAssigned = 0;
			int numEmpty = 0;
			//cout << "Before Random Initialization " << componentVariants[5] << endl;
			
			//NEW: Fix the unique EPSs if needed;
			int numMobileESPs = 0;
			vector <int> mobileESPs;
			/* Choice 1: Move all of them! */
			//for(int i = 0; i<numObservedESPs; i++){ mobileESPs.push_back(i); numMobileESPs++; }
			
			/* Choice 2: Only multiple mappings witll be moved */
			//Note: The mobile ESPs is a global variable! 
			for(int i = 0; i<numObservedESPs; i++){
				//This ESP is not movable; keep going;
				if(componentESPs[i].getPossibleMappings() == 1){ }
				else if(DO_UNIQUE == 1 ){
					//Do not allow this to be moved!
				}
				else{
					mobileESPs.push_back(i);
					numMobileESPs++;
				}
			}
			
			cout << "We are about to begin and have " << numMobileESPs  << "/" << numObservedESPs << " mobile ESPs.\n";
			
			randomInitialization(numMobileESPs,mobileESPs,componentESPs,componentVariants,assignedVariants,numAssigned,emptyVariants,numEmpty,p_err);
	
			
			cout << "After Random Initialization " << endl;
			
			//(0) Stay the same; (1) Naive; (2) Add; (3) Remove; (4) Swap
			//Needed Code
			vector<int> modifiedESPs; 
			modifiedESPs.clear();
			int numModifiedESPs = 0;
			vector<int> modifiedVariants; 
			modifiedVariants.clear();
			int numModifiedVariants = 0;
		
			int NUM_MOVES = 5;
			vector<int> moveTally(5,0);
			vector<int> moveAcceptanceTally(5,0);

			int total = SAMPLE+BURN_IN;
		
			double MAX_LIKELIHOOD = 1;
			stringstream MAX_LIKELIHOOD_STREAM;
			string MAX_LIKELIHOOD_STRING;
			MAX_LIKELIHOOD_STREAM.str("");
			MAX_LIKELIHOOD_STRING = MAX_LIKELIHOOD_STREAM.str();
			
			//Set the MLE Support Values;
			vector<int>MLE_SUPPORT;
			for(int i = 0; i<numObservedVariants;i++){ MLE_SUPPORT.push_back(0); }
			
			//outLIKE << "STEPS\tLOG_LIKELIHOOD\tNUM_ESPS\tESP_ASSIGNMENTS\n";
			
	//Note: No moves are possible; only unique entries so we will update all ESP mappings.
	int movesPossible = 0;
	if(numMobileESPs == 0){ movesPossible = 1; }		
	
	if(movesPossible == 1){
		int stepsToIncrement = SAMPLE-1;
		
		for(unsigned int j = 0; j<numObservedESPs; j++){ componentESPs[j].clearOccupancy(); }
		for(unsigned int j = 0; j<numObservedVariants; j++){ componentVariants[j].clearOccupancy(); }
		MAX_LIKELIHOOD = 1;
		MAX_LIKELIHOOD_STREAM << ""; //Need to empty a stream.
		MAX_LIKELIHOOD_STRING = MAX_LIKELIHOOD_STREAM.str();
	
		
		//Update the ESPs
		for(int i = 0; i<numObservedESPs; i++){
			componentESPs[i].updateOccupancy(stepsToIncrement);
		}
		
		//Update the Variants;
		for(int i = 0; i<numObservedVariants;i++){
			componentVariants[i].updateOccupancy(stepsToIncrement,p_err,COVERAGE,COVERAGE_SCALED,LAVG,LDIS); 
		}
		
		
		//Get the current likelihood:
		double currentLikelihood = -1;
		int moveToMake = 0;
		int accept = acceptMove(numMobileESPs,mobileESPs,moveToMake,currentLikelihood,componentESPs,componentVariants,assignedVariants,numAssigned,emptyVariants,numEmpty,modifiedESPs,numModifiedESPs,modifiedVariants,numModifiedVariants,p_err,COVERAGE,COVERAGE_SCALED,LAVG,LDIS);
		
		//outLIKE << 0 << "\t" << currentLikelihood << "\t" << componentESPs.size() << "\t";
		
		int MAX_FLAG = 0;
		if(currentLikelihood > MAX_LIKELIHOOD || MAX_LIKELIHOOD > 0){ 
			MAX_FLAG = 1;
			MAX_LIKELIHOOD = currentLikelihood;
			MAX_LIKELIHOOD_STREAM.str("");
			MAX_LIKELIHOOD_STRING = MAX_LIKELIHOOD_STREAM.str();
			for(int i = 0; i<numObservedVariants;i++){ MLE_SUPPORT[i] = 0; }
		}
		
		//if(MAX_FLAG == 1){ MAX_LIKELIHOOD_STREAM << steps << "\t" << currentLikelihood << "\t" << componentESPs.size() << "\t"; }
		
		for(unsigned int j = 0; j<componentESPs.size(); j++){
			//outLIKE << componentESPs[j].getVariant(componentESPs[j].getCurrentState()) << "\t";
			if(MAX_FLAG == 1){ 
				int ESPVariant = componentESPs[j].getVariant(componentESPs[j].getCurrentState());
				if(ESPVariant>=0){ MLE_SUPPORT[ESPVariant]++; }
				//We can print out whatever we want, I'm not sure why we would print it out though.
				//MAX_LIKELIHOOD_STREAM << componentESPs[j].getVariant(componentESPs[j].getCurrentState()) << "\t";
			}
		}
		//outLIKE << endl;
		if(MAX_FLAG == 1){ 
			for(int i = 0; i<numObservedVariants; i++){
				if(MLE_SUPPORT[i]>0){ 
					MAX_LIKELIHOOD_STREAM << componentVariants[i].getName() 
					<< "_" << componentVariants[i].getLikelihoodVariantGiven(p_err,COVERAGE,COVERAGE_SCALED,LAVG,LDIS,MLE_SUPPORT[i]) 
					<< "_" << componentVariants[i].getLikelihoodErrorGiven(p_err,COVERAGE,COVERAGE_SCALED,LAVG,LDIS,MLE_SUPPORT[i])
					<< "\t" << MLE_SUPPORT[i] << "\t" << componentVariants[i].getTheRest() << endl;
				}
			}
			
			MAX_LIKELIHOOD_STRING = MAX_LIKELIHOOD_STREAM.str(); 
			
		}
		
		
	
	}
	else{
	//	cout << "Starting the MCMC now\n" << flush;
			for(int steps  = 1; steps<total; steps++){
				if(steps%100000 == 0){ cout << "On Step " << steps << " of " << total << endl; }
				
				if(steps == BURN_IN){
					cout << "Clearing the occupancies after the burnin.\n";
					for(unsigned int j = 0; j<numObservedESPs; j++){ componentESPs[j].clearOccupancy(); }
					for(unsigned int j = 0; j<numObservedVariants; j++){ componentVariants[j].clearOccupancy(); }
					MAX_LIKELIHOOD = 1;
					MAX_LIKELIHOOD_STREAM << ""; //Need to empty a stream.
					MAX_LIKELIHOOD_STRING = MAX_LIKELIHOOD_STREAM.str();
				}
				
				double randomMoveType = ((double) rand()*1.0 / (RAND_MAX+1.0));
				int moveToMake = -1;

				if(randomMoveType < 0.5 || numMobileESPs == 0){
					moveTally[0]++;
					moveToMake = 0;
				}
				else{
					//moveToMake = 1;
							
					 moveToMake = -1;
					 while(moveToMake == -1){
						double whichMove = ((double) rand()*1.0 / (RAND_MAX+1.0));;
						 double numMoves = 3.0;
						 if(whichMove <= 1.0/numMoves){ moveToMake = 1;}
						else if(whichMove<= 2.0/numMoves){ moveToMake = 2;}
						else if(whichMove<=3.0/numMoves){ moveToMake = 3;}
						else{ moveToMake = 4;}
						
						if(moveToMake == 2 && numEmpty == 0){ moveToMake = -1; } //Cannot add a variant now;
						 if(moveToMake == 3 && numAssigned == 0){ moveToMake = -1; } //Cannot remove a variant now;
						//if(moveToMake == 4 && swapPossible == 0){ moveToMake = 1; } //Cannot swap a variant now;
						 if(moveToMake == 4){ moveToMake = - 1; } //Cannot swap a variant now;
					 }
				
				
					moveTally[moveToMake]++;
					
					if(moveToMake == 1){
						//cout << "MOVE " << steps << "--> NAIVE\n";
						//for(int i = 0; i<5; i++){ cout << i << "--> " << moveTally[i] << "; "; } cout << endl;
						proposeNaiveMove(numMobileESPs,mobileESPs,componentESPs,componentVariants,assignedVariants,numAssigned,emptyVariants,numEmpty,modifiedESPs,numModifiedESPs,modifiedVariants,numModifiedVariants,p_err);
					}
					else if(moveToMake == 2){
						//cout << "MOVE " << steps << "--> ADD\n";
						//for(int i = 0; i<5; i++){ cout << i << "--> " << moveTally[i] << "; "; } cout << endl;
						proposeAddMove(numMobileESPs,mobileESPs,componentESPs,componentVariants,assignedVariants,numAssigned,emptyVariants,numEmpty,modifiedESPs,numModifiedESPs,modifiedVariants,numModifiedVariants,p_err);
					}
					else if(moveToMake == 3){
						//cout << "MOVE " << steps << "--> REMOVE\n";
						//for(int i = 0; i<5; i++){ cout << i << "--> " << moveTally[i] << "; "; } cout << endl;
						proposeRemoveMove(numMobileESPs,mobileESPs,componentESPs,componentVariants,assignedVariants,numAssigned,emptyVariants,numEmpty,modifiedESPs,numModifiedESPs,modifiedVariants,numModifiedVariants,p_err);
					}
					else if(moveToMake == 4){
						//cout << "MOVE " << steps << "--> SWAP\n";
						//for(int i = 0; i<5; i++){ cout << i << "--> " << moveTally[i] << "; "; } cout << endl;
						//proposeSwapMove(numMobileESPs,mobileESPs,componentESPs, componentVariants, componentPairs,runningDependencies,TotalDependencies,visitedDependencies,numDependencies, steps,assignedVariants,numAssigned,emptyVariants,numEmpty,modifiedESPs,numModifiedESPs,modifiedVariants,numModifiedVariants,p_err);
					}
					else{
						cout << "Problem in random number generation?\n"; 
					}
				}
				
				double currentLikelihood = -1;
				int accept = acceptMove(numMobileESPs,mobileESPs,moveToMake,currentLikelihood,componentESPs,componentVariants,assignedVariants,numAssigned,emptyVariants,numEmpty,modifiedESPs,numModifiedESPs,modifiedVariants,numModifiedVariants,p_err,COVERAGE,COVERAGE_SCALED,LAVG,LDIS);
				
				/*cout << "Accepted --> " << accept << endl;
				for(int i = 0; i<numAssigned; i++){ 
				 cout << "Variant " << assignedVariants[i] << "\t";
				 
				 cout << componentVariants[assignedVariants[i]] << endl;
				 }
				cout << "******\n";
				*/
				
				//outLIKE << steps << "\t" << currentLikelihood << "\t" << componentESPs.size() << "\t";
				
				int MAX_FLAG = 0;
				if(currentLikelihood > MAX_LIKELIHOOD || MAX_LIKELIHOOD > 0){ 
					MAX_FLAG = 1;
					MAX_LIKELIHOOD = currentLikelihood;
					MAX_LIKELIHOOD_STREAM.str("");
					MAX_LIKELIHOOD_STRING = MAX_LIKELIHOOD_STREAM.str();
					for(int i = 0; i<numObservedVariants;i++){ MLE_SUPPORT[i] = 0; }
				}
				
				//if(MAX_FLAG == 1){ MAX_LIKELIHOOD_STREAM << steps << "\t" << currentLikelihood << "\t" << componentESPs.size() << "\t"; }
		
				for(unsigned int j = 0; j<componentESPs.size(); j++){
					//outLIKE << componentESPs[j].getVariant(componentESPs[j].getCurrentState()) << "\t";
					if(MAX_FLAG == 1){ 
						int ESPVariant = componentESPs[j].getVariant(componentESPs[j].getCurrentState());
						if(ESPVariant>=0){ MLE_SUPPORT[ESPVariant]++;}
						//MAX_LIKELIHOOD_STREAM << componentESPs[j].getVariant(componentESPs[j].getCurrentState()) << "\t";
					}
				}
				//outLIKE << endl;
				if(MAX_FLAG == 1){ 
					for(int i = 0; i<numObservedVariants; i++){
						if(MLE_SUPPORT[i]>0){ 
								MAX_LIKELIHOOD_STREAM << componentVariants[i].getName() 
										<< "_" << componentVariants[i].getLikelihoodVariantGiven(p_err,COVERAGE,COVERAGE_SCALED,LAVG,LDIS,MLE_SUPPORT[i]) 
										<< "_" << componentVariants[i].getLikelihoodErrorGiven(p_err,COVERAGE,COVERAGE_SCALED,LAVG,LDIS,MLE_SUPPORT[i])
										<< "\t" << MLE_SUPPORT[i] << "\t" << componentVariants[i].getTheRest() << endl;
						}
					}
					
					MAX_LIKELIHOOD_STRING = MAX_LIKELIHOOD_STREAM.str(); 
				
				}
				
				moveAcceptanceTally[moveToMake] += accept;
				//cout << "Finished accepting the move.\n";
				
				//If we have made a move, figure out if a swap is possible.
			//	if(moveToMake > 0){ 
			//		 swapPossible = isSwapPossible(assignedVariants,numAssigned,componentVariants);
			//	}
				
				//cout << "\tTotalESPs: " << componentESPs.size() << endl;
				int RUNNING = 0;
				int TOTAL_CHECK = componentESPs.size();
				int numErrors = 0;
				for(unsigned int j = 0; j<componentESPs.size(); j++){ if(componentESPs[j].getVariant(componentESPs[j].getCurrentState()) == -1) { numErrors++; }}
				RUNNING += numErrors;
				//cout << "\tErrors: " << numErrors << "\n";
				for(unsigned int j = 0; j<componentVariants.size(); j++){ 
					//cout << "\tVariant_" <<  j << ": " << componentVariants[j].getCurrentAssigned() << endl; 
					RUNNING+= componentVariants[j].getCurrentAssigned();
				}
											
				if(RUNNING!=TOTAL_CHECK){ cerr << "Error in MCMC after move: We are missing assignments: Running: " << RUNNING << " vs Total: " << TOTAL_CHECK << "\n"; exit(-1); }
				//else{ cout << "Move " << steps << " ran correctly; Running: " << RUNNING << " vs Total: " << TOTAL_CHECK << "\n"; }
 
			}
			
			
		} //END if there were moves required;
			
			cout << "Finished; outputting results of MCMC.\n";
					
			//OUTPUT ESP INFORMATION
			outESP << "ESP_INDEX\tESP_NAME\tERROR_FRAC\tNUM_MAPPINGS\tVARIANT\tOCCUPANCIES\n";
			for(unsigned int i = 0; i<numObservedESPs; i++){
				outESP << i << "\t" << componentESPs[i].getName() << "\t";
				outESP << componentESPs[i].getOccupancy(-1,SAMPLE)*1.0/SAMPLE*1.0 << "\t";
				outESP << componentESPs[i].getPossibleMappings() << "\t";
				for(int j = 0; j<componentESPs[i].getPossibleMappings(); j++){
					int var = componentESPs[i].getVariant(j);
					outESP << componentVariants[var].getName() << " ";
					outESP << componentESPs[i].getOccupancy(j,SAMPLE)*1.0/SAMPLE*1.0 << "\t";
				}
				outESP << endl;
			}
			
			//OUTPUT VARIANT INFORMATION
			outVAR << "VAR_INDEX\tVAR_NAME\tNUM_POSSIBLE_MAPPINGS\tFRACTIONAL_SUPPORT\n";
			for(unsigned int i = 0; i<numObservedVariants;i++){
				outVAR << i << "\t" << componentVariants[i].getName() << "\t"
					   << componentVariants[i].getPossibleAssigned() << "\t";
				for(int j = 0; j<=componentVariants[i].getPossibleAssigned(); j++){
						outVAR << componentVariants[i].getSupport(j)*1.0/SAMPLE << "\t";
				}
				outVAR << endl;
			}
			
			//Average Occupancy
			for(unsigned int i = 0; i<numObservedVariants;i++){
				//Were we ever occupied? If yes, then output!
				if(componentVariants[i].getTimeOccupied() > 0){
					//outAVERAGE << componentVariants[i].getName() << "_" 
					//			 << componentVariants[i].getRunningLogLikelihood()*1.0/(SAMPLE*1.0) << "_0"
					//			 << "\t" << componentVariants[i].getPossibleAssigned() 
					//			 << "\t" << componentVariants[i].getTheRest() << endl;
					
					//outAVERAGEVAR << componentVariants[i].getName() << "_" 
					//	<< componentVariants[i].getRunningLikelihoodVar()*1.0/(SAMPLE*1.0) << "_0"
					//	<< "\t" << componentVariants[i].getPossibleAssigned() 
					//	<< "\t" << componentVariants[i].getTheRest() << endl;
					
				}
			}
			
			//OUTPUT MAX_LIKELIHOOD_RESULT;
			//cout << "Max-Likelihood-Result" << endl;
			//outMLE << MAX_LIKELIHOOD_STRING;
			
			//BEGIN: Thresholding the results:

			//Step 3: Final Assignment of All ESPs to a cluster:
			//(0) Clear all current assignments of variants;
			for(int i = 0; i<numObservedVariants; i++){ componentVariants[i].clearCurrentState(); }
			
			//cout << "Successfully cleared each variant\n" << flush;
			
			//(1) Determine the TRUE final state for each ESP; Require that we stay in the state at least 80%.
			double threshold = 0.80;
			for(int i = 0; i<numObservedESPs; i++){				
				int variantToUpdate = componentESPs[i].setOptimalPosition(threshold*SAMPLE); 
				if(variantToUpdate>=0){ 
					//cout << "Adding ESP " << componentESPs[i].getName() << " to variant " << componentVariants[variantToUpdate].getName() << endl;
					componentVariants[variantToUpdate].addESPMapping(i);
				}
			}
			
			//(2) For each variant; output the cluster along with the variant likelihood:
			for(int i = 0; i<numObservedVariants; i++){ 
				int current = componentVariants[i].getCurrentAssigned();
				if(current > 0){
					double LLR = componentVariants[i].getLikelihoodVariantGiven(p_err,COVERAGE,COVERAGE_SCALED,LAVG,LDIS,current) - componentVariants[i].getLikelihoodErrorGiven(p_err,COVERAGE,COVERAGE_SCALED,LAVG,LDIS,current);
					if(LLR >= LRTHRESHOLD || PRINTALL ){
						outTHRESHOLD << componentVariants[i].getName() 
							<< "\t" << componentVariants[i].getCurrentAssigned() 
							<< "\t" << componentVariants[i].getTheRest() 
							<< "\t" << LLR << endl;
                        
                        /*
                        cout << "Name: " << componentVariants[i].getName() << endl;
                        cout << "Current Assigned: " << componentVariants[i].getCurrentAssigned() << endl;
                        cout << "The Rest: " << componentVariants[i].getTheRest() << endl;
                        cout << "LLR: " << LLR << endl;
                        */
					}
				}
			}

			outESP.close();
			outVAR.close();		
			outTHRESHOLD.close();
		
			subset_inp.close();
		}
	
	
	}
		
	return 0;
}
