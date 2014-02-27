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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <strings.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include "IntervalTree.h"  //INTERVALTREE Courtesy of Erik Garrison <https://github.com/ekg/intervaltree>
#include <vector>		   
#include <algorithm>

using namespace std;

struct uInt
{
	int start;
	int end;
	int length;
	
	int bases;
	float uniqueness;
};

struct CHelper
{
	int coverage;
	bool divergent;
};

struct point{
	int x;
	int y;
};	


int pointInPoly(int *STARTS, int *ENDS, int npoints, unsigned int xt, unsigned int yt);
int findLocation(int *STARTS, int NUM_FRAGMENTS, int START); 
double getUniqueness(int* START_UNIQUE,int* END_UNIQUE,int *VALUE_UNIQUE,int NUM_UNIQUE,int startpoint,int endpoint);
int getCoverage(int *STARTS, int *ENDS,int NUM_SEGMENTS, int startpoint, int endpoint, int BUFFER);
int getCoverage(int *STARTS, int *ENDS,int NUM_SEGMENTS, int point, int BUFFER);
int getCoverage(int start, int end, IntervalTree<CHelper*> *tree);
int getCoverage(int point, vector<pair<int,int> >* genome);
bool intervalCompare(pair<int,int> intervalA, pair<int,int> intervalB);
bool genomeCompare(pair<int,int> pointA, pair<int,int> pointB);
double getUniqueness(int start, int end, IntervalTree<uInt*>* tree, double MAX_UNIQUE_VAL);
double getUniquenessNEW(int start, int end, IntervalTree<uInt*> *tree, double MAX_UNIQUE_VAL);



//Goal:
// Compute the CDF (Prob(<=coverage)), reject cases that lie too far out in the tail.
int validCoverage(double mean, double coverage, double tolerance){
	int COV = (int) floor(coverage);
	int FACTOR;
	//Use normal approximation; Poisson is essentially normal w/ mean Lambda std Lambda;
	if(mean >= 100){
		if(tolerance>=0.05){ FACTOR = 2; } 
		else if(tolerance>=0.002){ FACTOR = 3;}
		else if(tolerance>=0.00006){ FACTOR = 4;}
		else if(tolerance>=0.0000005){ FACTOR = 5;}
		else{ FACTOR = 6; }
		if(COV >= FACTOR*mean){
			return 0; 
		}
		else{
			return 1;
		}
	}
	double rolling = 0;
	double VAL = exp(-mean);
	int keepGoing = 1;
	for(int i = 0; i<=COV && keepGoing == 1;i++){
		rolling += exp(i*log(mean) - lgamma(i+1));
		double test = rolling*VAL;
		if( (1.0 - test) < tolerance){ keepGoing = 0; }
	}
	if(keepGoing == 0){
		return 0;
	}
	else{ 
		return 1; 
	}
}

//Determines if possible breakpoint pair is on the boundary
int onBoundary(int* X_COORDS,int* Y_COORDS, int NUM_VALUES, int X, int Y){
	int match = 0;
	//Are we one of the corer points?
	for(int i = 0; i<NUM_VALUES;i++){
		if(X_COORDS[i] == X && Y_COORDS[i] == Y){
			match++;
		}
	}
	if(match>0){
		return match;
	}
	
	//Are we along one of the line segments?
	point c;
	c.x = X; c.y = Y;
	for(int i = 0; i<NUM_VALUES;i++){
		//Set up points on the line segment:
		point a, b;
		if(i == 0){ a.x = X_COORDS[NUM_VALUES-1]; a.y = Y_COORDS[NUM_VALUES-1]; }
		else{ a.x = X_COORDS[i-1]; a.y = Y_COORDS[i-1]; }
		b.x = X_COORDS[i]; b.y = Y_COORDS[i];
		
		int doWeMatch = 0;
		
		//Determine Relevant Values
		int crossproduct = (c.y - a.y) * (b.x - a.x) - (c.x - a.x) * (b.y - a.y);
		if(abs(crossproduct) > 0.0 ){ } //No match, skip to next edge
		else{
			int dotproduct = (c.x - a.x) * (b.x - a.x) + (c.y - a.y)*(b.y - a.y);
			if( dotproduct < 0){ } //No match, skip to next edge;
			else{
				int squaredlengthba = (b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y);
				if(dotproduct > squaredlengthba){ }//No match, skip to next edge;
				else{
					match++; //We pass all tests, we lie on this edge.
				}
			}
		}
	}
	
	return match;
}

//Log probability of coverage fragments mapped in error
double probError(double perr, double coverage){
	return log(perr)*coverage;
}

//Log Poisson Distribution with mean lambda*length evaluated at coverge
double probVariant(double lambda, double length, double coverage){
	double prob;
	
	double newLambda = lambda*length;
	
	int COV = (int) floor(coverage);
	long double logFact = 0;
	logFact = lgamma(COV+1);
	
	prob = COV*log(newLambda) - 1*(newLambda) - logFact;
	
	return prob;
}

/*
 To Dos:
 
 - Change integer/double arrays into vectors without limited size.
 */



/*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*/
/*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*/
/*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*/
           /****************/
           /* MAIN PROGRAM */
           /****************/
/*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*/
/*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*/
/*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*/

int main(int argc, char* argv[]){

 	 int PRINT_FLAG = 0;

    //argc = 2; Parameters file only (cluster file specified in parameters file);
    //argc = 3; Parameters file and clusters file;
    //argc = 4; Parameters file, clusters file, SR file;
    if( (argc != 2 && argc != 3 && argc!=4) ){
                //                 0        1              2                      3              4       5         6       7        8       9           // 10
				cout << "\nGASVPro: Geometric Analysis of Structural Variants, Probabilistic" << endl;
				cout << "Version: 1.2" << endl << endl;
				cout << "Usage: Pass parameter file of the following format AND a clusters file:" << endl;
                cout << "       \tConcordantFile: {path/to/file}\n\tUNIQUEFile: {path/to/file}\n\tLavg: {value} \n\tReadLen: {value}\n\tLambda: {value}\n\tPerr: {value}\n\tLimit: {value}\n\tTolerance: {value}\n\tVerbose: {value}\n\tMaxChrNumber: {Value}\n\tMaxUniqueValue: {value}\n\tMinScaledUniqueness: {Value}\n\tTranslocations: {Value}\n\tTransOnly: {value}\n\tAmbiguous: {value}\n\tLRThreshold: {value}\n\trunningMCMC: {value}\n\n";
                cout << "      ./exe {parametersfile} {clusterfile}\n " << endl;
                cout << "Arguments = " << argc << " (must be 2,3 or 4)\n";
         exit(-1);
	 }
    

	//SET DEFAULTS HERE
	string CLUSTERINPUT;
	string ESPINPUT;
	string UNIQUEINPUT;
    string SPLITREADINPUT;
	int Lavg = 0;
	int ReadLen = 0;
	double Lambda = 0;
	double Perr = 0.01;
	int Limit = 1000;
	double Tolerance = .0001;
	int MAX_CHR_NUMBER = 100;
	int numIgnored = 0;
	double MAX_UNIQUE_VALUE = -1;
	double MIN_SCALED_UNIQUE = 0.3;
	double LRTHRESHOLD = 0; 
	bool PRINTALL = false;
	bool IS_UNIQUE = true;
	bool TRANSLOCATIONS_ON = false;
	bool TRANS_ONLY = false;
	bool AMBIG_MODE = false;
    bool USE_SPLIT_READ = false;
    int READ_LEN_BUFFER = 50;
	int translocationCount = 0;
	bool MAXMODE= false;
	int numbadloc = 0;
	int DIVERGENTS = 0;
	bool MCMCMODE = false;
	

	//Step 0:
	cout << "Step 0: Setting up to run.\n";
	
	//Input Processing
	int *COUNTS_LEFT  = new int[100000];
	int *COUNTS_RIGHT = new int[100000];

	double *UNIQ_LEFT = new double[100000];
	double *UNIQ_RIGHT= new double[100000];
	
	int *STARTS = new int[120000000];
	int *ENDS   = new int[120000000];
	vector<Interval<CHelper*> > deletions;			//*** NEW deletion intervals for tree;
	IntervalTree<CHelper*> deletionTree;
	
	vector<pair<int,int> > inversionIntervals;			//for inversions
	vector<pair<int,int> > relevantGenome;				//for inversions + concordant covg.
	
		
	//int MAX_UNIQUE_VAL = 100000000;
	//int *START_UNIQUE = new int[10000000];
	//int *END_UNIQUE   = new int[10000000];
	//int *VALUE_UNIQUE = new int[10000000];
	vector<Interval<uInt*> > allIntervals;
	IntervalTree<uInt*> uniquenessTree;
	
	int MAX_UNIQUE = 550000;
	int *UNIQUE = new int[550000];
	int *TOTAL  = new int[550000];
	
	if(argc == 2 || argc == 3 || argc == 4)
	{
		cout << "|||INPUT PARAMETERS|||" << endl;
		string temp;
		ifstream p_file(argv[1], ios::in);
		if(p_file.is_open()){cout << "   Parameters file found." << endl;}
		else{cout << "ERROR: Parameters file \"" << argv[1] << "\" cannot be opened." << endl; exit(1);}
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
				
				if(term == "clusterFile:")
				{
						CLUSTERINPUT = value;
						cout << "   Cluster file: " << CLUSTERINPUT << endl;
						continue;
				}
				if(term == "ConcordantFile:") {
						ESPINPUT = value;
						cout << "   Concordant File: " << ESPINPUT << endl;
						continue;
				}
				if(term == "UNIQUEFile:"){
						UNIQUEINPUT = value;
						cout << "   Uniqueness File " << UNIQUEINPUT << endl;
						continue;
				}
				if(term == "Lavg:")
				{
						Lavg = atoi(value.c_str());
						cout << "   Lavg: " << Lavg << endl;
						continue;
				}
				if(term == "ReadLen:"){
						ReadLen = atoi(value.c_str());
						cout << "   Read Length: " << ReadLen << endl;
						continue;
				}
				if(term == "Lambda:"){
						Lambda = atof(value.c_str());
						cout << "   Lambda: " << Lambda << endl;
						continue;
				}
				if(term == "Perr:"){
						Perr = atof(value.c_str());
						cout << "   Perr: " << Perr << endl;
						continue;
				}
				if(term == "Limit:"){
						Limit = atoi(value.c_str());
						cout << "   Limit: " << Limit << endl;
						continue;
				}
				if(term == "Tolerance:"){
						Tolerance = atof(value.c_str());
						cout << "   Tolerance: " << Tolerance << endl;
						continue;
				}
				if(term == "Verbose:"){
						if(value == "Y" || value == "y" || value == "yes"){
							PRINT_FLAG = 1;
							cout << "   ***Verbose Mode Enabled***" << endl;}
						else
							PRINT_FLAG = 0;
				}	
				if(term == "MaxChrNumber:"){
					MAX_CHR_NUMBER = atoi(value.c_str());
					cout << "   Max Chromosome Number: " << MAX_CHR_NUMBER << endl;
				}
				if(term == "MaxUniqueValue:"){
					MAX_UNIQUE_VALUE = atof(value.c_str());
					cout << "   Max Unique Value: " << MAX_UNIQUE_VALUE << endl;
				}
				if(term == "MinScaledUniqueness:"){
					MIN_SCALED_UNIQUE = atof(value.c_str());
					cout << "   Min Scaled Uniqueness: " << MIN_SCALED_UNIQUE << endl;
				}
				if(term == "Translocations:"){
					if(value == "true" || value == "True"){
						TRANSLOCATIONS_ON = true;
						cout << "   Translocations Mode Activated." << endl;
					}
					else
						TRANSLOCATIONS_ON = false;
				}
				if(term == "TransOnly:"){
					if(value == "true" || value == "True"){
						TRANS_ONLY = true;
						cout << "***WARNING: ONLY TRANSLOCATIONS WILL BE DETECTED.***" << endl;
					}
				}
				if(term == "LRThreshold:"){
					if(value == "all" || value == "All"){
						PRINTALL = true;
						cout << "   LR Threshold: PRINT_ALL" << endl;
					}
					else{
						LRTHRESHOLD = atof(value.c_str());
						PRINTALL = false;
						cout << "   LR Threshold: " << LRTHRESHOLD << endl;
					}
				}
				if(term == "Ambiguous:"){
					if(value == "True" || value == "true")
						AMBIG_MODE = true;
				}
				if(term == "maxmode:"){
					if(value == "true")
						MAXMODE = true;
				}
				if(term == "runningMCMC:"){
					if(value == "true"){
						MCMCMODE = true;
						cout << "   MCMC Mode: true" << endl;
					}
				}
			}
		}
	}
	cout << "|||INPUT PARAMETERS|||" << endl;
	
    if(argc == 3 || argc == 4){
        CLUSTERINPUT = argv[2];
        if(argc==4){
            SPLITREADINPUT = argv[3];
            USE_SPLIT_READ= true;
            cout << "   SplitRead Mode = true\n";
        }
    }

	cout << "|||INPUT PARAMETERS|||" << endl << endl;


	bool exit_flag1 = false;
	bool exit_flag2 = false;
	cout << "Checking for input errors...";
	if(Lavg <= 0){ cerr << "\n\t\tERROR: Lavg must be positive\n"; exit(-1); }
	if(ReadLen <= 0){ cerr << "\n\t\tERROR: ReadLen must be positive\n"; exit(-1); }
	if(Lambda <= 0){ cerr << "\n\t\tERROR: Lambda must be positive\n"; exit(-1); }
	if(Perr <= 0){ cerr << "\n\t\tERROR: Perr must be positive\n"; exit(-1); }
	if(Tolerance <=0){ cerr << "\n\t\tERROR: Tolerance must be positive!\n"; exit(-1); }
	if(Limit <=0){ cerr << "Limit on the smallest deletion allowed must be positive!\n"; exit(-1); }
	if(MIN_SCALED_UNIQUE == 1) { cerr << "\nWARNING: A minimum scaled uniqueness value of 1 sets all regions to the same uniqueness value.\n\t"; }
	if(!(Lavg >= 2*ReadLen + READ_LEN_BUFFER)) {exit_flag1 = true;}
	if(!(Lavg >= 2*ReadLen)) {exit_flag2 = true;}
	cout << "OK." << endl;

	
	
	
	if(exit_flag1)
	{
		cout << "\n\t\tERROR: We require Lavg >= 2*ReadLen + (50bp). You may change this 50bp buffer (READ_LEN_BUFFER) in the source.\n";
	}
	if(exit_flag2)
	{
		cout << "\n\t\tERROR: We require that Lavg >= 2*ReadLength. Please check your calculations. If BAMToGASV generated ReadLength and Lavg, please recompute them manually.\n" << endl;
	}
	if(exit_flag1 || exit_flag2){exit(-1);}
	
	char *BPOUT1 = new char[200];
	
	if(MCMCMODE){
		PRINTALL = true;
	}
	
	sprintf(BPOUT1,"%s.GASVPro.coverage",CLUSTERINPUT.c_str());
	ofstream outFile1(BPOUT1,ios::out);
	
	
	char *BPOUT2 = new char[200];
	ofstream outFile2;
	ofstream graph;
	if(!AMBIG_MODE){
		sprintf(BPOUT2,"%s.GASVPro.clusters",CLUSTERINPUT.c_str());
		outFile2.open(BPOUT2,ios::out);
	}
	if(AMBIG_MODE){
		graph.open("graphArgs",ios::out);
	}
	
	outFile2 << "#Cluster_ID:\tNum PRS:\tLocalization:\tType:\tList of PRS:\tLeftChr:\tRightChr:\tBoundary Points:\tLog_Likelihood_Ratio:" << endl << flush;
	outFile1 << "#Cluster_ID:\tType:\tNumDiscordants:\tLog_Prob_Variant:\tLog_Prob_No_Variant:\tBestA:\tBestB:\tCovL:\tCovR:\tMapL:\tMapR:\tcode:" << endl << flush;
	
	
	
	//Max Cluster Size:
	//Note: Have a maximum line length of 15000000; should change this to a string; 
	char* Y1_INIT = new char[15000000];	
//	string Y1_INIT;
	char* CLUSTER_FOR_OUTPUT = new char[15000000];
	char* clusterName = new char[20];
	char* type = new char[20];
	int BUFFER = 1000;
	int index;
	int *X_COORDS = new int[100];
	int *Y_COORDS   = new int[100];
	
	int INT_MAX1 = 200000000;
	
	//cout << "\t(i)Done allocating memory.\n";
	cout << "Done with Step 0.\n" << flush;
	
	cout << "Step 1: Running each chromosome.\n" << flush; 
	
	int LASTSTREAMPOS = 0;
	
	for(int chrNumber = 1; chrNumber<=MAX_CHR_NUMBER && !TRANS_ONLY; chrNumber++){
	//for(int chrNumber = 16; chrNumber<=24; chrNumber++){
		cout << "\n\tStarting Chromosome " << chrNumber << "...." << flush << endl;
		cout << "\tInitializing memory...";
		int LARGESTPOS = 0;
		inversionIntervals.clear();
		relevantGenome.clear();
		allIntervals.clear();
		cout << "done. " << endl;
		
		int clustersSeen = 0;
		
		
		char *INPUT = new char[200];
		//sprintf(INPUT,"%s.chr%i.simple.sorted",ESPINPUT.c_str(),chrNumber);
		ifstream inFile(ESPINPUT.c_str(),ios::in);
		if(!inFile.is_open())
		{
			cout << "Error, cannot open concordant file " << INPUT << endl;
			continue;
		}
		
		char *UNIQUEFILE = new char[200];
		//sprintf(UNIQUEFILE,"%s/chr%i.unique.parsed.sorted",UNIQUEINPUT.c_str(),chrNumber);
		ifstream uFile(UNIQUEINPUT.c_str(),ios::in);
		if(!uFile.is_open())
		{
			cout << "NOTICE: No valid uniqueness file specified (" << UNIQUEINPUT << ")" << endl;
			IS_UNIQUE = false;
			//continue;
		}
		if(IS_UNIQUE && MAX_UNIQUE_VALUE == -1)
		{
			cout << "ERROR: Please input a maximimum uniqueness value and restart." << endl;
			exit(1);
		}
		
		for(int j = 0; j<MAX_UNIQUE;j++){ UNIQUE[j] = 0; TOTAL[j]=0; }
	
		ifstream clusterFile(CLUSTERINPUT.c_str(),ios::in);
        
        ifstream splitReadFile(SPLITREADINPUT.c_str(),ios::in);
        
		if(!clusterFile.is_open())
		{
			cout << "Error, cannot open clusterfile at chr " << chrNumber << "." << endl;
			continue;
		}
	
		//Begin File Processing
		int i = 0;
		int tmpChr,tmpStart,tmpEnd;
		//cout << "\t\tBegin processing ESP file:\n";
		//cout << "\t" << INPUT << "\n";
		/*while( inFile >> tmpStart >> tmpEnd){ 
			if(i%20000000==0){ cout << "\tProcessed " << i << " ESPs\n"; }
			STARTS[i] = tmpStart; ENDS[i] = tmpEnd; i++;
		}*/
		inFile.close();
		int NUM_SEGMENTS = i;
		//cout << "\t\t\tDone -->Total ESPs: " << NUM_SEGMENTS << endl << flush;

		//cout << "\t\tBegin processing Unique File\n";
		i = 0;
		float uniqueVal;
		/*while(uFile >> tmpStart >> tmpEnd >> uniqueVal){
			if(i%5000000==0){ cout << "\tProcessed " << i << " unique\n"; }
			START_UNIQUE[i] = tmpStart; END_UNIQUE[i] = tmpEnd;
			VALUE_UNIQUE[i] = (int) floor(uniqueVal);
			if(i >= MAX_UNIQUE_VAL){ cout << "Increase room for unique values\n"; exit(-1); }
			i++; 
		}*/
		uFile.close();
		int NUM_UNIQUE = i;	
		//cout << "\t\t\tDone --> Total Unique Values: " << NUM_UNIQUE << endl << flush;
	
		
		//BEGIN --> Variables needed for processing clusters:
		
		
		//Misc Variables;
		int numClustersProcessed = 0;
		int numBeyondTolerance;
		struct tm *current;
		time_t now;
		double discordantProbNoVariant, discordantProbHomoVariant, discordantProbHeteroVariant;
		int numDiscordants;
		int caseValue = 0;	
		
		int minA,minB,minCovLeft,minCovRight;
		minA = minB = minCovLeft = minCovRight = -1;

		double minUniqueLeft,minUniqueRight;
		int covLeft,covRight;	
		double uniqueLeft,uniqueRight, scaledCovLeft, scaledCovRight;
		int numLines = 0;
	
	
		//OLD VARIABLES:
		//Scaled LambdaD
		//lambda_d = lambda * (avgFragLength - 2*readLen)/avgFragLength;
		//**but then what is the length to consider? We observe this # of fragments in (essentially) a single position in the genome.
		//**My guess; length is simply Lmax; why? We have Lmax region to look at...
		//double Lambda_d = Lambda*(Lavg - 2*ReadLen)/Lavg; //= Lambda  * (200.0- 2*50.0)/200.0;
		while(	clusterFile.getline(Y1_INIT,1500000000) ){
            
			numLines++;
			strcpy(CLUSTER_FOR_OUTPUT,Y1_INIT);
			numBeyondTolerance = 0;
			caseValue++;	
			double local = 0;
			int localChr = 0;
			int localType = 0;
			int n = 0;
			char* token1 = strtok(Y1_INIT, "\t");
			int coordinate = 0;
			
			while(token1 != NULL){
				if(n == 0){ sprintf(clusterName,token1); }
				else if(n == 1){ numDiscordants = atoi(token1); }
				else if(n == 2){ local = atof(token1); }
				else if(n == 3){ sprintf(type,token1); }
				else if(n == 5){ localChr = atoi(token1); }
				else if(n == 7){
					index = 0;
					char *token2 = strtok(token1,",");
					while(token2 != NULL){
						if(coordinate%2 == 0){ X_COORDS[index] = atoi(token2); }
						else{ Y_COORDS[index] = atoi(token2); index++; }
						token2 = strtok(NULL,", ");
						coordinate++;
						if(index>=99){ cerr << "Need more room for coordinates!\n"; exit(-1);}
					}
				}
				n++;
				token1 = strtok(NULL,"\t");
			}
			
			if(localChr != chrNumber){continue;}

            		if(clusterName[0] == '#'){continue;}
            
			//Note: We can only support clusters that are deletions, inversions, or translocations!
			if(strcmp(type,"D") == 0){ localType = 0; }
			else if(strcmp(type,"I+") == 0 || strcmp(type,"I-") == 0 || strcmp(type,"IR") == 0){ localType = 1;}
			else if(type[0] == 'T'){ localType = 2; translocationCount++; }
			//else if(strcmp(type, "V") == 0 || strcmp(type, "DV") == 0) {localType = 3;}
			else{ 
				if(numIgnored == 0)
					cout << "WARNING: Found first cluster (cluster " << clusterName << ") of non-supported type \"" << type << "\", it will be ignored." << endl;
				cout << "ignored " << type << endl;				
				numIgnored++;
				localType = 4;
			}
			
			//cout << "Cluster:\t" << CLUSTER_FOR_OUTPUT << endl;
			//cout << "LocalType:\t" << localType << endl;
			

			
			if(localChr == chrNumber && (localType <=1 || localType == 3)){
				clustersSeen++;
				int minX, maxX, minY, maxY;
				minX = maxX = X_COORDS[0];
				minY = maxY = Y_COORDS[0];
				for(int k = 0; k<index; k++){
					if(X_COORDS[k]<minX){ minX = X_COORDS[k];}
					if(X_COORDS[k]>maxX){ maxX = X_COORDS[k];}
					if(Y_COORDS[k]<minY){ minY = Y_COORDS[k];}
					if(Y_COORDS[k]>maxY){ maxY = Y_COORDS[k];}
				}
			
					
				if(localType == 0 || localType == 3) 		//i.e. Deletion or Divergent
				{
					bool flag = false;
					if((minY - maxX) < Limit)	//SMALL DELETION, use beRD for uniqueness and coverage
					{
						flag = true;
						if(PRINT_FLAG == 1 && localType == 0) {cout << "\tNOTE: Small deletion, using BeRD model. " << endl << flush;}
						if(PRINT_FLAG == 1 && localType == 3) {cout << "\tNOTE: Small divergent, using BeRD model." << endl << flush;}
						inversionIntervals.push_back(make_pair(maxX, minY));
						inversionIntervals.push_back(make_pair(minX, maxX));
						inversionIntervals.push_back(make_pair(minY, maxY));
						
						uInt* temp = new uInt;
						temp->end = maxX+Lavg;
						temp->bases = 0;
						temp->uniqueness = 0;
						if(minX-Lavg >= 0){temp->start = minX-Lavg; temp->length = abs(maxX-minX)+(2*Lavg); allIntervals.push_back(Interval<uInt*>(minX-Lavg, maxX+Lavg, temp));}
						else{temp->start = minX; temp->length = abs(maxX-minX)+Lavg; allIntervals.push_back(Interval<uInt*>(minX, maxX+Lavg, temp));}
						
						temp = new uInt;
						temp->end = maxY+Lavg;
						temp->bases = 0;
						temp->uniqueness = 0;
						if(minY-Lavg >= 0){temp->start = minY-Lavg; temp->length = abs(maxY-minY)+(2*Lavg); allIntervals.push_back(Interval<uInt*>(minY-Lavg, maxY+Lavg, temp));}
						else{temp->start = minY; temp->length = abs(maxY-minY)+Lavg; allIntervals.push_back(Interval<uInt*>(minY, maxY+Lavg, temp));}
					
						if(minX > LARGESTPOS)
						    LARGESTPOS = minX;
						if(maxX > LARGESTPOS)
						    LARGESTPOS = maxX;
						if(minY > LARGESTPOS)
						    LARGESTPOS = minY;
						if(maxY > LARGESTPOS)
						    LARGESTPOS = maxY;
					}
					
					if(!flag)
					{
						CHelper* temp1 = new CHelper;
						temp1->coverage = 0;
						if(localType == 0){temp1->divergent = false; deletions.push_back(Interval<CHelper*>(maxX, minY, temp1));}
						if(localType == 3){
							temp1->divergent = true;
							deletions.push_back(Interval<CHelper*>(maxX, minY, temp1));
						}
					
						int int_start = min(maxX, minY);
						int int_stop = max(maxX, minY);
						
						uInt* temp = new uInt;
						if(int_start-Lavg >= 0){temp->start = int_start-Lavg; temp->length = abs(minY-maxX) + (2*Lavg);}
						else{temp->start = int_start; temp->length = abs(minY-maxX)+Lavg;}
						temp->end = int_stop+Lavg;	
						temp->bases = 0;
						temp->uniqueness = 0;

						allIntervals.push_back(Interval<uInt*>(int_start, int_stop, temp));
						if(maxX > LARGESTPOS)
							LARGESTPOS = maxX;
						if(minY > LARGESTPOS)
							LARGESTPOS = minY;					
					}
				}
					
				if(localType == 1)
				{				  
					inversionIntervals.push_back(pair<int,int>(minX,maxX));
					inversionIntervals.push_back(pair<int,int>(minY,maxY));
					uInt* temp = new uInt;
					temp->end = maxX+Lavg;
					temp->bases = 0;
					temp->uniqueness = 0;
					if(minX-Lavg >= 0){temp->start = minX-Lavg; temp->length = abs(maxX-minX)+(2*Lavg); allIntervals.push_back(Interval<uInt*>(minX-Lavg, maxX+Lavg, temp));}
					else{temp->start = minX; temp->length = abs(maxX-minX)+Lavg; allIntervals.push_back(Interval<uInt*>(minX, maxX+Lavg, temp));}
					
					temp = new uInt;
					temp->end = maxY+Lavg;
					temp->bases = 0;
					temp->uniqueness = 0;
					if(minY-Lavg >= 0){temp->start = minY-Lavg; temp->length = abs(maxY-minY)+(2*Lavg); allIntervals.push_back(Interval<uInt*>(minY-Lavg, maxY+Lavg, temp));}
					else{temp->start = minY; temp->length = abs(maxY-minY)+Lavg; allIntervals.push_back(Interval<uInt*>(minY, maxY+Lavg, temp));}
				
					if(minX > LARGESTPOS)
					    LARGESTPOS = minX;
					if(maxX > LARGESTPOS)
					    LARGESTPOS = maxX;
					if(minY > LARGESTPOS)
					    LARGESTPOS = minY;
					if(maxY > LARGESTPOS)
					    LARGESTPOS = maxY;
				}
					

			}
			
		}
		clusterFile.close();
        splitReadFile.close();
        
		if(clustersSeen == 0)
		{
			cout << "\tNOTE: No clusters at chr " << chrNumber << " (not including translocations)." << endl;
			continue;
		}
		
		cout << "\t\tSorting intervals for beRD model (" << inversionIntervals.size() << ")....";
		sort(inversionIntervals.begin(), inversionIntervals.end());
		cout << "done." << endl;
		
		
		
		cout << "\t\tProcessing intervals for beRD model....";
		vector<pair<int, int> > result;
		
		if(inversionIntervals.size() != 0)
		{
			vector<pair<int, int> >::iterator LARGER = inversionIntervals.begin();
			pair<int,int> SMALLER = *(LARGER)++;
			while (LARGER != inversionIntervals.end()){
			if (SMALLER.second >= LARGER->first){ // you might want to change it to >=
				SMALLER.second = std::max(SMALLER.second, LARGER->second); 
			} else {
				result.push_back(SMALLER);
				SMALLER = *(LARGER);
			}
			LARGER++;
			}
			result.push_back(SMALLER);
		}
   		cout << "done." << endl;
		

		
		cout << "\t\tBuilding genome vector...";	//QUESTIONABLE

		for(int i = 0; i < result.size(); i++)
		{
			int int_start = min(result[i].first, result[i].second);
			int int_end = max(result[i].first, result[i].second);			
			for(int j = int_start; j <= int_end; j++)
			{
				relevantGenome.push_back(make_pair(j,0));
				
			}
			
		}
		cout << "done. Size: " << relevantGenome.size() << endl;
		
		
		cout << "\t\tBuilding Deletion Interval tree....";
		deletionTree = IntervalTree<CHelper*>(deletions);
		cout << "done." << endl;
		clusterFile.close();
		
		cout << "\t\tBuilding uniqueness tree....";
		uniquenessTree = IntervalTree<uInt*>(allIntervals);
		cout << "done." << endl;

		
		if(IS_UNIQUE){uFile.open(UNIQUEINPUT.c_str(),ios::in);}
		if(IS_UNIQUE){cout << "\t\tAnalyzing Uniqueness...." << endl;}
		int numDone = 0;
		int count = 0;
		while(IS_UNIQUE && (uFile >> tmpChr >> tmpStart >> tmpEnd >> uniqueVal))
		{
			
			if(count%200000 == 0){cout << "\t\t\tProcessed " << count << " fragments." << endl;}
			count++;
			if(tmpChr > chrNumber)
				break;
			if(tmpChr != chrNumber)
				continue;
			vector<Interval<uInt*>* > overlap;
			
			uniquenessTree.findOverlapping(tmpStart, tmpEnd, overlap);
			//LASTSTREAMPOS = uFile.tellg();
			for(int a = 0; a < overlap.size(); a++)
			{
				//deal with overlap here.
				if(tmpStart >= overlap[a]->value->start && tmpEnd >= overlap[a]->value->end)
				{
					overlap[a]->value->bases += overlap[a]->value->end - tmpStart;
					overlap[a]->value->uniqueness += (uniqueVal*(overlap[a]->value->end - tmpStart));
					if(overlap[a]->value->bases == overlap[a]->value->length)
					{
						double preUnique = overlap[a]->value->uniqueness;
						int length = overlap[a]->value->length;
						overlap[a]->value->uniqueness = preUnique/length;
						numDone++;
						//cout << "completed1A" << endl;
					}
					continue;
				}
				if(tmpStart <= overlap[a]->value->start && tmpEnd >= overlap[a]->value->end)
				{
					int length = overlap[a]->value->length;
					overlap[a]->value->bases = length;
					overlap[a]->value->uniqueness = uniqueVal;
					numDone++;
					//cout << "Completed2" << endl;
					continue;
				}
				if(tmpStart <= overlap[a]->value->start && tmpEnd <= overlap[a]->value->end)
				{
					int overlapVal = tmpEnd - overlap[a]->value->start;
					overlap[a]->value->bases += overlapVal;
					overlap[a]->value->uniqueness += (uniqueVal*overlapVal);
					if(overlap[a]->value->bases == overlap[a]->value->length)
					{
						double preUnique = overlap[a]->value->uniqueness;
						int length = overlap[a]->value->length;
						overlap[a]->value->uniqueness = preUnique/length;
						numDone++;
						//cout << "completed3A" << endl;
					}
					continue;
				}	
				if(tmpStart >= overlap[a]->value->start && tmpEnd <= overlap[a]->value->end)
				{
					int overlapVal = tmpEnd-tmpStart;
					overlap[a]->value->bases += overlapVal;
					overlap[a]->value->uniqueness += (uniqueVal*overlapVal);
					if(overlap[a]->value->bases == overlap[a]->value->length)
					{
						double preUnique = overlap[a]->value->uniqueness;
						int length = overlap[a]->value->length;
						overlap[a]->value->uniqueness = preUnique/length;
						numDone++;
						//cout << "completed4A" << endl;
					}
				}
			}			
			i++; 
		}

		
	
		uFile.close();
		
		if(IS_UNIQUE && allIntervals.size() != numDone)
		{
			vector<Interval<uInt*>* > overlap;
			uniquenessTree.findOverlapping(0, LARGESTPOS, overlap);
			for(int i = 0; i < overlap.size(); i++)
			{
				if(overlap[i]->value->bases != overlap[i]->value->length)
				{	
					//cout << overlap[i]->value->bases << " " << overlap[i]->value->length << endl;
					double notCovered = overlap[i]->value->length - overlap[i]->value->bases;
					overlap[i]->value->uniqueness += (notCovered*MAX_UNIQUE_VALUE);
					double preUnique = overlap[i]->value->uniqueness;
					overlap[i]->value->uniqueness = preUnique/overlap[i]->value->length;
				}
			}	
		}
		int uniqErrors = allIntervals.size() - numDone;
		numDone = 0;
		
		
		
		
		inFile.open(ESPINPUT.c_str(), ios::in);
		count = 0;
		cout << "\t\tAnalyzing Concordance...." << endl;
		while( inFile >> tmpChr >> tmpStart >> tmpEnd)
		{
			inFile.ignore(255,'\n');
			if(tmpChr > chrNumber)
				break;
			if(tmpChr != chrNumber)
				continue;
			if(count%200000 == 0){cout << "\t\t\tProcessed " << count << " fragments." << endl;}
			count++;
			vector<Interval<CHelper*>* > overlap;
			deletionTree.findOverlapping(tmpStart,tmpEnd,overlap);
			//cout << overlap.size() << endl;
			for(int i = 0; i < overlap.size(); i++)
			{
				if(!(overlap[i]->value->divergent)){overlap[i]->value->coverage++;}
				else{
					if(tmpStart >= overlap[i]->start && tmpEnd <= overlap[i]->stop)
					{
						 overlap[i]->value->coverage++;
					}
				}
			}	
			overlap.clear();
			
			//If no inversions:
			if(relevantGenome.size() == 0)
			{ 
				continue;
			}
			for(int x = tmpStart; x <= tmpEnd; x++)
			{
				//cout << "going from " << tmpStart << " to " << tmpEnd << " AT " << x << endl;
				pair<int,int> testPair = make_pair(x,0);				
				vector<pair<int,int> >::iterator low = lower_bound(relevantGenome.begin(), relevantGenome.end(), testPair, genomeCompare);
				//cout << "\t found " << (*low).first << endl;
				//cout << "\t was " << (*low).second;
				//cout << "does " << (*low).first << " = " << x << "?" << endl;
				if((*low).first == x)
				{
					//cout << "INCREMENTING" << endl;
					(*low).second++;
				}	
			}	
		}
		inFile.close();

		
		cout << "\t\tRelevant Genome Size: " << relevantGenome.size() << ", # beRD Invervals: " << inversionIntervals.size() << endl;

		
		cout << "\t\tProcessing Clusters\n";
		//**********************************//
		//*** BEGIN READING CLUSTER FILE ***//
		//**********************************//
		//Example Line:
		//0     1     2     3          4            5   6            7
		//c1	1	204.2	D	SRR004856.7363154	1	1	746128, 748510, 746333, 748510, 746027, 748204, 746027, 748409
		clusterFile.open(CLUSTERINPUT.c_str(),ios::in);
        splitReadFile.open(SPLITREADINPUT.c_str(),ios::in);
        
		
		while(!TRANS_ONLY && clusterFile.getline(Y1_INIT,15000000) ){
			strcpy(CLUSTER_FOR_OUTPUT,Y1_INIT);
			numBeyondTolerance = 0;
			caseValue++;	
			double local = 0;
			int localChr = 0;
			int localType = 0;
			int n = 0;
			char* token1 = strtok(Y1_INIT, "\t");
			int coordinate = 0;
			while(token1 != NULL){
				if(n == 0){ sprintf(clusterName,token1); }
				else if(n == 1){ numDiscordants = atoi(token1); }
				else if(n == 2){ local = atof(token1); }
				else if(n == 3){ sprintf(type,token1); }
				else if(n == 5){ localChr = atoi(token1); }
				else if(n == 7){
					index = 0;
					char *token2 = strtok(token1,",");
					while(token2 != NULL){
						if(coordinate%2 == 0){ X_COORDS[index] = atoi(token2); }
						else{ Y_COORDS[index] = atoi(token2); index++; }
						token2 = strtok(NULL,", ");
						coordinate++;
						if(index>=99){ cerr << "Need more room for coordinates!\n"; exit(-1);}
					}
				}
				n++;
				token1 = strtok(NULL,"\t");
			}
			
			//Note: We can only support clusters that are deletions or inversions!
			if(strcmp(type,"D") == 0){ localType = 0; }
			else if(strcmp(type,"I+") == 0 || strcmp(type,"I-") == 0 || strcmp(type,"IR") == 0){ localType = 1;}
			else if(strcmp(type, "V") == 0) { DIVERGENTS++; continue;}//localType = 2; } IGNORING DIVERGENTS
            else if(type[0] == 'T'){ localType = 3; }
		
			//cout << "Cluster:\t" << CLUSTER_FOR_OUTPUT << endl;
			//cout << "LocalType:\t" << localType << endl;
			
            //Read in the same split-read line!
            char *clusterNameSR = new char[20];
            int numSplitReads = 0;
            if(USE_SPLIT_READ){
                char * SR_DATA = new char[100000];
                splitReadFile.getline(SR_DATA,100000);
                char* TK1 = strtok(SR_DATA, "\t");
                sprintf(clusterNameSR,TK1);
                TK1 = strtok(NULL,"\t");
                numSplitReads = atoi(TK1);
            }
            
			if(localChr == chrNumber && localType <=2 ){
				
                
                if(USE_SPLIT_READ){
                    if(strcmp(clusterName,clusterNameSR)==0){
                        //cout << "Hooray we match!\n";
                        //cout << "Clusters: " << clusterName << endl;
                        //cout << "SplitRead: " << clusterNameSR << endl;
                    }
                    else{
                        cout << "Error: Clusters names do not match in the regular and split read file.\n";
                        cout << "Clusters: " << clusterName << endl;
                        cout << "SplitRead: " << clusterNameSR << endl;
                        cout << "Try re-running without the split-read option\n";
                        exit(-1);
                    }
                }
                
                numClustersProcessed++;
				int minX, maxX, minY, maxY;
				minX = maxX = X_COORDS[0];
				minY = maxY = Y_COORDS[0];
				for(int k = 0; k<index; k++){
					if(X_COORDS[k]<minX){ minX = X_COORDS[k];}
					if(X_COORDS[k]>maxX){ maxX = X_COORDS[k];}
					if(Y_COORDS[k]<minY){ minY = Y_COORDS[k];}
					if(Y_COORDS[k]>maxY){ maxY = Y_COORDS[k];}
				}
			
				if(PRINT_FLAG == 1){
					cout << "\t\t\tProcessing Cluster " << numClustersProcessed << "\t"  << clusterName << endl;
					//<< "\t" << numDiscordants << "\t"; 
					//cout << minX << "\t" << maxX << "\t" << minY << "\t" << maxY << "\t" << endl;
					time(&now);
					current = localtime(&now);
					printf("\t\t\tTime:\t%i:%i:%i\n", current->tm_hour, current->tm_min, current->tm_sec);
				}
				//*********************************//
				//*** DONE READING CLUSTER LINE ***//
				//*********************************//
				
				int a,b;
				minA = minB = minCovLeft = minCovRight = -1;
				minUniqueLeft = minUniqueRight = scaledCovLeft = scaledCovRight = 0;
				
				//Values to output.
				double PROB_VARIANT;
				double PROB_NO_VARIANT;
				int code; //Code = 0; used the internal breakpoint location; = 1 used length of deletions;
				int variantCopyNumber;
				int bestA, bestB;
				int NUMCC_L, NUMCC_R;
				double MAP_L, MAP_R; 
				//For the default values if we do NOT find them!
				NUMCC_L = NUMCC_R = -1;
				MAP_L = MAP_R = -1;
				
				double bestHomozygousScore, bestHeterozygousScore, bestErrorScore;
				double bestOverall;
				int bestCopyNum = 0;
				bestOverall = -2*INT_MAX1;
				bestHomozygousScore = -2*INT_MAX1;
				bestHeterozygousScore = -2*INT_MAX1;
				bestErrorScore = -2*INT_MAX1;
						
				//NEW CODE: Want to scale the discordant fragments by uniqueness as well!
				//Begin --> Comment out to remove scaling of discordants by mapability;
				//double xLeftUnique  = getUniqueness(START_UNIQUE,END_UNIQUE,VALUE_UNIQUE,NUM_UNIQUE,minX-Lavg,minX);
				//double xRightUnique = getUniqueness(START_UNIQUE,END_UNIQUE,VALUE_UNIQUE,NUM_UNIQUE,maxX,maxX+Lavg);
				//double yLeftUnique  = getUniqueness(START_UNIQUE,END_UNIQUE,VALUE_UNIQUE,NUM_UNIQUE,minY-Lavg,minY);
				//double yRightUnique = getUniqueness(START_UNIQUE,END_UNIQUE,VALUE_UNIQUE,NUM_UNIQUE,maxY,maxY+Lavg);
				
				int scaledNumDiscordants;
				if(IS_UNIQUE){
					double xLeftUnique  = getUniqueness(minX-Lavg,int(minX), &uniquenessTree, MAX_UNIQUE_VALUE);
					double xRightUnique = getUniqueness(maxX,maxX+Lavg, &uniquenessTree, MAX_UNIQUE_VALUE);
					double yLeftUnique  = getUniqueness(minY-Lavg,minY, &uniquenessTree, MAX_UNIQUE_VALUE);
					double yRightUnique = getUniqueness(maxY,maxY+Lavg, &uniquenessTree, MAX_UNIQUE_VALUE);
					double discordantUniqueness = (xLeftUnique + xRightUnique + yLeftUnique + yRightUnique)/4.0;
					scaledNumDiscordants = (int) (floor)((1.0*numDiscordants)/((MIN_SCALED_UNIQUE)+(1-MIN_SCALED_UNIQUE)*discordantUniqueness));
				}
				else{
					scaledNumDiscordants = 1;
				}
				
				scaledNumDiscordants = numDiscordants;	//this line resets the above calculation.
				
                // Properly Accounting for the length for IR events 
                int discordantLength = Lavg - 2.0*ReadLen;
                if (strcmp(type,"I+") == 0 || strcmp(type,"I-") == 0 || strcmp(type,"D") == 0){
                    discordantLength = Lavg - 2.0*ReadLen;
                }
                else if (strcmp(type,"IR") == 0){
                    discordantLength = 2*(Lavg-2.0*ReadLen);
                }
                else{
                    //Should not be here:
                }
                
				discordantProbNoVariant     = probError(Perr,numDiscordants);
				discordantProbHomoVariant   = probVariant(Lambda,discordantLength, scaledNumDiscordants);
                discordantProbHeteroVariant = probVariant(Lambda/2,discordantLength, scaledNumDiscordants);
                
                if(USE_SPLIT_READ){
                    discordantProbNoVariant     = probError(Perr,numDiscordants+numSplitReads);
                    discordantProbHomoVariant   = probVariant(Lambda,discordantLength, scaledNumDiscordants+numSplitReads); //,covArray,MAX_COV);
                    discordantProbHeteroVariant = probVariant(Lambda/2,discordantLength, scaledNumDiscordants+numSplitReads); //,covArray,MAX_COV);
                }
                
				/*
				cout << "Lambda:\t" << Lambda << endl;
				cout << "Length:\t" << Lavg-2*ReadLen << endl;
				cout << "discordantProbNoVariant:\t" << discordantProbNoVariant << endl;
				cout << "discordantProbHomoVariant:\t" << discordantProbHomoVariant << endl;
				cout << "discordantProbHeteroVariant:\t" << discordantProbHeteroVariant << endl;
				*/
				
				long checking = 0;	
				
				int covBeyondTolerance = 0;
				
				//Q1: Positive Localization;
				//A1: YES
				if(local>0){
					//Q2: Are we below the minsize limit?	
					//A2: YES (or are we inversions!)
					//if( ( ( (minY - maxX) < Limit) && ( (minY-maxX) > 0) ) || localType == 1){
					
					
					if( ( (minY - maxX) < Limit ) || localType == 1)
					{
						if(localType == 2){ cout << ":\t\t\t\tNOTICE: SMALL divergent." << endl; code = 3;}
						if(localType == 0) {
							cout << "\t\t\t\tNOTICE: Small deletion at cluster " << clusterName << " is being processed using beRD algorithm." << endl;
							code = 0;
						}
						code = 0;
						//Store the coverage values;
						for(a = minX; a<=maxX; a++){
							//COUNTS_LEFT[a - minX] = getCoverage(STARTS,ENDS,NUM_SEGMENTS,a,BUFFER);
							//COUNTS_LEFT[a - minX] = chromosomeCoverage[a];
							COUNTS_LEFT[a - minX] = getCoverage(a,&relevantGenome);
							//cout << "cov at " << a << " " << COUNTS_LEFT[a-minX] << endl;
							int start_frag = min(a-Lavg, a+Lavg);
							int end_frag = max(a-Lavg, a+Lavg);
							if(IS_UNIQUE){ UNIQ_LEFT[a-minX] = getUniquenessNEW(start_frag,end_frag, &uniquenessTree, MAX_UNIQUE_VALUE);}//getUniqueness(START_UNIQUE,END_UNIQUE,VALUE_UNIQUE,NUM_UNIQUE,a-Lavg,a);}
							//UNIQ_LEFT[a-minX] = getUniqueness(a-Lavg,a,&uniquenessTree);
							//cout << clusterName << ": " << a << " in " << minX << " " << maxX << endl;
							//cout << "	OLD: " << getCoverage(STARTS,ENDS,NUM_SEGMENTS,a,BUFFER) << " _ NEW: " << COUNTS_LEFT[a - minX] << endl;
						}



						for(b = minY; b<=maxY; b++){		
							//COUNTS_RIGHT[b - minY] = getCoverage(STARTS,ENDS,NUM_SEGMENTS,b,BUFFER);	
							//COUNTS_RIGHT[b - minY] = chromosomeCoverage[b];
							COUNTS_RIGHT[b - minY] = getCoverage(b, &relevantGenome);
							//cout << "cov at " << b << " " << COUNTS_RIGHT[b-minY] << endl;
							int start_frag = min(b-Lavg, b+Lavg);
							int end_frag = max(b-Lavg, b+Lavg);
							if(IS_UNIQUE){ UNIQ_RIGHT[b-minY] = getUniquenessNEW(start_frag,end_frag,&uniquenessTree,MAX_UNIQUE_VALUE);}//getUniqueness(START_UNIQUE,END_UNIQUE,VALUE_UNIQUE,NUM_UNIQUE,b-Lavg,b); }
							//UNIQ_RIGHT[b-minY] = getUniqueness(b-Lavg,b,&uniquenessTree);
							//cout << clusterName << ": " << b << " in " << minY << " " << maxY << endl;
							//cout << "	OLD: " << getCoverage(STARTS,ENDS,NUM_SEGMENTS,b,BUFFER) << " _ NEW: " << COUNTS_RIGHT[b - minY] << endl;
						}

						for(a = minX; a<=maxX; a++){
							for(b = minY; b<=maxY; b++){
								checking++;
								int valid	 = pointInPoly(X_COORDS,Y_COORDS,index,a,b);
								int boundary = onBoundary(X_COORDS,Y_COORDS,index,a,b);

								if(valid == 1 || boundary > 0){
									//cout << "\tThey are valid!\n";
									double scaledUniqLeft = 1;
									double scaledUniqRight = 1;

									if(IS_UNIQUE)
									{
										scaledUniqLeft = ((1-MIN_SCALED_UNIQUE)*UNIQ_LEFT[a-minX]+(MIN_SCALED_UNIQUE));
										scaledUniqRight = ((1-MIN_SCALED_UNIQUE)*UNIQ_RIGHT[b-minY]+(MIN_SCALED_UNIQUE));
									}	
									else{
										scaledUniqLeft = scaledUniqRight = 1;
									}
									scaledCovLeft  = (int) floor((COUNTS_LEFT[a-minX]/scaledUniqLeft)+.5);
									scaledCovRight = (int) floor((COUNTS_RIGHT[b-minY]/scaledUniqRight)+.5); 
							
									/*if(scaledCovLeft + scaledCovRight >= 1000000){
										//cerr << "HUGE Concordant Coverage at:\t" << a << b << endl;
										numBeyondTolerance++;
									}*/
								
									covLeft = COUNTS_LEFT[a-minX]; covRight = COUNTS_RIGHT[b-minY];
									if(IS_UNIQUE){uniqueLeft = UNIQ_LEFT[a - minX]; uniqueRight = UNIQ_RIGHT[b-minY];}
									else{uniqueLeft = uniqueRight = 1;}

									scaledCovLeft = covLeft;   //
									scaledCovRight = covRight; //these two lines reset the above calculations
									
									//Here need to check for deletion, divergent inversion
									//if(divergent){ 
									
									//}
									//else if(deletion, inversions )
									
									double localHomozygous   = probError(Perr,scaledCovLeft) + probError(Perr,scaledCovRight);
									double localHeterozygous = probVariant(Lambda/2,1.0,scaledCovLeft) + probVariant(Lambda/2,1.0,scaledCovRight); 
									double localError        = probVariant(Lambda,1.0,scaledCovLeft) + probVariant(Lambda,1.0,scaledCovRight);
									
									
									double combinedDiffHomo   = (discordantProbHomoVariant + localHomozygous) - (localError + discordantProbNoVariant);
									double combinedDiffHetero = (discordantProbHeteroVariant + localHeterozygous) - (localError + discordantProbNoVariant);
									
									int validC = validCoverage(Lambda*2,scaledCovLeft+scaledCovRight,Tolerance);
									if(validC == 0){
										numBeyondTolerance++;
										covBeyondTolerance = scaledCovLeft+scaledCovRight;
									}
									else{
										if(combinedDiffHomo > bestOverall || combinedDiffHetero > bestOverall){
											bestA = a; bestB = b;
											NUMCC_L = covLeft;
											NUMCC_R = covRight;
											MAP_L = uniqueLeft;
											MAP_R = uniqueRight;
											bestHomozygousScore = localHomozygous;
											bestHeterozygousScore = localHeterozygous;
											bestErrorScore = localError;
											if(combinedDiffHomo > bestOverall){ bestOverall = combinedDiffHomo; bestCopyNum = 2;}
											if(combinedDiffHetero > bestOverall){ bestOverall = combinedDiffHetero; bestCopyNum = 1;}
										}
									}
								}//End if Valid							
							}//End For b
						}//End for a
						
						//Q3: Did we have any points outside of tolerance?
						//A3: Yes;
						if(numBeyondTolerance > 0){
							PROB_NO_VARIANT = 0;
							PROB_VARIANT = -2*INT_MAX1;
							variantCopyNumber = -1;
							bestA = maxX;
							bestB = minY;							
							NUMCC_L = (int) ceil(covBeyondTolerance/2.0);
							NUMCC_R = (int) ceil(covBeyondTolerance/2.0);
							MAP_L = MAP_R = 1.0;

						}
						//A4: No - Everything is fine!
						else{
							PROB_NO_VARIANT = bestErrorScore + discordantProbNoVariant;
							PROB_VARIANT = bestOverall + PROB_NO_VARIANT;
							variantCopyNumber = bestCopyNum;
							//outFile << clusterName << "\t" << numDiscordants << "\t" << type << "\t" << PROB_NO_VARIANT << "\t" << PROB_VARIANT << "\t" << variantCopyNumber << "\t" << bestA << "\t" << bestB << endl;
						}
					}
					//End if Bigger than Limit
					//A2: NO
					else if( (minY - maxX) >= Limit){
						code = 2;
						//Note: We should only have deletions and divergents at this point.
						if(localType == 2){						
							//Need to get coverage, only strict contaimnent;
							int combinedCoverage = getCoverage(maxX, minY, &deletionTree); 
							int lengthConcordant = minY - maxX; //Think!

							double combinedUniqueness = -1;
							double scaledUniqueness = 1;
							if(IS_UNIQUE){combinedUniqueness     = getUniqueness(maxX, minY, &uniquenessTree, MAX_UNIQUE_VALUE);}
							if(IS_UNIQUE){scaledUniqueness       = ((1-MIN_SCALED_UNIQUE)*combinedUniqueness+(MIN_SCALED_UNIQUE));}
							int combinedScaledCoverage           = (int) floor((combinedCoverage/scaledUniqueness)+.5); 
							
							int validC = validCoverage(Lambda*(lengthConcordant), combinedScaledCoverage, Tolerance);
							
							//Note: We should be rounding combinedScaledCoverage to the nearest int.
							//Round --> newNum = floor(num+0.5);
							
							double  prob0Copies, prob1Copy, prob2Copies, probGe2Copies;
							
							//NumDiscordants -> scaledNumDiscordants;
							
							//Prob0Copies --> Everything is an error;
							prob0Copies = probError(Perr,combinedScaledCoverage+scaledNumDiscordants);
							
							//Prob1Copy   --> Discordants are an error; concordants are half coverage
							prob1Copy = probError(Perr,scaledNumDiscordants) + probVariant(Lambda/2,lengthConcordant,combinedScaledCoverage);
							
							//Prob2Copies --> Discordants are error; concordants are at full coverage;
							prob2Copies =  probError(Perr,scaledNumDiscordants) + probVariant(Lambda,lengthConcordant,combinedScaledCoverage);
							
							//probGe2Copies = 1 - prob0Copies -prob1Copy - prob2Copies
							probGe2Copies = 1 - prob0Copies -prob1Copy - prob2Copies;
							
							
							if(validC == 0){
								numBeyondTolerance++;
								PROB_NO_VARIANT = 0;
								PROB_VARIANT = -2*INT_MAX1;
								variantCopyNumber = -1;
								bestA = maxX;
								bestB = minY;
								NUMCC_L = combinedCoverage;
								MAP_L = combinedUniqueness; //0.7*combinedUniqueness+0.3;
								NUMCC_R = 0;
								MAP_R = 0;
								//outFile << clusterName << "\t" << numDiscordants << "\t" << type << "\t" << 0				<< "\t" << -2*INT_MAX1	<< "\t" << -1				<< "\t" << -1		<< "\t" << -1 << endl;
							}
							else{
								
								//Need to define differently?
								PROB_VARIANT = probGe2Copies;
								PROB_NO_VARIANT = prob2Copies;
								
								double combinedDiff = PROB_VARIANT - PROB_NO_VARIANT;
								if(combinedDiff > bestOverall){
									bestA = maxX; bestB = minY;
									NUMCC_L = combinedCoverage;
									NUMCC_R = 0;
									MAP_L   = combinedUniqueness;
									MAP_R = 0;
									bestCopyNum = 3; //Really want to say >2
									bestOverall = combinedDiff;
								
								}
								variantCopyNumber = 3;
								PROB_VARIANT = probGe2Copies;
								PROB_NO_VARIANT = prob2Copies;
				
								
							}
							
							continue; 
						}				
						code = 1;
						double combinedUniqueness = 1;
						double scaledUniqueness = 1;
						//int combinedCoverage          = getCoverage(STARTS,ENDS,NUM_SEGMENTS,maxX,minY,BUFFER);
						int combinedCoverage = getCoverage(maxX, minY, &deletionTree);
						if(IS_UNIQUE){combinedUniqueness     = getUniquenessNEW(maxX-Lavg, minY+Lavg, &uniquenessTree, MAX_UNIQUE_VALUE);}
						
						if(IS_UNIQUE){scaledUniqueness = ((1-MIN_SCALED_UNIQUE)*combinedUniqueness+(MIN_SCALED_UNIQUE));}
						int combinedScaledCoverage    = (int) floor((combinedCoverage/scaledUniqueness)+.5);
						
						//int validC = validCoverage(Lambda*(minY-maxX-Lavg), combinedScaledCoverage, Tolerance);
						int lengthConcordant = minY - maxX + Lavg;
						int validC = validCoverage(Lambda*(lengthConcordant), combinedScaledCoverage, Tolerance);
						//Q3: Are we outside of tolerance? (i.e., too high concordant coverage?)
						//A3:YES;
						if(validC == 0){
							numBeyondTolerance++;
							PROB_NO_VARIANT = 0;
							PROB_VARIANT = -2*INT_MAX1;
							variantCopyNumber = -1;
							bestA = maxX;
							bestB = minY;
							NUMCC_L = combinedCoverage;
							MAP_L = combinedUniqueness; //0.7*combinedUniqueness+0.3;
							NUMCC_R = 0;
							MAP_R = 0;
							//outFile << clusterName << "\t" << numDiscordants << "\t" << type << "\t" << 0				<< "\t" << -2*INT_MAX1	<< "\t" << -1				<< "\t" << -1		<< "\t" << -1 << endl;
						}
						else{
							double localHomozygous		= probError(Perr,combinedScaledCoverage);
							double localError			= probVariant(Lambda,lengthConcordant,combinedScaledCoverage); //,covArray,MAX_COV);
							double localHeterozygous	= probVariant(Lambda/2,lengthConcordant,combinedScaledCoverage); //,covArray,MAX_COV);
							
							double combinedDiffHomo   = (discordantProbHomoVariant + localHomozygous) - (localError + discordantProbNoVariant);
							double combinedDiffHetero = (discordantProbHeteroVariant + localHeterozygous) - (localError + discordantProbNoVariant);
							
							if(combinedDiffHomo > bestOverall || combinedDiffHetero > bestOverall){
								bestA = maxX; bestB = minY;
								NUMCC_L = combinedCoverage;
								NUMCC_R = 0;
								MAP_L = combinedUniqueness;
								MAP_R = 0;
								bestHomozygousScore = localHomozygous;
								bestHeterozygousScore = localHeterozygous;
								bestErrorScore = localError;
								if(combinedDiffHomo > bestOverall){ bestOverall = combinedDiffHomo; bestCopyNum = 2;}
								if(combinedDiffHetero > bestOverall){ bestOverall = combinedDiffHetero; bestCopyNum = 1;}
							}
						
							PROB_NO_VARIANT = bestErrorScore + discordantProbNoVariant;
							PROB_VARIANT = bestOverall + PROB_NO_VARIANT;
							variantCopyNumber = bestCopyNum;
							
							//outFile << clusterName << "\t" << numDiscordants << "\t" << type << "\t" << PROB_NO_VARIANT << "\t" << PROB_VARIANT << "\t"<< variantCopyNumber << "\t" << bestA << "\t" << bestB << endl;	
						}
				
					}
					//A2: We have negative length? Or cross ourselves?
					//We cross ourselves, this is most unfortunate. We'll ignore CC in this case. (hopefully!)
					else{
						cout << "Erroneous cluster; skipping:\t" << "cluster Name: " << clusterName << "\t" << "Distance: " << minY << " - " << maxX << " = " << minY - maxX << endl;
						code = -1;
						PROB_NO_VARIANT = 0;
						PROB_VARIANT = -2*INT_MAX1;
						variantCopyNumber = -1;
						bestA = maxX;
						bestB = minY;
						NUMCC_L = NUMCC_R = -1;
						MAP_L = MAP_R = -1;
					}
				}
				//A1: NO
				// -1 localization;
				else{
					if(MAXMODE){ numbadloc++; }					
					code = -1;
					PROB_NO_VARIANT = 0;
					PROB_VARIANT = -2*INT_MAX1;
					variantCopyNumber = -1;
					bestA = -1;
					bestB = -1;
					NUMCC_L = NUMCC_R = -1;
					MAP_L = MAP_R = -1;

				}
				
				//Write to output files:
				//(1) Write to the GASVCC file;
				//Format: ClusterName Type NumDiscordants	Start	End		NumCC	Mapability	GC
				//		  c1			D	1				1347	218875	10079	0.0254701	0.329
				//		  c2			D	1				66312	557898	29317	0.0138073	0.29977

				/*
				For coverage file:
				    output the LRCLUSTER = PROB_VARIANT - PROB_NO_VARIANT
					
					ClusterName type NumDiscordants LRCLUSTER
				*/

				double LRCLUSTER = PROB_VARIANT - PROB_NO_VARIANT;
				//If LRCLUSTER >= LR then output
				//outFile1 is the coverage file
				//outFile2 is the clusters file
				
				if((LRCLUSTER >= LRTHRESHOLD || PRINTALL) && (code != -1)){
				outFile1 << clusterName << "\t" << type  << "\t" << numDiscordants << "\t" << PROB_VARIANT << "\t" << PROB_NO_VARIANT << "\t" << bestA << "\t" << bestB << "\t" << NUMCC_L << "\t" << NUMCC_R << "\t" << MAP_L << "\t" << MAP_R << "\t" << code << endl;
				if(!AMBIG_MODE){
					outFile2 << CLUSTER_FOR_OUTPUT << "\t" << LRCLUSTER;
					outFile2 << endl;
				}
				}
				
				//cout << clusterName << "\t" << variantCopyNumber << endl;
	
				//(2) Write to the clusters file; (more complicated!!)
				// 0    1        2      3                      4                        5       6        7
				//c1	1	189	D	chr17_41883_42105_0:0:0_2:0:0_15789d3	17	17	41642, 42623, 41821, 42623, 41532, 42334, 41532, 42513
				//token1 = strtok(CLUSTER_FOR_OUTPUT, "\t");
				//n = 0;
		
			/*
			LRCLUSTER = PROB_VARIANT - PROB_NO_VARIANT;
			
			LR - User Parameter (default is 0.5 or all
			
			if(LRCLUSTER>=LR){ //output cluster}
			else{ skip cluster}
		
			*/
		
				/*while(token1 != NULL){
					
					if(n == 0){ 
						sprintf(clusterName,token1); 
						outFile2 << clusterName << "_" << PROB_VARIANT << "_" << PROB_NO_VARIANT << "\t";
					}
					else{ 
						outFile2 << token1; 
						if(n<=6){ outFile2 << "\t"; }
					}
					n++;
					token1 = strtok(NULL,"\t");
				}*/
				//not output the LRCluster before this last endl.

							
			}//End if the right chromosome number
		
	}//End While
	inversionIntervals.clear();
	relevantGenome.clear();
	cout << "\nDONE w/ Chromosome; Processed --> " << numClustersProcessed << endl << flush;
	if(IS_UNIQUE && uniqErrors > 0){cout << "***NOTICE: There were " << uniqErrors << " uniqueness coverage errors. The missing values were replaced with the Max Unique Value parameter and scaled." << endl;}

	}//End For
	int numBeyondTolerance;
	int caseValue;
	int numDiscordants;
	
	if(TRANSLOCATIONS_ON){
		for(int i = 1; i <= MAX_CHR_NUMBER-1; i++)
		{
			for(int j = i+1; j <= MAX_CHR_NUMBER; j++)
			{
				int numclusters = 0;
				
				vector<pair<int,int> > translocationsLEFT;
				vector<pair<int,int> > relevantGenomeLEFT;
				vector<pair<int,int> > translocationsRIGHT;
				vector<pair<int,int> > relevantGenomeRIGHT;
				
				vector<Interval<uInt*> > allIntervalsLEFT;
				vector<Interval<uInt*> > allIntervalsRIGHT;
				
				int LARGESTPOSLEFT = 0;
				int LARGESTPOSRIGHT = 0;
				
				ifstream clusterFile(CLUSTERINPUT.c_str(), ios::in);
				while(	clusterFile.getline(Y1_INIT,15000000) ){
					strcpy(CLUSTER_FOR_OUTPUT,Y1_INIT);
					numBeyondTolerance = 0;
					caseValue++;	
					double local = 0;
					int LEFTChr = 0;
					int RIGHTChr = 0;
					int localType = 0;
					int n = 0;
					char* token1 = strtok(Y1_INIT, "\t");
					int coordinate = 0;
					while(token1 != NULL){
						if(n == 0){ sprintf(clusterName,token1); }
						else if(n == 1){ numDiscordants = atoi(token1); }
						else if(n == 2){ local = atof(token1); }
						else if(n == 3){ sprintf(type,token1); }
						else if(n == 5){ LEFTChr = atoi(token1); }
						else if(n == 6){ RIGHTChr = atoi(token1); }
						else if(n == 7){
							index = 0;
							char *token2 = strtok(token1,",");
							while(token2 != NULL){
								if(coordinate%2 == 0){ X_COORDS[index] = atoi(token2); }
								else{ Y_COORDS[index] = atoi(token2); index++; }
								token2 = strtok(NULL,", ");
								coordinate++;
								if(index>=99){ cerr << "Need more room for coordinates!\n"; exit(-1);}
							}
						}
						n++;
						token1 = strtok(NULL,"\t");
					}
					
					if(RIGHTChr > j)
						break;
					
					if(type[0] == 'T')
						localType = 2;
						
					
			
					if(localType == 2 && LEFTChr == i && RIGHTChr == j)
					{

						numclusters++;
						int minX, maxX, minY, maxY;
						minX = maxX = X_COORDS[0];
						minY = maxY = Y_COORDS[0];
						for(int k = 0; k<index; k++){
							if(X_COORDS[k]<minX){ minX = X_COORDS[k];}
							if(X_COORDS[k]>maxX){ maxX = X_COORDS[k];}
							if(Y_COORDS[k]<minY){ minY = Y_COORDS[k];}
							if(Y_COORDS[k]>maxY){ maxY = Y_COORDS[k];}
						}
						
						translocationsLEFT.push_back(make_pair(minX, maxX));
						translocationsRIGHT.push_back(make_pair(minY, maxY));
						
						
						uInt* temp = new uInt;
						temp->start = minX;
						temp->end = maxX;
						temp->length = abs(maxX-minX);
						temp->bases = 0;
						temp->uniqueness = 0;
						allIntervalsLEFT.push_back(Interval<uInt*>(minX, maxX, temp));
						
						temp = new uInt;
						temp->start = minY;
						temp->end = maxY;
						temp->length = abs(maxY-minY);
						temp->bases = 0;
						temp->uniqueness = 0;
						allIntervalsRIGHT.push_back(Interval<uInt*>(minY, maxY, temp));
						
						if(minX > LARGESTPOSLEFT)
							LARGESTPOSLEFT = minX;
						if(maxX > LARGESTPOSLEFT)
							LARGESTPOSLEFT = maxX;
						if(minY > LARGESTPOSRIGHT)
							LARGESTPOSRIGHT = minY;
						if(maxY > LARGESTPOSRIGHT)
							LARGESTPOSRIGHT = maxY;
						
					}//endif
				}//endclusterlinewhile
				clusterFile.close();
				
				if(translocationsLEFT.size() == 0)
					continue;				
				
				cout << "\t\tSorting translocation intervals on " << i << " " << j << " (" << translocationsLEFT.size() << "," << translocationsRIGHT.size() << ")....";
				sort(translocationsLEFT.begin(), translocationsLEFT.end());
				sort(translocationsRIGHT.begin(),translocationsRIGHT.end());
				cout << "done." << endl;
				
				cout << "\t\tProcessing translocation intervals....";
				vector<pair<int, int> > resultLEFT;
				vector<pair<int, int> > resultRIGHT;
				

				
				if(translocationsLEFT.size() != 0)
				{
					vector<pair<int, int> >::iterator LARGER = translocationsLEFT.begin();
					pair<int,int> SMALLER = *(LARGER)++;
					while (LARGER != translocationsLEFT.end()){
					if (SMALLER.second >= LARGER->first){ // you might want to change it to >=
						SMALLER.second = std::max(SMALLER.second, LARGER->second); 
					} else {
						resultLEFT.push_back(SMALLER);
						SMALLER = *(LARGER);
					}
					LARGER++;
					}
					resultLEFT.push_back(SMALLER);
				}
				
				if(translocationsRIGHT.size() != 0)
				{
					vector<pair<int, int> >::iterator LARGER = translocationsRIGHT.begin();
					pair<int,int> SMALLER = *(LARGER)++;
					while (LARGER != translocationsRIGHT.end()){
					if (SMALLER.second >= LARGER->first){ // you might want to change it to >=
						SMALLER.second = std::max(SMALLER.second, LARGER->second); 
					} else {
						resultRIGHT.push_back(SMALLER);
						SMALLER = *(LARGER);
					}
					LARGER++;
					}
					resultRIGHT.push_back(SMALLER);
				}
				cout << "done." << endl;
				
				cout << "\t\tBuilding relevant translocation genome vectors....";
				for(int x = 0; x < resultLEFT.size(); x++)
				{
					for(int y = resultLEFT[x].first; y <= resultLEFT[x].second; y++)
					{
						relevantGenomeLEFT.push_back(make_pair(y,0));
					}	
				}
				for(int x = 0; x < resultRIGHT.size(); x++)
				{
					for(int y = resultRIGHT[x].first; y <= resultRIGHT[x].second; y++)
					{
						relevantGenomeRIGHT.push_back(make_pair(y,0));
					}
				}
				cout << "done." << endl;
				
				
				ifstream inFile(ESPINPUT.c_str(), ios::in);
				int tmpChr, tmpStart, tmpEnd;
				int count = 0;
				cout << "\t\tAnalyzing Concordance...." << endl;
				if(!inFile.is_open())
				{
					cout << "Concordant file error (file does not exist). Program aborted." << endl;
					exit(22);
				}
				while(inFile >> tmpChr >> tmpStart >> tmpEnd)
				{
					inFile.ignore(255,'\n');
					if(tmpChr > j)
						break;
					if(tmpChr != j && tmpChr != i)
						continue;
					if(count%200000 == 0){cout << "\t\t\tProcessed " << count << " fragments." << endl;}
					count++;

					if(tmpChr == i){
						//cout << "\tFOUND ON LEFT" << endl;
						for(int x = tmpStart; x <= tmpEnd; x++)
						{
							pair<int,int> testPair = make_pair(x,0);				
							vector<pair<int,int> >::iterator low = lower_bound(relevantGenomeLEFT.begin(), relevantGenomeLEFT.end(), testPair, genomeCompare);
							
							if((*low).first == x)
							{
								(*low).second++;
								//cout << "\t" <<  (*low).second << endl;
							}			
						}
					}
					if(tmpChr == j){
						//cout << "\tFOUND ON RIGHT" << endl;
						for(int x = tmpStart; x <= tmpEnd; x++)
						{
							pair<int,int> testPair = make_pair(x,0);				
							vector<pair<int,int> >::iterator low = lower_bound(relevantGenomeRIGHT.begin(), relevantGenomeRIGHT.end(), testPair, genomeCompare);
							
							if((*low).first == x)
							{
								(*low).second++;
								//cout << "\t" <<  (*low).second << endl;
							}			
						}
					}
				}
				inFile.close();
				cout << "done" << endl;
				
				//****UNIQUENESS CALCULATION PROCESS****//
				//**Uncomment this code AND some code***//
				//***in cluster processing to enable****//
				//*****uniqueness calculation for ******//
				//**********translocations.*************//
				
				/*cout << "\t\tBuilding uniqueness tree....";
				IntervalTree<uInt*> uniquenessTreeLEFT = IntervalTree<uInt*>(allIntervalsLEFT);
				IntervalTree<uInt*> uniquenessTreeRIGHT = IntervalTree<uInt*>(allIntervalsRIGHT);
				cout << "done." << endl;
		
				ifstream uFile;
				if(IS_UNIQUE){uFile.open(UNIQUEINPUT.c_str(),ios::in);}
				if(IS_UNIQUE){cout << "\t\tAnalyzing Uniqueness...." << endl;}
				if(!uFile.is_open()) {
					cout << "wrong" << endl;}
				int numDoneLEFT = 0;
				int numDoneRIGHT = 0;
				count = 0;
				int uniqueVal;
				while(IS_UNIQUE && (uFile >> tmpChr >> tmpStart >> tmpEnd >> uniqueVal))
				{
					
					if(count%200000 == 0){cout << "\t\t\tProcessed " << count << " fragments." << endl;}
					count++;
					if(tmpChr > j)
						break;
					if(tmpChr != i && tmpChr != j)
						continue;
					
						
					vector<Interval<uInt*>* > overlap;
					
					
					if(tmpChr == i){
						uniquenessTreeLEFT.findOverlapping(tmpStart, tmpEnd, overlap);
						//LASTSTREAMPOS = uFile.tellg();
						for(int a = 0; a < overlap.size(); a++)
						{
							//deal with overlap here.
							if(tmpStart >= overlap[a]->value->start && tmpEnd >= overlap[a]->value->end)
							{
								overlap[a]->value->bases += overlap[a]->value->end - tmpStart;
								overlap[a]->value->uniqueness += (uniqueVal*(overlap[a]->value->end - tmpStart));
								if(overlap[a]->value->bases == overlap[a]->value->length)
								{
									double preUnique = overlap[a]->value->uniqueness;
									int length = overlap[a]->value->length;
									overlap[a]->value->uniqueness = preUnique/length;
									numDoneLEFT++;
									//cout << "completed1AL" << endl;
								}
								continue;
							}
							if(tmpStart <= overlap[a]->value->start && tmpEnd >= overlap[a]->value->end)
							{
								int length = overlap[a]->value->length;
								overlap[a]->value->bases = length;
								overlap[a]->value->uniqueness = uniqueVal;
								numDoneLEFT++;
								//cout << "Completed2L" << endl;
								continue;
							}
							if(tmpStart <= overlap[a]->value->start && tmpEnd <= overlap[a]->value->end)
							{
								int overlapVal = tmpEnd - overlap[a]->value->start;
								overlap[a]->value->bases += overlapVal;
								overlap[a]->value->uniqueness += (uniqueVal*overlapVal);
								if(overlap[a]->value->bases == overlap[a]->value->length)
								{
									double preUnique = overlap[a]->value->uniqueness;
									int length = overlap[a]->value->length;
									overlap[a]->value->uniqueness = preUnique/length;
									numDoneLEFT++;
									//cout << "completed3AL" << endl;
								}
								continue;
							}	
							if(tmpStart >= overlap[a]->value->start && tmpEnd <= overlap[a]->value->end)
							{
								int overlapVal = tmpEnd-tmpStart;
								overlap[a]->value->bases += overlapVal;
								overlap[a]->value->uniqueness += (uniqueVal*overlapVal);
								if(overlap[a]->value->bases == overlap[a]->value->length)
								{
									double preUnique = overlap[a]->value->uniqueness;
									int length = overlap[a]->value->length;
									overlap[a]->value->uniqueness = preUnique/length;
									numDoneLEFT++;
									//cout << "completed4AL" << endl;
								}
							}
						}
					}
					if(tmpChr == j){
						uniquenessTreeRIGHT.findOverlapping(tmpStart, tmpEnd, overlap);
						//LASTSTREAMPOS = uFile.tellg();
						for(int a = 0; a < overlap.size(); a++)
						{
							//deal with overlap here.
							if(tmpStart >= overlap[a]->value->start && tmpEnd >= overlap[a]->value->end)
							{
								overlap[a]->value->bases += overlap[a]->value->end - tmpStart;
								overlap[a]->value->uniqueness += (uniqueVal*(overlap[a]->value->end - tmpStart));
								if(overlap[a]->value->bases == overlap[a]->value->length)
								{
									double preUnique = overlap[a]->value->uniqueness;
									int length = overlap[a]->value->length;
									overlap[a]->value->uniqueness = preUnique/length;
									numDoneRIGHT++;
									//cout << "completed1AR" << endl;
								}
								continue;
							}
							if(tmpStart <= overlap[a]->value->start && tmpEnd >= overlap[a]->value->end)
							{
								int length = overlap[a]->value->length;
								overlap[a]->value->bases = length;
								overlap[a]->value->uniqueness = uniqueVal;
								numDoneRIGHT++;
								//cout << "Completed2R" << endl;
								continue;
							}
							if(tmpStart <= overlap[a]->value->start && tmpEnd <= overlap[a]->value->end)
							{
								int overlapVal = tmpEnd - overlap[a]->value->start;
								overlap[a]->value->bases += overlapVal;
								overlap[a]->value->uniqueness += (uniqueVal*overlapVal);
								if(overlap[a]->value->bases == overlap[a]->value->length)
								{
									double preUnique = overlap[a]->value->uniqueness;
									int length = overlap[a]->value->length;
									overlap[a]->value->uniqueness = preUnique/length;
									numDoneRIGHT++;
									//cout << "completed3AR" << endl;
								}
								continue;
							}	
							if(tmpStart >= overlap[a]->value->start && tmpEnd <= overlap[a]->value->end)
							{
								int overlapVal = tmpEnd-tmpStart;
								overlap[a]->value->bases += overlapVal;
								overlap[a]->value->uniqueness += (uniqueVal*overlapVal);
								if(overlap[a]->value->bases == overlap[a]->value->length)
								{
									double preUnique = overlap[a]->value->uniqueness;
									int length = overlap[a]->value->length;
									overlap[a]->value->uniqueness = preUnique/length;
									numDoneRIGHT++;
									//cout << "completed4AR" << endl;
								}
							}
						}
					}
						
					//i++; 
				}
				
				uFile.close();
				
				if(IS_UNIQUE && allIntervalsLEFT.size() != numDoneLEFT)
				{
					vector<Interval<uInt*>* > overlap;
					uniquenessTreeLEFT.findOverlapping(0, LARGESTPOSLEFT, overlap);
					for(int i = 0; i < overlap.size(); i++)
					{
						if(overlap[i]->value->bases != overlap[i]->value->length)
						{	
							//cout << overlap[i]->value->bases << " " << overlap[i]->value->length << endl;
							double notCovered = overlap[i]->value->length - overlap[i]->value->bases;
							overlap[i]->value->uniqueness += (notCovered*MAX_UNIQUE_VALUE);
							double preUnique = overlap[i]->value->uniqueness;
							overlap[i]->value->uniqueness = preUnique/overlap[i]->value->length;
						}
					}	
				}
				if(IS_UNIQUE && allIntervalsRIGHT.size() != numDoneLEFT)
				{
					vector<Interval<uInt*>* > overlap;
					uniquenessTreeRIGHT.findOverlapping(0, LARGESTPOSRIGHT, overlap);
					for(int i = 0; i < overlap.size(); i++)
					{
						if(overlap[i]->value->bases != overlap[i]->value->length)
						{	
							//cout << overlap[i]->value->bases << " " << overlap[i]->value->length << endl;
							double notCovered = overlap[i]->value->length - overlap[i]->value->bases;
							overlap[i]->value->uniqueness += (notCovered*MAX_UNIQUE_VALUE);
							double preUnique = overlap[i]->value->uniqueness;
							overlap[i]->value->uniqueness = preUnique/overlap[i]->value->length;
						}
					}	
				}*/
				
						
				clusterFile.open(CLUSTERINPUT.c_str(),ios::in);
                
				//BEGIN --> Variables needed for processing clusters:
				
				//Misc Variables;
				int numClustersProcessed = 0;
				int numBeyondTolerance;
				struct tm *current;
				time_t now;
				double discordantProbNoVariant, discordantProbHomoVariant, discordantProbHeteroVariant;
				int numDiscordants;
				int caseValue = 0;	
				
				int minA,minB,minCovLeft,minCovRight;
				minA = minB = minCovLeft = minCovRight = -1;
		
				double minUniqueLeft,minUniqueRight;
				int covLeft,covRight;	
				double uniqueLeft,uniqueRight, scaledCovLeft, scaledCovRight;
				int numLines = 0;
				
				
				
				cout << "\t\tProcessing translocation clusters..." << endl;
				while(clusterFile.getline(Y1_INIT,15000000) ){
					strcpy(CLUSTER_FOR_OUTPUT,Y1_INIT);
					numBeyondTolerance = 0;
					caseValue++;	
					double local = 0;
					int LEFTChr = 0;
					int RIGHTChr = 0;
					int localType = 0;
					int n = 0;
					char* token1 = strtok(Y1_INIT, "\t");
					int coordinate = 0;
					while(token1 != NULL){
						if(n == 0){ sprintf(clusterName,token1); }
						else if(n == 1){ numDiscordants = atoi(token1); }
						else if(n == 2){ local = atof(token1); }
						else if(n == 3){ sprintf(type,token1); }
						else if(n == 5){ LEFTChr = atoi(token1); }
						else if(n == 6){ RIGHTChr = atoi(token1); }
						else if(n == 7){
							index = 0;
							char *token2 = strtok(token1,",");
							while(token2 != NULL){
								if(coordinate%2 == 0){ X_COORDS[index] = atoi(token2); }
								else{ Y_COORDS[index] = atoi(token2); index++; }
								token2 = strtok(NULL,", ");
								coordinate++;
								if(index>=99){ cerr << "Need more room for coordinates!\n"; exit(-1);}
							}
						}
						n++;
						token1 = strtok(NULL,"\t");
					}
					
					//Note: We can only support clusters that are deletions or inversions!
					
					if(type[0] == 'T')
						localType = 2;
					
					
					//cout << "Cluster:\t" << CLUSTER_FOR_OUTPUT << endl;
					//cout << "LocalType:\t" << localType << endl;
					
					if(localType == 2 && LEFTChr == i && RIGHTChr == j){
						numClustersProcessed++;
						int minX, maxX, minY, maxY;
						minX = maxX = X_COORDS[0];
						minY = maxY = Y_COORDS[0];
						for(int k = 0; k<index; k++){
							if(X_COORDS[k]<minX){ minX = X_COORDS[k];}
							if(X_COORDS[k]>maxX){ maxX = X_COORDS[k];}
							if(Y_COORDS[k]<minY){ minY = Y_COORDS[k];}
							if(Y_COORDS[k]>maxY){ maxY = Y_COORDS[k];}
						}
						
						if(PRINT_FLAG == 0)
						{
							// 	double pct = (numClustersProcessed/numclusters)*100.0;
							cout << "\t\t\tProcessed " << numClustersProcessed << " of " << numclusters << " clusters. \r" << flush;
							//cout << "\t\t\t" << pct << "% of clusters processed. \r" << flush;
						}
						
						if(PRINT_FLAG == 1){
							cout << "\t\t\tProcessing Cluster " << numClustersProcessed << "\t"  << clusterName << endl;
							//<< "\t" << numDiscordants << "\t"; 
							//cout << minX << "\t" << maxX << "\t" << minY << "\t" << maxY << "\t" << endl;
							time(&now);
							current = localtime(&now);
							printf("\t\t\tTime:\t%i:%i:%i\n", current->tm_hour, current->tm_min, current->tm_sec);
						}
						//*********************************//
						//*** DONE READING CLUSTER LINE ***//
						//*********************************//
						
						int a,b;
						minA = minB = minCovLeft = minCovRight = -1;
						minUniqueLeft = minUniqueRight = scaledCovLeft = scaledCovRight = -1;
						
						//Values to output.
						double PROB_VARIANT;
						double PROB_NO_VARIANT;
						int code; //Code = 0; used the internal breakpoint location; = 1 used length of deletions;
						int variantCopyNumber;
						int bestA, bestB;
						int NUMCC_L, NUMCC_R;
						double MAP_L, MAP_R; 
						//For the default values if we do NOT find them!
						NUMCC_L = NUMCC_R = -1;
						MAP_L = MAP_R = -1;
						
						double bestHomozygousScore, bestHeterozygousScore, bestErrorScore;
						double bestOverall;
						int bestCopyNum = 0;
						bestOverall = -2*INT_MAX1;
						bestHomozygousScore = -2*INT_MAX1;
						bestHeterozygousScore = -2*INT_MAX1;
						bestErrorScore = -2*INT_MAX1;
						
						
						int scaledNumDiscordants;
						/*if(IS_UNIQUE){
							double xLeftUnique  = getUniqueness(minX-Lavg,int(minX), &uniquenessTreeLEFT, MAX_UNIQUE_VALUE);
							double xRightUnique = getUniqueness(maxX,maxX+Lavg, &uniquenessTreeLEFT, MAX_UNIQUE_VALUE);
							double yLeftUnique  = getUniqueness(minY-Lavg,minY, &uniquenessTreeRIGHT, MAX_UNIQUE_VALUE);
							double yRightUnique = getUniqueness(maxY,maxY+Lavg, &uniquenessTreeRIGHT, MAX_UNIQUE_VALUE);
							double discordantUniqueness = (xLeftUnique + xRightUnique + yLeftUnique + yRightUnique)/4.0;
							scaledNumDiscordants = (int) (floor)((1.0*numDiscordants)/((MIN_SCALED_UNIQUE)+(1-MIN_SCALED_UNIQUE)*discordantUniqueness));
						}
						else
							scaledNumDiscordants = 1;*/
							
						scaledNumDiscordants = numDiscordants;	//this line resets the above calculation.
						
                        int discordantLength = (Lavg-2.0*ReadLen);
                        if (strcmp(type,"TN+1") == 0 || strcmp(type,"TN+2") == 0 || strcmp(type,"TN-1") == 0 || strcmp(type,"TN-2") == 0){
                            discordantLength = (Lavg-2.0*ReadLen);
                        }
                        else if (strcmp(type,"TR+") == 0 || strcmp(type,"TR-") == 0){
                            discordantLength = 2*(Lavg-2.0*ReadLen);
                        }
                        else{
                            //Should not be here at this point.
                        }
                        
						discordantProbNoVariant     = probError(Perr,numDiscordants);
                        discordantProbHomoVariant   = probVariant(Lambda,discordantLength, scaledNumDiscordants);
                        discordantProbHeteroVariant = probVariant(Lambda/2,discordantLength, scaledNumDiscordants);
                        
                        //Need to add split-reads for translocations!!

						/*
						cout << "Lambda:\t" << Lambda << endl;
						cout << "Length:\t" << Lavg-2*ReadLen << endl;
						cout << "discordantProbNoVariant:\t" << discordantProbNoVariant << endl;
						cout << "discordantProbHomoVariant:\t" << discordantProbHomoVariant << endl;
						cout << "discordantProbHeteroVariant:\t" << discordantProbHeteroVariant << endl;
						*/
						
						long checking = 0;	
						
						int covBeyondTolerance = 0;
						
						//Q1: Positive Localization;
						//A1: YES
						if(local>0){
							//Q2: Are we below the minsize limit?	
							//A2: YES (or are we inversions!)
							//if( ( ( (minY - maxX) < Limit) && ( (minY-maxX) > 0) ) || localType == 1){
							
							if( ( (minY - maxX) < Limit ) || localType == 2)
							{
								if(localType == 0) {cout << "\t\t\t\tNOTICE: Small deletion at cluster " << clusterName << " is being processed using inversion algorithm." << endl;}
								code = 0;
								
								//Store the coverage values;
								for(a = minX; a<=maxX; a++){
									//COUNTS_LEFT[a - minX] = getCoverage(STARTS,ENDS,NUM_SEGMENTS,a,BUFFER);
									//COUNTS_LEFT[a - minX] = chromosomeCoverage[a];
									//COUNTS_LEFT[a - minX] = getCoverage(a,&relevantGenomeLEFT);
		
									//if(IS_UNIQUE){ UNIQ_LEFT[a-minX] = getUniqueness(a-Lavg,a, &uniquenessTreeLEFT, MAX_UNIQUE_VALUE);}
								}
		
		
		
								for(b = minY; b<=maxY; b++){		
									//COUNTS_RIGHT[b - minY] = getCoverage(STARTS,ENDS,NUM_SEGMENTS,b,BUFFER);	
									//COUNTS_RIGHT[b - minY] = chromosomeCoverage[b];
									//COUNTS_RIGHT[b - minY] = getCoverage(b, &relevantGenomeRIGHT);
									//if(IS_UNIQUE){ UNIQ_RIGHT[b-minY] = getUniqueness(b-Lavg,b,&uniquenessTreeRIGHT,MAX_UNIQUE_VALUE);}
								}
		
								for(a = minX; a<=maxX; a++){
									for(b = minY; b<=maxY; b++){
										checking++;
										int valid	 = pointInPoly(X_COORDS,Y_COORDS,index,a,b);
										int boundary = onBoundary(X_COORDS,Y_COORDS,index,a,b);
		
										if(valid == 1 || boundary > 0){
											//cout << "\tThey are valid!\n";
											double scaledUniqLeft = 1;
											double scaledUniqRight = 1;
		
											if(true)
											{
												scaledUniqLeft = ((1-MIN_SCALED_UNIQUE)*UNIQ_LEFT[a-minX]+(MIN_SCALED_UNIQUE));
												scaledUniqRight = ((1-MIN_SCALED_UNIQUE)*UNIQ_RIGHT[b-minY]+(MIN_SCALED_UNIQUE));
											}
											
											scaledCovLeft  = COUNTS_LEFT[a-minX]/scaledUniqLeft;
											scaledCovRight = COUNTS_RIGHT[b-minY]/scaledUniqRight; 
									
											/*if(scaledCovLeft + scaledCovRight >= 1000000){
												//cerr << "HUGE Concordant Coverage at:\t" << a << b << endl;
												numBeyondTolerance++;
											}*/
										
											covLeft = COUNTS_LEFT[a-minX]; covRight = COUNTS_RIGHT[b-minY];
											uniqueLeft = UNIQ_LEFT[a - minX]; uniqueRight = UNIQ_RIGHT[b-minY];
		
											scaledCovLeft = covLeft;   //
											scaledCovRight = covRight; //these two lines reset the above calculations
											
											double localHomozygous   = probError(Perr,scaledCovLeft) + probError(Perr,scaledCovRight);
											double localHeterozygous = probVariant(Lambda/2,1.0,scaledCovLeft) + probVariant(Lambda/2,1.0,scaledCovRight); 
											double localError        = probVariant(Lambda,1.0,scaledCovLeft) + probVariant(Lambda,1.0,scaledCovRight);
											
											
											double combinedDiffHomo   = (discordantProbHomoVariant + localHomozygous) - (localError + discordantProbNoVariant);
											double combinedDiffHetero = (discordantProbHeteroVariant + localHeterozygous) - (localError + discordantProbNoVariant);
											
											int validC = validCoverage(Lambda*2,scaledCovLeft+scaledCovRight,Tolerance);
											if(validC == 0){
												numBeyondTolerance++;
												covBeyondTolerance = scaledCovLeft+scaledCovRight;
											}
											else{
												if(combinedDiffHomo > bestOverall || combinedDiffHetero > bestOverall){
													bestA = a; bestB = b;
													NUMCC_L = covLeft;
													NUMCC_R = covRight;
													MAP_L = uniqueLeft;
													MAP_R = uniqueRight;
													bestHomozygousScore = localHomozygous;
													bestHeterozygousScore = localHeterozygous;
													bestErrorScore = localError;
													if(combinedDiffHomo > bestOverall){ bestOverall = combinedDiffHomo; bestCopyNum = 2;}
													if(combinedDiffHetero > bestOverall){ bestOverall = combinedDiffHetero; bestCopyNum = 1;}
												}
											}
										}//End if Valid							
									}//End For b
								}//End for a
								
								//Q3: Did we have any points outside of tolerance?
								//A3: Yes;
								if(numBeyondTolerance > 0){
									PROB_NO_VARIANT = 0;
									PROB_VARIANT = -2*INT_MAX1;
									variantCopyNumber = -1;
									bestA = maxX;
									bestB = minY;							
									NUMCC_L = (int) ceil(covBeyondTolerance/2.0);
									NUMCC_R = (int) ceil(covBeyondTolerance/2.0);
									MAP_L = MAP_R = 1.0;
		
								}
								//A4: No - Everything is fine!
								else{
									PROB_NO_VARIANT = bestErrorScore + discordantProbNoVariant;
									PROB_VARIANT = bestOverall + PROB_NO_VARIANT;
									variantCopyNumber = bestCopyNum;
									//outFile << clusterName << "\t" << numDiscordants << "\t" << type << "\t" << PROB_NO_VARIANT << "\t" << PROB_VARIANT << "\t" << variantCopyNumber << "\t" << bestA << "\t" << bestB << endl;
								}
							}
							//End if Bigger than Limit
							//A2: NO
							//A2: We have negative length? Or cross ourselves?
							//We cross ourselves, this is most unfortunate. We'll ignore CC in this case. (hopefully!)
							else{
								cout << "SHOULD NEVER BE HERE!!!!\t" << "cluster Name: " << clusterName << "\t" << "Distance: " << minY << " - " << maxX << " = " << minY - maxX << endl;
								code = -1;
								PROB_NO_VARIANT = 0;
								PROB_VARIANT = -2*INT_MAX1;
								variantCopyNumber = -1;
								bestA = maxX;
								bestB = minY;
								NUMCC_L = NUMCC_R = -1;
								MAP_L = MAP_R = -1;
							}
						}
						//A1: NO
						// -1 localization;
						else{
							code = -1;
							PROB_NO_VARIANT = 0;
							PROB_VARIANT = -2*INT_MAX1;
							variantCopyNumber = -1;
							bestA = -1;
							bestB = -1;
							NUMCC_L = NUMCC_R = -1;
							MAP_L = MAP_R = -1;
		
						}
						
						//Write to output files:
						//(1) Write to the GASVCC file;
						//Format: ClusterName Type NumDiscordants LogProbVar LogProbNoVar	Start	End		NumCCLeft NumCCR	MapabilityL MapabilityR	Code
						//		  c1			D	1				1347	218875	10079	0.0254701	0.329
						//		  c2			D	1				66312	557898	29317	0.0138073	0.29977
						//If code = 0 or 1 depending on inversion or deletion model
						/*
						For coverage file:
							output the LRCLUSTER = PROB_VARIANT - PROB_NO_VARIANT
							
							ClusterName type NumDiscordants LRCLUSTER
						*/
		
						double LRCLUSTER = PROB_VARIANT - PROB_NO_VARIANT;
						//If LRCLUSTER >= LR then output
						//outFile1 is the coverage file
						//outFile2 is the clusters file
						
						if((LRCLUSTER >= LRTHRESHOLD || PRINTALL) && (code != -1)){
						outFile1 << clusterName << "\t" << type  << "\t" << numDiscordants << "\t" << PROB_VARIANT << "\t" << PROB_NO_VARIANT << "\t" << bestA << "\t" << bestB << "\t" << NUMCC_L << "\t" << NUMCC_R << "\t" << MAP_L << "\t" << MAP_R << "\t" << code << endl;
						if(!AMBIG_MODE){
							outFile2 << CLUSTER_FOR_OUTPUT << "\t" << LRCLUSTER;
							outFile2 << endl;
						}
						}
						//cout << clusterName << "\t" << variantCopyNumber << endl;
			
						//(2) Write to the clusters file; (more complicated!!)
						// 0    1        2      3                      4                        5       6        7
						//c1	1	189	D	chr17_41883_42105_0:0:0_2:0:0_15789d3	17	17	41642, 42623, 41821, 42623, 41532, 42334, 41532, 42513
						//token1 = strtok(CLUSTER_FOR_OUTPUT, "\t");
						//n = 0;
				
					/*
					LRCLUSTER = PROB_VARIANT - PROB_NO_VARIANT;
					
					LR - User Parameter (default is 0.5 or all
					
					if(LRCLUSTER>=LR){ //output cluster}
					else{ skip cluster}
				
					*/
				
						/*while(token1 != NULL){
							
							if(n == 0){ 
								sprintf(clusterName,token1); 
								outFile2 << clusterName << "_" << PROB_VARIANT << "_" << PROB_NO_VARIANT << "\t";
							}
							else{ 
								outFile2 << token1; 
								if(n<=6){ outFile2 << "\t"; }
							}
							n++;
							token1 = strtok(NULL,"\t");
						}*/
						//not output the LRCluster before this last endl.
		
									
					}//End if the right chromosome number
				
				}//End While
				
				
			}//endjfor
			allIntervals.clear();
		}//endifor
	} 
	
	if(AMBIG_MODE){
		string dirname = CLUSTERINPUT + ".GASVPro.graphdir";
		graph << CLUSTERINPUT << endl;
		graph << BPOUT1 << endl;
		graph << dirname << endl;
		graph.close();
	}
	
	outFile1.close();
	if(!AMBIG_MODE){ outFile2.close(); }
	if(MAXMODE) {cout << "\n***NOTICE: There were " << numbadloc << " clusters with bad localization (-1). " << endl;}
	cout << "***NOTICE: Ignored " << numIgnored << " clusters of non-supported types. See GASVPro documentation. " << endl;
	if(!TRANSLOCATIONS_ON){cout << "***NOTICE: Found " << translocationCount << " translocations. Enable translocation mode or TransOnly mode to process them. " << endl;}
	if(DIVERGENTS > 0){cout << "***NOTICE: Found " << DIVERGENTS << " divergent variants that were not processed (see documentation)." << endl;}	
	return 0;	
}

/******************************/
/* BEGIN FUNCTION DEFINITIONS */
/******************************/

//Point in Poly

/***************************************************************************
 *   INPOLY.C                                                              *
 *   Copyright (c) 1995-1996 Galacticomm, Inc.  Freeware source code.      *
 *                                                                         *
 *   Please feel free to use this source code for any purpose, commercial  *
 *   or otherwise, as long as you don't restrict anyone else's use of      *
 *   this source code.  Please give credit where credit is due.            *
 *                                                                         *
 *   Point-in-polygon algorithm, created especially for World-Wide Web     *
 *   servers to process image maps with mouse-clickable regions.           *
 *                                                                         *
 *	 http://www.visibone.com/inpoly/inpoly.c.txt                           *
 *   http://www.visibone.com/inpoly/inpoly.c                               *
 *                                                                         *
 *                                       6/19/95 - Bob Stein & Craig Yap   *
 *                                       stein@visibone.com                *
 *                                       craig@cse.fau.edu                 *
 ***************************************************************************/
/*int									1=inside, 0=outside                */
/*inpoly(								is target point inside a 2D polygon*/
/*	   unsigned int poly[][2],			polygon points, [0]=x, [1]=y       */
/*	   int npoints						number of points in polygon        */
/*	   unsigned int xt					x (horizontal) of target point     */
/*	   unsigned int yt					y (vertical) of target point       */
/***************************************************************************/
int pointInPoly(int *STARTS, int *ENDS, int npoints, unsigned int xt, unsigned int yt){
	unsigned int xnew,ynew;
	unsigned int xold,yold;
	unsigned int x1,y1;
	unsigned int x2,y2;
	int i;
	int inside=0;
	
	if (npoints < 3) {
		return(0);
	}
	xold=STARTS[npoints-1];
	yold=ENDS[npoints-1];
	for (i=0 ; i < npoints ; i++) {
		xnew=STARTS[i]; //poly[i][0];
		ynew=ENDS[i]; //poly[i][1];
		if (xnew > xold) {
			x1=xold;
			x2=xnew;
			y1=yold;
			y2=ynew;
		}
		else {
			x1=xnew;
			x2=xold;
			y1=ynew;
			y2=yold;
		}
		if ((xnew < xt) == (xt <= xold)         /* edge "open" at left end */
			&& ((long)yt-(long)y1)*(long)(x2-x1)
            < ((long)y2-(long)y1)*(long)(xt-x1)) {
			inside=!inside;
		}
		xold=xnew;
		yold=ynew;
	}
	return(inside);
}


//Find Location
int findLocation(int* sortedArray, int NUM_SEGMENTS, int key){
	int first = 0;
	int last = NUM_SEGMENTS - 1;
	int mid;
	int FOUND = 0;
 
   	while (first <= last && FOUND == 0) {
       		mid = (first + last) / 2;  // compute mid point.
       		if (key > sortedArray[mid]) 
           		first = mid + 1;  // repeat search in top half.
       		else if (key < sortedArray[mid]) 
           		last = mid - 1; // repeat search in bottom half.
       		else
           		FOUND = 1;     // found it. return position /////
   	}
	if(FOUND == 0){
  		//return -(first + 1);    // failed to find key
		return first; //This is where we should insert it.
	}
	else{ 
		while( key == sortedArray[mid] && mid>=0){ mid--; }
		return mid;
	}
}

//Get Uniqueness:
double getUniqueness(int* START_UNIQUE,int* END_UNIQUE,int *VALUE_UNIQUE,int NUM_UNIQUE,int startpoint,int endpoint){
	int totalSum = 0;
	double UNIQUE_VAL;
	int locationStart, locationEnd,tmpStart,tmpEnd,uniqueOffBeginning,uniqueOffEnd;
	uniqueOffBeginning = uniqueOffEnd = 0;
	int location = findLocation(START_UNIQUE,NUM_UNIQUE,startpoint);
	if(location > 0 && START_UNIQUE[location] != startpoint){ location--; }
	locationStart = location;
	int length = endpoint - startpoint + 1;
	while( (endpoint)>=START_UNIQUE[location] && location<NUM_UNIQUE){
		tmpStart = START_UNIQUE[location] - startpoint;
		if(tmpStart<0){ tmpStart = 0;}
		tmpEnd = END_UNIQUE[location] - startpoint;
		if(tmpEnd > length){ tmpEnd = (length-1); }
		if(tmpEnd <0){ tmpEnd = 0;}
		for(int i = tmpStart;i<tmpEnd;i++){
			totalSum += VALUE_UNIQUE[location];
		}
		location++;
	}
	
	locationEnd = location;			
	
	if( (locationStart == 0) && (locationEnd == 0) ){
		uniqueOffBeginning = 1;
		UNIQUE_VAL = -1;
	}
	else if( (locationStart == (NUM_UNIQUE-1) ) && (locationEnd == (NUM_UNIQUE-1) ) ){
		uniqueOffEnd = 1;
		UNIQUE_VAL = -1;
	}
	else{
		UNIQUE_VAL = (totalSum*1.0)/(length*1.0)/(37.0);
	}
	
	return UNIQUE_VAL;
}

//Currently Inefficient: Should update to be an interval tree calculation.
//Get Coverge for a point.
int getCoverage(int *STARTS, int *ENDS,int NUM_SEGMENTS, int point, int BUFFER){
	//Get Coverage:
	int location = findLocation(STARTS,NUM_SEGMENTS,point-BUFFER);
	//cerr << "\t\tFound ESP location: " << STARTS[location] << endl;
	int total = 0;
	while( (point+BUFFER)>=STARTS[location] && location>=0){ //just need point? not point+buffer
		//NEW COVERAGE -- Just has to overlap
		if(STARTS[location]<=point && ENDS[location]>=point){
			total++;		
		}
		location++;
	}
	
	return total;
}

//Get Coverge for an interval;
//Currently Inefficient: Should update to be an interval tree calculation
int getCoverage(int *STARTS, int *ENDS,int NUM_SEGMENTS, int startpoint, int endpoint, int BUFFER){
	//Get Coverage:
	int location = findLocation(STARTS,NUM_SEGMENTS,startpoint-BUFFER);
	
	int locationStart = findLocation(STARTS,NUM_SEGMENTS,startpoint);
	int locationEnd   = findLocation(STARTS,NUM_SEGMENTS,endpoint);
	
	int total = locationEnd - locationStart;
	
	//cerr << "\t\tFound ESP location: " << STARTS[location] << endl;
	//int total = 0;
	//while( (endpoint+BUFFER)>=STARTS[location] && location>=0){
	
	//Check the ones that are before!
	while(location < locationStart){
		//NEW COVERAGE -- Just has to overlap
		if(STARTS[location]>=endpoint || ENDS[location]<=startpoint){
			//No overlap START>=end or END<=start
		}
		else{
			total++;
		}
		location++;
	}
	
	return total;
}

int getCoverage(int start, int end, IntervalTree<CHelper*> *tree)
{
	vector<Interval<CHelper*> > results;
	tree->findContained(start, end, results);
	if(results.size() == 0){return 010101010;}
	for(int i = 0; i < results.size(); i++)
	{
		if(results[i].start == start && results[i].stop == end)
			return results[i].value->coverage;
	}
	return 0101010101;
}

bool intervalCompare(pair<int,int> intervalA, pair<int,int> intervalB)
{
	if(intervalA.first < intervalB.first)
		return true;
		
	else 
		return false;
}

bool genomeCompare(pair<int,int> pointA, pair<int,int> pointB)
{
	if(pointA.first < pointB.first)
		return true;
	else
		return false;
}

int getCoverage(int point, vector<pair<int,int> >* genome)
{
	vector<pair<int,int> >::iterator low;
	pair<int,int> testPair = make_pair(point,0);
	low = lower_bound(genome->begin(), genome->end(), testPair, genomeCompare);
	
	return (*low).second;
}

double getUniqueness(int start, int end, IntervalTree<uInt*> *tree, double MAX_UNIQUE_VAL)
{
	vector<Interval<uInt*> > results;
	tree->findContained(start, end, results);
	if(results.size() == 0){return 200.2;}
	
	for(int i = 0; i < results.size(); i++)
	{
		cout << "\t\tfor " << start << " " << end << " FOUND " << results.size() << endl; 
		cout << "\t\t" << results[i].value->start << " "  << results[i].value->end << endl;
		if(results[i].value->start >= start && results[i].value->end <= end)
		{
			return results[i].value->uniqueness/(MAX_UNIQUE_VAL);
		}
	}
	

	return -200.1;
}

double getUniquenessNEW(int start, int end, IntervalTree<uInt*> *tree, double MAX_UNIQUE_VAL)
{
	vector<Interval<uInt*> > results;
	tree->findOverlapping(start, end, results);
	if(results.size() == 0){return MAX_UNIQUE_VAL;} //NOTE: should be changed to MAX_UNIQUE_VAL 
	double sum = 0;
	double goodones = 0.0;
	for(int i = 0; i < results.size(); i++)
	{
	    if(start >= results[i].value->start && end <= results[i].value->end) { goodones++; sum += (results[i].value->uniqueness)/MAX_UNIQUE_VAL;}
	}
	return sum/goodones;
}
