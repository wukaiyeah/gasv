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

#include <vector>
#include <limits.h>
#include <algorithm>
#include <list>

//STD Library Linking
#include<sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
using namespace std;

//int numMobileESPs;
//vector <int> mobileESPs;

#ifndef PROBLEMINSTANCE_H
#define PROBLEMINSTANCE_H

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


int findLocationLocal(vector<int>&locations,int size, int key){
	int first = 0;
	int last = size-1;
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
	if(FOUND == 0){
		//cout << "\t NO: did not find it.\n";
		return -(first + 1);    // failed to find key
	}
	else{ 
		//cout << "\tYES: Found it at index " << mid << "\n";
		return mid;
	}
}

/*?*?*?*?*?*?*?*?*?*?*?*?*/
// Probability Functions //
/*?*?*?*?*?*?*?*?*?*?*?*?*/

double probVariant(int dis_n, int length, double coverage, int num){
	double retVal;
	
	double mean = length*coverage;
	double mean_2 = mean/2.0;
		
	//gamma(n) = (n-1)!
	
	if(num == 2){
		retVal = -1*mean + dis_n *log(mean) - lgamma(dis_n+1); 
	}
	else if(num == 1){
		retVal = -1*mean_2 + dis_n*log(mean_2) - lgamma(dis_n+1);
	}
	else{
		cerr << "Invalid Copy Number --> " << num << endl; exit(-1);
	}
	return retVal;
}

class ESP{
	public:
		//ESP(): possibleVariants(100) {state = -1; possibleMappings = 0; }
		ESP(){ possibleMappings = 0; currentState = proposedState = -2; } //The -2 state lets us know we are uninitialized to anything.
		~ESP(){};
	
	int setInitialMapping(int map){
		int variant;
		if(map<possibleMappings){
			currentState = map;
			variant = currentState; //Will be equal to -1 if unchanged;
			if(currentState>=0){ variant = possibleVariants[currentState]; }
		}
		else{
			cerr << "Error in ESP setMapping\n";
			cerr << "Only " << possibleMappings << " possible, but trying to access mapping " << map << endl << flush;
			exit(-1);
		}
		return variant;
	}
	
		int setUpForMCMC(){
			//(1) Sort the variants; 
			sort( possibleVariants.begin(), possibleVariants.end() );
			
			//(2) Set all the occupancy values as 0;
			std::vector<int>::iterator ip;
			variantOccupancy.resize(possibleVariants.size());
			for(ip = variantOccupancy.begin(); ip != variantOccupancy.end(); ip++){  *ip = 0;} 
			return 0;
		}
		
		int getCurrentState(){ return currentState; }
		int currentlyUnAssigned(){ if(currentState == -1){ return 1;} else{ return 0;} }
		int getProposedState(){ return proposedState; }
		int proposedUnAssigned(){ if(proposedState == -1){ return 1;} else{ return 0;} }

		int proposeRandomMove(double PERR){
			int MAP = chooseRandomMapping(PERR);
			proposedState = MAP;
			return 0;
		}
	
		//Returns 1 if ESP supports variant V and 0 otherwise.
		int canESPSupport(int V){
			int loc = findLocationLocal(possibleVariants,possibleMappings,V);
			if(loc >= 0){ return 1; }
			else{ return 0; }
		}
		int proposeSwapMove(int V){
			int loc = findLocationLocal(possibleVariants,possibleMappings,V);
			if(loc >= 0){ proposedState = loc; }
			return 0;
		}
	
		//Propose a move to variant A; return the current variant assigned to or -1 if nothing.
		int proposeAddMove(int A){
			int loc = findLocationLocal(possibleVariants,possibleMappings,A);
			if(loc<0){ cerr << "Impossible proposed add move, check ESP for consistency!\n"; exit(-1); }
			else{
				//We need to propose this move;
				proposedState = loc;
				if(currentState>=0){ return getVariant(currentState); }
				else{ return -1; }
			}
		}
	
		//Propose a remove from variant A; return to another variant or -1 with prob. perr:
		//Return either the new proposed variant or -1 for an error;
		int proposeRemoveMove(int A, double perr){
			//Set to error w/ prob PERR
			double p = (1.0*rand())/(RAND_MAX*1.0);
			
			//If only one possible mapping, we move this to an error;
			if(possibleMappings == 1){ proposedState = -1; return -1; }
			
			//Otherwise, move to an error with prob perr and otherwise select another mapping.
			if(p < perr){ proposedState = -1; return -1; }
			else{
				int randomMapping = rand() % (possibleMappings);
				proposedState = randomMapping;
				while(proposedState == currentState){
					randomMapping = rand() % (possibleMappings);
					proposedState = randomMapping;
				}
				return getVariant(proposedState);
			}
		}
	
		int getVariant(int index){ 
			if(index == -1){ return -1; } 
			else if(index<possibleMappings){ return possibleVariants[index]; } 
			else{ cerr << "Invalid variant index in ESP named " << name << "!\n"; exit(-1); return -2; }
		}
		
		//Note: We have to report the total possible mappings:
		int getOccupancy(int i,int SAMPLE){
			int retVal = -1;
			if(i == -1){
				retVal = SAMPLE;
				for(int j = 0; j<possibleMappings; j++){ 
					//cout << "Mapping " << j << " had occupancy " << variantOccupancy[j] << "\t";
					retVal -= variantOccupancy[j]; 
					//cout << "RetVal = " << retVal << endl;
				}
				return retVal;
			}
			else{
				retVal = variantOccupancy[i];
			}
			return retVal;
		}
	
		int getPossibleMappings(){ return possibleMappings;}
		int getMapping(){ return currentState;}
		
		int setProposedMapping(){ currentState = proposedState; return 0; } //Set a proposed mapping;
	
		int setMapping(int MAP){
			if(MAP == -1 || MAP<possibleMappings){ currentState = MAP; return 0;}
			else{ 	cerr << "Error in SetMapping! --> Attempting to set a variant to an impossible mapping\n"; exit(-1); }
		}
	
		int chooseRandomMapping(double PERR){
			//Set to error w/ prob PERR
			double p = (1.0*rand())/(RAND_MAX*1.0);
			if(p < PERR){ return -1; }
		
			//Set to a location with the remaining probability; [from 0 to (possibleMappings-1)]
			int randomMapping =  rand() % (possibleMappings);
			if(randomMapping>=0 && randomMapping<possibleMappings){ return randomMapping; }
			else{ cerr << "Error in choosing random mapping in class ESP.\n"; exit(-1); }
		}
	
	int initializeRandomMapping(double perr){
		currentState =  chooseRandomMapping(perr);
		int variant = currentState; //Will be equal to -1 if unchanged;
		if(currentState>=0){ variant = possibleVariants[currentState]; }
		return variant;
	}
	
	int matchProposedCurrent(){
		proposedState = currentState;
		for(int i = 0; i<possibleMappings; i++){ variantOccupancy[i] = 0; }
		if(currentState>=0){ variantOccupancy[currentState]++; }
		
		return 0;
	}
	
	//Initialize the occupancy;
	int clearOccupancy(){
		for(int i = 0; i<possibleMappings; i++){ variantOccupancy[i] = 0; }
		return 0;
	}
	
	void outputOccupancy(int total,ofstream & outESP){
		int errorPosition = total;
		int nonErrorOccupancy = 0;
		for(int i = 0; i<possibleMappings;i++){ 
			outESP << getVariant(i) << "\t" << variantOccupancy[i] << endl; 
			nonErrorOccupancy+=variantOccupancy[i];
			errorPosition -= variantOccupancy[i];
		}
		outESP << "Occupancy not in errors -> " << nonErrorOccupancy << endl;
		outESP << -1 << "\t" << errorPosition << endl;
	}
	
	int updateOccupancy(){
		if(currentState>=0) {variantOccupancy[currentState]++; }
		return 0;
	}
	
	int updateOccupancy(int STEPS){
		if(currentState>=0){ variantOccupancy[currentState]+=STEPS; }
		return 0;
	}

	int acceptMove(){
		currentState = proposedState;
		return 0;
	}
	
	int setOptimalPosition(double THRESHOLD){
		int TRUE_MAPPING = -1;
		int VARIANT = -1;
		for(int i = 0; i<possibleMappings; i++){
			if(variantOccupancy[i]>=THRESHOLD){ TRUE_MAPPING = i; }
		}
		
		if(TRUE_MAPPING >= 0){
			VARIANT = getVariant(TRUE_MAPPING);
		}
		
		return VARIANT;
	}
	
	int rejectMove(){ proposedState = currentState; return 0; }
	
	int addVariant(int variant){ possibleVariants.push_back(variant); possibleMappings++; return 0;}
	
	int setName(int i){ name = i; return 0;}
	int getName(){ return name;}
	private:
		int name; //This is simply its value from being read in the db.
		int currentState; //Either -1 or an integer with the cluster it is mapped to!
		int proposedState;
		int possibleMappings;
		std::vector<int> possibleVariants;
		std::vector<int> variantOccupancy; //Note if this is an error, we don't need it.
};

class variant{
	public:
	variant(){ 
		currentAssigned = 0; possibleAssigned = 0; numDependencies = 0; runningLogLikelihood = 0; runningLikelihoodVar = 0;
	}
	
	~variant(){};
		
		int checkMobility(int numMobileESPs, vector<int> &mobileESPs){
			int currentMobile = 0; //Assume we are NOT mobile;
			for(int i = 0; i<currentAssigned;i++){
				int loc = findLocationLocal(mobileESPs,numMobileESPs,currentESPs[i]);
				if(loc>=0){ currentMobile++; }
			}
			//All are mobile!
			if(currentMobile == currentAssigned){ 
				return 1;
			}
			//This variant contains some immobile ESPs, can not move it.
			else{
				return 0;
			}
		}
	
		int getProposedAssigned(){ return proposedAssigned;}
		int getCurrentAssigned(){ return currentAssigned; }
		int getCurrentESP(int index){ return currentESPs[index]; }
		
		int getPossibleAssigned(){ return possibleAssigned; }
		int getPossibleESP(int index){ return possibleESPs[index]; }
	
		string getName(){return name; }
		int setName(string N){ name = N; return 0;}
        int setType(string T){ type = T; return 0; }
	
		int setTheRest(string R){ theRest = R; return 0; }
		string getTheRest(){ return theRest;}
	
		int checkESPAssignment(int ESP){
			int possible = findLocationLocal(possibleESPs,possibleAssigned,ESP);
			int current = findLocationLocal(currentESPs,currentAssigned,ESP);
			if(possible<0){ 
				cout << "This ESP is not even possible in variant " << name << "\n"; 
				cout << "ESP: " << ESP << "\n";
				for(int i = 0; i<possibleAssigned; i++){ cout << possibleESPs[i] << " "; } cout << endl;
				exit(-1);
			}
			return current;
		}
	
		int addESP(int ESP){ 
			possibleESPs.size();
			possibleESPs.push_back(ESP); 
			possibleAssigned++; 
			return 0;
		}
	
		int sortPossibleESPs(){
			sort( possibleESPs.begin(), possibleESPs.begin() + possibleAssigned );
			return 0;
		}
	
		int addESPMapping(int ESP){
			//Check if the ESP is in the proposed list:
			int loc = findLocationLocal(possibleESPs,possibleAssigned,ESP);
			if(loc>=0){
				int tmpLoc = findLocationLocal(currentESPs,currentAssigned,ESP);
				int trueLoc = -1*tmpLoc - 1; 
				if(trueLoc>=0){ currentESPs.insert(currentESPs.begin() + trueLoc,ESP); currentAssigned++;}
			}
			else{
				cout << "Fatal Problem: We are trying to add an ESP that can not support variant " << name << "!\n"; exit(-1);
			}
			return 0;
		}
	
		int rejectMove(){
			proposedESPs = currentESPs;
			proposedAssigned = currentAssigned;
			return 0;
		}
		
		int clearCurrentState(){
			currentAssigned = 0;	
			currentESPs.clear(); 
			return 0;
		}
	
		int matchProposedCurrent(){
			proposedESPs = currentESPs;
			proposedAssigned = currentAssigned;
			for(int i = 0; i<=possibleAssigned; i++){ ESPoccupancy.push_back(0); }
			ESPoccupancy[currentAssigned]++;
			return 0;
		}
	
		int setValues(int CL, int CR, int NCL, int NCR, double ML, double MR, int CODE){
			coordLeft = CL; coordRight = CR;
			numConcordantsLeft = NCL; numConcordantsRight = NCR;
			mapLeft = ML; mapRight = MR; code = CODE;
					
			length = CR-CL;
			return 0;
		}
	
		int proposeAddingESP(int ESP){
			//Note need to make more efficient:
			int loc = findLocationLocal(proposedESPs,proposedAssigned,ESP);
			int trueLoc = -1*loc - 1; //loc = -(first+1); -1*loc - 1 = -1*-1*(first+1) - 1 = first+1 - 1 = first;
			if(trueLoc>=0){ proposedESPs.insert(proposedESPs.begin() + trueLoc,ESP); proposedAssigned++;}
			else{ cerr << "Error in adding ESP in proposal!\n"; exit(-1); }
			return proposedAssigned; 
		}
	
		int proposeRemovingESP(int ESP){
			int loc = findLocationLocal(proposedESPs,proposedAssigned,ESP);
			if(loc>=0){ proposedESPs.erase(proposedESPs.begin() + loc); proposedAssigned--; return proposedAssigned;}
			else{ 
				cerr << "Error in removing ESP in proposal!\n"; 
				cerr << "Trying to remove esp " << ESP << " from the following list\n";
				for(int i = 0; i<proposedAssigned; i++){
					cerr << i << "\t" << proposedESPs[i] << endl;
				}
				cerr << "*********\n";
				exit(-1); 
			}
		}

	
		double getLikelihoodErrorGiven(double perr, double COVERAGE, double COVERAGE_SCALED, int LAVG, int LDIS, int DISCORDANT){
			int numDiscordants = DISCORDANT;
			double probDisError = log(perr)*numDiscordants;
			
			double probConHomo = 0.0;
			int validCoverageFlag = 1;
			//beRD Model:
			//Do not scale, consider left and right endpoints separately;
			if(code == 0){ 
				probConHomo  = probVariant(numConcordantsLeft,LAVG,COVERAGE,2) + probVariant(numConcordantsRight,LAVG,COVERAGE,2);
				validCoverageFlag = validCoverage(COVERAGE*LAVG*2,numConcordantsLeft+numConcordantsRight,0.000100);
			}
			//Interval Model:
			//Scale concordant coverage appropriately;
			else if(code == 1){ 
				double scaledCovLeft  = (numConcordantsLeft*1.0)/(0.3+0.7*mapLeft);
				probConHomo  = probVariant(scaledCovLeft,length+LAVG,COVERAGE,2); 
				validCoverageFlag = validCoverage(COVERAGE*(length+LAVG),scaledCovLeft,0.000100);
				//return 0; 
			}
			else{
				cerr << "Invalid Variant Code: Check variant codes from GASVCC Output\n";
				exit(-1);
			}
			
			//We have exceeded the valid coverage; this MUST be an error;
			if(validCoverageFlag == 0){ return 0; }
			else{ return probDisError + probConHomo; }
			
			//int validC = validCoverage(Lambda*Lavg*2,scaledCovLeft+scaledCovRight,Tolerance);
			//If the valid coverage is TOO HIGH; then this must be an error; prob 0.
			
			
		}
	
		double getLikelihoodVariantGiven(double perr, double COVERAGE, double COVERAGE_SCALED, int LAVG, int LDIS, int DISCORDANT){
			int numDiscordants = DISCORDANT;
            
            int LDIS_LOCAL = LDIS;
            if(type.compare("IR")==0 || type.compare("TR+")==0 || type.compare("TR-") == 0){ LDIS_LOCAL = 2*LDIS; }
            double probDisHomo   = probVariant(numDiscordants,LDIS_LOCAL,COVERAGE,2);
			double probDisHetero = probVariant(numDiscordants,LDIS_LOCAL,COVERAGE,1);
			
			double probConError,probConHetero;
			probConError = probConHetero = 0.0;
			
			//REMOVE CONCORDANT COVERAGE!
			double numConcordantsOrig = -1;
			double numConcordants = -1;
			
			int validCoverageFlag = 1;
						
			//beRD Model:
			//Do not scale, consider left and right endpoints separately;
			if(code == 0){ 
				probConError   = log(perr)*numConcordantsLeft + log(perr)*numConcordantsRight;
				probConHetero  = probVariant(numConcordantsLeft,LAVG,COVERAGE,1) + probVariant(numConcordantsRight,LAVG,COVERAGE,1);
				validCoverageFlag = validCoverage(COVERAGE*LAVG*2,numConcordantsLeft+numConcordantsRight,0.000100);
				//return 0; 
			}
			//Interval Model:
			//Scale concordant coverage appropriately;
			else if(code == 1){ 
				double scaledCovLeft  = (numConcordantsLeft*1.0)/(0.3+0.7*mapLeft);
				numConcordantsOrig = numConcordantsLeft;
				numConcordants = scaledCovLeft;
				probConError   = log(perr)*scaledCovLeft;
				probConHetero  = probVariant(scaledCovLeft,length+LAVG,COVERAGE,1); 
				validCoverageFlag = validCoverage(COVERAGE*(length+LAVG),scaledCovLeft,0.000100);

				//return 0; 
			}
			else{
				cerr << "Invalid Variant Code: Check variant codes from GASVCC Output\n";
				exit(-1);
			}
			
			double pHomo   = probDisHomo + probConError;
			double pHetero = probDisHetero + probConHetero;
			
			/*cout << "CURRENT--> Name:\t" << name << "\tCode:\t" << code << "\t" 
			 << "ProbDiscordant:\t" << probDisHomo << "\tProbConcordantErrors:\t " << probConError 
			 << "\tNumConcordants Orig:\t" << numConcordantsOrig
			 << "\tNumConcordants Scaled:\t" << numConcordants << endl;
			 */
			
			if(validCoverageFlag == 0){
				double errorState = -2*200000000*1.0;
				return errorState;
			}
			else{
				if(pHomo>pHetero){ return pHomo;}
				else{ return pHetero; }
			}
			
		}
	
		//Note: COVERAGE = Lambda Parameter; COVERAGE_SCALED = Lambda_D, Lavg = Average Fragment Length; Ldis = Discordant Length
		double probCurrentlyAssigned(double perr, double COVERAGE, double COVERAGE_SCALED, int LAVG, int LDIS){
			//Gamma(n) = (n-1)!
			int numDiscordants = currentAssigned;
            
            int LDIS_LOCAL = LDIS;
            if(type.compare("IR")==0 || type.compare("TR+")==0 || type.compare("TR-") == 0){ LDIS_LOCAL = 2*LDIS; }
            double probDisHomo   = probVariant(numDiscordants,LDIS_LOCAL,COVERAGE,2);
			double probDisHetero = probVariant(numDiscordants,LDIS_LOCAL,COVERAGE,1);

			
			double probConError,probConHetero;
			probConError = probConHetero = 0.0;
			
			//REMOVE CONCORDANT COVERAGE!
			double numConcordantsOrig = -1;
			double numConcordants = -1;
			
			int validCoverageFlag = 1;

			
			//beRD Model:
			//Do not scale, consider left and right endpoints separately;
			if(code == 0){ 
				probConError   = log(perr)*numConcordantsLeft + log(perr)*numConcordantsRight;
				probConHetero  = probVariant(numConcordantsLeft,LAVG,COVERAGE,1) + probVariant(numConcordantsRight,LAVG,COVERAGE,1);
				validCoverageFlag = validCoverage(COVERAGE*LAVG*2,numConcordantsLeft+numConcordantsRight,0.000100);
				//return 0; 
			}
			//Interval Model:
			//Scale concordant coverage appropriately;
			else if(code == 1){ 
				double scaledCovLeft  = (numConcordantsLeft*1.0)/(0.3+0.7*mapLeft);
				numConcordantsOrig = numConcordantsLeft;
				numConcordants = scaledCovLeft;
				probConError   = log(perr)*scaledCovLeft;
				probConHetero  = probVariant(scaledCovLeft,length+LAVG,COVERAGE,1); 
				validCoverageFlag = validCoverage(COVERAGE*(length+LAVG),scaledCovLeft,0.000100);

				//return 0; 
			}
			else{
				cerr << "Invalid Variant Code: Check variant codes from GASVCC Output\n";
				exit(-1);
			}
			
			double pHomo   = probDisHomo + probConError;
			double pHetero = probDisHetero + probConHetero;
						
			/*cout << "CURRENT--> Name:\t" << name << "\tCode:\t" << code << "\t" 
			<< "ProbDiscordant:\t" << probDisHomo << "\tProbConcordantErrors:\t " << probConError 
			<< "\tNumConcordants Orig:\t" << numConcordantsOrig
			<< "\tNumConcordants Scaled:\t" << numConcordants << endl;
			*/
			if(validCoverageFlag == 0){ 
				double errorState = -2*200000000*1.0;
				return errorState;
			}
			else{
				if(pHomo>pHetero){ return pHomo;}
				else{ return pHetero; }
				//return pHomo;
			}
		}
	
		double probProposedAssigned(double perr, double COVERAGE, double COVERAGE_SCALED, int LAVG, int LDIS){ 
			//Gamma(n) = (n-1)!
			int numDiscordants = proposedAssigned;
            
            int LDIS_LOCAL = LDIS;
            if(type.compare("IR")==0 || type.compare("TR+")==0 || type.compare("TR-") == 0){ LDIS_LOCAL = 2*LDIS; }
            double probDisHomo   = probVariant(numDiscordants,LDIS_LOCAL,COVERAGE,2);
			double probDisHetero = probVariant(numDiscordants,LDIS_LOCAL,COVERAGE,1);
			
			double probConError,probConHetero;
			probConError = probConHetero = 0.0;
			//REMOVE CONCORDANT COVERAGE!
			
			double numConcordantsOrig = -1;
			double numConcordants = -1;
			
			int validCoverageFlag = 1;
			
			//beRD Model:
			//Do not scale, consider left and right endpoints separately;
			if(code == 0){ 
				probConError = log(perr)*numConcordantsLeft + log(perr)*numConcordantsRight;
				probConHetero  = probVariant(numConcordantsLeft,LAVG,COVERAGE,1) + probVariant(numConcordantsRight,LAVG,COVERAGE,1);
				validCoverageFlag = validCoverage(COVERAGE*LAVG*2,numConcordantsLeft+numConcordantsRight,0.000100);
			}
			//Interval Model:
			//Scale concordant coverage appropriately;
			else if(code == 1){ 
				double scaledCovLeft  = (numConcordantsLeft*1.0)/(0.3+0.7*mapLeft);
				numConcordantsOrig = numConcordantsLeft;
				numConcordants = scaledCovLeft;
				probConError   = log(perr)*scaledCovLeft;
				probConHetero  = probVariant(scaledCovLeft,length+LAVG,COVERAGE,1); 
				validCoverageFlag = validCoverage(COVERAGE*(length+LAVG),scaledCovLeft,0.000100);
			}
			else{
				cerr << "Invalid Variant Code: Check variant codes from GASVCC Output\n";
				exit(-1);
			}
			
			double pHomo   = probDisHomo + probConError;
			double pHetero = probDisHetero + probConHetero;
		
			if(validCoverageFlag == 0){ 
				double errorState = -2*200000000*1.0;
				return errorState;
			}
			//return pHomo;
			if(pHomo>pHetero){ return pHomo;}
			else{ return pHetero; }		
		}
	
		int getSupport(int numFrags){ 
			if(numFrags<=possibleAssigned){ return ESPoccupancy[numFrags]; }
			else{ cout << "The number of fragments " << numFrags << " is not possible for variant " << name << endl; exit(-1); return -1;}
		}
		
	int incrementLogLikelihood(double value){ runningLogLikelihood+= value; return 0; }
	double getRunningLogLikelihood(){ return runningLogLikelihood; }
	
	int incrementLikelihoodVar(double value){ runningLikelihoodVar += value; return 0; }
	double getRunningLikelihoodVar(){ return runningLikelihoodVar; }
	
		int updateOccupancy(){	
			ESPoccupancy[currentAssigned]++; 
			return 0; 
		}
		int updateOccupancy(int STEPS, double perr, double COVERAGE, double COVERAGE_SCALED, int LAVG, int LDIS){ 
			ESPoccupancy[currentAssigned]+=STEPS; 
			if(currentAssigned > 0){
				incrementLogLikelihood(STEPS*(getLikelihoodVariantGiven(perr,COVERAGE,COVERAGE_SCALED,LAVG,LDIS,currentAssigned) - getLikelihoodErrorGiven(perr,COVERAGE,COVERAGE_SCALED,LAVG,LDIS,currentAssigned))); 
				incrementLikelihoodVar(STEPS*(getLikelihoodVariantGiven(perr,COVERAGE,COVERAGE_SCALED,LAVG,LDIS,currentAssigned))); 

			}
			return 0; 
		
		}
	

		
		int clearOccupancy(){ 
			for(int i = 0; i<=possibleAssigned; i++){ ESPoccupancy[i] = 0; }
			runningLogLikelihood = runningLikelihoodVar= 0;
			return 0;
		}
	
		int acceptMove(){			
			int retVal = 0;
			
			//cout << "Before accepting variant move:\n";
			//cout << "\tCurrent Support = " << currentAssigned << endl;
			//cout << "\tProposed Support = " << proposedAssigned << endl;
		
			if(currentAssigned == 0 && proposedAssigned > 0){ retVal = 1; } //Move from 0 to supported;
			else if(currentAssigned > 0 && proposedAssigned == 0){ retVal = -1; } //Move from supported to 0;			
			currentAssigned = proposedAssigned;
			currentESPs = proposedESPs;
					
			//cout << "After accepting variant move:\n";
			//cout << "\tCurrent Support = " << currentAssigned << endl;
			//cout << "\tProposed Support = " << proposedAssigned << endl;
			
			return retVal;
		}
	
		friend ostream& operator<<(ostream& output, const variant& V) {
			output << "Name: " << V.name << "\n";
			output << "Possible Assigned ESPs: " << V.possibleAssigned << "\n";
			for(int i = 0; i<V.possibleAssigned; i++){
				output << V.possibleESPs[i] << " ";
			}
			output << endl;
			
			output << "Current Assigned --> " << V.currentAssigned << "\n";
			for(int i = 0; i<V.currentAssigned; i++){
				output << V.currentESPs[i] << " ";
			}
			output << endl;
			
			output << "Proposed Assigned --> " << V.proposedAssigned << "\n";
			for(int i = 0; i<V.proposedAssigned; i++){
				output << V.proposedESPs[i] << " ";
			}
			output << endl;
			
			return output;
		}
	
		int addDependency(int ESP){
			int valid = findLocationLocal(possibleESPs,possibleAssigned,ESP);
			if(valid<0){ cerr << "Can not add this ESP " << ESP << " to the variant " << name << "\n"; exit(-1); }
			int loc = findLocationLocal(possibleDependentESPs,numDependencies,ESP);
			if(loc<0){ 
				int trueLoc = -1*loc-1;  //-1*-1(first+1) = first + 1; so -1
				possibleDependentESPs.insert(possibleDependentESPs.begin()+trueLoc,ESP);
				numDependencies++;
			}
			return 0;
		}
	
		// We return a 1 if this ESP is currently assigned and 0 otherwise;
		int currentlyAssigned(int ESP){
			int loc = findLocationLocal(currentESPs,currentAssigned,ESP);
			if(loc < 0){ return 0; }
			else{ return 1; }
		}
	
		//We now define a variant as "swapable" if there is at least one ESP assigned that this variant
		// shares with another one.
		//Return 0 if we cannot swap and 1 otherwise;
		int canWeSwap(){
			int retVal = 0;
			if(numDependencies == 0){ return retVal; } //If no dependencies, then we can never swap!
			for(int i = 0; retVal == 0 && i<currentAssigned; i++){
				int loc = findLocationLocal(possibleDependentESPs,numDependencies,currentESPs[i]);
				if(loc>=0){ retVal++; }
			}
			return retVal;
		}
		
	
		int getTimeOccupied(){	
			int retVal = 0;
			for(int i = 1; i<=possibleAssigned; i++){ retVal+= ESPoccupancy[i]; }
			return retVal;
		}
	
	private:
		string name;
        string type;
		string theRest; //This is needed for outputting the resulting SVs

		std::vector<int> possibleDependentESPs; //The value is 0 if no dependencies and 1 otherwise.
		std::vector<int> possibleESPs;
		std::vector<int> currentESPs;
		std::vector<int> proposedESPs;
		std::vector<int> ESPoccupancy; //A tally of the number of times X ESPs have occupied this variant;
	
		double runningLogLikelihood;
		double runningLikelihoodVar;
	
		int currentAssigned;
		int proposedAssigned;
		int possibleAssigned;
		int numDependencies;
	
		int coordLeft;
		int coordRight;
		int numConcordantsLeft;
		int numConcordantsRight;
		double mapLeft;
		double mapRight;
		int code;
		int length;
	
};

class pairedDependency{
	
	public:
		pairedDependency(){ variant1 = variant2 = -1; numESPsInCommon = 0;};
		pairedDependency(int var1, int var2, int E){ variant1 = var1; variant2 = var2; ESPsInCommon.push_back(E); numESPsInCommon = 1; }
		~pairedDependency(){};
	
		friend bool operator<( const pairedDependency&, const pairedDependency &);
		friend bool operator>( const pairedDependency&, const pairedDependency &);

		int getVariant1(){ return variant1;}
		int getVariant2(){ return variant2;}
	
		int addDependency(int ESP){ ESPsInCommon.push_back(ESP); numESPsInCommon++; return 0;}
	
		int getNumESPsInCommon(){ return numESPsInCommon; }
		int getESPByIndex(int index){ return ESPsInCommon[index];}

	private:
		int variant1;
		int variant2;
		int numESPsInCommon;
		std::vector<int> ESPsInCommon;
};



// post: return true iff lhs < rhs
bool operator< ( const pairedDependency& As, const pairedDependency& Bs ) {
	if(As.variant1 == Bs.variant1){ 
		if(As.variant2 < Bs.variant2){ return true; }
		else{ return false; }
	}
	else{ 
		if(As.variant1 < Bs.variant1){ return true; }
		else{ return false; }
	}
}

// post: return true iff lhs > rhs
bool operator> (const pairedDependency& As, const pairedDependency& Bs){
		if(As.variant1 == Bs.variant1){ 
			if(As.variant2 > Bs.variant2){ return true; }
			else{ return false; }
		}
		else{ 
			if(As.variant1 > Bs.variant1){ return true; }
			else{ return false; }
		}
}

int randomInitialization(int numMobileESPs, vector<int> mobileESPs, vector<ESP> & ESPVector, vector<variant> &VariantVector, vector<int> &AssignedVariants, int &numAssigned, vector<int> &EmptyVariants, int &numEmpty, double perr){
	//Goal: For each ESP vector, with probability perr, assign ESP to an error; w/ probability (1-perr) assign to one of the
	//      random locations; Update the variants;
	
	//(1) Set a random assignment for each ESP;
	for(unsigned int i = 0; i<ESPVector.size(); i++){
		int numVariants = VariantVector.size();
		int chosenVariant = ESPVector[i].initializeRandomMapping(perr);
		
		//Fix if Needed:
		int movable = findLocationLocal(mobileESPs,numMobileESPs,i);
		//We can not move this; put it in its unique or location;
		if(movable < 0){ 
			//Set to unique mapping;
			if(ESPVector[i].getPossibleMappings() == 1){
				chosenVariant = ESPVector[i].setInitialMapping(0);
			}
			//Set to error;
			else{
				chosenVariant = ESPVector[i].setInitialMapping(-1);
			}
			
		}
		else{ /* Nothing, we already picked the move! */ }
	
		if(chosenVariant == -1){
			//cout << "ESP " << i << " is unassigned.\n";
		}
		else if(chosenVariant >= 0 && chosenVariant<numVariants){
			VariantVector[chosenVariant].addESPMapping(i);
		}
		else{
			//cerr << "Possible error with MCMC initialization. Investigate.\n";
		}
	}
	
	//(2) Determine all assigned and empty variants;

	std::vector<int>::iterator front;
	for(unsigned int i = 0; i<VariantVector.size(); i++){
		if(VariantVector[i].getCurrentAssigned() > 0){
			if(VariantVector[i].checkMobility(numMobileESPs,mobileESPs) == 1){
				front = AssignedVariants.begin(); 
				AssignedVariants.insert(front,i); 
				numAssigned++; 
			}
		}
		else{ 
			//If empty, then we KNOW that each possible ESP is a mobile one!
			front = EmptyVariants.begin(); 
			EmptyVariants.insert(front,i); 
			numEmpty++; 
		}
	}
	sort( AssignedVariants.begin(), AssignedVariants.begin() + numAssigned );
	sort( EmptyVariants.begin(), EmptyVariants.begin() + numEmpty); 
		

	//(3) For all variants and ESPs set the proposed move to the current move:
	for(unsigned int i  = 0; i<ESPVector.size(); i++){ ESPVector[i].matchProposedCurrent(); }
	for(unsigned int i = 0; i<VariantVector.size(); i++){ VariantVector[i].matchProposedCurrent(); }
	
	return 0;
}

int proposeNaiveMove(int numMobileESPs, vector<int> mobileESPs, vector<ESP> & ESPVector, vector<variant> &VariantVector, vector<int> &AssignedVariants, int &numAssigned, vector<int> &EmptyVariants, int &numEmpty, vector<int> &modifiedESPs, int &numModifiedESPs, vector<int> &modifiedVariants, int &numModifiedVariants, double perr){
	//Goal: Pick an ESP at random and change its assignment. NOTE: Must actually change the assignment;
	
	//Step 1: Choose a random ESP to modify;
	int randomESP =  rand() % (ESPVector.size());
	//cout << "Choosing to move " << randomESP <<"\n";
	
	//Fix as random!
	int randomIndex = rand() % (numMobileESPs);
	randomESP = mobileESPs[randomIndex];
	
	modifiedESPs.push_back(randomESP);
	numModifiedESPs++;
	
	int i = randomESP;
	ESPVector[i].proposeRandomMove(perr);
	while( ESPVector[i].getProposedState() == ESPVector[i].getCurrentState() ){ 
		ESPVector[i].proposeRandomMove(perr); 
	}
	
	int currentVariant = ESPVector[i].getVariant(ESPVector[i].getCurrentState());
	int proposedVariant = ESPVector[i].getVariant(ESPVector[i].getProposedState());

/*	
	cout << "Trying to move ESP " << i << endl;
	cout << "\tCurrent State:\t" << ESPVector[i].getCurrentState() << endl;
	cout << "\tCurrent Variant:\t" << currentVariant << endl;
	cout << "\tProposed State:\t" << ESPVector[i].getProposedState() << endl;
	cout << "\tProposed Variant:\t" << proposedVariant << endl;
*/
	
	if(currentVariant == -1){ 
		//cout << "\tMoved from error\n";
	} //Do nothing, the ESP was an error;
	else{
		VariantVector[currentVariant].proposeRemovingESP(i);
		modifiedVariants.push_back(currentVariant); numModifiedVariants++;
	}
	if(proposedVariant == -1){ 
		//cout << "\tMoved to error.\n";
	} //Do nothing, the ESP was moved to an error;
	else{
		VariantVector[proposedVariant].proposeAddingESP(i);
		modifiedVariants.push_back(proposedVariant); numModifiedVariants++;
	}
	
//	cout << "Finished proposing the move!\n" << flush;
	
	return 0;
}

int proposeAddMove(int numMobileESPs, vector<int> &mobileESPs, vector<ESP> & ESPVector, vector<variant> &VariantVector, vector<int> &AssignedVariants, int &numAssigned, vector<int> &EmptyVariants, int &numEmpty, vector<int> &modifiedESPs, int &numModifiedESPs, vector<int> &modifiedVariants, int &numModifiedVariants, double perr){

	//(1) Select a variant to add proportional to the number of empty variants;
	int randomEmptyVariantIndex =  rand() % (numEmpty);
	
	int A = EmptyVariants[randomEmptyVariantIndex];
	modifiedVariants.push_back(A); numModifiedVariants++;

	/*cout << "Trying to add variant " << A << endl;
	cout << VariantVector[A] << endl;
	cout << "***************\n";
	 */
	
	//(2) Move ESPs to this variant with probability 1/2, guarantee at least 1 ESP is moved.
	while(numModifiedESPs == 0){
		//If we go all the way through and don't change anything; we just repeat. 
		for(int i = 0; i<VariantVector[A].getPossibleAssigned(); i++){
			int ESPtoChange = VariantVector[A].getPossibleESP(i);
			double prob = ((double) rand()*1.0 / (RAND_MAX+1.0));
			//With Prob 1/2 add this to A:
			if(prob<=0.50){
				//cout << "\tAdding ESP " << ESPtoChange << endl;
				modifiedESPs.push_back(ESPtoChange); numModifiedESPs++;
				int otherVariant = ESPVector[ESPtoChange].proposeAddMove(A);
				//cout << "\tRemoving ESP from variant " << otherVariant << endl;
				if(otherVariant >=0){
					//(1) Remove ESP from this Variant:
					VariantVector[otherVariant].proposeRemovingESP(ESPtoChange);
					
					//(2) Indicate this variant also needs to be processed
					int loc = findLocationLocal(modifiedVariants,numModifiedVariants,otherVariant);
					//Do we need to add this to the modified variant list?
					//No: We already have seen it;
					if(loc >= 0){ } 			
					//Yes: We haven't seen it yet;
					else{ 
						int trueLoc = -1*loc - 1;
						modifiedVariants.insert(modifiedVariants.begin() + trueLoc, otherVariant); 
						numModifiedVariants++;
					}
				}
			}
		}
	}
	
	//cout << "Propose Adding variant " << A << endl;
	//cout << "We have proposed adding " << numModifiedESPs << " to variant " << A << endl;

	//Now propose adding these ESPs to the variant;
	for(int i = 0; i<numModifiedESPs; i++){
		VariantVector[A].proposeAddingESP(modifiedESPs[i]);
	}
	
		
	return 0;
}

int proposeRemoveMove(int numMobileESPs, vector<int> &mobileESPs, vector<ESP> & ESPVector, vector<variant> &VariantVector, vector<int> &AssignedVariants, int &numAssigned, vector<int> &EmptyVariants, int &numEmpty, vector<int> &modifiedESPs, int &numModifiedESPs, vector<int> &modifiedVariants, int &numModifiedVariants, double perr){

	//(1) Select a variant to add proportional to the number of empty variants;
	int randomAssignedVariantIndex =  rand() % (numAssigned);
	
	int A = AssignedVariants[randomAssignedVariantIndex];
	modifiedVariants.push_back(A); numModifiedVariants++;
	
	/*
	cout << "Trying to remove variant " << A << " out of " << numAssigned << endl;
	cout << VariantVector[A] << endl;
	cout << "Currently Assigned:\t" << endl;
	for(int i = 0; i<VariantVector[A].getCurrentAssigned(); i++){
		cout << i << " = " << VariantVector[A].getCurrentESP(i) << endl;
	}
	cout << "***************\n";
	*/
	
	if(numModifiedESPs >0){ cout << "Problem here, right?" << endl; }
	
	numModifiedESPs = 0;
	//(2) Remove ALL ESPs currently assigned!
	for(int i = 0; i<VariantVector[A].getCurrentAssigned(); i++){
		int ESPtoChange = VariantVector[A].getCurrentESP(i);
		//cout << "\tAdding ESP " << ESPtoChange << endl;
		modifiedESPs.push_back(ESPtoChange); numModifiedESPs++;
		int otherVariant = ESPVector[ESPtoChange].proposeRemoveMove(A,perr);
		if(otherVariant >=0){
			//(1) Add ESP from this Variant:
			//cout << "\tMoving ESP " << ESPtoChange << " from variant (" << A << "): " << VariantVector[A].getName() << " to (" << otherVariant<< "):" << VariantVector[otherVariant].getName() << endl;
			VariantVector[otherVariant].proposeAddingESP(ESPtoChange);
					
			//(2) Indicate this variant also needs to be processed
			int loc = findLocationLocal(modifiedVariants,numModifiedVariants,otherVariant);
			//Do we need to add this to the modified variant list?
			//No: We already have seen it;
			if(loc >= 0){ } 			
			//Yes: We haven't seen it yet;
			else{ 
				int trueLoc = -1*loc - 1; //loc = -(first+1) + 1
				//cout << "\tAdding variant " << otherVariant << " at offset " << trueLoc << endl;
				modifiedVariants.insert(modifiedVariants.begin() + trueLoc, otherVariant); 
				numModifiedVariants++;
			}
		}
		else{
			//cout << "\tMoving ESP " << ESPtoChange << " from variant (A): " << VariantVector[A].getName() << " to error" <<  endl;			
		}
		//cout << "\t\tDone processing ESP " << ESPtoChange << endl;
	}
	
	//cout << "In total we are adding esps to the following variants:\n";
	//for(int i = 0; i<numModifiedVariants; i++){
	//	cout << i << " " << modifiedVariants[i] << endl;
	//}
	
	//cout << "Proposing Removing variant " << A << endl;
	
	//Now propose adding these ESPs to the variant;
	for(int i = 0; i<numModifiedESPs; i++){
		VariantVector[A].proposeRemovingESP(modifiedESPs[i]);
	}
	
	return 0;
}

int isSwapPossible(vector<int> &assignedVariants, int numAssigned, vector<variant> &VariantVector){
	int retVal = 0;
	int v = 0;
	while(retVal == 0 and v<numAssigned){
		if(VariantVector[assignedVariants[v]].canWeSwap()>0){ retVal++; }
		v++;
	}
	return retVal;
}

int proposeSwapMove(int numMobileESPs, vector<int> &mobileESPs, vector<ESP> &ESPVector, vector<variant> &VariantVector, vector<pairedDependency> &Dependencies, vector<int> &runningDependency, int TotalDependency, vector<int> &visitedDependencies, int numDependencies, int step, vector<int> &AssignedVariants, int &numAssigned, vector<int> &EmptyVariants, int &numEmpty, vector<int> &modifiedESPs, int &numModifiedESPs, vector<int> &modifiedVariants, int &numModifiedVariants, double perr){
	
	//0: Select a valid swap move:
	int variant1 = -1;
	int variant2 = -1;
	int swapMove = -1;
	while(variant1 == -1 && variant2 == -1){
		//Pick a random index out of all dependencies;
		int randomSwapValue =  rand() % (TotalDependency);
		//Find the swapable location we should visit;
		 swapMove = findLocationLocal(runningDependency,numDependencies,randomSwapValue);
		if(swapMove < 0){ swapMove = -1*(swapMove + 1); }
		
		//If we haven't attempted this before; consider this swap.
		if(visitedDependencies[swapMove] != step){
			//Mark as visited;
			visitedDependencies[swapMove] = step;
			
			//Which variants are involved?
			variant1 = Dependencies[swapMove].getVariant1();
			variant2 = Dependencies[swapMove].getVariant2();
			
			//Is this swap move even valid?
			int numOccupied = 0;
			for(int i = 0; i<Dependencies[swapMove].getNumESPsInCommon() && numOccupied == 0 ; i++){
				int ESP = Dependencies[swapMove].getESPByIndex(i);
				numOccupied += VariantVector[variant1].currentlyAssigned(ESP) + VariantVector[variant2].currentlyAssigned(ESP);
			}
			
			if(numOccupied == 0){ variant1 = variant2 = -1; }
		}
	}
	
	cout << "Chose a valid swap move --> " << swapMove << " between " << variant1 << " and " << variant2 << endl;
	
	//Note: We insert from the front and know that variant1 < variant2;
	modifiedVariants.insert(modifiedVariants.begin(), variant2); 
	numModifiedVariants++;
	modifiedVariants.insert(modifiedVariants.begin(), variant1); 
	numModifiedVariants++;
	
	
	//1: Propose the swap move by; going through all indices of variant1 and variant2 and moving them to the other with probability (1/2):
	while(numModifiedESPs == 0){ 
		cout << "Processing variant " << variant1 << " with " << VariantVector[variant1].getCurrentAssigned() << endl;
		
		for(int i = 0; i<VariantVector[variant1].getCurrentAssigned(); i++){
			int ESP = VariantVector[variant1].getCurrentESP(i); 
			//Q: Could this be moved to variant2?
			int possibleToSwap = ESPVector[ESP].canESPSupport(variant2);
			if(possibleToSwap>0){
				double prob = ((double) rand()*1.0 / (RAND_MAX+1.0));
				if(prob<=0.5){ 
					ESPVector[ESP].proposeSwapMove(variant2); 
					modifiedESPs.insert(modifiedESPs.begin(), ESP); 
					numModifiedESPs++;
					
					//Remove ESP from variant1 add to variant2;
					VariantVector[variant1].proposeRemovingESP(ESP);
					VariantVector[variant2].proposeAddingESP(ESP);
				}
			}
		}
		
		cout << "Finished going through " << variant1 << " and moved " << numModifiedESPs << "\n";
		
		cout << "Next processing variant " << variant2 << " with " << VariantVector[variant2].getCurrentAssigned() << endl;
		
		for(int i = 0; i<VariantVector[variant2].getCurrentAssigned(); i++){
			int ESP = VariantVector[variant2].getCurrentESP(i); 
			//Q: Could this be moved to variant2?
			int possibleToSwap = ESPVector[ESP].canESPSupport(variant1);
			if(possibleToSwap>0){
				double prob = ((double) rand()*1.0 / (RAND_MAX+1.0));
				if(prob<=0.5){ 
					ESPVector[ESP].proposeSwapMove(variant1); 
					modifiedESPs.insert(modifiedESPs.begin(), ESP); 
					numModifiedESPs++;
					
					VariantVector[variant2].proposeRemovingESP(ESP);
					VariantVector[variant1].proposeAddingESP(ESP);
				}
			}
		}
		
		cout << "Finished going through " << variant2 << " and now have " << numModifiedESPs << "\n";
		
	}

	
	return 0;
}


//(n choose k) = n!/(k!(n-k)!)
//gamma(n) = (n-1)!
double logNChooseK(int n, int k){
	double retVal = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
	return retVal;
}


int acceptMove(int numMobileESPs, vector<int> &mobilesESPs, int moveToMake, double& currentLikelihood, vector<ESP> & ESPVector, vector<variant> &VariantVector, vector<int> &AssignedVariants, int &numAssigned, vector<int> &EmptyVariants, int &numEmpty, vector<int> &modifiedESPs, int &numModifiedESPs, vector<int> &modifiedVariants, int &numModifiedVariants, double perr, double COVERAGE,double COVERAGE_SCALED,int LAVG,int LDIS){

//	cout << "Trying to accept move.\n" << flush;
	
	int retVal = -1;
	
	//Note: moveToMake --> 0 = selfedge; 1 = Naive; 2 = Add; 3 = Remove;
	
	//The issue: For each move type you need to figure out if it could have been made by a different move:
	
	//Compute the proposal distribution --> Our chosen move is A to Aprime;
	
	//Is this an add?
	// --> All ESPs are moved to the same variant; had 0 support originally;
	
	//Is this a remove?
	// --> All ESPs are removed from the same variant; now has 0 support.
	
	//If there is only 1 ESP moved, this has to be a naive move.
	
	
	
	double proposeAtoAPrime = 0;
	double proposeAPrimetoA = 0;
	
	int considerAddMove = 0;
	double proposeAtoAPrimeAdd = 0;
	double proposeAPrimetoAAdd = 0;

	int considerRemoveMove = 0;
	double proposeAtoAPrimeRemove = 0;
	double proposeAPrimetoARemove = 0;
	
	//Claim 1:
	//   If A --> A' by an add/remove move, then A' --> A by a remove/add move.
	// YES: Add take from anywhere to a single variant; Remove take from a single variant to anywhere. 
	
	//Check for Add move:
	vector <int> currentVariants;
	vector <int> proposedVariants;
	int numCurrentErrors = 0;
	int numProposedErrors = 0;
	int numCurrentVariants = 0;
	int numProposedVariants = 0;
	for(int i = 0; i<numModifiedESPs; i++){
		int ESPInQuestion = modifiedESPs[i];
		int currentState = ESPVector[ESPInQuestion].getCurrentState();
		int proposedState = ESPVector[ESPInQuestion].getProposedState();
				
		if(currentState == -1){ numCurrentErrors++; }
		else{
			int currentVariant = ESPVector[ESPInQuestion].getVariant(currentState);
			int loc = findLocationLocal(currentVariants,numCurrentVariants,currentVariant);
			if(loc < 0){ 
				int trueLoc = -1*loc -1;
				currentVariants.insert(currentVariants.begin()+trueLoc,currentVariant);
				numCurrentVariants++;
			}
		}
		
		if(proposedState == -1){ numProposedErrors++; }
		else{
			int proposedVariant = ESPVector[ESPInQuestion].getVariant(proposedState);
			int loc = findLocationLocal(proposedVariants,numProposedVariants,proposedVariant);
			if(loc < 0){ 
				int trueLoc = -1*loc -1;
				proposedVariants.insert(proposedVariants.begin()+trueLoc,proposedVariant);
				numProposedVariants++;
			}
		}
	}
		
	//Remove Move: (Only moved ESPs from 1 variant; none were set to errors)
	if(numCurrentErrors == 0 && numCurrentVariants == 1){
		int v = currentVariants[0];
		//cout << "Trying for remove\n";
		if(VariantVector[v].getCurrentAssigned() > 0 && VariantVector[v].getProposedAssigned() == 0){
			considerRemoveMove = 1;
			//if(moveToMake == 3){ cout << "CORRECT Remove correctly!\n"; }			
		}
	}
	
	//Add Move: (Do not propose any errors; moving everything to 1 variant)
	if(numProposedErrors == 0 && numProposedVariants == 1){
		int v = proposedVariants[0];
		if(VariantVector[v].getCurrentAssigned() == 0 && VariantVector[v].getProposedAssigned() > 0){
			considerAddMove = 1;
			//if(moveToMake == 2){ cout << "CORRECT Add correctly!\n"; }
		}
	}
	
	if(moveToMake == 2 && considerAddMove == 0){ 
		int v = proposedVariants[0];
		cerr << "We missed an add move flag!\n"; 
		cerr << "Variant " << v << endl;
		cout << VariantVector[v] << endl;
		exit(-1); 
	}
	if(moveToMake == 3 && considerRemoveMove == 0){ 
		int v = currentVariants[0];
		cerr << "We missed a remove move flag!\n"; 
		cerr << "Variant " << v << endl;
		cout << VariantVector[v] << endl;
		exit(-1); 
	}
	
	//From A to A' we add a variant; and from A' to A we remove a variant;
	if(considerAddMove == 1){
		int v = proposedVariants[0];
		
		int numFullVariantsAPrime = 0;
		for(int i = 0; i<VariantVector.size(); i++){ if(VariantVector[v].getProposedAssigned()>0){ numFullVariantsAPrime++; } }

		//Picking the right variant.
		proposeAtoAPrimeAdd = log(1.0) - log(numEmpty);
		proposeAPrimetoAAdd = log(1.0) - log(numFullVariantsAPrime); //Guaranteed this is non-zero;
		
		//For add move; we need to look at the probability of picking this many.
		//We had to choose this move! (<-- Only 1 ESP that could be moved, we know exactly where it went!)
		if(VariantVector[v].getPossibleAssigned()==1){
			proposeAtoAPrimeAdd = log(1.0); //We had to add this particular case here.
		}
		else{
			//(1) Probability of selecting this variant: Prob(k successes)/Prob(k>0)
			//(2) Need to divide/subtract by the probability of choosing AT LEAST 1; --> (1 - (1/2)^possible)
			proposeAtoAPrimeAdd += logNChooseK(VariantVector[v].getPossibleAssigned(),numModifiedESPs)-log(2.0)*numModifiedESPs - log(1 - pow(0.5,VariantVector[v].getPossibleAssigned()));
		}
		
		//For each ESP to move, we need to consider each ESP;
		for(int i = 0; i<numModifiedESPs; i++){
			int E = modifiedESPs[i];
			int currentMapping  = ESPVector[E].getCurrentState();
			int proposedMapping = ESPVector[E].getProposedState();
			int numMappings     = ESPVector[E].getPossibleMappings();
				
			//Then in the remove move we set this to an error automatically.;
			if(numMappings == 1){
				proposeAPrimetoAAdd += log(1.0);
			}
			//More than 1 mapping
			else{
				//We are removing this from 1 variant to another; so have to pick another mapping;
				if(currentMapping >= 0){
					proposeAPrimetoAAdd += log(1.0/(numMappings*1.0-1.0)) + log(1-perr);
				}
				//We are moving this from 1 variant to an error; to had to pick an error;
				else if(currentMapping == -1){
					proposeAPrimetoAAdd += log(perr);
				}
			}
		}
	}
	
	
	if(considerRemoveMove == 1){
		int v = currentVariants[0];
		
		int numEmptyVariantsAPrime = 0;
		for(int i = 0; i<VariantVector.size(); i++){ if(VariantVector[v].getCurrentAssigned()>0){ numEmptyVariantsAPrime++; } }

		//Picking the right variant
		proposeAtoAPrimeRemove = log(1) - log(numAssigned);
		proposeAPrimetoARemove = log(1) - log(numEmptyVariantsAPrime); //Guaranteed this is non-zero;

		//Ok, when we remove things we need to consider each place we could have moved them too;
		//For each ESP to move, we need to consider each ESP;
		for(int i = 0; i<numModifiedESPs; i++){
			int E = modifiedESPs[i];
			int currentMapping  = ESPVector[E].getCurrentState();
			int proposedMapping = ESPVector[E].getProposedState();
			int numMappings     = ESPVector[E].getPossibleMappings();
			
			//Then in the remove move we set this to an error automatically.;
			if(numMappings == 1){
				proposeAtoAPrimeRemove += log(1.0);
			}
			//More than 1 mapping
			else{
				//We are removing this from 1 variant to another; so have to pick another mapping;
				if(proposedMapping >= 0){
					proposeAtoAPrimeRemove += log(1.0/(numMappings*1.0-1.0)) + log(1-perr);
				}
				//We are moving this from 1 variant to an error; to had to pick an error;
				else if(proposedMapping == -1){
					proposeAtoAPrimeRemove += log(perr);
				}
			}
		}
		
		//Now we consider the reverse, the add move;
		if(VariantVector[v].getPossibleAssigned()==1){
			proposeAPrimetoARemove = log(1.0); //We had to add this particular case here.
		}
		else{
			//(1) Probability of selecting this variant: Prob(k successes)/Prob(k>0)
			//(2) Need to divide/subtract by the probability of choosing AT LEAST 1; --> (1 - (1/2)^possible)
			proposeAPrimetoARemove += logNChooseK(VariantVector[v].getPossibleAssigned(),numModifiedESPs)-log(2.0)*numModifiedESPs - log(1 - pow(0.5,VariantVector[v].getPossibleAssigned()));
		}
	}
	
	
	int considerNaiveMove = 0;
	double proposeAtoAPrimeNaive = 0;
	double proposeAPrimetoANaive = 0;
	
	if(numModifiedESPs == 1){ 
		considerNaiveMove = 1;
		int ESPInQuestion = modifiedESPs[0];
		int numMappings = ESPVector[ESPInQuestion].getPossibleMappings();
		if(numMappings == 1){
			proposeAtoAPrimeNaive = log(1) - log(numMobileESPs);
			proposeAPrimetoANaive = log(1) - log(numMobileESPs);
		}
		else{
			int currentMapping = ESPVector[ESPInQuestion].getCurrentState();
			int proposedMapping = ESPVector[ESPInQuestion].getProposedState();
			
			if(currentMapping >= 0 && proposedMapping >= 0){ 
				proposeAtoAPrimeNaive = log(1) - log(numMobileESPs) + log(1.0/(numMappings*1.0-1.0)) + log(1-perr);
				proposeAPrimetoANaive = log(1) - log(numMobileESPs) + log(1.0/(numMappings*1.0-1.0)) + log(1-perr);
			}
			else if(currentMapping == -1 && proposedMapping >= 0){
				proposeAtoAPrimeNaive = log(1) - log(numMobileESPs) + log(1.0/(numMappings*1.0));
				proposeAPrimetoANaive = log(1) - log(numMobileESPs) + log(perr);
			}
			else if(currentMapping >= 0 && proposedMapping == -1){
				proposeAtoAPrimeNaive = log(1) - log(numMobileESPs) + log(perr);
				proposeAPrimetoANaive = log(1) - log(numMobileESPs) + log(1.0/(numMappings*1.0));
			}
			else{
				cerr << "Accept Move --> Should NEVER be here!\n"; exit(-1);
			}
		}
	}
	//else if(numModifiedESPs > 1){
	//	cerr << "Only stable for naive only moves. Remove this halt to do a more general chain.\n"; exit(-1);
	//}
	
	double probMoveCurrent = 0;
	double probMoveProposed = 0;
	
	for(unsigned int i = 0; i<ESPVector.size(); i++){
		
		probMoveCurrent  += log(perr)*ESPVector[i].currentlyUnAssigned();
		probMoveProposed += log(perr)*ESPVector[i].proposedUnAssigned();
		
		/*
		if(ESPVector[i].currentlyUnAssigned() == 1 || ESPVector[i].currentlyUnAssigned() == 1){ 
			cout << "Unassigned move considered!\n"; 
			cout << "Current --> " << probMoveCurrent << endl;
			cout << "Proposed --> " << probMoveProposed << endl;
			exit(-1); 
		}
		*/
		
	}
		
	for(unsigned int i = 0; i<VariantVector.size(); i++){
		if(VariantVector[i].getCurrentAssigned() > 0){
			probMoveCurrent  += VariantVector[i].probCurrentlyAssigned(perr,COVERAGE,COVERAGE_SCALED,LAVG,LDIS);
		}
		if(VariantVector[i].getProposedAssigned() > 0){
			probMoveProposed += VariantVector[i].probProposedAssigned(perr,COVERAGE,COVERAGE_SCALED,LAVG,LDIS); 
		}
	}
		
	//Proposal Distribution: 
	// Numerator:   --> log(P(A'|D) q(A|A')) = probMoveProposed + proposeAPrimetoA;
	// Denominator: --> log(P(A |D) q(A'|A)) = probMoveCurrent + proposeAtoAPrime;
	
	if(considerNaiveMove == 1){
		proposeAPrimetoA = log(exp(proposeAPrimetoA) + exp(proposeAPrimetoANaive));
		proposeAtoAPrime = log(exp(proposeAtoAPrime) + exp(proposeAtoAPrimeNaive));
	}
	
	if(considerAddMove == 1){
		proposeAPrimetoA = log(exp(proposeAPrimetoA) + exp(proposeAPrimetoAAdd));
		proposeAtoAPrime = log(exp(proposeAtoAPrime) + exp(proposeAtoAPrimeAdd));
	}
	
	if(considerRemoveMove == 2){
		proposeAPrimetoA = log(exp(proposeAPrimetoA) + exp(proposeAPrimetoARemove));
		proposeAtoAPrime = log(exp(proposeAtoAPrime) + exp(proposeAtoAPrimeRemove));
	}
	
	if(considerNaiveMove == 0 && considerAddMove == 0 && considerRemoveMove == 0 && moveToMake != 0){
		cerr << "Our proposal probability looks like a self-loop, but it isn't check proposal distribution.\n"; exit(-1); 
	}
		
	double RATIO = exp((probMoveProposed + proposeAPrimetoA) - (probMoveCurrent + proposeAtoAPrime));
	
	//alpha(M,M') = min{1,RATIO}; if it's better we always accept. Otherwise we do not.
	double alpha = 1;
	if(RATIO < 1){ alpha = RATIO; }
	
	if(moveToMake == 0){ alpha = 1; } //Automatically accept a self-edge
	
	double prob = ((double) rand()*1.0 / (RAND_MAX+1.0));
	if(prob < alpha){
		retVal = 1; //The return value is 1, we accepted the move!
		currentLikelihood = probMoveProposed;

		
		//cout << "Accept Move" << endl << flush; 	
		//For each case accept the move and then update the assignments.
		for(int i = 0; i<numModifiedESPs; i++){ ESPVector[modifiedESPs[i]].acceptMove(); }
		
		for(int i = 0; i<numModifiedVariants; i++){ 
			//0 = Nothing;
			//-1 = Move from Assigned to Empty;
			// 1 = Move from Empty to Assigned;
		
			int variant = modifiedVariants[i];
			int result = VariantVector[variant].acceptMove(); 
			
			//VariantVector[i].checkMobility(numMobileESPs,mobileESPs) == 1
			//HERE!
			if(result == 1){
				int location = findLocationLocal(EmptyVariants,numEmpty,variant);
				EmptyVariants.erase(EmptyVariants.begin() + location);
				numEmpty--;
				AssignedVariants.insert(AssignedVariants.begin(),variant); 
				numAssigned++;
				sort( AssignedVariants.begin(), AssignedVariants.begin() + numAssigned);
			}
			else if(result == -1){
				int location = findLocationLocal(AssignedVariants,numAssigned,variant);
				AssignedVariants.erase(AssignedVariants.begin() + location);
				numAssigned--;

				//Note: We CAN NOT pushback, because the sorted order will be screwed up!!!
				EmptyVariants.insert(EmptyVariants.begin(),variant); 
				numEmpty++;
				sort( EmptyVariants.begin(), EmptyVariants.begin() + numEmpty );			
			}
		}
	}
	else{
		//cout << "Do not accept move!\n";
		
		currentLikelihood = probMoveCurrent;
		
		retVal = 0; //We did not accept the move, the return value is 0,
		for(int i = 0; i<numModifiedESPs; i++){ ESPVector[modifiedESPs[i]].rejectMove(); }
		for(int i = 0; i<numModifiedVariants; i++){ VariantVector[modifiedVariants[i]].rejectMove(); }
	}
	//cout << "Accepted each of the Variant moves\n" << flush;
	
	//We need to update the occupancy of each ESP and each variant:
	for(unsigned int i = 0; i<ESPVector.size(); i++){ ESPVector[i].updateOccupancy(); }
	for(unsigned int i = 0; i<VariantVector.size(); i++){ 
		VariantVector[i].updateOccupancy(); 
		int D = VariantVector[i].getCurrentAssigned();
		//Note: Here is where we could rightly impose a threshold of some kind.
		if(D > 0){
			VariantVector[i].incrementLogLikelihood(VariantVector[i].getLikelihoodVariantGiven(perr,COVERAGE,COVERAGE_SCALED,LAVG,LDIS,D) - VariantVector[i].getLikelihoodErrorGiven(perr,COVERAGE,COVERAGE_SCALED,LAVG,LDIS,D)); 
			VariantVector[i].incrementLikelihoodVar(VariantVector[i].getLikelihoodVariantGiven(perr,COVERAGE,COVERAGE_SCALED,LAVG,LDIS,D)); 
		}
	}
	
	modifiedESPs.clear(); numModifiedESPs = 0;
	modifiedVariants.clear(); numModifiedVariants = 0;
	
	return retVal;
}


#endif

