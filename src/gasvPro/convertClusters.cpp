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
#include <string>
#include <fstream>
#include <time.h>
#include <vector>		   
#include <cstdlib>
#include <cstring>

using namespace std;

//itoa: Convert an integer to a string;
string itoa(int number)
{
    if (number == 0)
        return "0";
    string temp="";
    string returnvalue="";
    while (number>0)
    {
        temp+=number%10+48;
        number/=10;
    }
    for (int i=0;i<temp.length();i++)
        returnvalue+=temp[temp.length()-i-1];
    return returnvalue;
}

//String Tokenizer: Creates an array of strings from tokenized string
int split( vector<string> & theStringVector, const  string  & theString, const  string  & theDelimiter){
	int numFields = 0;
	size_t  start = 0, end = 0;
	
    while ( end != string::npos){
        end = theString.find( theDelimiter, start);
		
        // If at end, use length=maxLength.  Else use length=end-start.
        theStringVector.push_back( theString.substr( start, (end == string::npos) ? string::npos : end - start));
		
        // If at end, use start=maxSize.  Else use start=end+delimiter.
        start = (   ( end > (string::npos - theDelimiter.size()) )?  string::npos  :  end + theDelimiter.size());
		numFields++;
    }
	return numFields;
}

//String Tokenizer: Creates an array of integers from tokenized string
int splitToIntegers( vector<int> & theIntegerVector, const  string  & theString, const  string  & theDelimiter){
	int numFields = 0;
	size_t  start = 0, end = 0;
	
    while ( end != string::npos){
        end = theString.find( theDelimiter, start);
		
        // If at end, use length=maxLength.  Else use length=end-start.
		string tmpString = theString.substr( start, (end == string::npos) ? string::npos : end - start);
        theIntegerVector.push_back( atoi(tmpString.c_str()) );
		
        // If at end, use start=maxSize.  Else use start=end+delimiter.
        start = (   ( end > (string::npos - theDelimiter.size()) )?  string::npos  :  end + theDelimiter.size());
		numFields++;
    }
	return numFields;
}

/****************/
/* MAIN PROGRAM */
/****************/

int main(int argc, char* argv[]){


     if(argc != 3 && argc != 2){
                //                 0        1              2                      3              4       5         6       7        8       9           // 10
				cout << "\nconvertClusters: Convert a set of GASV clusters (regions format) to intervals or reads format" << endl;
				cout << "Version: 1.0" << endl << endl;
				cout << "Usage: ./exe {clusterfile} {type (optional)} " << endl;
				cout << "\n";
                cout << "Parameters Considered:\n";
				cout << "           type:  {GASV Format}    (Default: standard)\n";
				cout << "\n";
				cout << "Output: {clusterfile}\n";
                cout << "cluterfile is converted from regions to desired output";
				exit(-1);
	 }

	
	//SET DEFAULTS HERE
	int fileType = 0; //Default file type is intervals!
	//Step 0 Command Line Argument
	//cout << "Step 0: Processing Command Line Arguments.\n";
		
	if(argc == 3){
		string temp = argv[2];
        if(temp == "intervals" || temp == "INTERVALS" || temp == "Intervals"){
            fileType = 0;
	    cout << "Outputting to standard cluster format.\n";
        }
        else if(temp == "standard" || temp == "STANDARD" || temp == "Standard" ){
            fileType = 0;
	    cout << "Outputting to standard cluster format.\n";
        }
        else if(temp == "reads" || temp == "READS" || temp == "Reads"){
            fileType = 1;
	    cout << "Outputting to reads cluster format.\n";
        }
        else if(temp == "regions" || temp == "Regions" || temp == "REGIONS"){
            fileType = 2;
	    cout << "Outputting to regions cluster format.\n";
        }
        else{
            cout << "GASV Output Type \"" << argv[2] << "\" is invalid." << endl;
            cout << "Assuming standard output\n";
            fileType = 0;
        }
        
    }

    //Step 1: Process Clusters File;
	cout << "Formatting Final Clusters File.\n";
	
	int translocationCount = 0;
	int divergentCount = 0;	
	int inversionCount = 0;
	int deletionCount = 0;
	int numIgnored = 0;
	int numConsidered = 0;
	int numNegativeLocalization = 0;
	int numCorrectlyProcessed = 0;
	
	//Input File:
	string CLUSTERFILE = argv[1];
    string TMP_FILE = CLUSTERFILE + ".tmp";
    string CMD = "cp " + CLUSTERFILE + " " + TMP_FILE;
    system(CMD.c_str());
        
	ifstream clusterFile(TMP_FILE.c_str(),ios::in);
			
	ofstream outFile(CLUSTERFILE.c_str(),ios::out);
			
	string clusterLine;
	string clusterID;
    string prList;
        
    int support;
	int chrL;
	int chrR;
	int start; 
	int end;
	double localization;
	string type;
    double logLikeRatio;
	int localType; //Tracks which model
	
	bool tooManyFieldsWarning = false;
	bool processCluster;
	int numLines = 0;
	
	int numCasesBeforePrint = 1000; //print every numCasesBeforePrint cluster ids
	
	//Example Line:
	//0     1     2     3          4            5   6            7
	//c1	1	204.2	D	SRR004856.7363154	1	1	746128, 748510, 746333, 748510, 746027, 748204, 746027, 748409
	while(getline(clusterFile,clusterLine) ){
        
		/////////////////////////
		// Read and Parse Line //
		/////////////////////////
		
		vector<string> v;

		int numFields = split( v, clusterLine.c_str(), "\t" );	
		if(numFields > 9 && !tooManyFieldsWarning){
			cout << "\tError: There are too many fields in your clusters file check output.\n";
            cout << "We expect 8 or 9 fields and there were " << numFields << endl;
            for(int i = 0; i<numFields; i++){
                cout << "Field " << i << ": " << v[i] << endl;
            }
            cout << clusterLine << endl;
			tooManyFieldsWarning = true;
            exit(-1);
		}
		if(numFields < 8){
			cout << "\tError: We require clusters are in regions format. Clusters file has only " << numFields << " columns instead of 8.\n";
			exit(-1);
		}
        
		//Note: We allow 8 or 9 fields!
        
		if(v[0].substr(0,1) == "#"){
            //Note: Headerlines taken from ClusterESP.java for GASV (lines 35-37).
            if(fileType == 0){
                outFile << "#Cluster_ID:\tLeftChr:\tLeftBreakPoint:\tRightChr:\tRightBreakPoint:\tNum PRS:\tLocalization:\tType:";
                if(numFields == 9){ outFile << "\tLogLikelihoodRatio:"; }
                outFile << endl;

            }
            else if(fileType == 1){
                outFile << "#Cluster_ID:\tLeftChr:\tLeftBreakPoint:\tRightChr:\tRightBreakPoint:\tNum PRS:\tLocalization:\tType:\tList of PRS:";
                if(numFields == 9){ outFile << "\tLogLikelihoodRatio:"; }
                outFile << endl;
            }
            else if(fileType == 2){
                //Regions format, just output as written;
                outFile << "#Cluster_ID:\tNum PRS:\tLocalization:\tType:\tList of PRS:\t LeftChr:\tRightChr:\tBoundary Points:";
                if(numFields == 9){ outFile << "\tLogLikelihoodRatio:"; }
                outFile << endl;
            }
            else{
                cout << "Invalid File Type Set: Exiting. Repeat with no file type argument for intervals default format\n";
                exit(-1);
            }
        }
        else{
            if(fileType == 2){
                outFile << clusterLine << endl;
            }
            else{
            
                clusterID    = v[0];
                support      = atoi(v[1].c_str());
                localization = atof(v[2].c_str());
                type         = v[3];
                prList       = v[4];
            
                chrL         = atoi(v[5].c_str());
                chrR         = atoi(v[6].c_str());
                if(numFields == 9){
                    logLikeRatio = atof(v[8].c_str());
                }
            

                ////////////////////////
                // Processing Cluster //
                ////////////////////////
                
                numCorrectlyProcessed++;
                
                vector <int> bothCoordinates;
                int totalCoords = splitToIntegers(bothCoordinates,v[7],",");
                int minX, maxX, minY, maxY;
                int numCoords = totalCoords/2;
                minX = maxX = bothCoordinates[0];
                minY = maxY = bothCoordinates[1];
                for(int k = 2; k<totalCoords; k++){
                    if(k%2 == 0){
                        if(bothCoordinates[k]<minX){ minX = bothCoordinates[k];}
                        if(bothCoordinates[k]>maxX){ maxX = bothCoordinates[k];}
                    }
                    else{
                        if(bothCoordinates[k]<minY){ minY = bothCoordinates[k];}
                        if(bothCoordinates[k]>maxY){ maxY = bothCoordinates[k];}
                    }
                }
                
                if(fileType == 0){
                    outFile << clusterID << "\t"
                            << chrL << "\t" << minX << "," << maxX << "\t"
                            << chrR << "\t" << minY << "," << maxY << "\t"
                            << support << "\t" << localization << "\t" << type;
                    if(numFields == 9){
                        outFile << "\t" <<  logLikeRatio;
                    }
                    outFile << endl;
                }
                else if(fileType == 1){
                    outFile << clusterID << "\t"
                    << chrL << "\t" << minX << "," << maxX << "\t"
                    << chrR << "\t" << minY << "," << maxY << "\t"
                    << support << "\t" << localization << "\t" << type
                    << "\t" << prList;
                    if(numFields == 9){
                        outFile << "\t" <<  logLikeRatio;
                    }
                    outFile << endl;
                }
                else{
                    cout << "Invalid File Type Set: Exiting. Repeat with no file type argument for intervals default format\n";
                    exit(-1);
                }
            }
                        
            numLines++;
            
            if(numLines%numCasesBeforePrint==0){
                cout << "\tProcessing cluster " << numLines << " with clusterID " << clusterID << endl << flush;
            }
        }
	}
	
	clusterFile.close();
	outFile.close();
    
    string CLEAN = "rm " + TMP_FILE;
    system(CLEAN.c_str());
    
	cout << "\tSuccessfully Finished.\n";
		
	return 0;	
}
