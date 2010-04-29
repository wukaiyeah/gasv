/**
 * Copyright 2010 Benjamin Raphael, Suzanne Sindi, Hsin-Ta Wu, Anna Ritz, Luke Peng
 *
 *  This file is part of gasv.
 * 
 *  gasv is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  gasv is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with gasv.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <math.h>
#include <ctype.h>
using namespace std;
void process_getLminLmax(string fn, int mapq, string outfile, string prefix, string location, int cutoff, int truncate_length);
void outputResult(string o, string p, int min, int max);
vector<string> parsingSD_PERC(string q);


int main(int argc, char** argv) {

	if(argc != 8){
		cerr << "Please check parameters!" << argc << endl;
		cerr << "Paramters passed were: " << argv[1] << ", " << atoi(argv[2]) << ", " << argv[3] << ", " << argv[4] << ", " << argv[5] << ", " << atoi(argv[6]) << ", " << atoi(argv[7]);
	}
	else{
		process_getLminLmax(argv[1], atoi(argv[2]), argv[3], argv[4], argv[5], atoi(argv[6]), atoi(argv[7]));
	}
}

vector<string> parsingSD_PERC(string quantile){

	vector<string> qt_tokens;
        string de_qt = "=";
	string::size_type lastPos = quantile.find_first_not_of(de_qt, 0);
	string::size_type pos = quantile.find_first_of(de_qt, lastPos);
	while(string::npos != pos || string::npos != lastPos){
		qt_tokens.push_back(quantile.substr(lastPos, pos - lastPos));
		lastPos = quantile.find_first_not_of(de_qt, pos);
		pos = quantile.find_first_of(de_qt, lastPos);
	}

	if(qt_tokens[0] == "PCT"){
		// remove %
		string perc_tmp = qt_tokens[1].substr(0, qt_tokens[1].length() - 1);
		qt_tokens[1] = perc_tmp;
	}
	return qt_tokens;
}

string prefix_parser(string prefix){
	vector<string> tmp_tokens;
	string de = "/";
	string::size_type lastPos = prefix.find_first_not_of(de, 0);
	string::size_type pos = prefix.find_first_of(de, lastPos);
	while(string::npos != pos || string::npos != lastPos){
		tmp_tokens.push_back(prefix.substr(lastPos, pos - lastPos));
		lastPos = prefix.find_first_not_of(de, pos);
		pos = prefix.find_first_of(de, lastPos);
	}

	return tmp_tokens.back();
}

void process_getLminLmax(string fn, int mapq, string outfile, string prefix, string quantile, int cutoff, int truncate_length){

	int sd_int = 0;
	float perc = 0;
	float Lmin_q = 0;
	float Lmax_q = 0;
	int lmin_out = 0;
	int lmax_out = 0;

	int allcount = 0;
	int counter = 0;
	long int total_length = 0;
	map<int, int> fragl;
	map<string, int> reads;
	string line;
	string outprefix = prefix_parser(prefix);

	int minReadLength = 100000000; // a very high number to start off
	ifstream r_file(fn.c_str());
	if(r_file.is_open()){
		while(!r_file.eof()){
			getline(r_file, line);
			if(line != "" && allcount < cutoff){
				string delimiters = "	";
				vector<string> tokens;
				string::size_type lastPos = line.find_first_not_of(delimiters, 0);
				string::size_type pos = line.find_first_of(delimiters, lastPos);
				while(string::npos != pos || string::npos != lastPos){
					tokens.push_back(line.substr(lastPos, pos - lastPos));
					lastPos = line.find_first_not_of(delimiters, pos);
					pos = line.find_first_of(delimiters, lastPos);
				}
	
				if(atoi(tokens[4].c_str()) > mapq){
					if(reads[tokens[0]] == 0){
						reads[tokens[0]]++;
					}
					else{
						if (tokens[9].length() < minReadLength) {
							minReadLength = tokens[9].length();
						}
						counter++;
						int flength = abs(atoi(tokens[8].c_str()));
						if(truncate_length == 0){
							fragl[flength]++;
							total_length += flength;
							allcount++;
						}
						else{
							if(flength <= truncate_length){
								fragl[flength]++;
								total_length += flength;
								allcount++;
							}
						}
					}
				}
			}
			else if (allcount >= cutoff){
				break;
			}
		}
	}
	r_file.close();

	vector<string> qt_tokens = parsingSD_PERC(quantile);
	if(qt_tokens[0] == "SD"){
		sd_int = atoi(qt_tokens[1].c_str());
		double mean_length = (double)total_length/(double)allcount;
		// get sd
		double sd_up = 0;
		map<int, int>::iterator iter;
		for(iter = fragl.begin(); iter != fragl.end();++iter){
			sd_up += (iter->second)*pow(iter->first-mean_length,2);
		}
		double sd = sqrt(sd_up/allcount);
		lmin_out = (int)(mean_length - sd*sd_int);
		lmax_out = (int)(mean_length + sd*sd_int);
	}

	else if(qt_tokens[0] == "PCT"){
		perc = atof(qt_tokens[1].c_str());
		if (perc < 50) { 
			cerr << "WARNING: PCT=" << perc << "% would result in an Lmin LARGER than Lmax!!  Proceeding under assumption user meant to input ";
			perc = 100 - perc;
			cerr << "PCT=" << perc << "% instead!" << endl;
		}

		Lmin_q = (100 - perc)/100;
		Lmax_q = perc/100;

		int lmin_index = allcount*Lmin_q;
		int lmax_index = allcount*Lmax_q;
		int current_counting = 0;
		map<int, int>::iterator iter;
		for(iter = fragl.begin(); iter != fragl.end();++iter){
			if(current_counting <= lmin_index && current_counting+iter->second >= lmin_index){
				//cerr << "Lmin: " << iter->first << endl;
				lmin_out = iter->first;
			}
			if(current_counting <= lmax_index && current_counting+iter->second >= lmax_index){
				//cerr << "Lmax: " << iter->first << endl;
				lmax_out = iter->first;
			}
			current_counting += iter->second;
			//cout << iter->first << "      " << iter->second << endl;
		}
	}
	else{
		cerr << "Wrong Format of Quantile Tag: SD=(int), PCT=(float)%, EXACT=(int),(int)" << endl;
		exit(0);
	}

	if (lmin_out < (2 * minReadLength)) {
		cerr << "Calculated Lmin value of " << lmin_out << " is less than (2 * minReadLength) = " << (2*minReadLength) 
			<< " so automatically increasing Lmin value to " << (2*minReadLength) << endl;
		lmin_out = 2 * minReadLength;

	}
	outputResult(outfile, outprefix, lmin_out, lmax_out);

	// output distribution
	//map<int, int>::iterator iter;
	//for(iter = fragl.begin(); iter != fragl.end(); ++iter){
		//cerr << iter->first << "	" << iter->second << endl;
	//}
}

void outputResult(string outfile, string outprefix, int lmin_out, int lmax_out){
	//ofstream outESP(outfile.c_str());
	//outESP << outprefix << "	" << "ESP" << "	" << lmin_out << "	" << lmax_out << endl; 
	//outESP.close();
	cout << lmin_out << "\t" << lmax_out << endl;
}

