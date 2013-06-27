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
#include <fstream>
#include <string>
#include <stdio.h>
#include <cmath>
#include <sstream>
#include <map>
#include <vector>
#include <sys/stat.h>


//Additional dependencies.
#include "espmapping.h"
#include "strutils.h"
#include "LinkedList.h"

using namespace std;

/* Functions Needed to Access Array of Mappings */
int compareByIndex(const void *a, const void *b){ 
	const espmapping *ia = (const espmapping *)a;
	const espmapping *ib = (const espmapping *)b;
	return (ia->index - ib->index);
}
int compareByName(const void *a, const void* b){ 
	const espmapping *ia = (const espmapping *)a;
	const espmapping *ib = (const espmapping *)b;
	
	return ((ia->fragName).compare(ib->fragName));
	return 0; 
}

int findLocation(string value, espmapping *espMappingIndex, int maxValues){
	 int first = 0;
	 int last = maxValues - 1;
	 int mid = (first + last)/2; //Mid point
	 
	 while( first <= last){
		 mid = (first + last)/2;
		 if(value.compare(espMappingIndex[mid].getFrag()) > 0){ first = mid + 1; }
		 else if ( value.compare(espMappingIndex[mid].getFrag()) < 0){ last = mid - 1;}
		 else { return mid; }
	 }
	 return -(first + 1);
}

int findLocationIndex(int index, espmapping *espMappingIndex, int maxValues){
	int first = 0;
	int last = maxValues - 1;
	int mid = (first + last)/2; //Mid point
	
	while( first <= last){
		mid = (first + last)/2;
		if(index < espMappingIndex[mid].getIndex() ){ first = mid + 1; }
		else if (index > espMappingIndex[mid].getIndex() ){ last = mid - 1;}
		else { return mid; }
	}
	return -(first + 1);
}
/**************/

/* Functions Needed to Extract Information from Fragment Names */
string pruneFragment(string s){ 
	if(s.at(s.length()-1)==','){ s = s.substr(0,s.length()-1);}
	return s;
}

//(1) Determine the string following the last ".".
//(2) Determine if there are .INT_INT_INT
//--> Only if following the last period there are: (1) 2 dashes and (2) Only integers is this ambiguous;
string get_frg_name(string  s){
	//This is what ambiguous fragments look like: 
	//chr17_17261_17474_2:1:0_1:0:0_2ea8b32.0_2_1
	//Step 1: Get rid of comma if necessary;
	if(s.at(s.length()-1)==','){ s = s.substr(0,s.length()-1);}
	
	int slen = s.length();
	int c = slen-1;
	int maxLen = slen-1;
	while(c>=0 && s.at(c) != '.'){c--; }
	if(c < 0){ return s; }// cout << "Unique_No_Period\t" << s << endl; }
	else{
		int dashCount, numCount, expectedNumCount;
		dashCount = numCount = expectedNumCount = 0;
		for(int i = c+1; i<=maxLen; i++){
			if(s.at(i) == '_'){ dashCount++; }	else{ expectedNumCount++; }
			if(s.at(i) == '0' || s.at(i) == '1' || s.at(i) == '2' || s.at(i) == '3' ||
			   s.at(i) == '4' || s.at(i) == '5' || s.at(i) == '6' || s.at(i) == '7' ||
			   s.at(i) == '8' || s.at(i) == '9'){ numCount++; }
		}
		//Are we ambiguous??
		//YES:
		if(dashCount == 2 && numCount == expectedNumCount){
			//cout << "Ambiguous\t" << s << "\t" << s.substr(0,c) << endl;
			return s.substr(0,c);
		}
		//NO: Return the full string;
		else{
			//cout << "Unique_Failed_Count\t" << s << endl;
			return s; 
		}
	}	
}

int get_dist(string s){
	int c= 0;
	while( c<s.length() && s.at(c) != '_'){ c++; }
	return atoi(s.substr(0,c));	
}

/**************/

/* Selim's Previous Functions*/
int which_class(int* arr, int index){
	if(arr[index] < 0){ return -1 * arr[index]; }
	else{ return which_class(arr,arr[index]); } 
}

long int choose2(int ne) { return ne*(ne-1)/2;	}

struct int_pair{
	int i1;
	int i2;
	int_pair(int ii1, int ii2)
	:i1(ii1),i2(ii2)
	{}
};
/**************/


/*
 This program parses the clusters and coverage file several times
 
 (1) File Read 1: Removes all cases with -1 localization; (makes *new suffix)
 (2) File Read 2: Indexes the fragments with appropriate named prefix: "a_b_c" --> numberOfFragment, numberOfFragmentMapping, totalMappingIndex;
 (3) File Read 3: Determines (for each fragment how many clusters it supports)
 (4) Determines the connected components.
 (5) Outputs the connected components.
*/

int main(int argn, char* argv []){
	
	/********************************/
	/* Step 0: Processing Arguments */
	/********************************/

	if(argn < 4 ){ 
		cerr << "GASVPro-Graph: Separates Ambiguous Clusters in to Dependent Components\n";
		cerr << "Version:       1.2.1\n\n";
		
		cerr << "Usage: ./GASVPro-Graph {ClustersFile} {CoverageFile} {OutputDir} {MinimimSupport}\n";
		
		cerr << "\n";
		cerr << "\t\t*gasv.in.clusters         \\\\GASV ClustersFile\n";
		cerr << "\t\t*gasvPro-CC.in.coverage   \\\\GASVPro-CC Coverage File\n";
		cerr << "\t\t<dir>                     \\\\Desired Output Directory\n";
		cerr << "\t\t<int>                     \\\\Minimum Support for a Cluster to be Considered (Default: 1)\n";
		return -1;
	}
	
	string clusters_file, coverage_file, output_path, summary_out, strtemp;
	string esp_file_list;
//	string infos_fn = argv[1];
//	ifstream infos;

	int num_of_frg_threshold = 4; //Setting to the same default value as MIN_CLUSTER_SIZE in GASV.
	
	cerr << "Processing Command line arguments....\n";

	clusters_file = argv[1];
	coverage_file = argv[2];
	output_path   = argv[3];
	if(argn == 5){
		strtemp = argv[4];
		num_of_frg_threshold = atoi(strtemp.c_str());
		if(num_of_frg_threshold<=0){ cerr << "Error: Minimim Threshold must be positive.\n"; }
	}
	
	struct stat st;
	if(stat(output_path.c_str(),&st) == 0){
		cerr << "\tOutput directory " << output_path << " exists.\n";
	}	
	else{
		cerr << "\tCreating directory " << output_path << endl;
		string rem_com = "mkdir " + output_path;
		system(rem_com.c_str());
	}
	cerr << "\tDONE\n";
	
	/********************************************/
	/* Step 1: Read/Store All Deletion Clusters */
	/*                                          */
	/* Goal: Remove all clusters w/ -1 loc      */
	/********************************************/
	
	cerr << "Removing all clusters with -1 localization....";
	
	ifstream clus_inp_del1;
	clus_inp_del1.open(clusters_file.c_str());
	if(clus_inp_del1.fail()){ cerr<<"Input failed - 1"<<endl; }
	
	string cluster, inp_f_n = clusters_file + ".new";
	string inp_f_n_orig = inp_f_n;
	double loc;
	ofstream output_inp_del1;
	output_inp_del1.open(inp_f_n.c_str());
	ofstream output_for_cov;

	string coverage_file_new = coverage_file + ".new";
	output_for_cov.open(coverage_file_new.c_str());
	
	//Concordants == coverageFile
	ifstream concordants ;
	concordants.open(coverage_file.c_str());
	if(concordants.fail()) { cerr<<"Concordants failed"<<endl; }
	
	//cout << "Ready to process clusters.\n";
	while(getline(clus_inp_del1,cluster)){
		vector<string> tokens;		
		if(cluster[0] == '#'){getline(clus_inp_del1,cluster);}		
		istringstream line1(cluster);
		while(line1.good())
		{
			string temp;
			line1>>temp;
			tokens.push_back(temp);
		}
		loc = atof(tokens[2].c_str());
		getline(concordants,strtemp); //Read in the concordants; output them.
		if(strtemp[0] == '#'){getline(concordants, strtemp);}
		if(loc>=0){
			string finalout;
            //Remove the final column for Log-Likelihood!
			for(int i = 0; i < (tokens.size()-1); i++)
			{
				if(tokens[i][tokens[i].length()-1] != ',')
					finalout+=tokens[i]+"\t";
				else
					finalout+=tokens[i]+" ";
			}
			output_inp_del1<<finalout<<endl; //Remove the last two columns (LLR and #Copies) because we want what is output to be the original GASVclusters mode.
			output_for_cov<<strtemp<<endl;
		}
	}
	output_for_cov.close();
	output_inp_del1.close();
	clus_inp_del1.close();
	concordants.clear();
	concordants.close();
	concordants.open(coverage_file_new.c_str());
	if(concordants.fail()){ cerr<<"Concordants new failed"<<endl; }
	
	cerr << "DONE\n";

	
	/**********************************************************/
	/* Step 2: Count Fragments and Mappings output to renamed */
	/**********************************************************/

	cerr << "Idexing all discordant fragments.....";
	
	clus_inp_del1.open(inp_f_n.c_str());
	if(clus_inp_del1.fail()) { cerr<<"Input failed - 2"<<endl; }
	else { /*cout<<"Input - pruned clusters : OK"<<endl;*/ }

	//cout << "Reading file --> " << inp_f_n.c_str() << endl;

	ofstream output_renamed_frgs;
	inp_f_n += ".renamed";
	output_renamed_frgs.open(inp_f_n.c_str()); //The renamed file contains the new fragment names; 
	map<string,int_pair> amb_frgs_from_same_frg;
	map<string,int_pair>::iterator map_itr;
	string cn;
	int num_of_frg;
	string frg;
	string frg_name;
	int dist_frg_cntr = 0;
	int frg_amb_count = 0;
	int total_frg_ctr = 0;
	string type;
	
	//Read in name;
	while(clus_inp_del1>>cn){
		//Read in num fragments, localizaiton and type DEL
		clus_inp_del1>>num_of_frg>>loc>>type;
		//Output to RENAMED file, cluster name, number of fragments and localization
		output_renamed_frgs<<cn<<"\t"<<num_of_frg<<"\t"<<loc<<"\t"<<type<<"\t";
		
		//Compute the indexing for each fragment based on:
		//(1) number of times we've seen it; 
		//(2) number of mappings we've seen of it so far;
		//(3) total number of mappings in total so far
		//map_itr->i1 -- Index of Read (1)
		//map_itr->i2 -- number of mappings so far (2) 
		for(int i =1; i<=num_of_frg; i++){
			clus_inp_del1>>frg;	
			//int frgLocation = findLocation(pruneFragment(frg),espMappingIndex,numMappings);
			//if(frgLocation<0){ cerr << "Invalid: A fragment name we do not have in the ESP file:\t" << pruneFragment(frg) << ".\n"; exit(-1); }
			
			frg_name = get_frg_name(frg);
			map_itr = amb_frgs_from_same_frg.find(frg_name);
			
			int valA, valB, valC;
			
			//Q: Have we seen this fragment before?
			//NO: Create a new entry for it.
			if(map_itr == amb_frgs_from_same_frg.end()){
				dist_frg_cntr++;
				total_frg_ctr++;
				frg_amb_count = 1;
				amb_frgs_from_same_frg.insert(pair<string,int_pair>(frg_name,int_pair(dist_frg_cntr,frg_amb_count)));
				valA = dist_frg_cntr;
				valB = frg_amb_count;
				valC = total_frg_ctr;
				output_renamed_frgs<<dist_frg_cntr<<"_"<<frg_amb_count<<"_"<<total_frg_ctr<<"_"<<frg<<" ";
			}
			//YES: Increment the count and output information;
			else{	
				(map_itr->second).i2++;
				total_frg_ctr++;
				valA = (map_itr->second).i1;
				valB = map_itr->second.i2;
				valC = total_frg_ctr;
				output_renamed_frgs<<(map_itr->second).i1<<"_"<<map_itr->second.i2<<"_"<<total_frg_ctr<<"_"<<frg<<" ";
			}

		}
		getline(clus_inp_del1,strtemp);
		//output_renamed_frgs<<"\t"<<strtemp<<endl;
		output_renamed_frgs<<strtemp<<endl;
		//Commented Out: cout << "\tFinished Cluster name:\t" << cn << endl;
	}
	output_renamed_frgs.close();
	
	cerr << "DONE\n";
	
	/************************************************************************/
	/* Step 3: Determine the Number of Clusters Supported by Each Fragment  */
	/************************************************************************/
	
	cerr << "Determining clusters supported by each fragment....";
	
	List<int>* frgs_supporting_sv; // indexed by frg
	frgs_supporting_sv = new List<int> [dist_frg_cntr+1]; //Ok - each fragment gets a list of the SVs it supports.
	ListItr<int> itr_for_frgs_sup_sv; //And we make an iterator for this (so we can go next?)
	int* frg_sup_sv_cntr = new int[dist_frg_cntr+1];

	//Initialize all fragments NOW support zero variants.
	for(int i= 0; i<= dist_frg_cntr; i++){ frg_sup_sv_cntr[i]=0; }

	//Q: Why are these two things indexed from 0 to dist_frg_cntr and 1 to dist_frg_cntr respectively?
	int* which_frg_sup_which_sv = new int[dist_frg_cntr+1]; //Now an interger list indexed by fragment; // new int[total_frg_cntr+1];
	for(int t=1; t<= dist_frg_cntr; t++){ which_frg_sup_which_sv[t] = 0; }

	//Making the new_sv.clusters file;
	ofstream new_sv;
	new_sv.open("new_sv.clusters");

	string conc_str;
	ofstream new_ver_con_inp;
	new_ver_con_inp.open("new_con.coverage");
	string all_cluster;
	int num_of_sup_frg;
	int temp_int;
	int sv_cntr = 0;
	double localization;
	string con_info_line;

	clus_inp_del1.clear();
	clus_inp_del1.close();
	clus_inp_del1.open(inp_f_n.c_str());

	while(getline(clus_inp_del1,all_cluster)){
		//Make two input string streams from the cluster
		istringstream clus_inp1(all_cluster);
		clus_inp1>>strtemp>>num_of_sup_frg>>localization>>type;
		
		//Get the concordant file information;
		getline(concordants,con_info_line);

		if(num_of_sup_frg>=num_of_frg_threshold){
			sv_cntr++;
			new_ver_con_inp<<con_info_line<<endl; //Simply write out the SAME information as in the concordant file.
			new_sv<<all_cluster<<endl;
			for(int s= 1; s<= num_of_sup_frg; s++){
				//Separate the Unique Identifier:
				//strtemp:	1_1_1_IL17_297:2:56:439:948.8_4_2
				//tmp_int:	1
				clus_inp1>>strtemp;		
				temp_int = get_dist(strtemp);
				which_frg_sup_which_sv[temp_int] = sv_cntr; //This stores the current cluster for a fragment
															//Note if we had mappings in the same cluster
				itr_for_frgs_sup_sv = frgs_supporting_sv[temp_int].zeroth(); //Set iterator to the 0th position
				frgs_supporting_sv[temp_int].insert(sv_cntr,itr_for_frgs_sup_sv); //Insert the cluster at the head
				frg_sup_sv_cntr[temp_int]++; //Increment the number of clusters this fragment supports.
			}			
		}
	}
	concordants.close();
	new_ver_con_inp.close();
	new_sv.close();
	clus_inp_del1.close();

	cerr << "DONE.\n";
	
	/*****************************************/
	/* Step 4: Determine Connected Components*/
	/*****************************************/
	
	cerr << "Determining Connected Components....";
	
	int* related_svs = new int[sv_cntr+1]; //These are the categories of related variants.
	for(int i=0; i<= sv_cntr; i++){ related_svs[i] = -1*i; }
	
	int size, min;
	int* classes, *svs_for_class;

	for(int f= 1; f<= dist_frg_cntr; f++){
		size= frg_sup_sv_cntr[f];  //The number of SVs a cluster belongs to;
		classes = new int[size+1]; //
		svs_for_class = new int[size+1];
		itr_for_frgs_sup_sv = frgs_supporting_sv[f].first();
		min = sv_cntr;
		for(int s=1; s<= size; s++){
			svs_for_class[s]=itr_for_frgs_sup_sv.retrieve();
			classes[s] = which_class(related_svs,itr_for_frgs_sup_sv.retrieve());
			if(min > classes[s]){ min = classes[s];}
			itr_for_frgs_sup_sv.advance();
		}

		for(int s=1; s<= size; s++){
			if(classes[s]!=min){ related_svs[classes[s]] = min; }
		}
		
		delete [] classes;
		delete [] svs_for_class;
		
	}
	int* new_related_svs = new int[sv_cntr+1];
	int* already_taken = new int[sv_cntr+1];
	for(int s=1 ; s<= sv_cntr; s++){ already_taken[s] = -1; }
	
	//Actually counting the SV clusters
	int num_of_sv_clusters = 0;
	for(int s=1 ; s<= sv_cntr; s++){
		temp_int= which_class(related_svs,s);
		if(0>already_taken[temp_int]){
			num_of_sv_clusters++;
			already_taken[temp_int]= num_of_sv_clusters;
		}
		new_related_svs[s] = already_taken[temp_int];
	}

	delete[] related_svs;
	related_svs = new_related_svs;
	new_related_svs = NULL;
					
	cerr << "DONE.\n";
	
	/******************************************/
	/* Step 5: Output the Connected Components*/
	/******************************************/
	
	cerr << "Outputting Connected Components....\n";
	
	int* sv_clus_counter = new int[num_of_sv_clusters+1]; 
	for(int nsvc = 1 ; nsvc <= num_of_sv_clusters; nsvc++){ sv_clus_counter[nsvc] = 0; }
	for(int s=1 ; s<= sv_cntr; s++){
		sv_clus_counter[related_svs[s]]++;
	}
	ofstream* outputs;
	outputs = new ofstream[num_of_sv_clusters+1];

	int num_of_singletons = 0;
	for(int s=1; s<= num_of_sv_clusters; s++){
		if(sv_clus_counter[s] == 1){ //Only singletons;
			num_of_singletons += sv_clus_counter[s];
		}
	}
	
	sv_clus_counter[0] = num_of_singletons; 
	string filename;
	int max_num_of_svs = 0;
	int numOpenedFiles = 0;
	ofstream outputMaster;
	string tmpFile = output_path + "/svFileList.summary"; 
	outputMaster.open(tmpFile.c_str());
	
	int sv_c = 0; // counter
	int max_num_of_ele= 0;
	int fileRound = 0;
	int sizePerRound = 200;
	int total_sv_c = 0;
	//cerr << "Printing in rounds because a limit on the number of file streams open!\n";
	
	for(int k = 0; k<=num_of_sv_clusters; k=k+sizePerRound){
		fileRound++;
		if(k%500 == 0){
			cerr << "\tOutputting components from " << k << " to " << k+sizePerRound-1 << endl;
			cerr << "\t";
			system("date");
		}
		
		ifstream new_ver_sv2;
		new_ver_sv2.open("new_sv.clusters");
		if(new_ver_sv2.fail()){cerr<<"new-sv failed"<<endl;}
		
		ifstream new_ver_con;
		new_ver_con.open("new_con.coverage");
		if(new_ver_con.fail()){ cerr<<"conc-failed"<<endl; }		
		
		for(int i= k; i<= (k+sizePerRound-1) && i<=num_of_sv_clusters; i++){  // outputs[0] is for svs supported by uniqely mapped  or - svs not in contact with any other
			filename = output_path + "/sv_" + itoa(i) + ".sv" ;
			if( i ==  0){
				outputMaster << filename << "\t" << sv_clus_counter[i] << endl;
				outputs[i].open(filename.c_str());
				outputs[i]<<sv_clus_counter[0]<<endl;
				numOpenedFiles++;
			}
			else if(sv_clus_counter[i] > 1){
				outputMaster << filename << "\t" << sv_clus_counter[i] << endl;
				numOpenedFiles++;
				outputs[i].open(filename.c_str());
				outputs[i]<<sv_clus_counter[i]<<endl;
				if(max_num_of_svs<sv_clus_counter[i]) { max_num_of_svs = sv_clus_counter[i];}
			}	
		}
				
		sv_c = 0;
		while(getline(new_ver_sv2,strtemp)){
            
			sv_c++;
			int tmpVal = related_svs[sv_c];
			getline(new_ver_con,con_info_line);
			if( tmpVal>=k && tmpVal<=(k+sizePerRound-1)){
				if(sv_clus_counter[related_svs[sv_c]] <= 1){
					outputs[0]<<strtemp<<endl;
					outputs[0]<<con_info_line<<endl;
				}
				else{
					if(max_num_of_ele <sv_clus_counter[related_svs[sv_c]]){ max_num_of_ele = sv_clus_counter[related_svs[sv_c]];}
					outputs[related_svs[sv_c]]<<strtemp<<endl;
					outputs[related_svs[sv_c]]<<con_info_line<<endl; 	
				}
				total_sv_c++;
			}
		}
		
		for(int i= k; i<=(k+sizePerRound-1) && i<=num_of_sv_clusters; i++){ 
			if(i != 0){ outputs[i].close(); }
		}
		
		new_ver_sv2.close();
		new_ver_con.close();
	}
	outputs[0].close();
	outputMaster.close();
	
	cerr << "\tDONE\n"; 
	
	if(total_sv_c != sv_cntr){ cerr<<"Check Possible error: total_sv_c("<<total_sv_c<<") does not equal sv_cntr("<<sv_cntr << ")" <<endl; }
	
	/**************************************/
	/* Step 6: Output Summary Information */
	/**************************************/
	
	cerr << "Outputting summary statistics....";
	
	ofstream sum;
	summary_out = output_path + "/"  + "p_star.summary";
	sum.open(summary_out.c_str());
	sum<<num_of_sv_clusters<<" "<<total_frg_ctr<<" "<<dist_frg_cntr<<endl;
	sum.close();
	
	cerr << "DONE\n";

	/********************/
	/* Step 7: Clean-Up */
	/********************/
	delete [] outputs;
	delete [] related_svs;
	delete [] which_frg_sup_which_sv;
	delete [] frg_sup_sv_cntr;
	delete [] frgs_supporting_sv;

	//Remove temporary files:
    string clean0 = "rm -rf " + inp_f_n_orig;
	string clean1 = "rm -rf " + inp_f_n;
	string clean2 = "rm -rf " + coverage_file_new;
	string clean3 = "rm -rf new_sv.clusters";
	string clean4 = "rm -rf new_con.coverage";
	system(clean0.c_str());
	system(clean1.c_str());
	system(clean2.c_str());
	system(clean3.c_str());
	system(clean4.c_str());
    
	
	
	cerr << "GASVPro-graph finished successfully.\n";
	
	return 0;
}
