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

#ifndef ESPMAPPING_H
#define ESPMAPPING_H

#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <strings.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <limits.h>


#include <iostream>
#include <sstream>
#include <string>
using namespace std;



class espmapping{
	
	public:
		espmapping(){
			fragName = "ZZZZZtmpName";
			index = mismatchL = mismatchR = positionL = positionR = INT_MAX;
			fragNum = fragMapNum = INT_MAX;
		}
	
		int setEqual(espmapping rhs){
			setIndex(rhs.getIndex()); 
			setFrag(rhs.getFrag());
			setFragNum(rhs.getFragNum());
			setFragMapNum(rhs.getFragMapNum());
			setMismatchL(rhs.getMismatchL()); 
			setMismatchR(rhs.getMismatchR()); 
			setPositionL(rhs.getPositionL()); 
			setPositionR(rhs.getPositionR()); 
			return 0;
		}
	
		const espmapping& operator=(const espmapping &rhs){
			if( &rhs != this){
				setIndex(rhs.index); 
				setFrag(rhs.fragName);
				setFragNum(rhs.fragNum);
				setFragMapNum(rhs.fragMapNum);
				setMismatchL(rhs.mismatchL); 
				setMismatchR(rhs.mismatchR); 
				setPositionL(rhs.positionL); 
				setPositionR(rhs.positionR); 
			}
			return *this;
		}
	
		friend istream&operator>>(istream &input, espmapping& e){
			input >> e.fragNum >> e.fragMapNum >> e.index >> e.fragName >> e.mismatchL >> e.mismatchR >> e.positionL >> e.positionR;
			return input;
		}
	
		friend ostream& operator<<(ostream & output, const espmapping& e){
			output << e.fragNum << " " << e.fragMapNum << " " << e.index << " " << e.fragName << " " << e.mismatchL << " " << e.mismatchR << " " << e.positionL << " " << e.positionR;
			return output;
		}
	
		bool operator==(const espmapping &rhs) const{
			if(fragNum!=rhs.fragNum ||fragMapNum!=rhs.fragMapNum || index!=rhs.index || 
				fragName!=rhs.fragName || mismatchL!=rhs.mismatchL || mismatchR!=rhs.mismatchR ||
				positionL!=rhs.positionL || positionR!=rhs.positionR){ return false; }
			else{ return true; }
		}
		bool operator!=(const espmapping&rhs) const{ return !(*this == rhs);}
	
		int getIndex(){return index;}
		int setIndex(int val){ index = val; return 0; }

		string getFrag(){ return fragName;}
		int setFrag(string s){ fragName = s; return 0;}
	
		int getFragNum(){ return fragNum; }
		int setFragNum(int val){ fragNum = val; return 0; }
		int getFragMapNum(){ return fragMapNum; }
		int setFragMapNum(int val){ fragMapNum = val; return 0; }
		
		int getMismatchL(){ return mismatchL; }
		int setMismatchL(int val){ mismatchL = val; return 0;}
		int getMismatchR(){ return mismatchR; }
		int setMismatchR(int val){ mismatchR = val; return 0;}
		int getPositionL(){ return positionL; }
		int setPositionL(int val){ positionL = val; return 0;}
		int getPositionR(){ return positionR; }
		int setPositionR(int val){ positionR = val; return 0;}
	
		//Example:
		//chr17_41883_42105_0:0:0_2:0:0_15789d3	17	41532	41582	+	17	42573	42623	-
		int defineFromMappingFile(string s){
			istringstream iss(s);
			string tmpString;
			int tmp;
			iss >> tmpString;
			defineFromFragName(tmpString);
			//      17	   41532	  41582	  +	     17	     42573	    42623	-
			iss >> tmp >> positionL >> tmp >> tmp >> tmp >> positionR >> tmp >> tmp;
			return 0;
		}
	
		//BEGIN --> Define from fragName;
		int defineFromFragName(string s){
			int retVal = 0;
			
			fragName = s;
			
			if(s.at(s.length()-1)==','){ s = s.substr(0,s.length()-1); }
			int slen = s.length();
			int c = slen-1;
			int maxLen = slen-1;
			while(c>=0 && s.at(c) != '.'){ c--; }
			if(c < 0){ mismatchL = mismatchR = 0; return retVal; }
			else{
				int valL, valR;
				valL = valR = 0;
				int dashCount,numCount,expectedNumCount,index;
				dashCount = numCount = expectedNumCount = 0;
				index = 1;
				for(int i = maxLen; i>=c+1; i--){
					if(s.at(i) == '_'){ dashCount++; index = 1;} else{ expectedNumCount++; }
					if(s.at(i) == '0' || s.at(i) == '1' || s.at(i) == '2' || s.at(i) == '3' || s.at(i) == '4' || s.at(i) == '5' || s.at(i) == '6' || s.at(i) == '7' || s.at(i) == '8' || s.at(i) == '9'){ 
						if(dashCount == 0){	valR += atoi(s.substr(i,i+1).c_str())*index; }
						else if(dashCount == 1){ valL += atoi(s.substr(i,i+1).c_str())*index; }
						numCount++; 
						index=index*10;
					}
				}
				//Are we ambiguous??
				//YES:
				if(dashCount == 2 && numCount == expectedNumCount){ mismatchL = valL; mismatchR = valR; return 1;}
				//NO: Return the full string;
				else{  mismatchL = mismatchR = 0; return retVal; }
			}	
		}
		//END --> Define from fragName;
	
		string fragName; //Must be unique
		int fragNum;     //Corresponds to the output in the clusters file;
		int fragMapNum;  //Corresponds to the output in the clusters file;
		int index;       //Must be unique; corresponds to the output in the clusters file;
	
		
	
	private:

		int mismatchL;
		int mismatchR;
		int positionL;
		int positionR;
};

//#include "espmapping.cpp"
    
#endif 
