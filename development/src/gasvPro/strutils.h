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

#ifndef _STRUTILS_H
#define _STRUTILS_H

#include <iostream>
#include <string>
using namespace std;

void ToLower(string & s) {
	int len = s.length();
	for(int k=0; k < len; k++){
		s[k] = tolower(s[k]);
	}
}

void ToUpper(string & s){
    int len = s.length();
    for(int k=0; k < len; k++){
        s[k] = toupper(s[k]);
    }
}

void StripPunc(string & word){
    int first = 0;
    int len = word.length();
    while (first < len && ispunct(word[first])){ first++;}
    // assert: first indexes either '\0' or non-punctuation
    
    // now find last non-nonpunctuation character
    int last = len - 1;  // last char in s
	while(last >= 0 && ispunct(word[last])){last--;}
    word = word.substr(first,last-first+1);
}

void StripWhite(string & word){
    int first = 0;
    int len = word.length();
	while (first < len && isspace(word[first])){ first++;}
    // assert: first indexes either '\0' or non-punctuation
	// now find last non-nonpunctuation character
    int last = len - 1;  // last char in s
	while(last >= 0 && isspace(word[last])){last--;}
    word = word.substr(first,last-first+1);
}

// Return lowercase equivalent of s
string UpperString(const string & s){
    string word = s;
    ToUpper(word);
    return word;
}

//Returns lowercase equivalent of s
string LowerString(const string & s){
	string word = s;
    ToLower(word);
    return word;
}

int atoi(const string & s){ return atoi(s.c_str()); } // returns int equivalent
double atof(const string & s){return atof(s.c_str());} // returns double equivalent

// returns string equivalent
string itoa(int n){
    ostringstream output;
    output << n;   
    return output.str();
}

//convert int to string
string tostring(int n){ return itoa(n);}

//convert double to string
string tostring(double n){
    ostringstream output;
    output << n;   
    return output.str();
}
    
#endif 
