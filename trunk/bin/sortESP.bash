#!/bin/bash
# Copyright 2010 Benjamin Raphael, Suzanne Sindi, Hsin-Ta Wu, Anna Ritz, Luke Peng
# 
#   This file is part of gasv.
#  
#   gasv is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#  
#   gasv is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with gasv.  If not, see <http://www.gnu.org/licenses/>.
#

#Usage:
# ./sortESP.bash someFile [output_directory]
#Output:
# someFile.sorted
if [[ -e $2 ]]
then
	echo "Use $2 as sorting temporary directory"
fi
if [[ -e $1 ]]
then
	echo "Sorting $1"
else
	echo "File $1 does not exist"
	echo "usage: ./sortESP.bash espFileToSort [output_directory]"
	exit 
fi


# need to "flip" start and end coordinates if a read is negative for sorting purposes
cat $1 |  awk '{r1Start=$3; r1End=$4; r2Start=$7; r2End=$8; if ($5 == "-" || $5 == "MINUS" || $5 == "Minus") $3=r1End; if ($5 == "-" || $5 == "MINUS" || $5 == "Minus") $4=r1Start; if ($9 == "-" || $9 == "MINUS" || $9 == "Minus") $7=r2End; if ($9 == "-" || $9 == "MINUS" || $9 == "Minus") $8=r2Start; print}' >$1.tmp.toSortFlipped

# First change any chr2 or chrx into just "2" and "x", respectively.
# Next change any X or Y to 23 or 24, respectively
# Add a 0 in front of any single digits for sorting purposes
# Finally re-order so that chr and pos columns are first
cat $1.tmp.toSortFlipped |  awk '{if (length($2)==4) $2=substr($2,4,1);if (length($2)==5) $2=substr($2,4,2);if (length($6)==4) $6=substr($6,4,1);if (length($6)==5) $6=substr($6,4,2); if ($2=="X" || $2=="x") $2=23; if ($2=="Y" || $2=="y") $2=24; if ($6=="X" || $6=="x") $6=23; if ($6=="Y" || $6=="y") $6=24; if (length($2)==1) $2=0$2; if (length($6)==1) $6=0$6; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | awk '{print $2"\t"$6"\t"$3"\t"$7"\t"$4"\t"$8"\t"$5"\t"$9"\t"$1}'> $1.tmp.toSort

# Do a sort
#sort -n -k 1,1 -k 2,2 -k 3,3 -k 4,4 tmp.toSort > tmp.toSort.sorted
if [[ -e $2 ]]
then
	sort -T $2 -n -k 1,1 -k 2,2 -k 3,3 -k 4,4 $1.tmp.toSort > $1.tmp.toSort.sorted
else
	sort -n -k 1,1 -k 2,2 -k 3,3 -k 4,4 $1.tmp.toSort > $1.tmp.toSort.sorted
fi


# Remove the 0 in front of any single digit and put columns back in correct order
cat $1.tmp.toSort.sorted | awk '{if (substr($1,1,1)==0) $1=substr($1,2,1); if (substr($2,1,1)==0) $2=substr($2,2,1); print $9"\t"$1"\t"$3"\t"$5"\t"$7"\t"$2"\t"$4"\t"$6"\t"$8}' > $1.tmp.toSortFlipped

# need to "flip" start and end coordinates for negative reads back to their previous order 
cat $1.tmp.toSortFlipped | awk '{r1Start=$3; r1End=$4; r2Start=$7; r2End=$8; if ($5 == "-" || $5 == "MINUS" || $5 == "Minus") $3=r1End; if ($5 == "-" || $5 == "MINUS" || $5 == "Minus") $4=r1Start; if ($9 == "-" || $9 == "MINUS" || $9 == "Minus") $7=r2End; if ($9 == "-" || $9 == "MINUS" || $9 == "Minus") $8=r2Start; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' > $1.sorted
rm -f $1.tmp.toSort
rm -f $1.tmp.toSort.sorted 
rm -f $1.tmp.toSortFlipped

