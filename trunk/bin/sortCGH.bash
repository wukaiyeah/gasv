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
# ./sortCGH.bash someFile
#Output:
# someFile.sorted
if [[ -e $1 ]]
then
	echo "Sorting $1"
else
	echo "File $1 does not exist"
	echo "usage: ./sortCGH.bash cghFileToSort"
	exit 
fi

# First change any chr2 or chrx into just "2" and "x", respectively.
# Next change any X or Y to 23 or 24, respectively
# Add a 0 in front of any single digits for sorting purposes
# Finally re-order so that chr and pos columns are first
cat $1 |  awk '{if (length($3)==4) $3=substr($3,4,1);if (length($3)==5) $3=substr($3,4,2);if ($3=="X" || $3=="x") $3=23; if ($3=="Y" || $3=="y") $3=24; if (length($3)==1) $3=0$3; print $3"\t"$4"\t"$5"\t"$1"\t"$2}' > $1.tmp.toSort

# Do a sort
sort -n -k 1,1 -k 2,2 -k 3,3 $1.tmp.toSort > $1.tmp.toSort.sorted

# Remove the 0 in front of any single digit and put columns back in correct order
cat $1.tmp.toSort.sorted | awk '{if (substr($1,1,1)==0) $1=substr($1,2,1); print $4"\t"$5"\t"$1"\t"$2"\t"$3}'> $1.sorted
rm -f $1.tmp.toSort
rm -f $1.tmp.toSort.sorted 

