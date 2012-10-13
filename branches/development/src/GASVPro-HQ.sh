#!/bin/sh

####
###
## Copyright 2012 Benjamin Raphael, Suzanne Sindi, Anthony Cannistra, Hsin-Ta Wu, Luke Peng, Selim Onal
##
##  This file is part of the GASVPro code distribution.
##
##  GASVPro is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
## 
##  GASVPro is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public License
##  along with gasv.  If not, see <http://www.gnu.org/licenses/>.
###
####

#### GASVPro High Quality #####

###############################
##### SET PARAMETERS HERE #####
###############################

#REQUIRED:
BAMFILE= ##BAMFILE
GASVDIR= ##GASVDIRECTORY               

#OPTIONAL (set to NULL if not being used):
UNIQUEFILE= ##UNIQUENESSFILE
MAXUNIQUEVAL=NULL               #MUST SPECIFY IF UNIQUEFILE IS GIVEN
MINSCALEDUNIQUE=NULL            #MUST SPECIFY IF UNIQUEFILE IS GIVEN
LRTHRESHOLD=NULL                #default 0
MINCLUSTER=NULL                 #default 4
MAXIMAL=FALSE			#use GASV's --maximal flag. (use TRUE or FALSE)


###############################
#DO NOT MODIFY BELOW THIS POINT
###############################
DATEPREFIX=$BAMFIILE_$(date +%F_%R)
### Check Input ###

echo "\n** GASVProHQ Running... **\n" 
echo "====Input===="
if [ -d $GASVDIR ]; then
	echo "      GASVDIR...ok" 
else
	echo "!! ERROR: Please configure this script with the correct GASV directory."
	exit 1
fi
if [ -r $BAMFILE ]; then
    echo "      Valid .bam file: $BAMFILE"
else
    echo "      !! WARNING: Cannot open .BAM FILE. Aborting.\n"
    exit 1
fi
if [ -r $UNIQUEFILE ]; then
    echo "      Unique file: $UNIQUEFILE"
    if [ "$MAXUNIQUEVAL" != "NULL" -a "$MINSCALEDUNIQUE" != "NULL" ]; then
    	echo "          Max Uniqueness Value: $MAXUNIQUEVAL"
    	echo "          Min Scaled Uniqueness: $MINSCALEDUNIQUE"
    else
    	echo "      !! ERROR: Max uniqueness value and minimum scaled uniqueness value must be provided. Please edit your parameters and restart.\n"
    	exit 1
    fi
else
    echo "      No unique file Provided."
fi





if [ "$LRTHRESHOLD" = "NULL" ]; then
    echo "      No LR Threshold Provided"
else
    echo "      LR Threshold: " $LRTHRESHOLD
fi
if [ "$MINCLUSTER" = "NULL" ]; then
    echo "      No Min Cluster Size provided"
else
    echo "      Min Cluster Size: " $MINCLUSTER
fi
echo "=============\n"

echo $S1 $S2



### Run BAMtoGASV ###

echo "===================================\n\n *** Running BAMToGASV....*** \n\n===================================\n"
java -jar -Xms512m -Xmx2g $GASVDIR/bin/BAMToGASV.jar $BAMFILE -GASVPRO true -LIBRARY_SEPARATED all -OUTPUT_PREFIX $DATEPREFIX
OUT=$?
if [ "$OUT" != 0 ]; then
	echo "Warning: BAMToGASV aborted. Stopping pipeline."
	exit  1
fi 

### Run GASV ###

if [ -r $DATEPREFIX.gasv.in ]; then
    echo "====================\n\n *** Running GASV....*** \n\n===================="
else
    echo "\n\n!! ERROR: necessary parameters file \"$DATEPREFIX.gasv.in\" does not exist. Ensure BAMToGASV ran correctly and restart.\n"
    exit 1
fi

if [ "$MINCLUSTER" = "NULL" ]; then
    if [ "$MAXIMAL" = "TRUE" ]; then
	java -jar -Xms512m -Xmx2g $GASVDIR/bin/GASV.jar  --output regions --maximal --batch $DATEPREFIX.gasv.in
    else
	java -jar -Xms512m -Xmx2g $GASVDIR/bin/GASV.jar --output regions --batch $DATEPREFIX.gasv.in
    fi
else
    if [ "$MAXIMAL" = "TRUE" ]; then
	java -jar -Xms512m -Xmx2g $GASVDIR/bin/GASV.jar --output regions --maximal --minClusterSize $MINCLUSTER --batch $DATEPREFIX.gasv.in
    else
    	java -jar -Xms512m -Xmx2g $GASVDIR/bin/GASV.jar  --output regions --minClusterSize $MINCLUSTER --batch $DATEPREFIX.gasv.in
    fi
fi

OUT=$?
if [ "$OUT" != 0 ]; then
	echo "Warning: GASV aborted. Stopping pipeline."
	exit  1
fi 

### Run GASV-CC ###

if [ -r $DATEPREFIX.gasvpro.in ]; then
    echo "====================\n\n *** Running GASVPro-CC with the following parameters...***"
else
    echo "\n\n!! ERROR: necessary parameters file \"$DATEPREFIX.gasvpro.in\" does not exist. Ensure BAMToGASV ran correctly and restart.\n"
    exit 1
fi

if [ -r $UNIQUEFILE ]; then
    echo "UNIQUEFile: $UNIQUEFILE" >> $DATEPREFIX.gasvpro.in
    echo "MaxUniqueValue: $MAXUNIQUEVAL" >> $DATEPREFIX.gasvpro.in
    echo "MinScaledUniqueness: $MINSCALEDUNIQUE" >> $DATEPREFIX.gasvpro.in
fi

if [ "$LRTHRESHOLD" != "NULL" ]; then
    echo "LRThreshold: $LRTHRESHOLD" >> $DATEPREFIX.gasvpro.in
fi

if [ "$MAXIMAL" = "FALSE" ]; then
    echo "maxmode: true" >> $DATEPREFIX.gasvpro.in
fi

cat $DATEPREFIX.gasvpro.in
echo "====================\n\n"

$GASVDIR/bin/GASVPro-CC $DATEPREFIX.gasvpro.in $DATEPREFIX.gasv.in.clusters

OUT=$?
if [ "$OUT" != 0 ]; then
	echo "Warning: GASVPro-cc aborted. Stopping pipeline."
	exit  1
fi 

###Prune Clusters###

echo "\n===================================\n\n *** Pruning Clusters... *** \n\n===================================\n"

$GASVDIR/scripts/GASVPruneClusters.pl $DATEPREFIX.gasv.in.clusters.GASVPro.clusters


echo "===================\n\n GASVPro-HQ Completed Successfully \n\n===================\n"








