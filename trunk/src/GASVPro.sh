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

########### GASVPro ###########

###############################
##### SET PARAMETERS HERE #####
###############################

#REQUIRED:
BAMFILEHQ=/home/tonyc/gasv/gasvprotest/wgsim.bam
BAMFILELQ=/home/tonyc/gasv/gasvprotest/wgsim.bam.novo.bam
WORKINGDIR=/home/tonyc/gasv/gasvprotest/outputdir  			#GIVE FULL PATH!
GASVDIR=/home/tonyc/gasv
                  

#OPTIONAL (set to NULL if not being used):
UNIQUEFILE=NULL
MAXUNIQUEVAL=NULL               #MUST SPECIFY IF UNIQUEFILE IS GIVEN
MINSCALEDUNIQUE=NULL            #MUST SPECIFY IF UNIQUEFILE IS GIVEN
LRTHRESHOLD=NULL                #default 0
MINCLUSTER=1                 #default 4
MAXIMAL=TRUE			#use GASV's --maximal flag. (use TRUE or FALSE)


###############################
#DO NOT MODIFY BELOW THIS POINT
###############################

### Check Input ###

echo "\n** GASVPro Running... **\n" 
echo "====Input===="
if [ -d $GASVDIR ]; then
	echo "      GASVDIR...ok" 
else
	echo "!! ERROR: Please configure this script with the correct GASV directory."
	exit 1
fi

if [ -r $BAMFILEHQ ]; then
	echo "      $BAMFILEHQ...OK"
else
	echo "!! ERROR: first .bam file is not readable. Aborting."
	exit 1
fi

if [ -r $BAMFILELQ ]; then
	echo "      $BAMFILELQ...OK"
else
	echo "!! ERROR: second .bam file is not readable. Aborting."
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
    echo "      No LR Threshold Provided."
else
    echo "      LR Threshold: " $LRTHRESHOLD
fi
if [ "$MINCLUSTER" = "NULL" ]; then
    echo "      No Min Cluster Size provided."
else
    echo "      Min Cluster Size: " $MINCLUSTER
fi
echo "=============\n"

### Run BAMToGASV and BAMToGASV_AMBIG ###

echo "===================================\n\n *** Running BAMToGASV....*** \n\n===================================\n"
java -Xms512m -Xmx2048m -jar $GASVDIR/bin/BAMToGASV.jar $BAMFILEHQ -GASVPRO true -LIBRARY_SEPARATED all -OUTPUT_PREFIX BAMToGASV

echo "===================================\n\n *** Running BAMToGASV_AMBIG....*** \n\n===================================\n"
java -jar $GASVDIR/bin/BAMToGASV_AMBIG.jar $BAMFILELQ BAMToGASV.gasv.in -OUTPUT_PREFIX BAMToGASV_AMBIG

### Run GASV ### 

if [ -r BAMToGASV_AMBIG.gasv.combined.in ]; then
    echo "====================\n\n *** Running GASV....*** \n\n===================="
else
    echo "\n\n!! ERROR: necessary parameters file \"BAMToGASV_AMBIG.gasv.combined.in\" does not exist. Ensure BAMToGASV ran correctly and restart.\n"
    exit 1
fi

if [ "$MINCLUSTER" = "NULL" ]; then
    if [ "$MAXIMAL" = "TRUE" ]; then
	java -jar -Xms512m -Xmx2048m $GASVDIR/bin/GASV.jar --nohead --output regions --maximal --batch BAMToGASV_AMBIG.gasv.combined.in
    else
	java -jar -Xms512m -Xmx2048m $GASVDIR/bin/GASV.jar --nohead --output regions --batch BAMToGASV_AMBIG.gasv.combined.in
    fi
else
    if [ "$MAXIMAL" = "TRUE" ]; then
	java -jar -Xms512m -Xmx2048m $GASVDIR/bin/GASV.jar --nohead --output regions --maximal --minClusterSize $MINCLUSTER --batch BAMToGASV_AMBIG.gasv.combined.in
    else
    	java -jar -Xms512m -Xmx2048m $GASVDIR/bin/GASV.jar --nohead --output regions --minClusterSize $MINCLUSTER --batch BAMToGASV_AMBIG.gasv.combined.in
    fi
fi

### Run GASVPro-CC ###

if [ -r BAMToGASV.gasvpro.in ]; then
    echo "====================\n\n *** Running GASV-CC with the following parameters...***"
else
    echo "\n\n!! ERROR: necessary parameters file \"BAMToGASV.gasvpro.in\" does not exist. Ensure BAMToGASV ran correctly and restart.\n"
    exit 1
fi

if [ -r $UNIQUEFILE ]; then
    echo "UNIQUEFile: $UNIQUEFILE" >> BAMToGASV.gasvpro.in
    echo "MaxUniqueValue: $MAXUNIQUEVAL" >> BAMToGASV.gasvpro.in
    echo "MinScaledUniqueness: $MINSCALEDUNIQUE" >> BAMToGASV.gasvpro.in
fi

if [ "$LRTHRESHOLD" != "NULL" ]; then
    echo "LRThreshold: $LRTHRESHOLD" >> BAMToGASV.gasvpro.in
fi

if [ "$MAXIMAL" = "FALSE" ]; then
    echo "maxmode: true" >> BAMToGASV.gasvpro.in
fi

cat BAMToGASV.gasvpro.in
echo "====================\n\n"

$GASVDIR/bin/GASVPro-CC BAMToGASV.gasvpro.in BAMToGASV_AMBIG.gasv.combined.in.clusters


### Run GASVPro-graph ###

if [ -r BAMToGASV_AMBIG.gasv.combined.in.clusters.GASVPro.clusters -a -r BAMToGASV_AMBIG.gasv.combined.in.clusters.GASVPro.coverage ]; then
	echo "====================\n\n *** Running GASVPro-Graph....*** \n\n===================="
else
	echo "!! ERROR: Check that GASVPro-CC completed successfully and produced BAMToGASV_AMBIG.gasv.combined.in.clusters.GASVPro.clusters and BAMToGASV_AMBIG.gasv.combined.in.clusters.GASVPro.coverage files."
	exit 1
fi



$GASVDIR/bin/GASVPro-graph BAMToGASV_AMBIG.gasv.combined.in.clusters.GASVPro.clusters BAMToGASV_AMBIG.gasv.combined.in.clusters.GASVPro.coverage $WORKINGDIR

### Run GASVPro-mcmc ###

echo "====================\n\n *** Running GASVPro-mcmc....*** \n\n===================="

$GASVDIR/bin/GASVPro-mcmc BAMToGASV.gasvpro.in $WORKINGDIR

echo "====================\n\n *** GASVPro complete *** \n\n===================="


