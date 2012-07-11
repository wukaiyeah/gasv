#!/bin/sh

########### GASVPro ###########

###############################
##### SET PARAMETERS HERE #####
###############################

#REQUIRED:
BAMFILEHQ=/home/tonyc/gasv/example.bam
BAMFILELQ=/home/tonyc/gasv/protest_smallbam_venterbwa/VenterBWA.AMBI.short.bam
WORKINGDIR=outputdir  
GASVDIR=/home/tonyc/gasv
                  

#OPTIONAL (set to NULL if not being used):
UNIQUEFILE=NULL
MAXUNIQUEVAL=NULL               #MUST SPECIFY IF UNIQUEFILE IS GIVEN
MINSCALEDUNIQUE=NULL            #MUST SPECIFY IF UNIQUEFILE IS GIVEN
LRTHRESHOLD=NULL                #default 0
MINCLUSTER=NULL                 #default 4
MAXIMAL=FALSE			#use GASV's --maximal flag. (use TRUE or FALSE)


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
	echo "      $BAMFILELQ-Xms512m -Xmx2048m...OK"
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
java -Xms512m -Xmx2048m -jar $GASVDIR/bin/BAMToGASV.jar $BAMFILEHQ -GASVPRO true -LIBRARY_SEPARATED all

echo "===================================\n\n *** Running BAMToGASV_AMBIG....*** \n\n===================================\n"
java -jar $GASVDIR/bin/BAMToGASV_AMBIG.jar $BAMFILELQ $BAMFILEHQ.gasv.in

### Run GASV ### 

if [ -r $BAMFILEHQ.gasv.in ]; then
    echo "====================\n\n *** Running GASV....*** \n\n===================="
else
    echo "\n\n!! ERROR: necessary parameters file \"$BAMFILE.gasv.in\" does not exist. Ensure BAMToGASV ran correctly and restart.\n"
    exit 1
fi

if [ "$MINCLUSTER" = "NULL" ]; then
    if [ "$MAXIMAL" = "TRUE" ]; then
	java -jar -Xms512m -Xmx2048m $GASVDIR/bin/GASV.jar --nohead --output regions --maximal --batch $BAMFILELQ.gasv.combined.in
    else
	java -jar -Xms512m -Xmx2048m $GASVDIR/bin/GASV.jar --nohead --output regions --batch $BAMFILELQ.gasv.combined.in
    fi
else
    if [ "$MAXIMAL" = "TRUE" ]; then
	java -jar -Xms512m -Xmx2048m $GASVDIR/bin/GASV.jar --nohead --output regions --maximal --minClusterSize $MINCLUSTER --batch $BAMFILELQ.gasv.combined.in
    else
    	java -jar -Xms512m -Xmx2048m $GASVDIR/bin/GASV.jar --nohead --output regions --minClusterSize $MINCLUSTER --batch $BAMFILELQ.gasv.combined.in
    fi
fi

### Run GASVPro-CC ###

if [ -r $BAMFILEHQ.gasvpro.in ]; then
    echo "====================\n\n *** Running GASV-CC with the following parameters...***"
else
    echo "\n\n!! ERROR: necessary parameters file \"$BAMFILE.gasvpro.in\" does not exist. Ensure BAMToGASV ran correctly and restart.\n"
    exit 1
fi

if [ -r $UNIQUEFILE ]; then
    echo "UNIQUEFile: $UNIQUEFILE" >> $BAMFILEHQ.gasvpro.in
    echo "MaxUniqueValue: $MAXUNIQUEVAL" >> $BAMFILEHQ.gasvpro.in
    echo "MinScaledUniqueness: $MINSCALEDUNIQUE" >> $BAMFILEHQ.gasvpro.in
fi

if [ "$LRTHRESHOLD" != "NULL" ]; then
    echo "LRThreshold: $LRTHRESHOLD" >> $BAMFILEHQ.gasvpro.in
fi

if [ "$MAXIMAL" = "FALSE" ]; then
    echo "maxmode: true" >> $BAMFILEHQ.gasvpro.in
fi

cat $BAMFILEHQ.gasvpro.in
echo "====================\n\n"

$GASVDIR/bin/GASVPro-CC $BAMFILEHQ.gasvpro.in $BAMFILEHQ.gasv.in.clusters


### Run GASVPro-graph ###

$GASVDIR/bin/GASVPro-graph $BAMFILEHQ.gasv.in.clusters $BAMFILEHQ.gasv.in.clusters.GASVPro.coverage $WORKINGDIR



