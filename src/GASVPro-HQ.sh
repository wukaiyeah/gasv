#!/bin/sh

#### GASVPro High Quality #####

###############################
##### SET PARAMETERS HERE #####
###############################

#REQUIRED:
BAMFILE=/home/tonyc/gasv-read-only/example.bam   
GASVDIR=/home/tonyc/gasv                 

#OPTIONAL (set to NULL if not being used):
UNIQUEFILE=NULL
MAXUNIQUEVAL=NULL               #MUST SPECIFY IF UNIQUEFILE IS GIVEN
MINSCALEDUNIQUE=NULL            #MUST SPECIFY IF UNIQUEFILE IS GIVEN
LRTHRESHOLD=NULL                #default 0
MINCLUSTER=NULL                 #default 4
MAXIMAL=TRUE			#use GASV's --maximal flag. (use TRUE or FALSE)


###############################
#DO NOT MODIFY BELOW THIS POINT
###############################

### Check Input ###

echo "\n** GASVProHQ Running... **\n" 
echo "====Input===="
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
java -jar -Xms512m -Xmx4g $GASVDIR/bin/BAMToGASV.jar $BAMFILE -GASVPRO true -LIBRARY_SEPARATED all

### Run GASV ###

if [ -r $BAMFILE.gasv.in ]; then
    echo "====================\n\n *** Running GASV....*** \n\n===================="
else
    echo "\n\n!! ERROR: necessary parameters file \"$BAMFILE.gasv.in\" does not exist. Ensure BAMToGASV ran correctly and restart.\n"
    exit 1
fi

if [ "$MINCLUSTER" = "NULL" ]; then
    if [ "$MAXIMAL" = "TRUE" ]; then
	java -jar -Xms512m -Xmx4g $GASVDIR/bin/GASV.jar  --output regions --maximal --batch $BAMFILE.gasv.in
    else
	java -jar -Xms512m -Xmx4g $GASVDIR/bin/GASV.jar --output regions --batch $BAMFILE.gasv.in
    fi
else
    if [ "$MAXIMAL" = "TRUE" ]; then
	java -jar -Xms512m -Xmx4g $GASVDIR/bin/GASV.jar --output regions --maximal --minClusterSize $MINCLUSTER --batch $BAMFILE.gasv.in
    else
    	java -jar -Xms512m -Xmx4g $GASVDIR/bin/GASV.jar  --output regions --minClusterSize $MINCLUSTER --batch $BAMFILE.gasv.in
    fi
fi

### Run GASV-CC ###

if [ -r $BAMFILE.gasvpro.in ]; then
    echo "====================\n\n *** Running GASV-CC with the following parameters...***"
else
    echo "\n\n!! ERROR: necessary parameters file \"$BAMFILE.gasvpro.in\" does not exist. Ensure BAMToGASV ran correctly and restart.\n"
    exit 1
fi

if [ -r $UNIQUEFILE ]; then
    echo "UNIQUEFile: $UNIQUEFILE" >> $BAMFILE.gasvpro.in
    echo "MaxUniqueValue: $MAXUNIQUEVAL" >> $BAMFILE.gasvpro.in
    echo "MinScaledUniqueness: $MINSCALEDUNIQUE" >> $BAMFILE.gasvpro.in
fi

if [ "$LRTHRESHOLD" != "NULL" ]; then
    echo "LRThreshold: $LRTHRESHOLD" >> $BAMFILE.gasvpro.in
fi

if [ "$MAXIMAL" = "FALSE" ]; then
    echo "maxmode: true" >> $BAMFILE.gasvpro.in
fi

cat $BAMFILE.gasvpro.in
echo "====================\n\n"

$GASVDIR/bin/GASVPro-CC $BAMFILE.gasvpro.in $BAMFILE.gasv.in.clusters

###Prune Clusters###

echo "\n===================================\n\n *** Pruning Clusters... *** \n\n===================================\n"

$GASVDIR/scripts/GASVPruneClusters.pl $BAMFILE.gasv.in.clusters.GASVPro.clusters


echo "===================\n\n GASVPro-HQ Completed Successfully \n\n===================\n"









