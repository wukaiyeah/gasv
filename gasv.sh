#!/bin/bash
##this script is aimed to call strcuture viaration with delly software

DATA_DIR=./

GASVProHQ=./GASVPro-HQ.multifiles.sh


for ALL in $DATA_DIR/*.rmdup.bam
do
	SAMPLE_NAME=${ALL%%.rmdup.bam}
	SAMPLE_NAME=${SAMPLE_NAME##*/}
	echo $SAMPLE_NAME
##--usage sh $GASVProHQ {bam file direction} {sample name}
	sh $GASVProHQ $ALL $SAMPLE_NAME

done

