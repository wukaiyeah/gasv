#!/usr/bin/perl
use strict;

#Goal: *Prune a set of clusters in the BreakDancer Likelihood Metric.
#      *Breakdancer Metric --> The intersection of two variants is > 50% of the union of the two variants.
#      
#      What about pruning a set of clusters AND concordant coverage information? (Wait until later!)

my $SEED = time^$$;
srand($SEED);


if($#ARGV<0){
	die "Usage: ./GASVPruneClusters.pl {GASVClusters}\n";
}

#This is the number of tab separated fields in a original Gasv clusters file.
my $NUM_NORMAL_ARGS = 8;
#This the position of the likelihood value in the clusters file;
my $LIKELIHOOD_ARG = 8;

open(F,$ARGV[0]) or die "Could not open clusters file $ARGV[0]\n";

###Set $DEBUG = "true" to print usage information;
my $DEBUG = "true";
$DEBUG = "false";

my $METHOD = $ARGV[1];
if(length($METHOD)<=0){ 
	$METHOD = "LIKE";
}

print "The selected method is:\t$METHOD\n";

if($METHOD eq "ESP" or $METHOD eq "LIKE"){ }
else{ die "PRUNE_METHOD must be either ESP or LIKE\n";}

my $OUTPUT_FILE; 
if($ARGV[1] eq "ESP"){ $OUTPUT_FILE = $ARGV[0].".pruned.clusters"; }
else{ $OUTPUT_FILE = $ARGV[0].".pruned.clusters"; }
open(O,">".$OUTPUT_FILE);

#Example Line;
my @printFlag;
my @variants;
my $numVariants = 0;
my @type;
my @names;

my $HEADER = "";
my $totalVariants = 0;
while(my $line=<F>){
	chomp($line);
	
	my $firstChar = substr($line,0,1);
	if($firstChar eq "#"){
		$HEADER = $line;
	}
	else{
		my @args=split(/\t/,$line); 
		my $loc = $args[2];
	
		my $flag = 1;
		if($METHOD eq "ESP"){ $flag = 1;}
		else{
			$flag = 1;
			#if(getLike($line)>=0){ $flag = 1;}
			#else{ $flag = 0; }
		}
 
		if($loc>0 and $flag == 1){
			$names[$numVariants] = $args[0];
			$type[$numVariants] = substr($args[3],0,1);
			$variants[$numVariants] = $line;
			$printFlag[$numVariants] = 0;
			$numVariants++;
		}
		$totalVariants++;
	}
}
close(F);

my $varToRemove;
for(my $i = 0; $i<$numVariants; $i++){
	if($DEBUG eq "true" or $i%10 == 0){
		print "****************\n";
		print "Processing Variant $i $names[$i]\n";
	}
	if($printFlag[$i] == 0){
		my $origType = $type[$i];
		for(my $j = $i+1; $j<$numVariants; $j++){
			my $newType = $type[$j];
			if($printFlag[$j] == 0 and $printFlag[$i] == 0 and ($origType eq $newType) ){
				my $result = overlapVariant($variants[$i],$variants[$j]);
				#The variants are dependent; keep the one with most support;
				if($result>0){
					if($DEBUG eq "long"){ 
						print "Variant $names[$i] overlaps $names[$j]\n";
						print "Print Flag:\n";
						print "\t$names[$i] $printFlag[$i]\n";
						print "\t$names[$j] $printFlag[$j]\n";
					}
					my $supportResult;
					if($METHOD eq "ESP"){
						$supportResult = biggestSupport($variants[$i],$variants[$j]);
					}
					else{ 
						$supportResult = biggestLike($variants[$i],$variants[$j]);
					}
					#Same support; need to keep one at random!
					if($supportResult == -1){ 
						my $range = 1;
						$supportResult = rand($range);
					}
					if($supportResult <= 0.5){ $varToRemove = $j; }
					else{ $varToRemove = $i; }
			
					my $tmpVar = $variants[$varToRemove];
					if($DEBUG eq "long"){ print "\tRemoving $names[$varToRemove]\n";}
				
					$printFlag[$varToRemove]++;
					if($DEBUG eq "long"){ print "\tPrint Flag: $printFlag[$varToRemove]\n";}
				}
			}
		}
	}
}

if(length($HEADER)>0){ print O "$HEADER\n"; }
for(my $i = 0; $i<$numVariants; $i++){
	if($printFlag[$i] == 0){
		print O "$variants[$i]\n";
	}
}

sub getLike{
	my $val = $_[0];
	my @args = split(/\t/,$val);
	my @lvals = split(/\_/,$args[0]);
	my $real  = split(/\_/,$args[0]);
	
	if($real<=1){
		die "File does not seem to contain likelihoods:\n$val\n";
	}

	my $LLR = $lvals[1] - $lvals[2];
	
	return $LLR;
}

#but prob(variant) and prob(no variant) attached to cluster names
sub biggestLike{
	my $val1 = $_[0];
	my $val2 = $_[1];
	
	my @args1 = split(/\t/,$val1);
	my @args2 = split(/\t/,$val2);
	
	my $numArgs1 = split(/\t/,$val1);
	my $numArgs2 = split(/\t/,$val2);
	
	if($numArgs1 <= $NUM_NORMAL_ARGS or $numArgs2<=$NUM_NORMAL_ARGS){
		die "Clusters file does not seem to have likelihoods.\n";
	}
	
	my $llr1 = $args1[$LIKELIHOOD_ARG];
	my $llr2 = $args2[$LIKELIHOOD_ARG];
	
	if($llr1>$llr2){ return 0; }
	elsif($llr1<$llr2){ return 1;}
	else{ return -1; } 
}

sub biggestSupport{
	my $var1 = $_[0];
	my $var2 = $_[1];
	
	my @args1 = split(/\t/,$var1);
	my @args2 = split(/\t/,$var2);
	
	my $sup1 = $args1[1];
	my $sup2 = $args2[1];
	
	if($sup1>$sup2){ return 0; }
	elsif($sup1<$sup2){ return 1; }
	else{ return -1; }
}


#Example: 
#0       1        2      3                         4                      5       6            7
#c1      1       189     D       chr17_41883_42105_0:0:0_2:0:0_15789d3   17      17      41642, 42623, 41821, 42623, 41532, 42334, 41532, 42513

sub overlapVariant{
	my $var1 = $_[0];
	my $var2 = $_[1];

	my $retVal = 0;
	
	my @args1 = split(/\t/,$var1);
	my @args2 = split(/\t/,$var2);
	
	my $chr1L = $args1[5];
	my $chr1R = $args1[6];
	my $chr2L = $args2[5];
	my $chr2R = $args2[6];	
	
	##Note: The below is ONLY valid in the case of chrL = chrR
	##      Need a different criteria for translocations!
	
	if( ($chr1L eq $chr2L) and ($chr1R eq $chr2R)){
		#Step 1: Get Intervals;
		my $coordList1 = $args1[7];
		$coordList1 =~s/\,//g;
		my $coordList2 = $args2[7];
		$coordList2 =~s/\,//g;
	
		my @coords1    = split(/\s+/,$coordList1);
		my $numCoords1 = split(/\s+/,$coordList1);
		my @coords2    = split(/\s+/,$coordList2);
		my $numCoords2 = split(/\s+/,$coordList2);

        my $minX1 = $coords1[0];
		my $maxX1 = $coords1[0];
		my $minY1 = $coords1[1];
        my $maxY1 = $coords1[1];
		for(my $i = 2; $i<$numCoords1; $i = $i+2){
            if($coords1[$i]  < $minX1){ $minX1 = $coords1[$i];   }
			if($coords1[$i]  > $maxX1){ $maxX1 = $coords1[$i];   }
			if($coords1[$i+1]< $minY1){ $minY1 = $coords1[$i+1]; }
            if($coords1[$i+1]> $maxY1){ $maxY1 = $coords1[$i+1]; }
		}
        
        my $minX2 = $coords2[0];
		my $maxX2 = $coords2[0];
		my $minY2 = $coords2[1];
        my $maxY2 = $coords2[1];
		for(my $i = 2; $i<$numCoords2; $i = $i+2){
            if($coords2[$i]  < $minX2){ $minX2 = $coords2[$i];   }
			if($coords2[$i]  > $maxX2){ $maxX2 = $coords2[$i];   }
			if($coords2[$i+1]< $minY2){ $minY2 = $coords2[$i+1]; }
            if($coords2[$i+1]> $maxY2){ $maxY2 = $coords2[$i+1]; }
		}
		
		#Step 2: Do the intervals overlap?
        #We are Inversion/Deletions/Divergent (On the same chromosome)
        if($chr1L eq $chr1R){
            #NO
            if($minY2 < $maxX1 or $minY1 < $maxX2){ }
            #YES
            else{
                #Do we have containment?
                if($maxX1 < $maxX2 and $minY1 > $minY2){ $retVal = 1; }
                elsif($maxX2 < $maxX1 and $minY2 < $minY1){ $retVal = 1; }
			
                #THEN FIGURE OUT THE UNION AND INTERSECTION!
                my $INT_X;
                my $UNION_X;
                if($maxX1 < $maxX2){ $UNION_X = $maxX1; $INT_X = $maxX2;}
                else{ $UNION_X = $maxX2; $INT_X = $maxX1; }
			
                my $INT_Y;
                my $UNION_Y;
                if($minY1 < $minY2){ $UNION_Y = $minY2; $INT_Y = $minY1; }
                else{ $UNION_Y = $minY1; $INT_Y = $minY2; }
			
                my $UNION_LEN = $UNION_Y - $UNION_X + 1;
                my $INT_LEN   = $INT_Y - $INT_X + 1;
			
                if($INT_LEN >=0.5*$UNION_LEN){ $retVal = 1; }
            }
        }
        #We are Translocaitons (different chromosomes)
        else{
            #Do we overlap in X?
            if( ($maxX1 < $minX2) || ($maxX2 < $minX1) ){
                #NO:
            }
            else{
                #YES:
                #Do we overlap in Y?
                if( ($maxY1 < $minY2) || ($maxY2 < $minY1) ){
                    #NO
                }
                else{
                    #YES: Overlap in BOTH;
                    $retVal = 1;
                }
            }
        }
	}

	
	return $retVal;
	
}
