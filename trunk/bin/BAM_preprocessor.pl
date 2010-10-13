#!/usr/bin/perl

#use strict;
use warnings;
use Getopt::Long;
use LWP::Simple qw();

my $BAMFILE = $ARGV[0];
my $SAMPATH = "samtools";
my $LIBSEP = "separated-library";
my $OUTPUTPREFIX = "NO_SPE";
my $MAPPINGQ = 10;
my $CUTOFF = "PCT=99%";
my $FIRST_READ = 500000;
my $NAMING_INDEXING_FILE = "DEFAULT";
my $PROPER_L = 10000;
my $LMINLMAX_FILE = "";
my $PLATFORM = "";

GetOptions("SAMTOOLS_PATH=s", \$SAMPATH, "LIBRARY_SEPARATED=s", \$LIBSEP, "OUTPUT_PREFIX=s", \$OUTPUTPREFIX, "MAPPING_QUALITY=s", \$MAPPINGQ, "CUTOFF_LMINLMAX=s", \$CUTOFF, "USE_NUMBER_READS=s", \$FIRST_READ, "CHROMOSOME_NAMING=s", \$NAMING_INDEXING_FILE, "PROPER_LENGTH=s", \$PROPER_L, "PLATFORM=s", \$PLATFORM);
die("
	Usage: perl BAM_preprocessor.pl <bam file>
	Options:
	-SAMTOOLS_PATH     STRING  A path of samtools
	-LIBRARY_SEPARATED STRING  Extract fragments from separated libraries or whole BAM file
	-OUTPUT_PREFIX     STRING  The prefix for output files
	-MAPPING_QUALITY   INT     Minimum mapping quality
	-CUTOFF_LMINLMAX   STRING  A lower and upper bounds on the fragment distribution
	-USE_NUMBER_READS  INT     The number of fragments (designated as a proper pair) use to computing Lmax and Lmin
	-CHROMOSOME_NAMING STRING  Specify chromosome naming in this file
	-PROPER_LENGTH     INT     A threshold to consider fragments as a proper pair
	-PLATFORM	   STRING  Specify Illumina or SOLiD as platform if missing PL tag in \@RG line
	\n
") unless (scalar(@ARGV) == 1);

if($OUTPUTPREFIX eq "NO_SPE"){
	$OUTPUTPREFIX = getDefaultPrefix($BAMFILE);
}

CheckRequiredFile($SAMPATH, $BAMFILE, $MAPPINGQ, $CUTOFF);
my $ERRLOG = $OUTPUTPREFIX.".bampreprocessor.err.log".$$.time();
my $debug_mode = 0;
my %LB;
my %LB_LMIN;
my %LB_LMAX;
my %ID_LB = ();
my %PLAT_LB = ();
my %LB_FH = ();
my $SEP_LIB = 0;
my $ESPfile = $OUTPUTPREFIX.".info";
my %LB_CUTOFF = ();
my $MISSING_PL_TAG = 0;
my @MISSING_LB_AY;
my $RG_TAG = 1;
`echo "LibraryName      Lmin    Lmax" > $ESPfile`;

print "\n===================================================================================================\n";
print "PATH OF SAMTOOLS:     $SAMPATH\n";
print "BAM FILE:             $BAMFILE\n";
print "USING DATASETS:       $LIBSEP\n";
print "OUTPUT PREFIX:        $OUTPUTPREFIX\n";
print "MAPPING QUALITY:      $MAPPINGQ\n";
print "CUTOFF SETTING:       $CUTOFF\n";
print "NUMBER OF READS:      $FIRST_READ\n";
if($NAMING_INDEXING_FILE ne "DEFAULT"){
	print "NAMING FILE:          $NAMING_INDEXING_FILE\n";
}
if($CUTOFF =~ /SD/){
	print "MAXIMUM PROPER LENGTH:$PROPER_L\n";
}

if($PLATFORM ne ""){
	if($PLATFORM !~ /^(SOLiD|Illumina)$/i){
		print "ERROR: Please specify your platform as Illumina or SOLiD.\n";
		exit(1);
	}
	print "PLATFORM:             ".$PLATFORM."\n";
}
else{
	$PLATFORM = ((readHeaderInfo() eq "illumina")?"Illumina":"SOLiD");
	print "PLATFORM:             ".$PLATFORM."\n";
}
print "===================================================================================================\n\n";

$FIRST_READ = $FIRST_READ * 2;
if ($MISSING_PL_TAG == 1){
	print STDERR "WARNING: Missing PL tag in the header. Using 'Illumina' as platform instead. \n";
	print STDERR "         Or you can specify the platform as Illumina or SOLiD using -PLATFORM. \n";
}
if (scalar(@MISSING_LB_AY) != 0){
	print join("\n", @MISSING_LB_AY)."\n\n";
}
if (scalar(keys %ID_LB) == 0) {
	print STDERR "Warning: no libraries found in BAM header information. Using all mode instead.\n";
	$LB{"all"} = "all";
	if($PLATFORM =~ /solid/i){
		$PLAT_LB{"all"} = "solid";
	}
	elsif($PLATFORM =~ /illumina/i){
		$PLAT_LB{"all"} = "Illumina";
	}
	$SEP_LIB = 0;
}
if ($CUTOFF =~ /FILE/){
	my $cut_file = substr $CUTOFF, 5;
	# read lower and upper bounds for each library
	open(IN, $cut_file);
	while(<IN>){
		chomp;
		(my $l, my $v) = split(/\t+/, $_);
		$LB_CUTOFF{$l} = $v;
	}
	close(IN);
}

my @lib_delete;
foreach my $library_id (keys %LB){
	my $USE_EXACT = 0;
	my $LMIN=0;
	my $LMAX=0;
	my $CHECKING_FILE = $OUTPUTPREFIX."_".$library_id.".bampreprocessor.checking.txt".$$.time();

	if(!CheckPairInfo($CHECKING_FILE, $SAMPATH, $BAMFILE, $library_id, $SEP_LIB, $FIRST_READ)){
		if($library_id eq "all"){
			print STDERR "\nBAM file you're using appears to lack pairing information!\n";
			print STDERR "Try to run 'samtools sort -n' and then 'samtools fixmate' to generate the pairing information! See http://samtools.sourceforge.net/samtools.shtml for details.\n";
			print STDERR "Or set more reads to check pairing information using option -USE_NUMBER_READS <number>.\n";
			system("rm -f $CHECKING_FILE");
			exit 0;
		}
		else {
			print STDERR "\nLibrary $library_id in BAM file you're using appears to lack pairing information!\n";
			print STDERR "Try to run 'samtools sort -n' and then 'samtools fixmate' to generate the pairing information! See http://samtools.sourceforge.net/samtools.shtml for details.\n";
			print STDERR "Or set more reads to check pairing information using option -USE_NUMBER_READS <number>.\n";
			print STDERR "Skip computing $library_id!\n\n";
			system("rm -f $CHECKING_FILE");
			my @id_delete;
			foreach my $id (keys %ID_LB){
				if($ID_LB{$id} eq $library_id){
					push(@id_delete, $id);
				}
			}
			foreach my $id (@id_delete){
				delete $ID_LB{$id};
			}
			push(@lib_delete, $library_id);
		}
		next if ($SEP_LIB == 1);
	}

	if ($CUTOFF =~ /EXACT=/) {
		$USE_EXACT = 1;
		my $LMINLMAX = substr $CUTOFF, 6;
		my @tmpLMINLMAX = split(/,/, $LMINLMAX);
		$LMIN = $tmpLMINLMAX[0];
		$LMAX = $tmpLMINLMAX[1];
		print "EXACT flag: $CUTOFF used, so using Lmin of ".$LMIN." and Lmax of ".$LMAX.".  Will NOT check distribution of concordant pairs!!\n"; 
	}

	elsif ($CUTOFF =~ /SD=/ || $CUTOFF =~/PCT=/){
		print "Computing Lmin and Lmax ...\n";
		($LMIN, $LMAX) = GetLminLmax($CHECKING_FILE, $library_id, $MAPPINGQ, $CUTOFF, $FIRST_READ, $ESPfile, $PROPER_L);
	}
	elsif ($CUTOFF =~ /FILE=/) {
		if(exists $LB_CUTOFF{$library_id}){
			if ($LB_CUTOFF{$library_id} =~ /EXACT=/) {
				$USE_EXACT = 1;
				my $LMINLMAX = substr $LB_CUTOFF{$library_id}, 6;
				my @tmpLMINLMAX = split(/,/, $LMINLMAX);
				$LMIN = $tmpLMINLMAX[0];
				$LMAX = $tmpLMINLMAX[1];
				print "EXACT flag: $LB_CUTOFF{$library_id} used, so using Lmin of ".$LMIN." and Lmax of ".$LMAX.".  Will NOT check distribution of concordant pairs!!\n";
			}
			elsif ($LB_CUTOFF{$library_id} =~ /SD=/ || $LB_CUTOFF{$library_id} =~ /PCT=/){
				print "Computing Lmin and Lmax ...\n";
				($LMIN, $LMAX) = GetLminLmax($CHECKING_FILE, $library_id, $MAPPINGQ, $LB_CUTOFF{$library_id}, $FIRST_READ, $ESPfile, $PROPER_L);
			}
			else{
				print STDERR "A wrong CUTOFF tag in the file.\n";
				$LMIN = 0;
				$LMAX = 0;
			}
		}
		else{
			print STDERR "A library $library_id does not record in the file.\n";
			print STDERR "Use default setting PCT=99% as the cutoff of Lmin and Lmax!\n";
			print "Computing Lmin and Lmax ...\n";
			my $dfile_cutoff = "PCT=99%";
			($LMIN, $LMAX) = GetLminLmax($CHECKING_FILE, $library_id, $MAPPINGQ, $dfile_cutoff, $FIRST_READ, $ESPfile, $PROPER_L);
		}
	}
	else{
		print STDERR "Check your CUTOFF_LMINLMAX tag, $CUTOFF!\n";
		exit 0;
	}

	CheckAllTypes($CHECKING_FILE, $LMIN, $LMAX);

	$LB_LMIN{$library_id} = $LMIN;
	$LB_LMAX{$library_id} = $LMAX;

	print "Lmax and Lmin are recorded in $ESPfile\n";

	my $filesize = -s $CHECKING_FILE;
	if($filesize == 0){
		print STDERR "Size of Concordant pairs is ZERO.\n";
	}
	if ((!defined($LMIN)) || (!defined($LMIN)) || $LMIN < 1 || $LMAX < 1) {
		print STDERR "WARNING: Lmin: $LMAX, and Lmax: $LMIN for library $library_id are either undefined or set to invalid values!\n";
		#exit;
	}
	system("rm -f $CHECKING_FILE");
}

foreach my $library_id(@lib_delete){
	delete $LB{$library_id};
}

if ($SEP_LIB == 1) {
	print "Parsing BAM file with different library ids ...\n";
}
else {
	print "Parsing BAM file ...\n";
}

my $cmd = "$SAMPATH view -F 0x000C -f 0x0001 $BAMFILE 2>$ERRLOG|";

# open files
my $fh_c = 0;
foreach my $library_id(keys %LB){
	my $outdisfile = $OUTPUTPREFIX."_".$library_id.".read.discordant";
	my $fh = "FH".$fh_c;
	$fh_c++;
	$LB_FH{$library_id} = $fh;
	open($fh, ">$outdisfile");
}

open(BAM, $cmd) || die "unable to open $BAMFILE\n";
while(my $line = <BAM>){
	chomp $line;
	my @temp = split(/\t+/, $line);
	my $flag = $temp[1];
	my $id;
	if($RG_TAG == 1){ # Header contains Reading Group information
		($id)=($line=~/RG\:Z\:(\S+)\s/);
	}
	else{ # no RG info. Assign RG as all.
		$id = "all";
	}

	if($SEP_LIB == 0 || exists $ID_LB{$id}){
		my ($TAG,$LMIN,$LMAX);
		if ($SEP_LIB == 0) {
			$TAG = $LB_FH{"all"};
			$LMIN = $LB_LMIN{"all"};
			$LMAX = $LB_LMAX{"all"};
			$ID_LB{$id} = "all";
		}
		else {
			$TAG = $LB_FH{$ID_LB{$id}};
			$LMIN = $LB_LMIN{$ID_LB{$id}};
			$LMAX = $LB_LMAX{$ID_LB{$id}};
		}	    

		if($temp[6] eq "="){
			my $q_ori = ($flag & 0x0010)?'-':'+';
			my $m_ori = ($flag & 0x0020)?'-':'+';

			if($PLAT_LB{$ID_LB{$id}} =~ /SOLID/i) { # flip the orientation of the second read
				if ($flag & 0x0040){
					$m_ori = ($m_ori eq '-')?'+':'-';
				}
				elsif ($flag & 0x0080){
					$q_ori = ($q_ori eq '-')?'+':'-';
				}
				else{
					print STDERR "WARNING: Read ID $temp[0] does not contain either the flag 0x0040 or 0x0080! This read will be skipped.\n";
					next;
				}
			}

			if($q_ori eq $m_ori){ # discordant file, putative inversion pairs
				print $TAG $line."\n";
			}
			else{
				if($q_ori eq "+" && $m_ori eq "-" && $temp[8] > 0){ # concordant or delete

					my $left_start = $temp[3];
					my $right_end = $temp[7]+length($temp[9]);
					my $flength = $right_end-$left_start;

					# 20101012, based on Coordinations rather than Insert Size tag.
					if ($flength <= $LMAX && $flength >= $LMIN){ 
					#if ($temp[8] <= $LMAX && $temp[8] >= $LMIN){ # output as concordant   Lmin <= L <= Lmax
						if($debug_mode == 1){ # continually write into concordant only debug mode was opened
							print CON $line."\n";
						}
					}
					else{ # output as discordant, putative deletion pairs
						print $TAG $line."\n";
					}
				}
				elsif($q_ori eq "-" && $m_ori eq "+" && $temp[8] < 0){ # concordant or delete
					my $left_start = $temp[7];
					my $right_end = $temp[3]+length($temp[9]);
					my $flength = $right_end-$left_start;

					if (abs($flength) <= $LMAX && abs($flength) >= $LMIN){
					#if (abs($temp[8]) <= $LMAX && abs($temp[8]) >= $LMIN){ # output as concordant
						if($debug_mode == 1){
							print CON $line."\n";
						}
					}
					else{ # output as discordant
						print $TAG $line."\n";
					}
				}
				else{ # discordant, putative divergent pairs
					print $TAG $line."\n";
				}
			}
		}	
		else{ # putative translocation pairs
			print $TAG $line."\n";
		}
	}
	else{
		print STDERR "Warning: Missing reading group $id in the header! Check header by samtools view -H! \n";
	}
}
close(BAM);

foreach my $library_id(keys %LB){
	my $prefix_esp = $OUTPUTPREFIX."_".$library_id;
	my $outdisfile = $OUTPUTPREFIX."_".$library_id.".read.discordant";

	# close files
	close($LB_FH{$library_id});

	my $DELTHRE = $LB_LMAX{$library_id};
	my $disfilesize = -s $outdisfile;
	if($disfilesize == 0){
		print STDERR "Size of Discordant pairs is ZERO.\n";
		HandleSamtoolsErrors();
	}
	else{
		sleep(3);
		HandleSamtoolsErrors();
		print "Generating GASV inputs... \n";
		print "$prefix_esp\n";
		# add naming parameter $NAMING_INDEXING_FILE (2010.01.27)
		# add Solid platform tag (2010.07.22)
		if (CheckFileInPATH("generate_GASV.pl")) {
			system("./generate_GASV.pl ALL $outdisfile $MAPPINGQ $DELTHRE $prefix_esp $NAMING_INDEXING_FILE $PLAT_LB{$library_id}");
		} else {
			system("./generate_GASV.pl ALL $outdisfile $MAPPINGQ $DELTHRE $prefix_esp $NAMING_INDEXING_FILE $PLAT_LB{$library_id}");
		}
		print "GASV input generation complete. Sorting Files... \n";
		system("bash sortESP.bash $prefix_esp.deletion");
		system("bash sortESP.bash $prefix_esp.divergent");
		system("bash sortESP.bash $prefix_esp.inversion");
		system("bash sortESP.bash $prefix_esp.translocation");
		system("mv -f $prefix_esp.deletion.sorted $prefix_esp.deletion");
		system("mv -f $prefix_esp.divergent.sorted $prefix_esp.divergent");
		system("mv -f $prefix_esp.inversion.sorted $prefix_esp.inversion");
		system("mv -f $prefix_esp.translocation.sorted $prefix_esp.translocation");
		print "Sorting Complete. \n";
		print "GASV inputs: \n";
		print "- $prefix_esp.deletion\n";
		print "- $prefix_esp.divergent\n";
		print "- $prefix_esp.inversion\n";
		print "- $prefix_esp.translocation\n";
		print "- lmin: ".$LB_LMIN{$library_id}.", lmax: ".$LB_LMAX{$library_id}."\n";
		MissingChromosomeErrors($prefix_esp);
		if($debug_mode == 0){
			system("rm -f $outdisfile");
		}
	}
}

if(-e $ERRLOG){
	system("rm -f $ERRLOG");
}

sub MissingChromosomeErrors {
	my $missingfile = $_[0].".missing_chromosome";
	my $missingfilesize = -s $missingfile;
	print "\n\n";

	if(-e $missingfile){
		print STDERR "The following chromosome(s) are not defined in chromosome naming file nor by default settings!\n";
		open(IN, $missingfile);
		while(<IN>){
			print STDERR $_;
		}
		close(IN);
	}
	system("rm -f $missingfile");
}

sub HandleSamtoolsErrors {
	open(ERR_IN, $ERRLOG) or die "cannot open $ERRLOG: $!\n";;
	while (my $line = <ERR_IN>) {
		chomp $line;
		if($line =~ /bam_header_read(.+?)EOF marker is absent/){
			#print "EOF line found: ".$line."\n";
			#ignore this error message, as it seems to come up with some 1000G bam files, but processing of header appears to be otherwise ok.

		} else {
			print $line."\n";
		}

	}
	close(ERR_IN);
	#system("rm -f $ERRLOG");
}

sub CheckPairInfo {

	my $in_checkf = $_[0];
	my $in_sampath = $_[1];
	my $in_bamfile = $_[2];
	my $in_library = $_[3];
	my $lib_sep = $_[4];
	my $in_fread = $_[5] * 4;
	my $cmd;
	my $pairing_flag = 0;

	print "Checking pairing information on library $in_library... ";
	if ($lib_sep == 1) {
		$cmd = "$in_sampath view -l $in_library $in_bamfile 2>$ERRLOG| head -n $in_fread > $in_checkf";
	}
	else {
		$cmd = "$in_sampath view $in_bamfile 2>$ERRLOG| head -n $in_fread > $in_checkf";
	}
	system($cmd);

	open(CHECK, $in_checkf) || die "unable to open $in_checkf\n";
	while(my $line = <CHECK>){
		chomp $line;

		my @temp = split(/\t+/, $line);

		if($temp[1] & 0x0001 && $temp[6] ne "*"){
			$pairing_flag = 1;
		}

		last if($pairing_flag == 1);
	}
	close(CHECK);

	if($pairing_flag == 1){
		print "\n";
	}

	return $pairing_flag;

}

sub GetLminLmax{

	my $in_checkf = $_[0];
	my $in_library = $_[1];
	my $in_mappingq = $_[2];
	my $in_cutoff = $_[3];
	my $in_num_reads = $_[4];
	my $in_espfile = $_[5];
	my $in_truncate_length = $_[6];

	if (CheckFileInPATH("LminLmaxProcessor")) {
		print STDERR "$in_checkf $in_mappingq $in_library $in_library $in_cutoff $in_num_reads $in_truncate_length\n";
		open(PROC_IN, "LminLmaxProcessor $in_checkf $in_mappingq $in_library $in_library $in_cutoff $in_num_reads $in_truncate_length |");
	} else {
		open(PROC_IN, "./LminLmaxProcessor $in_checkf $in_mappingq $in_library $in_library $in_cutoff $in_num_reads $in_truncate_length |");
	}

	my $line = <PROC_IN>;
	if (! $line) {
		print "LminLmaxProcessor failed - Could not find Lmin/Lmax!!\n";
		exit 0;
	}
	chomp $line;
	my @tmpESP = split(/\t/, $line);
	my $in_LMIN = $tmpESP[0];
	my $in_LMAX = $tmpESP[1];
	printf "\n%20s	Lmin	Lmax\n", " ";
	printf "%20s	%s	%s\n\n", $in_library,$in_LMIN,$in_LMAX;

	`echo "$in_library      $in_LMIN   $in_LMAX" >> $in_espfile`;
	close(PROC_IN);

	return ($in_LMIN, $in_LMAX);
}

sub CheckAllTypes {

	my $in_checkf = $_[0];
	my $in_LMIN = $_[1];
	my $in_LMAX = $_[2];
	# check four type data in this BAM file
	my $tl = 0;
	my $di = 0;
	my $in = 0;
	my $de = 0;

	print "Checking four types of structural variants in BAM file ... ";

	open(CHECK, $in_checkf);
	while(my $line = <CHECK>){
		chomp $line;
		my $id;
		my @temp = split(/\t+/, $line);
		my $flag = $temp[1];
		if($RG_TAG == 1){
			($id)=($line=~/RG\:Z\:(\S+)\s/);
		}
		else{
			$id = "all";
		}

		if ($SEP_LIB == 0) {
			$ID_LB{$id} = "all";
		}
		if($temp[6] eq "="){
			my $q_ori = ($temp[1] & 0x0010)?'-':'+';
			my $m_ori = ($temp[1] & 0x0020)?'-':'+';

			if($PLAT_LB{$ID_LB{$id}} =~ /SOLID/i) { # flip the orientation of the second read
				if ($flag & 0x0040){
					$m_ori = ($m_ori eq '-')?'+':'-';
				}
				elsif ($flag & 0x0080){
					$q_ori = ($q_ori eq '-')?'+':'-';
				}
				else{
					print STDERR "WARNING: Read ID $temp[0] does not contain either the flag 0x0040 or 0x0080! This read will be skipped.\n";
					next;
				}
			}


			if($q_ori eq $m_ori){ # discordant file, putative inversion pairs
				$in = 1;
			}
			else{
				if($q_ori eq "+" && $m_ori eq "-" && $temp[8] > 0){ # concordant or delete
					my $left_start = $temp[3];
					my $right_end = $temp[7]+length($temp[9]);
					my $flength = $right_end-$left_start;
					# 20101012, based on Coordinations rather than Insert Size tag.
					if ($flength <= $in_LMAX && $flength >= $in_LMIN){
					}
					else{ # output as discordant, putative deletion pairs
						$de = 1;
					}
				}
				elsif($q_ori eq "-" && $m_ori eq "+" && $temp[8] < 0){ # concordant or delete
					my $left_start = $temp[7];
					my $right_end = $temp[3]+length($temp[9]);
					my $flength = $right_end-$left_start;

					if (abs($flength) <= $in_LMAX && abs($flength) >= $in_LMIN){
					}
					else{ # output as discordant
						$de = 1;
					}
				}
				else{ # discordant, putative divergent pairs
					$di = 1;
				}
			}
		}
		else{
			if ($temp[6] ne "*") {
				$tl = 1;
			}
		}

		last if($tl == 1 && $di == 1 && $in == 1 && $de == 1)
	}
	close(CHECK);

	if($tl == 1 && $di == 1 && $in == 1 && $de == 1){
		print "\n";
	}

	if($tl == 0){
		print STDERR "\nWARNING: BAM file you're using appears to lack translocation information\n";
	}
	if($di == 0){
		print STDERR "\nWARNING: BAM file you're using appears to lack divergent information\n";
	}
	if($in == 0){
		print STDERR "\nWARNING: BAM file you're using appears to lack inversion information\n";
	}
	if($de == 0){
		print STDERR "\nWARNING: BAM file you're using appears to lack deletion/insertion information\n";
	}
}
sub CheckFileInPATH{
	my $file = $_[0];
	my $foundInPath = 0;
	my @PATH = split(/:/, $ENV{PATH});
	foreach my $pathDir (@PATH) {
		if (-e $pathDir."/".$file) {
			$foundInPath = 1;
			last;
		}
	}
	return $foundInPath;
}
sub CheckRequiredFile{
	my $Spath = $_[0];
	my $Bfile = $_[1];
	my $Mq = $_[2];
	my $Co = $_[3];
	my $ans = "";

	if (!(CheckFileInPATH("generate_GASV.pl") ) && !(-e "generate_GASV.pl")) {
		print STDERR "ERROR: generate_GASV.pl script not found!  BAM_Preprocessor.pl must be run from same directory in which generate_GASV.pl is located, or generate_GASV.pl must be in PATH.\n";
		exit 0;
	}

	if ($Spath eq "samtools") {
		#case where no path was set, so assumes samtools exists in the current PATH!
		if (!(CheckFileInPATH("samtools"))) {
			print STDERR "ERROR: path to samtools, currently set to $Spath , is not in PATH! Please either use -SAMTOOLS_PATH <path> to set it or add it to the PATH.\n";
			exit 0;
		}
	}
	elsif(!(-e $Spath)){
		print STDERR "ERROR: The path to samtools, currently set to $Spath , does not exist! Please use -SAMTOOLS_PATH <path> to set it.\n";
		exit 0;
	}
	elsif(-d $Spath){
		print STDERR "ERROR: The path to samtools $Spath is set as a directory, but it must actually indicate the samtools executable. Please correct -SAMTOOLS_PATH <path> to include the executable name, not just the containing directory.\n";
		exit 0;
	}
	elsif (!($Spath =~ /samtools/)) {
		print STDERR "ERROR: The path to samtools, currently set to $Spath , points to a file which exists, but does not actually point to the samtools command (the path does not include any files named 'samtools')! Please use -SAMTOOLS_PATH <path> to set it correctly.\n";
		exit 0;
	}

	# **** check BAM file format ****
	if ($Bfile =~ /^ftp:\/\//i || $Bfile =~ /^http:\/\//i){ # remote file check
		if(!LWP::Simple::head($Bfile)){
			print STDERR "The BAM file $Bfile does not exist!\n\n";
			exit 0;
		}
	}
	elsif (!(-e $Bfile)) { # local file check

		print STDERR "ERROR: BAM file $Bfile does not exist!\n\n";
		exit 0;
	}

	if ($Mq < 0) {
		print STDERR "ERROR: Mapping Quality must be larger than zero!\n";
		exit 0;
	}

	if ($Co !~ /^EXACT=/ && $Co !~ /^SD=/ && $Co !~ /^PCT=/ && $Co !~ /^FILE=/){
		print STDERR "Check your CUTOFF_LMINLMAX tag, $CUTOFF!\n";
		exit 0;
	}
	else{
		if (!($Co =~ /EXACT/)) {
			#check LminLmaxProcessor
			if (!(CheckFileInPATH("LminLmaxProcessor") ) && !(-e "LminLmaxProcessor")) {
				print STDERR "ERROR: LminLmaxProcessor does not exist!  LminLmaxProcessor must be built first and in PATH or current directory.  Run install.sh to build.\n";
				exit 0;
			}
		}
	}

}

sub readHeaderInfo{

	my $display_pl;
	my $readinggp_tag = 0;
	if($LIBSEP eq "all"){
		$LB{"all"} = "all";
		open(SAM_IN, $SAMPATH." view -H $BAMFILE 2>$ERRLOG|"); # pipeline in
		while(my $line = <SAM_IN>){
			chomp $line;
			if($line =~ /^\@RG/){
				my ($platform) = ($line =~ /PL\:(\S+)/);
				my ($id)=($line=~/ID\:(\S+)/);
				if($line !~ /PL\:(\S+)/ && $PLATFORM eq ""){
					$MISSING_PL_TAG = 1;
					$platform = "Illumina";
				}

				if(!exists $PLAT_LB{"all"}){
					if($line =~ /PL\:(\S+)/){
						$PLAT_LB{"all"} = $platform;
					}
					else{
						$PLAT_LB{"all"} = $PLATFORM;
					}
				}
				else{
					if($line =~ /PL\:(\S+)/){
						if($PLAT_LB{"all"} ne $platform){
							print STDERR "ERROR: Platforms in BAM file are not uniform. Check your SAM/BAM files!\n";
							exit(1);
						}
					}
				}
				$readinggp_tag = 1;
			}
		}
		close(SAM_IN);
		$display_pl = $PLAT_LB{"all"};
	}

	elsif($LIBSEP eq "separated-library"){
		$SEP_LIB = 1;
		my %plset;
		open(SAM_IN, $SAMPATH." view -H $BAMFILE 2>$ERRLOG|"); # pipeline in
		while(my $line = <SAM_IN>){
			chomp $line;
			if($line =~ /^\@RG/){
				my ($id)=($line=~/ID\:(\S+)/);
				my ($lib)=($line=~/LB\:(\S+)/);
				my ($platform)=($line=~/PL\:(\S+)/);
				my ($sample)=($line=~/SM\:(\S+)/);
				my ($insertsize)=($line=~/PI\:(\d+)/);

				if($line !~ /PL\:(\S+)/ && $PLATFORM eq ""){
					$MISSING_PL_TAG = 1;
					$platform = "Illumina";
				}
				if($line !~ /LB\:(\S+)/){
					#print STDERR "WARNING: Missing LB tag in \@RG $id lines. \n\n";
					push (@MISSING_LB_AY, "WARNING: Missing LB tag in \@RG $id lines.");
				}
				else{
					if($PLATFORM ne ""){
						$PLAT_LB{$lib} = $PLATFORM;
						$plset{lc $PLATFORM} = 1;
					}
					else{
						$PLAT_LB{$lib} = $platform;
						$plset{lc $platform} = 1;
					}
					$ID_LB{$id} = $lib;
					$LB{$lib} = $lib;
				}
				$readinggp_tag = 1;
			}
		}
		close(SAM_IN);
		#HandleSamtoolsErrors();
		my @tmp = keys %plset;
		if(scalar(@tmp) == 1){
			$display_pl = $tmp[0];
		}
		else{
			if($PLATFORM ne ""){
				$display_pl = $PLATFORM;
			}
			else{
				if(scalar(@tmp) > 1){
					print STDERR "ERROR: Platforms in BAM file are not uniform. Check your SAM/BAM files!\n";
					exit(1);
				}
				else{
					$display_pl = "illumina";
				}
			}
		}
		 if($readinggp_tag == 0){
			 $MISSING_PL_TAG = 1;
			 $display_pl = "illumina";
			 $RG_TAG = 0;
		 }
	 }
	 else{
		 print STDERR "Bad LIBRARY_SEPARATED_FLAG: ".$ARGV[0]."\n";
		 $display_pl = "";
		 exit;
	}

	return lc($display_pl);
}

sub getDefaultPrefix{
	my $in_bam = $_[0];
	my @temp = split(/\//,$in_bam);
	my $bamname = pop(@temp);

	if($bamname =~ /^(.+?)\.bam$/i){
		return $1;
	}
	else{
		return $bamname;
	}
}
