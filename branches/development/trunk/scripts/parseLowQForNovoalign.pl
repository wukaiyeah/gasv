#!/usr/bin/perl

use strict;

# Decoding Bam File FLAGS
# Bit 7: 1 --> The first read in the template.

if($#ARGV<0){
	die "Usage: ./exe.pl {Samfile}\n";
}

my $SAM_FILE = $ARGV[0];
open(BAM,$SAM_FILE) or die "Could not open samfile $SAM_FILE\n";

open(READ1,">".$SAM_FILE.".read1.fastq");
open(READ2,">".$SAM_FILE.".read2.fastq");


my $bothMappedHighQual   = 0;
my $bothMappedLowQual    = 0;
my $bothUnmapped         = 0;
my $oneMapped            = 0;
my $other                = 0;

my %otherFlags;

my $numLines = 0;

my $read1 = 0;
my $read2 = 0;

while(my $line = <BAM>){
	$numLines++;
	if($numLines%500000 == 0){ print "Processed $numLines\n"; }
	
	chomp($line);

	my $firstChar = substr($line,0,1);
	if($firstChar eq "@"){ 
		print "Processing Header:\t$line\n";
	}
	else {
		my @temp = split(/\t+/, $line);

		#Name				Flag	Chr	Pos 	MapQ
		#0                               1        2     3        4      5
		#CHR17_1_450_1:0:0_3:0:0_8b4ab4	163	CHR17	1	29	50M	=	401	450	AAGCTTCTCACCCTGTTCCTGCATAGATAATTGCATGACAATTGCGTTGT	22222222222222222222222222222222222222222222222222	XT:A:U	NM:i:1	SM:i:29	AM:i:29	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:45C4
	
		my $FLAG  = $temp[1];
		my $name  = $temp[0];
		my $pos   = $temp[3];
		my $map   = $temp[4];
		my $CIGAR = $temp[5];
		my $seq   = $temp[9];
		my $qual  = $temp[10];

		my $q_ori = ($FLAG & 0x0010)?'-':'+';
		my $m_ori = ($FLAG & 0x0020)?'-':'+';

                my $finalSeq = "";
                if($q_ori eq "+"){$finalSeq = $seq; }
                else{ $finalSeq = revComp($seq); }

		#@CHR17_4205961_4206161_0:0:0_1:0:0_0/1
		#TTTACCAAATCAAGGAAATTTACTTCTATTCCTAGTTTCTTTCTTTTTTA
		#+
		#22222222222222222222222222222222222222222222222222

		if($FLAG & 0x0040){
			print READ1 "\@$name/1\n$seq\n+\n$qual\n";
			$read1++;
		}
		elsif($FLAG & 0x0080){
			print READ2 "\@$name/2\n$seq\n+\n$qual\n";
			$read2++;
		}
		else{
			die "First or second read flag is not correctly set in BAM file --> $line\n";
		}	
	}
}
close(BAM);

print "Finished processing reads; we found $read1 first reads and $read2 second reads.\n";

if($read1==$read2){
	print "Finished processing $SAM_FILE successfully...\n";
}
else{ 
	print "The numbers of read1 and read2 do not match up.\n";
	print "This is usually not reasonable output. Check file for errors.\n";
}

#Take the reverse complement of a sequence.
sub revComp {
  my $dna = $_[0];
  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}


