#!/usr/bin/perl
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

%first_end;
%naming_chrom;
%missing_naming_chrom;

open(FHD, $ARGV[1]);
my $tag = $ARGV[0];
my $MAPQ = $ARGV[2];
my $output_prefix = $ARGV[4];
my $chrom_flag = "df"; # df: default ns: non-standard

if($tag eq "ALL"){
	my $TRLFILE = $ARGV[4].".translocation";
        open(TRL, "> $TRLFILE");
	my $INVFILE = $ARGV[4].".inversion";
	my $DELFILE = $ARGV[4].".deletion";
	my $DISVFILE = $ARGV[4].".divergent";
	open(INV, "> $INVFILE");
	open(DEL, "> $DELFILE");
	open(DISV, "> $DISVFILE");
	$T_LENGTH = $ARGV[3];
}

# handle non-standard chromosome naming problem, read chromosome naming table
if($ARGV[5] ne "DEFAULT"){
	open(IN, $ARGV[5]) or die $!;
	$chrom_flag = "ns";
	while(my $line = <IN>){
		chomp $line;
		my @temp = split(/\t/, $line);
		$naming_chrom{$temp[0]} = $temp[1];
	}
	close(IN);
}

while(my $line = <FHD>){
	chomp $line;

	my $chrom_tmp;
	my $chrom_exist = 0;
	my @temp_info = split(/\t/, $line);
	my $flag = $temp_info[1];

	my $q_ori = ($temp[1] & 0x0010)?'-':'+';
	($chrom_tmp, $chrom_exist) = checkChromosomeNaming($chrom_flag, $temp_info[2]);

	if(!($flag & 0x0400) && ($flag & 0x0001) && $temp_info[4] > $MAPQ && ($chrom_flag eq "df" && $chrom_exist == 1) || ($chrom_flag eq "ns" && $chrom_exist == 1)){
		if(!exists $first_end{$temp_info[0]}){
			$first_end{$temp_info[0]}++;
		}
		elsif(++$first_end{$temp_info[0]} == 2){
			delete $first_end{$temp_info[0]};
			my $left_coord_start;
			my $right_coord_start;
			my $left_coord_end;
			my $right_coord_end;
			my $left_ori;
			my $right_ori;
			my $left_chrom;
			my $right_chrom;
			my $current_chrom = $chrom_tmp;
			my $mate_chrom;

			if($temp_info[6] eq "="){
				$mate_chrom = $current_chrom;
			}
			else{
				($mate_chrom, $chrom_exist) = checkChromosomeNaming($chrom_flag, $temp_info[6]);
			}

			if($chrom_exist == 1){
				#if($temp_info[1] eq "1" || $temp_info[1] eq "65" || $temp_info[1] eq "129"){
				if(!($flag & 0x0010) && !($flag & 0x0020)){
					if($temp_info[3] < $temp_info[7]){ # this read as left mapped read
						$left_coord_start = $temp_info[3];
						$left_coord_end = $temp_info[3] + length($temp_info[9]);
						$left_chrom = $current_chrom;
						$left_ori = '+';
						$right_coord_end = $temp_info[7] + length($temp_info[9]);
						$right_coord_start = $temp_info[7];
						$right_chrom = $mate_chrom;
						$right_ori = '+';

					}
					else{ # this read as right mapped read
						$left_coord_start = $temp_info[7];
						$left_coord_end = $temp_info[7] + length($temp_info[9]);
						$left_chrom = $mate_chrom;
						$left_ori = '+';
						$right_coord_end = $temp_info[3] + length($temp_info[9]);
						$right_coord_start = $temp_info[3];
						$right_chrom = $current_chrom;
						$right_ori = '+';
					}
				}
				#elsif($temp_info[1] eq "17" || $temp_info[1] eq "81" || $temp_info[1] eq "145"){ # all with flag 0x0010 - strand of the query is reverse
				elsif(($flag & 0x0010) && !($flag & 0x0020)){					
					if($temp_info[3] < $temp_info[7]){ # this read as left mapped read
						$left_coord_start = $temp_info[3];
						$left_coord_end = $temp_info[3] + length($temp_info[9]);
						$left_chrom = $current_chrom;
						$left_ori = '-';
						$right_coord_end = $temp_info[7] + length($temp_info[9]);
						$right_coord_start = $temp_info[7];
						$right_chrom = $mate_chrom;
						$right_ori = '+';
					}
					else{
						$left_coord_start = $temp_info[7];
						$left_coord_end = $temp_info[7] + length($temp_info[9]);
						$left_chrom = $mate_chrom;
						$left_ori = '+';
						$right_coord_end = $temp_info[3] + length($temp_info[9]);
						$right_coord_start = $temp_info[3];
						$right_chrom = $current_chrom;
						$right_ori = '-';
					}
				}
				#elsif($temp_info[1] eq "33" || $temp_info[1] eq "161" || $temp_info[1] eq "97"){ # all with flag 0x0020 - strand of the mate is reverse
				elsif(!($flag & 0x0010) && ($flag & 0x0020)){					
					if($temp_info[3] < $temp_info[7]){ # this read as left mapped read
						$left_coord_start = $temp_info[3];
						$left_coord_end = $temp_info[3] + length($temp_info[9]);
						$left_chrom = $current_chrom;
						$left_ori = '+';
						$right_coord_end = $temp_info[7] + length($temp_info[9]);
						$right_coord_start = $temp_info[7];
						$right_chrom = $mate_chrom;
						$right_ori = '-';
					}
					else{
						$left_coord_start = $temp_info[7];
						$left_coord_end = $temp_info[7] + length($temp_info[9]);
						$left_chrom = $mate_chrom;
						$left_ori = '-';
						$right_coord_end = $temp_info[3] + length($temp_info[9]);
						$right_coord_start = $temp_info[3];
						$right_chrom = $current_chrom;
						$right_ori = '+';
					}
				}
				#elsif($temp_info[1] eq "49" || $temp_info[1] eq "113" || $temp_info[1] eq "177"){
				#if(($flag & 0x0010) && ($flag & 0x0020)){
				else{
					if($temp_info[3] < $temp_info[7]){ # this read as left mapped read
						$left_coord_start = $temp_info[3];
						$left_coord_end = $temp_info[3] + length($temp_info[9]);
						$left_chrom = $current_chrom;
						$left_ori = '-';
						$right_coord_end = $temp_info[7] + length($temp_info[9]);
						$right_coord_start = $temp_info[7];
						$right_chrom = $mate_chrom;
						$right_ori = '-';
					}
					else{ # this read as right mapped read
						$left_coord_start = $temp_info[7];
						$left_coord_end = $temp_info[7] + length($temp_info[9]);
						$left_chrom = $mate_chrom;
						$left_ori = '-';
						$right_coord_end = $temp_info[3] + length($temp_info[9]);
						$right_coord_start = $temp_info[3];
						$right_chrom = $current_chrom;
						$right_ori = '-';
					}
				}

				if ((($right_chrom eq "X") || ($right_chrom eq "x")) && $chrom_flag eq "df") {
					$right_chrom = 23;
				}
				if ((($right_chrom eq "Y") || ($right_chrom eq "y")) && $chrom_flag eq "df") {
					$right_chrom = 24;
				}
				#if ((($right_chrom eq "MT") || ($right_chrom eq "mt")) && $chrom_flag eq "df") {
				#	$right_chrom = 25;
				#}

				if ((($left_chrom eq "X") || ($left_chrom eq "x")) && $chrom_flag eq "df") {
					$left_chrom = 23;
				}
				if ((($left_chrom eq "Y") || ($left_chrom eq "y")) && $chrom_flag eq "df") {
					$left_chrom = 24;
				}
				#if ((($left_chrom eq "MT") || ($left_chrom eq "mt")) && $chrom_flag eq "df") {
				#	$left_chrom = 25;
				#}

				if($right_chrom ne $left_chrom && $tag eq "ALL"){
					# translocation
					if ($left_chrom < $right_chrom) {
						print TRL $temp_info[0]."	".$left_chrom."	".$left_coord_start."	".$left_coord_end."	".$left_ori."	".$right_chrom."	".$right_coord_start."	".$right_coord_end."	".$right_ori."\n";
					}
					else{
						print TRL $temp_info[0]."	".$right_chrom."	".$right_coord_start."	".$right_coord_end."	".$right_ori."	".$left_chrom."	".$left_coord_start."	".$left_coord_end."	".$left_ori."\n";
					}

				}
				elsif($left_chrom eq $right_chrom && $tag eq "ALL"){

					if(($left_ori eq '+' && $right_ori eq '+') || ($left_ori eq '-' && $right_ori eq '-')){
						# inversion
						if($right_coord_end - $left_coord_start > 0){
							print INV $temp_info[0]."	".$left_chrom."	".$left_coord_start."	".$left_coord_end."	".$left_ori."	".$right_chrom."	".$right_coord_start."	".$right_coord_end."	".$right_ori."\n";
						}
						else{
							print INV $temp_info[0]."	".$right_chrom."	".$right_coord_start."	".$right_coord_end."	".$right_ori."	".$left_chrom."	".$left_coord_start."	".$left_coord_end."	".$left_ori."\n";
						}

					}
					else{
						# deletion
						if($left_ori eq '+'){ #&& $right_coord_end - $left_coord_start > 0){
							my $frag_length = $right_coord_end - $left_coord_start;
							if($frag_length > $T_LENGTH){
								print DEL $temp_info[0]."	".$left_chrom."	".$left_coord_start."	".$left_coord_end."	".$left_ori."	".$right_chrom."	".$right_coord_start."	".$right_coord_end."	".$right_ori."\n";
							}
						}
						# divergent
						elsif($right_ori eq '+'){ # && $left_coord_end - $right_coord_start > 0){
							my $frag_length = $right_coord_end - $left_coord_start;
							if($frag_length > $T_LENGTH){
								print DISV $temp_info[0]."	".$left_chrom."	".$left_coord_start."	".$left_coord_end."	".$left_ori."	".$right_chrom."	".$right_coord_start."	".$right_coord_end."	".$right_ori."\n";
							}
						}
					}
				}
			}	

		}

	}
}

close(FHD);
if($tag eq "ALL"){
	close(TRL);
	close(INV);
	close(DEL);
	close(DISV);
}

if(scalar(keys %missing_naming_chrom) != 0){
	my $missing_file = $ARGV[4].".missing_chromosome";
	open(OUT, ">$missing_file");
	foreach (keys %missing_naming_chrom){
		print OUT $missing_naming_chrom{$_}."\n";
	}
	close(OUT);
}


sub checkChromosomeNaming{
	my $c_flag = $_[0];
	my $chrom = $_[1];
	my $r_chrom;
	my $r_chrom_exist;

	if($c_flag eq "df" && $chrom =~ /^chr(.+)$/i){
		$r_chrom = $1;
		if(($r_chrom >= 1 && $r_chrom <= 24) || $r_chrom eq 'X' || $r_chrom eq 'x' || $r_chrom eq 'Y' || $r_chrom eq 'y' ){#|| $r_chrom eq "MT" || $r_chrom eq "mt"){
			$r_chrom_exist = 1;
		}
		else{
			$missing_naming_chrom{$chrom} = $chrom." missed in DEFAULT NAMING!";
			$r_chrom_exist = 0;
		}
	}
	elsif($c_flag eq "df"){
		$r_chrom = $chrom;
		if(($r_chrom >= 1 && $r_chrom <= 24) || $r_chrom eq 'X' || $r_chrom eq 'x' || $r_chrom eq 'Y' || $r_chrom eq 'y' ){ #|| $r_chrom eq "MT" || $r_chrom eq "mt"){
			$r_chrom_exist = 1;
		}
		else{
			$missing_naming_chrom{$chrom} = $chrom." missed in DEFAULT NAMING!";
			$r_chrom_exist = 0;
		}

	}
	elsif($c_flag eq "ns"){
		if(exists $naming_chrom{$chrom}){
			$r_chrom = $naming_chrom{$chrom};
			$r_chrom_exist = 1;
		}
		else{
			$missing_naming_chrom{$chrom} = $chrom." missed in NAMING TABLE!";
			$r_chrom_exist = 0;
		}
	}
	else{
		$r_chrom_exist = 0;
	}

	return ($r_chrom, $r_chrom_exist);

}
