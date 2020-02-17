#!usr/bin/perl
# sam_jointmapsplit.pl
# reads are aligned to a combined fasta of two genomes.
# this script splits the reads into their respective genomes.

use strict;
use warnings;
use POSIX;
use File::Basename;

my $fai1 = $ARGV[0]; # index of first genome
my $fai2 = $ARGV[1]; # index of second genome
my $input1 = $ARGV[2]; # sam file

my $out1 = substr basename($input1), 0, -4;
my $out2 = substr basename($input1), 0, -4;
open OUT1, ">$out1.split1.sam";
open OUT2, ">$out2.split2.sam";

open FAI1, $fai1 or die;
open FAI2, $fai2 or die;
my %fa1;
my %fa2;
while (<FAI1>) {
	s/[\r\n]+$//;
	my @line = split "\t", $_;
	$fa1{$line[0]} ++;
}
while (<FAI2>) {
	s/[\r\n]+$//;
	my @line = split "\t", $_;
	$fa2{$line[0]} ++;
	if ($fa1{$line[0]}) {
		die "same contig names found in both references\n";
	}
}


open SAM, $input1 or die;

my $better1 = 0;
my $better2 = 0;
my $missing = 0;
my $PE_diff = 0;

while (<SAM>) {
	s/[\r\n]+$//;
	my @line = split "\t", $_;
#	print "$_\n";
	if ($_ =~ /^@/ ) {
#		print "$_\n";
		if ($line[0] =~ /^\@SQ/) {
			my $contig = substr $line[1], 3;
#			print "$contig\n";
#			<STDIN>;
			if ($fa1{$contig}) {
				print OUT1 "$_\n";
			} elsif ($fa2{$contig}) {
				print OUT2 "$_\n"; 
			} else {
#				print "nothing";
#				<STDIN>;
			}
		} else {
			print OUT1 "$_\n";
			print OUT2 "$_\n";
		}
		next;
	}
	
	if ($line[2] eq "*") {
		$missing ++;
	} else {
		if ($fa1{$line[2]}) {
			if ($fa2{$line[6]}) {
				$PE_diff ++;
			} else {
				$better1 ++;
				print OUT1 "$_\n";
			}
		} elsif ($fa2{$line[2]}) {
			if ($fa1{$line[6]}) {
				$PE_diff ++;
			} else {
				$better2 ++;
				print OUT2 "$_\n";
			}
			
		} else {
#			print "can't find contig\n";
#			<STDIN>;
		}
	}
}
print basename($ARGV[0]), "\t", basename($ARGV[1]), "\tunmapped_reads\tpairs_mapped_to_different_genomes\ttotal_reads\n";
print "$better1\t$better2\t$missing\t$PE_diff\t", $better1 + $better2 + $missing + $PE_diff, "\n";

