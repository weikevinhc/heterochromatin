#!usr/bin/perl
# samreads_TE.insert.pl
# identify read pairs where one maps to unique sequence in genome and the other to the repeat library
# paired end reads are separately mapped to either the repeat library or the genome
#

use strict;
use warnings;
use POSIX;
use File::Basename;

my $TE1 = $ARGV[0]; # pair1 mapping to TE
my $TE2 = $ARGV[1]; # pair2 mapping to TE
my $GE1 = $ARGV[2]; # pair1 mapping to the genome
my $GE2 = $ARGV[3]; # pair2 mapping to the genome

my %reads;

open T1, $TE1 or die;
while (<T1>) {
	next if ($_ =~ /^@/);
	s/[\r\n]+$//;
	my @line = split "\t", $_;
#	if (exists $reads{$line[0]}{"TE1"}) {
#		print "$_\n";
#	}
	if ($line[2] eq "*") {
		$reads{$line[0]}{"TE1"} = 0;
	} elsif ($line[5] =~ /H/) {
		$reads{$line[0]}{"TE1"} = 0;
#		print "$_\n";
	} else {
		$reads{$line[0]}{"TE1"} = 1;
	}
}
close T1;

open T2, $TE2 or die;
while (<T2>) {
	next if ($_ =~ /^@/);
	s/[\r\n]+$//;
	my @line = split "\t", $_;
	if ($line[2] eq "*") {
		$reads{$line[0]}{"TE2"} = 0;
	} elsif ($line[5] =~ /H/) {
#		print "$_\n";
		$reads{$line[0]}{"TE2"} = 0;
	} else {
		$reads{$line[0]}{"TE2"} = 1;
	}
}
close T2;

open G1, $GE1 or die;
while (<G1>) {
	next if ($_ =~ /^@/);
	s/[\r\n]+$//;
	my @line = split "\t", $_;
	if ($line[2] eq "*") {
		$reads{$line[0]}{"GE1"} = 0;
	} elsif ($line[5] =~ /H/) {
		$reads{$line[0]}{"GE1"} = 0;
	} else {
		$reads{$line[0]}{"GE1"} = $line[4];
	}
}
close G1;

open G2, $GE2 or die;
while (<G2>) {
	next if ($_ =~ /^@/);
	s/[\r\n]+$//;
	my @line = split "\t", $_;
	if ($line[2] eq "*") {
		$reads{$line[0]}{"GE2"} = 0;
	} elsif ($line[5] =~ /H/) {
		$reads{$line[0]}{"GE2"} = 0;
	} else {
		$reads{$line[0]}{"GE2"} = $line[4];
	}
}
close G2;
warn "number of starting pairs: ", scalar(keys %reads), "\n";
foreach my $r (keys %reads) {
	if (scalar(keys %{$reads{$r}}) != 4) {
		delete $reads{$r};
		next;
	}
	if ($reads{$r}{"TE1"} == $reads{$r}{"TE2"}) {
		delete $reads{$r} ;
		next;
	}
	if ($reads{$r}{"TE1"} == 1 && $reads{$r}{"GE2"} < 4) {
		delete $reads{$r};
		next;
	}
	if ($reads{$r}{"TE2"} == 1 && $reads{$r}{"GE1"} < 4) {
		delete $reads{$r};
		next;
	}
}
warn "number of remaining pairs: ", scalar(keys %reads), "\n";
#<STDIN>;

open TE, $TE1 or die;
my $out = $TE1 . ".junc";
open OUT, ">$out" or die;
while (<TE>) {
	s/[\r\n]+$//;
	if ($_ =~ /^@/) {
		print OUT "$_\n";
		next;
	}
	my @line = split "\t", $_;
	if (exists $reads{$line[0]}) {
		if ($reads{$line[0]}{"TE1"} == 1) {
			print OUT "$_\n";
		}
	}
}
close TE;
close OUT;

open GE, $GE2 or die;
$out = $GE2 . ".junc";
open OUT, ">$out" or die;
while (<GE>) {
	s/[\r\n]+$//;
	if ($_ =~ /^@/) {
		print OUT "$_\n";
		next;
	}
	my @line = split "\t", $_;
	if (exists $reads{$line[0]}) {
		if ($reads{$line[0]}{"TE1"} == 1) {
			print OUT "$_\n";
		}
	}
}
close GE;
close OUT;

open TE, $TE2 or die;
$out = $TE2 . ".junc";
open OUT, ">$out" or die;
while (<TE>) {
	s/[\r\n]+$//;
	if ($_ =~ /^@/) {
		print OUT "$_\n";
		next;
	}
	my @line = split "\t", $_;
	if (exists $reads{$line[0]}) {
		if ($reads{$line[0]}{"TE2"} == 1) {
			print OUT "$_\n";
		}
	}
}
close TE;
close OUT;

open GE, $GE1 or die;
$out = $GE1 . ".junc";
open OUT, ">$out" or die;
while (<GE>) {
	s/[\r\n]+$//;
	if ($_ =~ /^@/) {
		print OUT "$_\n";
		next;
	}
	my @line = split "\t", $_;
	if (exists $reads{$line[0]}) {
		if ($reads{$line[0]}{"TE2"} == 1) {
			print OUT "$_\n";
		}
	}
}
close TE;
close OUT;
