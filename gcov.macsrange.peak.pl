#!usr/bin/perl
# gcov.macsrange.peak.pl
# take enrichment files and for each peak in the macs peak caller provide single basepair resolution of enrichment depending on window size

use warnings;
use strict;
use POSIX;

my $macs = $ARGV[0];
my $enr = $ARGV[1];
my $win = $ARGV[2];


my %peaks;
open MAC, "$macs" or die;
while (<MAC>) {
	s/[\r\n]+$//;
	my @line = split "\t", $_;
	my $peak = $line[1] + $line[9];
	my $peakright = $peak + $win;
	$peaks{$line[0]}{$peakright} = $peak;
}
close MAC;

my %enrich_prior;
my $chr = "undef";
open EN, "$enr" or die;
my @cord;
while (<EN>) {
	s/[\r\n]+$//;
	#	print "$_\n";
	my @line = split "\t", $_;
	
	if ($line[0] eq $chr) {
		push @cord, $line[1];
		$enrich_prior{$line[0]}{$line[1]} = $line[2];
		if (scalar(keys %{$enrich_prior{$chr}}) > $win*3) {
			
			delete $enrich_prior{$chr}{$cord[0]};
			shift @cord;
		}
		
		if ($peaks{$line[0]}{$line[1]}) {
			my $peakright = $line[1];
			my $peak = $peaks{$line[0]}{$line[1]};
			my $peakleft = $peak - $win;
			
			print "$chr\t$peak";

			for my $i ($peakleft .. ($peak - 1)) {
				if (exists $enrich_prior{$line[0]}{$i}) {
					print "\t", $enrich_prior{$line[0]}{$i};
					
				} else {
					print "\t", "NA";
				}
			}
			print "\t", $enrich_prior{$line[0]}{$peak};
			for my $i (($peak + 1) .. $peakright) {
				if (exists $enrich_prior{$line[0]}{$i}) {
					print "\t", $enrich_prior{$line[0]}{$i};
				} else {
					print "\t", "NA";
				}
			}
#			print "$chr\t$peak\t", $peakleft, "\t", $peakright, "\n";
			print "\n";

#			print $cord[0],"\t", $cord[-1], "\n";
			delete $peaks{$chr}{$peakright};
#			<STDIN>;
		}
	} else { ## new chromosome. need to purge all hashes and arrays from previous chromosome if they exist
		foreach my $r (keys %{$peaks{$chr}}) {
			my $peakright = $r;
			my $peak = $peaks{$chr}{$r};
			my $peakleft = $peak - $win;
			if ($enrich_prior{$chr}{$peakleft}) {
				
				print "$chr\t$peak";
				
				for my $i ($peakleft .. ($peak - 1)) {
					if (exists $enrich_prior{$chr}{$i}) {
						print "\t", $enrich_prior{$chr}{$i};
					} else {
						print "\t", "NA";
					}
				}
				print "\t", $enrich_prior{$chr}{$peak};
				for my $i (($peak + 1) .. $peakright) {
					if (exists $enrich_prior{$chr}{$i}) {
						print "\t", $enrich_prior{$chr}{$i};
					} else {
						print "\t", "NA";
					}
				}
				print "\n";
#				print "old $chr\t$peak\t", $peakleft, "\t", $peakright, "\n";
#				print "old ", $cord[0],"\t", $cord[-1], "\n";
				delete $peaks{$chr}{$peakright};
#				<STDIN>;
			}
		}
		$chr = $line[0];
		undef @cord;
		undef %enrich_prior;

		$enrich_prior{$line[0]}{$line[1]} = $line[2];
		push @cord, $line[1];
	}
}
close EN;


