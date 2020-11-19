#!usr/bin/perl
# euTEjunc.pl
# take junction reads from samread_TE.insert.pl to find euchromatic insertion
# insertions where both left and right junctions are identified are label 2
# isenrtions where only on junction is identified are label 1

use strict;
use warnings;
use POSIX;
use File::Basename;

my $TE1 = $ARGV[0]; # pair1 mapping to TE junction
my $TE2 = $ARGV[1]; # pair2 mapping to TE junction
my $GE1 = $ARGV[2]; # pair1 mapping to the genome junction
my $GE2 = $ARGV[3]; # pair2 mapping to the genome junction

my %reads;

open T1, $TE1 or die;
while (<T1>) {
	next if ($_ =~ /^@/);
	s/[\r\n]+$//;
	my @line = split "\t", $_;
	my @TEname = split "#", $line[2];
	$reads{$line[0]}{"TE"} = $TEname[0];
}
close T1;

open G2, $GE2 or die;
while (<G2>) {
	next if ($_ =~ /^@/);
	s/[\r\n]+$//;
	my @line = split "\t", $_;
	
	$reads{$line[0]}{"c"} = $line[2];
	$reads{$line[0]}{"p"} = $line[3];
	$reads{$line[0]}{"C"} = $line[4];
	my @cig = cigar_parse($line[5]);
	my $end = $line[3];
	foreach my $c (@cig){
		if (substr($c, -1) eq "M") {
			$end += substr($c, 0, length($c) - 1);
		} elsif (substr($c, -1) eq "D") {
			$end -= substr($c, 0, length($c) - 1);
		}
	}
	$reads{$line[0]}{"e"} = $end;
	
	if ($line[1] & 0x10){ ##find reverse mapping by sam flag
		$reads{$line[0]}{"s"} = "-";	
	} else {
		$reads{$line[0]}{"s"} = "+";
		
	}
	
}
close G2;
#<STDIN>;

open T2, $TE2 or die;
while (<T2>) {
	next if ($_ =~ /^@/);
	s/[\r\n]+$//;
	my @line = split "\t", $_;
	my @TEname = split "#", $line[2];
	$reads{$line[0]}{"TE"} = $TEname[0];
}
close T2;

open G1, $GE1 or die;
while (<G1>) {
	next if ($_ =~ /^@/);
	s/[\r\n]+$//;
	my @line = split "\t", $_;
	
	$reads{$line[0]}{"c"} = $line[2];
	$reads{$line[0]}{"p"} = $line[3];
	$reads{$line[0]}{"C"} = $line[4];
	my @cig = cigar_parse($line[5]);
	my $end = $line[3];
	foreach my $c (@cig){
		if (substr($c, -1) eq "M") {
			$end += substr($c, 0, length($c) - 1);
		} elsif (substr($c, -1) eq "D") {
			$end -= substr($c, 0, length($c) - 1);
		}
	}
	$reads{$line[0]}{"e"} = $end;
	
	if ($line[1] & 0x10){ ##find reverse mapping by sam flag
		$reads{$line[0]}{"s"} = "-";
		
	} else {
		$reads{$line[0]}{"s"} = "+";
	}
}
close G1;
#<STDIN>;

#print scalar keys %reads, "\n";
my %chrpos;
foreach my $r (keys %reads) {
	unless (exists $reads{$r}{"p"} && exists $reads{$r}{"TE"}) {
		next;
	}
	my $chr = $reads{$r}{"c"};
	foreach my $p ($reads{$r}{"p"} .. ($reads{$r}{"e"})) {
#		print join("", $reads{$r}{"s"}, $reads{$r}{"TE"}), "\n";
		${$chrpos{$chr}{$p}}{$reads{$r}{"TE"}}{$reads{$r}{"s"}} ++;
		
	}
}

my %five;
my %three;
foreach my $chr (sort keys %chrpos) {
	my @ar5 = (0);
	my %h5TEs = ();
	my @ar3 = (0);
	my %h3TEs = ();

	foreach my $pos (sort {$a <=> $b} keys %{$chrpos{$chr}}) {
		
		foreach my $te (sort keys %{$chrpos{$chr}{$pos}}) {
#		print "$chr\t$pos\t$te\t";
			if (exists $chrpos{$chr}{$pos}{$te}{"+"}) {
#				print "+\t", $chrpos{$chr}{$pos}{$te}{"+"}, "\t";
				if ($pos > ($ar5[-1] + 100)) {
					if ($ar5[0] != 0) {
						$five{$chr}{$ar5[0]}{3} = $ar5[-1];
						my @tejoin;
						foreach my $t (sort keys %h5TEs) {
							my @f = sort {$b <=> $a} @{$h5TEs{$t}};
							push @tejoin, join("=", $t, $f[0]);
						}
						$five{$chr}{$ar5[0]}{"TE"} = join(";", @tejoin);
#						print "\n$ar5[0]\t$ar5[-1]\t", join(";", @tejoin), "\n";
#						<STDIN>;
						
					}
					@ar5 = ($pos);
					undef %h5TEs;
					push @{$h5TEs{$te}}, $chrpos{$chr}{$pos}{$te}{"+"};
				} else {
					push @ar5, $pos;
					push @{$h5TEs{$te}}, $chrpos{$chr}{$pos}{$te}{"+"};
				}
				
			}
			if (exists $chrpos{$chr}{$pos}{$te}{"-"}) {
#				print "-\t", $chrpos{$chr}{$pos}{$te}{"-"};
				if ($pos > ($ar3[-1] + 100)) {
					if ($ar3[0] != 0) {
						$three{$chr}{$ar3[0]}{3} = $ar3[-1];
						my @tejoin;
						foreach my $t (sort keys %h3TEs) {
							my @f = sort {$b <=> $a} @{$h3TEs{$t}};
							push @tejoin, join("=", $t, $f[0]);
						}
						$three{$chr}{$ar3[0]}{"TE"} = join(";", @tejoin);
#						print "\n$ar3[0]\t$ar3[-1]\t", join(";", @tejoin), "\n";
					}
					@ar3 = ($pos);
					undef %h3TEs;
					push @{$h3TEs{$te}},  $chrpos{$chr}{$pos}{$te}{"-"};
				} else {
					push @ar3, $pos;
					push @{$h3TEs{$te}},  $chrpos{$chr}{$pos}{$te}{"-"};
				}
				
				
				
#				<STDIN>;
				
			}
#			print "\n";
		}
		
	}
}
my $insert = 100;
my %in;
my %expos3;
foreach my $chr (sort keys %five) {
	foreach my $pos5 (sort {$a <=> $b} keys %{$five{$chr}}) {
		my %diff;
		foreach my $pos3 (keys %{$three{$chr}}) {
			next if (($pos3 + 10) < $five{$chr}{$pos5}{3});
			my $dist = abs($five{$chr}{$pos5}{3} - $pos3);
			$diff{$pos3} = $dist;
		}
		foreach my $pos3 (sort {$diff{$a} <=> $diff{$b}} keys %diff) {
			if ($diff{$pos3} <= $insert) {
				$expos3{$chr}{$pos3} ++;
				$in{$chr}{$pos5} = join "\t", $chr, $pos5."-".$five{$chr}{$pos5}{3}, $five{$chr}{$pos5}{"TE"}, $pos3."-".$three{$chr}{$pos3}{3}, $three{$chr}{$pos3}{"TE"};
				last;
			}
		}
	}
}
foreach my $chr (sort keys %five) {
	foreach my $pos (sort {$a <=> $b} keys %{$five{$chr}}) {
		if (exists $in{$chr}{$pos}) {
		
		} else {
			$in{$chr}{$pos} = join "\t", $chr, $pos."-".$five{$chr}{$pos}{3}, $five{$chr}{$pos}{"TE"}, "-", "na";
		}
	}
}
foreach my $chr (sort keys %three) {
	foreach my $pos (sort {$a <=> $b} keys %{$three{$chr}}) {
		if (exists $expos3{$chr}{$pos}) {
		
		} else {
			$in{$chr}{$pos} = join "\t", $chr, "-", "na", $pos."-".$three{$chr}{$pos}{3}, $three{$chr}{$pos}{"TE"};
		}		
	}
}
foreach my $chr (sort keys %in) {
	foreach my $pos (sort {$a <=> $b} keys %{$in{$chr}}) {
		print $in{$chr}{$pos}, "\n";
	}
}


sub cigar_parse {
	my @CIG = split/(?<=[A-Z])/, $_[0];
	return @CIG;
}
