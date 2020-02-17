#!usr/bin/perl
# peakcountbywindow.pl
# determine the number of peak within window along chromosome


use warnings;
use strict;
use POSIX;
use File::Basename;

my $win = 100;

my @macs = @ARGV;

my %peaks;
my %presence;
my %prename;
my %peakrange;
my %peakcount;

foreach my $file (@macs) {
	open PEAK, "$file" or die;
	my %chr_win2;
	while (<PEAK>) {
		s/[\r\n]+$//;
		my @line = split "\t", $_;
		
		my $peak = $line[1] + $line[9];
		if ($presence{$line[0]}{$peak}) {
		  my $peakname = $prename{$line[0]}{$peak};
		  my @peakarray = split(";", $peakname);
		  $peakcount{$line[0]}{$peakarray[2]}{$file} ++;
		  
		  foreach my $i (($line[1] - $win) .. ($line[2] + $win)) {
		    $presence{$line[0]}{$i} ++;
		    $prename{$line[0]}{$i} = $peakname;
		    $peakrange{$line[0]}{$peakname}{$i} ++;
		  }
		  next;
		  
		} else {
		  my $fname = basename $file;
		  $peaks{$line[0]}{$peak} = join(";", $fname, $line[0], $peak, $line[1], $line[2]);
		  $peakcount{$line[0]}{$peak}{$file} ++;
		  my $peakname = join(";", $fname, $line[0], $peak);
		  foreach my $i (($line[1] - $win) .. ($line[2] + $win/2)) {
		    $presence{$line[0]}{$i} ++;
		    $prename{$line[0]}{$i} = $peakname;
		    $peakrange{$line[0]}{$peakname}{$i} ++;
		  }
		  
		}
		
	}
	close PEAK, or die;
}

foreach my $chr (sort %peaks) {
  foreach my $p (sort {$a <=> $b} keys %{$peaks{$chr}}) {
    my @peakarray = split(";", $peaks{$chr}{$p});
    my $peakname = join(";", @peakarray[0..2]);
    my @range = sort {$a <=> $b} keys %{$peakrange{$chr}{$peakname}};
    
    
    print "$chr\t", $range[0] + $win, "\t", $range[-1] - $win, "\t", $peakname, "\t";
    print scalar keys %{$peakcount{$chr}{$p}}, "\t";
    
    
    print "-\t-\t-\t-\t", $p - ($range[0] + $win), "\n";
  }
}

