#!/usr/bin/perl
# A script to isolate position that are contain in a given interval
# The input files are a list of position in the format Scaffold\tStart\tEnd\tSample
# such as Hmel215006	125	1259	Sample1
# such as Hmel215006	1889	3859	Sample1
# and the intervals we are interested in for each sample, in the format Sample\tScaffold\tStart\t
# such as Sample1 Scaffold1	123	158900 
# such as Sample3 Scaffold6	5899	58987 
#
use strict;
use warnings;
use Getopt::Std;

my %options=();

getopts("i:o:p:", \%options);

open(POS, "<", $options{p}) || die; # List of interval
open(IN, "<", $options{i}) || die; #List of position
open(OUT, ">", $options{o}) || die; # OUTPUT


my %position=();

while(<POS>)
{
	my @Read=split(/[\b\n\t\s]+/, $_);
	$position{$Read[0]}{$Read[1]}=[$Read[2],$Read[3]];
}

while(<IN>)
{
	my @Read=split(/[\b\n\t\s]+/, $_);
	if (exists $position{$Read[3]})
	{
		if (exists $position{$Read[3]}{$Read[2]})
		{
			if ($Read[0]>$position{$Read[3]}{$Read[2]}[0] && $Read[1]<$position{$Read[3]}{$Read[2]}[1])
			{
				print OUT $_;
			}
		}
	}
}




