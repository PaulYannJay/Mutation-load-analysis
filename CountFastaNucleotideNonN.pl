#!/usr/bin/perl
### Just a small script to compute the total size of the sequences contained in a fasta file (not including N)
use strict;
use warnings;
use Getopt::Std;
my %options=();
getopts("i:", \%options);
#
open(IN, "<",$options{i}) || die ; 
my $seq_size;
my $totSize=0;

while(<IN>)
{
	if ( $_ =~ m/^>/)
	{
	}
	else
	{
		$seq_size= $_ =~ tr/A|T|C|G//;
		$totSize=$totSize + $seq_size;
	}
}
print $totSize;

