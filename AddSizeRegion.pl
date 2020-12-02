#!/usr/bin/perl
#strict used to add a "Region size" column in the repeat masker out table, for ploting with R for instance;
use strict;
use warnings;
use Getopt::Std;
my %options=();
getopts("i:o:s:", \%options);
open(IN, "<",$options{i}) || die ; 
open(SIZE, "<", $options{s}) || die;
open(OUT, ">",$options{o}) || die ; 

my %HashSize;
while(<SIZE>) #Keep the regions size in a Hashtable
{
	my @Read=split(/[\s\t\b\n]+/, $_);
	$HashSize{$Read[0]}=$Read[1];
}

while(<IN>) #Print Header
{
	chomp $_;
	print OUT $_, "\t", "RegionSize", "\n";
	last;
}

while(<IN>)
{
chomp $_;
my @TE=split(/[\t\b\n]+/, $_);
my $pos=$TE[17];
print OUT $_, "\t", $HashSize{$pos}, "\n"; #Read the region size in the hashTable and print it
}
