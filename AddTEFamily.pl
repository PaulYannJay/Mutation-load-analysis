#!/usr/bin/perl
#strict used to add a "family" column in the Repeat Masker out table, just based on TE name, for ploting with R for instance;
use strict;
use warnings;
use Getopt::Std;
my %options=();
getopts("i:o:", \%options);
#
open(IN, "<",$options{i}) || die ; 
open(OUT, ">",$options{o}) || die ; 

while(<IN>)
{
	chomp $_;
	print OUT $_, "\t", "Family", "\n";
	last;
}

while(<IN>)
{
	chomp $_;
	if ($_ =~ m/DNA/)
	{
		print OUT $_, "\t", "DNA\n";
	}
	elsif ($_ =~ m/LINE/)
	{
		print OUT $_, "\t", "LINE\n";
	}
	elsif ($_ =~ m/SINE/)
	{
		print OUT $_, "\t", "SINE\n";
	}
	elsif ($_ =~ m/LTR/)
	{
		print OUT $_, "\t", "LTR\n";
	}
	elsif ($_ =~ m/RC/)
	{
		print OUT $_, "\t", "RC\n";
	}
	else
	{
		print OUT $_, "\t", "Other\n";
	}
	
}
