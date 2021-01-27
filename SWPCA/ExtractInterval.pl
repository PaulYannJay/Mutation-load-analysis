#!/usr/bin/perl
### A code to extract the scaffold list, end and start from a vcf file
use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts(":i:o:", \%options);

if ($options{i} =~ /.gz$/) { open(INPUT, "gunzip -c $options{i} |") || die "can’t open pipe to $options{i}"; } 
else { open(INPUT, "<", $options{i}) || die "no input file"; }

open(OUTPUT, ">", $options{o}) || die "impossible open file";

my %Scaf;
my @Read;
my $deb;
my $fi;
my $chrom;

while (<INPUT>) #Only First line
{
	@Read=split(/\s/, $_);		#Put the line into an array
	$Scaf{$Read[0]}=$Read[1]; #put the scaffold in the hash table, with the value equal to the first position
	$chrom=$Read[0];
	last;
}

while (<INPUT>) # Start a the second line
{
	@Read=split(/\s/, $_);		#Put the line into an array
	if (exists( $Scaf{$Read[0]})) #If the scaffold is already defined in the hash, that mean that it is not the first position
		{
		$fi=$Read[1]; #Store the position in case it is the last one
		next;
		}
	else # if not, that mean that the previous line was the last position of the previous scaffold
		{
		print OUTPUT $chrom, "\t", $Scaf{$chrom}, "\t", $fi, "\n"; #Print the name of the scaffold, its first position and its last position
		$Scaf{$Read[0]}=$Read[1]; #put the scaffold in the hash table, with the value equal to the first position
		$fi=$Read[1]; #Store the position in case it is the last one
		$chrom=$Read[0];#Store the scaffold name
		}
}		
