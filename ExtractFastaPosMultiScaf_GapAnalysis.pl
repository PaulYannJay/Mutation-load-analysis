#!/usr/bin/perl

##################
##
##	Allow to extract fasta position using position in a file. In order to allow to extract protein fasta based on a CDS file, I have had a $order option describid the sens of the protein (reverse complement). In any other situation, order should be "+".
#
##################


use warnings;
use strict;

use Getopt::Std;

my %options=();
getopts("f:i:o:", \%options);

my $help="
\n\n
\ti = input fasta file 
\tf = file containing the position 
\to = output 
 ";

	  if (defined $options{o}) {print "\n\noutput file : $options{o}\n";} else {die("\nthe 'o' option was not defined. \n $help \n");}
	  if (defined $options{i}) {print "Input file : $options{i}\n";} else {die("\nthe 'i' option was not defined. \n $help \n");}
	  if (defined $options{f}) {print "End position : $options{f} \n\n\n";} else {die("\nthe 'f' option was not defined. \n $help \n");}


open(FASTA, '<', $options{i}) || die "can't open fasta file";
open(OUTPUT, '>', $options{o}) || die "can't write in output fasta file";
open(POS, '<', $options{f}) || die "can't open the position file";

my %hash;
my $key;
my $sc;
my @Pos;
my @Sc;
my $scaf;

while(<POS>) # Store the position to extract in a hashTab
{
	my @Read=split(/[\b\t\s\n]+/, $_);
	$hash{$Read[2]}{$Read[0]}=$Read[1]-$Read[0]+1;
	#$hash{$Read[1]}{$Read[2]}=$Read[3]-$Read[2]+1;
}


while (<FASTA>) # Read the (one line) fasta file
{
	if ($_ =~ m/^>/)
	{
		$scaf=$_;
		$scaf =~ s/>//;
		$scaf =~ s/ .*$//;
	}
	else
	{
		my $seq=$_;
		@Sc=keys %hash;
		foreach $sc (@Sc)
		{
			if ($sc == $scaf)
			{
				foreach $key (sort keys %{$hash{$sc}})
				{
					my $sub = substr $seq, ($key-1), $hash{$sc}{$key};
					print OUTPUT ">",$sc, "_",$key , "-", $key+$hash{$sc}{$key},"\n";
					print OUTPUT $sub, "\n";
				}
			}
		}
	}
}

