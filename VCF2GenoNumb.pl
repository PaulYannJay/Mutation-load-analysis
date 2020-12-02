#!/usr/bin/perl
use strict;
use warnings;
##############
# Reformate a vcf into a custom genotype file
## ####################
use Getopt::Std; 
use Scalar::Util qw(looks_like_number);

my %options=();
getopts("i:o:", \%options);
my $file=$options{i};
if ($file =~ /.gz$/) {
	open(VCF,"gunzip -c $file|") or die "can't open $file";
	}
	else {
	open(VCF,'<', $file) or die "can't open $file";
	}

open(OUT,'>', $options{o}) or die "can't open file"; 

my @Geno;
my $PIL;
my @Read;
while(<VCF>)
{
	if ($_ =~ m/^##/) # Skip vcf top header
	{
		next;
	}
	elsif ($_ =~ m/^#CHR/)  # Read and write header (name of sample)
	{
		@Read=split(/[\t\n\b]+/, $_);
		$PIL=@Read - 1;
		print OUT $Read[0], "\t", $Read[1],"\t", "Impact";
		for (my $i=9; $i<=$PIL; $i++)
		{
			print OUT "\t",$Read[$i];
		}
		print OUT "\n";
	}
	else
	{
		@Read=split(/[\t\n\b]+/, $_); 
		if ($Read[4] =~ m/,/) # Skip multiallelic snp
		{
			next;
		}
		print OUT $Read[0], "\t", $Read[1];
		if ($_ =~ m/HIGH/ || $_ =~ m/MODERATE/) #Base on SnpEff annotation, determine the effect of the snp
		{
			print OUT "\t", "missense"; #HIGH and MODERATE SnpEFF annotation correspond to missense variant
		}
		elsif ($_ =~ m/LOW/ )
		{
			print OUT "\t", "synonymous"; # LOW SnpEff annotation correspond to synonymous variant
		}
		elsif ( $_ =~ m/MODIFIER/)
		{
			print OUT "\t", "modifier"; #Variant in non coding region
		}
		else
		{
			die;
		}

		$PIL=@Read-1;
		for (my $i=9; $i<=$PIL; $i++) #For all individuals, write their genotype
		{
			@Geno=split(/[:\t\n\n]+/,$Read[$i]);
			if ($Geno[0] eq "./.")
			{
				print OUT "\t", ".";
			}
			elsif ($Geno[0] eq "0/0")
			{
				print OUT "\t", "0";
			}
			elsif ($Geno[0] eq "1/1")
			{
				print OUT "\t", "1";
			}
			elsif ($Geno[0] eq "2/2") #Since we do not consider biallelic snp, this is not usefull here, but I kept that here in case I add a biallelic option
			{
				print OUT "\t", "2";
			}
			elsif ($Geno[0] eq "3/3")
			{
				print OUT "\t", "3";
			}
			elsif ($Geno[0] eq "0/1")
			{
				print OUT "\t", "het01";
			}
			elsif ($Geno[0] eq "0/2")
			{
				print OUT "\t", "het02";
			}
			elsif ($Geno[0] eq "0/3")
			{
				print OUT "\t", "het03";
			}
			else
			{
				print OUT "\t", "HetSNP";
			}
		}
		print OUT "\n";
	}	
	
}
