#!/usr/bin/perl
use strict;
use warnings;
##############
#Script to compute MKT stat on a set of gene. It take in input the position of the gene and a geno file. It compute the number of synonymous and non synonymous fixed and polymorphic position for each gene. it output the mid position of gene, the number of each variant type, the dnds, the pnps and the pnps/dnds ratio (MKTest)
## ####################
use Getopt::Std; 
use Scalar::Util qw(looks_like_number);

my %options=();
getopts("i:o:g:", \%options);
my $file=$options{i};

open(GENE,'<', $options{g}) or die "can't open gene file"; #file containing the gene you are interested in. In  the format Scaffold\tStart\tEnd
open(OUT,'>', $options{o}) or die "can't open file"; 
if ($file =~ /.gz$/) { open(VCF,"gunzip -c $file|") or die "can't open $file"; }
else { open(VCF,'<', $file) or die "can't open $file"; }

my @Read;
my $size;
my @Pos;

print OUT "Scaf\t", "PosDeb","\t", "PosFin","\t", "FixMissense","\t", "FixSynonymous","\t", "HeteroMissense","\t", "HeteroSynonymous","\t","TotPos","\t","DoS", "\t", "dN/dS","\t","pN/pS","\t","MKTRatio", "\n";

my %ArVal;
my %ArPos;
my $GeneNum=1;
while(<GENE>) #Store the gene to analyse in a hashtab. Each gene is associated with a given number
{
	@Pos=split(/[\t\b\n]+/, $_);
	$ArVal{$GeneNum}=[0,0,0,0,0];#[FixNonSyn, FixSyn, HetNonSyn, HetSyn, TotPos] # To store the number of mutation of each type
	$ArPos{$GeneNum}=[$Pos[0],$Pos[1],$Pos[2]]; # To store gene info: Scaffold, Start, End
	$GeneNum++;
}

my $FocGene=1;
my $bool=0;
while(<VCF>) #Only for the header
	{
		@Read=split(/[\t\n\b]+/,$_); 
		$size=@Read; # Size fo the header *---> $size -3 = number of sample
		last;
	} 

while (<VCF>) # We slide along the genome, crossing gene position. for each SNP? we check if it is inside gene, and if it is the case, we incrememnt the mutation counter of this gene 
{
	$bool=0;
	@Read=split(/[\t\n\b]+/,$_); #Split the line
	if ( $Read[0] eq @{$ArPos{$FocGene}}[0] && $Read[1] > @{$ArPos{$FocGene}}[1] && $Read[1] < @{$ArPos{$FocGene}}[2]) #Check if the line is in the gene we are looking at
	{
		$bool=1;	#useless
	}
	else # The line is in a new gene, so we can print the data for the previous
	{
		HASHLOOP: for (keys %ArPos) #Find it which gene the SNP is
		{
			if ( $Read[0] eq @{$ArPos{$_}}[0] && $Read[1] > @{$ArPos{$_}}[1] && $Read[1] < @{$ArPos{$_}}[2])
			{
				$FocGene=$_; #Store the gene number
				$bool=1; #useless
				last HASHLOOP;
			}
		}
	}
	if ($bool == 1) #useles, always equal to 1...
	{
		my $indiv= join('', @Read[3..$#Read]); # To count the number of 1, 0 and heterozygous individual, concatenate the genotype of all sample
		if ( $_ =~ m/het/) #if the line contain an heterozygous individual, then this position is polymorph -> pNpS Counter
		{
			if ($Read[2] eq "missense") # Check the effect of the snp and increment the right variable
			{
				@{$ArVal{$FocGene}}[2]++; #pN counter
				@{$ArVal{$FocGene}}[4]++; #Total mutation nimber counter
			}
			elsif ($Read[2] eq "synonymous")
			{
				@{$ArVal{$FocGene}}[3]++; #pS counter
				@{$ArVal{$FocGene}}[4]++;
			}
			else  # the position is a Modifier
			{
				# Do nothing # 
			}
			next; #If one individual is hetero, the snp is not fixed. Go to next position ! 
		}
		my $N1 = () = $indiv =~ /1/g; # Number indiv homo 1/1
		my $N0 = () = $indiv =~ /0/g; # Number indiv homo 2/2
		my $Nno = () = $indiv =~ /\./g; # Number indiv non genotyped
		if ($Nno < 0.1*$size) # if more than 90% of individual genotyped
		{
			@{$ArVal{$FocGene}}[4]++;
			if ( ($N1 + $Nno) == ($size-3)) # Pop fixed for the alternative allele
			{
				if ($Read[2] eq "missense") # SNP missense
				{
					@{$ArVal{$FocGene}}[0]++; #dN counter
				}
				elsif ($Read[2] eq "synonymous") # SNP synonymous
				{
					@{$ArVal{$FocGene}}[1]++; #dS counter
				}
				else #SNP annotated has "modifier : intergenic or intron"
				{
					#Do nothing
				}
			}
			elsif ( ($N0 + $Nno) == ($size-3)) # If pop fixed for ancestal allele, do nothing
			{
					next;
			}
			else # Pop containing at least one alternate allele but not fixed for this one # If this occurs, that mean that at least one individual is homozygous for the alternative allele while the others are homozygous for the ancetral allele # We are facing a polymorphic position
			{
				if ($Read[2] eq "missense")
				{
					@{$ArVal{$FocGene}}[2]++; #pN counter
				}
				elsif ($Read[2] eq "synonymous")
				{
					@{$ArVal{$FocGene}}[3]++; #pS counter
				}
				else #SNP annotated has "modifier" : intergenic or intron, do Nothing
				{
					#next;
				}
			}

		}
	}
}

my $dnds;
my $pnps;
my $mktRatio;
my $dos;
for (keys %ArVal) # For each gene, print the data
{	
	if ((@{$ArVal{$_}}[0]+ @{$ArVal{$_}}[1]) > 0 &&  (@{$ArVal{$_}}[2]+ @{$ArVal{$_}}[3]) > 0)
	{
		$dos=@{$ArVal{$_}}[0]/(@{$ArVal{$_}}[0] + @{$ArVal{$_}}[1]) - @{$ArVal{$_}}[2]/(@{$ArVal{$_}}[2]+@{$ArVal{$_}}[3]); #DoS
	}
	else
	{
		$dos="NA";
	}
	if (@{$ArVal{$_}}[0] > 0 &&  @{$ArVal{$_}}[1]> 0 && @{$ArVal{$_}}[2]> 0 && @{$ArVal{$_}}[3] > 0 )
	{
		$dnds=@{$ArVal{$_}}[0]/@{$ArVal{$_}}[1];#dN/dS
		$pnps=@{$ArVal{$_}}[2]/@{$ArVal{$_}}[3];#pN/pS
		$mktRatio=$pnps/$dnds; #MKT ration
	}
	else
	{
		$dnds="NA";
		$pnps="NA";
		$mktRatio="NA";
	}
	print OUT @{$ArPos{$_}}[0], "\t",@{$ArPos{$_}}[1], "\t",@{$ArPos{$_}}[2], "\t",@{$ArVal{$_}}[0],  "\t",@{$ArVal{$_}}[1], "\t",@{$ArVal{$_}}[2], "\t",@{$ArVal{$_}}[3],"\t",@{$ArVal{$_}}[4],"\t", $dos, "\t", $dnds,"\t", $pnps,"\t", $mktRatio,"\n";
}
		
		
