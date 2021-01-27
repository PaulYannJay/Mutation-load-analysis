#!/bin/bash
##########
##
##
##	This workflow allow to create data table to produce the Sliding window PCA plot 
## 	Then, you can just plot the first axis of the PCA along the position to produce a Figure similar de Fig. 1B
## 
##
#######################################

while [ $# -gt 0 ] ; do
  case $1 in
    -vcf ) VCF="$2" ;echo "the vcf file : $VCF" >&2;;
    -s | --sample) Sample="$2" ;echo "the file containing the list of samples to analyse is $Sample" >&2;;
    -r | --region) Scaf="$2" ;echo "the file containing the list of scaffolds and regions to analyse is $Scaf" >&2;;
    -o | --output) OUT="$2" ;echo "the output prefix is $OUT" >&2;;
    -w | --window) Wind="$2" ;echo "the window size is $Wind" >&2;;
	-h | --help) echo -e "Option required:
-vcf \t the vcf file to analyse
-o/--output \t Output Prefix
-w/--window \t the window size to use for sliding window
Optional: 
-s/--sample \t set of samples to analyse. If not defined, all sample are used
-r/--region \t set of regions to analyse in format : Scaffold\tDebut\tFin. If not defined, all regions are analysed " >&2;exit 1;; 
  esac
  shift
done

if [ -z "$OUT" ] || [ -z "$VCF" ] || [ -z "$Wind" ]; then
	echo >&2 "Fatal error: Ouput, Vcf or Window size are not defined"
exit 2
fi

SEQ=$VCF

echo -e "\n"

mkdir $OUT

if [ -n "$Scaf" ] && [ -n "$Sample" ]; then # If Region and sample specified
	echo -e "Filtering vcf file using only samples present in $Sample and regions present in $Scaf"
	tabix -p vcf $SEQ 
	bcftools view -R $Scaf -S $Sample -O z -o $OUT/FilteredDataset.vcf.gz $VCF
	tabix -p vcf $OUT/FilteredDataset.vcf.gz
	SEQ=$OUT/FilteredDataset.vcf.gz
elif [ -n "$Scaf" ] && [ -z "$Sample" ] ; then # If Region specified
	echo -e "Filtering vcf file using only regions present in $Scaf "
	tabix -p vcf $SEQ 
	bcftools view -R $Scaf -O z -o $OUT/FilteredDataset.vcf.gz $VCF
	tabix -p vcf $OUT/FilteredDataset.vcf.gz
	SEQ=$OUT/FilteredDataset.vcf.gz
elif [ -z "$Scaf" ] && [ -n "$Sample" ] ; then # If Sample specified
	echo -e "Filtering vcf file using only samples present in $Sample "
	bcftools view -S $Sample -O z -o $OUT/FilteredDataset.vcf.gz $VCF
	tabix -p vcf $OUT/FilteredDataset.vcf.gz
	SEQ=$OUT/FilteredDataset.vcf.gz 
else #If no sample and region specified
	echo -e "Not filtering the vcf : no region and sample file defined "
	tabix -p vcf $SEQ 
fi

if [ -z "$Scaf" ]; then #If no region are specified, create a region file containing all positions !
	bcftools query -f '%CHROM %POS\n' $SEQ > Position.rm
	./ExtractInterval.pl -i Position.rm -o AllRegion.rm
	Scaf="AllRegion.rm"
fi

[ -e $OUT.PCA.txt ] && rm $OUT.PCA.txt
	
while read line ; do 
	stringarray=($line) #Put the content of $Scaf in a array, line by line
	Chr=${stringarray[0]} #Chromosome
	deb=${stringarray[1]} #start positionn
	fin=${stringarray[2]} #End position
	nbWindow=$((($fin-$deb)/$Wind)) #Allow to compute the number of entire window
	for i in `seq 1 $nbWindow`; do #For each window
		posdeb=$(((($i-1)*$Wind)+1))		#Define start and end position
		posfin=$(($i*$Wind))		
		bcftools view -r $Chr:$posdeb-$posfin -O z -o $OUT/$Chr:$posdeb-$posfin.FilteredDataset.vcf.gz $SEQ
		Rscript PCAsnpSA.R	$OUT/$Chr:$posdeb-$posfin.FilteredDataset.vcf.gz $Chr $posdeb $posfin $OUT.PCA.txt >/dev/null 2>&1 
		echo "$Chr	$posdeb	$posfin"
	done		
done < $Scaf
