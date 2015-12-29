#!/usr/bin/perl -w

#Margaret Antonio 15.12.26

#DESCRIPTION: Filter to identify which "essential" regions are actually "cold spots", 
# here defined as regions with 2000 bp before and after the designated region that also
# have low insertion rate (i.e. they are also defined as "essential" with significant 
# p-values)

#use strict;
use Getopt::Long;
use warnings;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR

#ASSIGN INPUTS TO VARIABLES
our ($infile,$h,$size);
GetOptions(
'in:s' => \$infile,
'usage' => \$h,
);

if (!$size) { $size=500 };   #set the default sliding window size to 500

#Two hashes: one where key:value is start:pval and other where key:value is end:pval
my @outArray;

open (IN, $infile) or die "can't open $infile: $!\n";

while (<IN>){
	chomp; 							#get rid of white spaces
	my @fields=split ",";
    push (@outArray,\@fields); #should be ref array
}

close IN;

my $lastLine=scalar @outArray;

for (my $i=0;$i<$lastLine;$i++){
	my @roi=$outArray[$i];
	print "$i ";
}
	
#Step1: Do backwards and forwards checks
# $roiStart=$start or $keys in %startHash
# $roiEnd=$roiStart+$size

#For backwards check all keys where $roiStart-2000 < $key < $roiStart
#For forwards check all keys where $roiStart+$steo

#Step2: Append result of both checks to %outHash value array for keys $start and $end corresponding
	


