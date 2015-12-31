#!/usr/bin/perl -w

#Margaret Antonio 15.12.26

#DESCRIPTION: Filter to identify which "essential" regions are actually "cold spots", 
# here defined as regions with 2000 bp before and after the designated region that also
# have low insertion rate (i.e. they are also defined as "essential" with significant 
# p-values)

use strict;
use Getopt::Long;
use warnings;

#USAGE from /8-essFilters
#perl ../Blueberries/test.pl --in slidingWindowF.csv

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR

#ASSIGN INPUTS TO VARIABLES
our ($infile,$h,$size, $bracket,$step,$defEss);
GetOptions(
'in:s' => \$infile,
'usage' => \$h,
'b'=>\$bracket,
'size'=>\$size,
'step'=>\$step,
'ess'=>\$defEss,
);

if (!$size) { $size=500 };   #set the default sliding window size to 500
if (!$step) {$step=10};
if (!$bracket) { $bracket=2000 };
if (!$defEss) {$defEss=.05} #define essentiality as meeting .05 pval cutoff
#GET DATA OUT OF FILE AND INTO 2D ARRAY
my @outArray;
open(DATA, '<', $infile) or die "Could not open '$infile' \n";
my $dummy=<DATA>;
while (my $entry = <DATA>) {
	chomp $entry;
	my @fields=split(",",$entry);
	push (@outArray,\@fields); #should be ref array
}
close DATA;

#Quick and dirty way: check back and forward 2000/step=2000/10=200 windows
#Better way to do it in case file is not from sliding windows: rec. until 2000 bef and aft

#Forward check function (recursively go through 2000 bp ahead to check if is "Essential")
function backRes($outArray,$roiIndex,$backStop){
	#dereference @outarray
	#roiIndex is going to be the index of the roi in @outArray
	#Need to go back ($bracket/$size)-1 windows/indices to cover bracket
	#OR Keep going back until reach forwardStop
	#Default status: assume  outside region is essential
	#If outside regions not ess, then not coldspot so stop looking and return status=1
	my $status=0;
	function isEssentialBack($outArray,$currentIndex, $status){
	    my @outArray=@{$outArray};
		my @currentRegion=@{$outArray[$currentIndex]};
		my $regionStart=$currentRegion[0];
		my $regionEss=$currentRegion[7]; #p-value for essentiality
		if ($regionStart<$backStop or $regionStart<=1){
			return $status;
			#Problem: for first window where roiStart=1, will return status=0
		}

		if ($regionEss>$defEss){
			$status=1; #because not cold spot----there is a non essential region
			return $status;
		}
		isEssentialBack(\@outArray,$currentIndex-1,$status);
	}
	my $resultBack=isEssentialBack(\@outArray,$roiIndex,$status);
	return $resultBack;
}

	
my $lastLine=scalar @outArray;
for (my $i=0;$i<$lastLine;$i++){
	my @roi=$outArray[$i];
	print "\n Start\tpvalue\tbackCold\n";
	
	#only check windows that are "essential" i.e. meet $defEss
	if ($roi[7]<$defEss){
	
		my $backStop=$roi[0]-2000;
    	my $backRes=backRec(\@outArray,$i,$backStop);
    	print $roi[0],"\t",$roi[7],"\t",$backRes,"\n";
		my $forwardStop=$roi[0]+2000+$size;
		#my $forwardRes=forwardRec(\@outArray,$i,$forwardStop);
		
	}
	
	
}
		
	

=begin comment

Step1: Do backwards and forwards checks
$roiStart=$start or $keys in %startHash
$roiEnd=$roiStart+$size

For backwards check all keys where $roiStart-2000 < $key < $roiStart
For forwards check all keys where $roiStart+$steo

Step2: Append result of both checks to %outHash value array for keys $start and $end corresponding
	
#Print tester	
	foreach (@roi) {
		foreach (@$_){
			print $_, "\t";
		}
		print "\n";
	}
=end comment
=cut

