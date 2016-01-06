#!/usr/bin/perl -w

#Margaret Antonio 15.12.26

#DESCRIPTION: Filter to identify which "essential" regions are actually "cold spots", 
# here defined as regions with 2000 bp before and after the designated region that also
# have low insertion rate (i.e. they are also defined as "essential" with significant 
# p-values)

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
no warnings 'recursion';

#USAGE from /8-essFilters
#perl ../Blueberries/test.pl --in slidingWindowF.csv

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR

#ASSIGN INPUTS TO VARIABLES
our ($infile,$h,$size, $bracket,$step,$defEss,$out);
GetOptions(
'in:s' => \$infile,
'out:s' =>\$out,
'usage' => \$h,
'b'=>\$bracket,
'size'=>\$size,
'step'=>\$step,
'ess'=>\$defEss,
);

if (!$size) { $size=500 };   #set the default sliding window size to 500
if (!$step) {$step=10};
if (!$bracket) { $bracket=2000 };
if (!$defEss) {$defEss=.05}; #define essentiality as meeting .05 pval cutoff
if (!$out) {$out="out.csv"};
#GET DATA OUT OF FILE AND INTO 2D ARRAY
my @outArray;
open(DATA, '<', $infile) or die "Could not open '$infile' \n";

my $line=<DATA>;
#print "This is the line: ",$line,"stoooooop";
chomp($line);
my @header=split(",",$line);
push (@header,'csTest');

while (my $entry = <DATA>) {
	chomp ($entry);
	my @fields=split(",",$entry);
	push (@outArray,\@fields); #should be ref array
}
close DATA;

#Quick and dirty way: check back and forward 2000/step=2000/10=200 windows
#Better way to do it in case file is not from sliding windows: rec. until 2000 bef and aft

#Forward check function (recursively go through 2000 bp ahead to check if is "Essential")



	#dereference @outarray
	#roiIndex is going to be the index of the roi in @outArray
	#Need to go back ($bracket/$size)-1 windows/indices to cover bracket
	#OR Keep going back until reach forwardStop
	#Default status: assume  outside region is essential
	#If outside regions not ess, then not coldspot so stop looking and return status=1


sub isEssentialBack{
	my ($wholeArray,$currentIndex,$status,$stop)=@_;
	if ($currentIndex==0){
		return $status;
		}
	my @outArray=@{$wholeArray};
	my @currentRegion=@{$outArray[$currentIndex]};
	my $regionStart=$currentRegion[0];
	my $regionEss=$currentRegion[7]; #p-value for essentiality
	if ($regionStart<$stop){
		return $status; #Problem: for first window where roiStart=2, will return status=1
	}
	elsif ($regionEss>$defEss){return 2}
	else{isEssentialBack(\@outArray,$currentIndex-1,$status,$stop)}
}

sub isEssentialForward{
	my ($wholeArray,$currentIndex,$status,$stop)=@_;
	my @outArray=@{$wholeArray};
	my $last=scalar @outArray-1;
	if ($currentIndex==$last){return $status}
	my @currentRegion=@{$outArray[$currentIndex]};
	my $regionStart=$currentRegion[0];
	my $regionEss=$currentRegion[7]; #p-value for essentiality
	if ($regionEnd>$stop){
		return $status; #Problem: for first window where roiStart=2, will return status=1
	}
	if ($regionEss>$defEss){   #Does this ever execute? Was error that wasn't caught
		return 2;
	}
	else{
		isEssentialForward(\@outArray,$currentIndex+1,$status,$stop);
	}
}

my @finalOutput;
my $lastLine=scalar @outArray;

for (my $roiIndex=0;$roiIndex<$lastLine;$roiIndex++){
	my @roi=@{$outArray[$roiIndex]};
	my $backResult=0; #Default no test, then backResult=0
	my $forwardResult=0;
	#only check windows that are "essential" i.e. satisfy $defEss
	if ($roi[7]<$defEss){
		my $backStop=$roi[0]-2000;
		my $forwardStop=$roi[1]+2000;
    	$backResult=isEssentialBack(\@outArray,$roiIndex,1,$backStop);
    	$forwardResult=isEssentialForward(\@outArray,$roiIndex,1,$forwardStop);
	}
	push (@roi,$backResult);
	push (@finalOutput,\@roi);
}

my $csv = Text::CSV->new({ binary => 1, auto_diag => 1, eol => "\n"}) or die "Cannot use CSV: " . Text::CSV->error_diag(); 
open (my $OUT, ">$out"); 
foreach (@header){
	print $_,"\t";
	}
	print "\n";
$csv->print($OUT, \@header); #header


foreach (@finalOutput){
	$csv->print ($OUT, $_);
	
}
close $OUT;

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

