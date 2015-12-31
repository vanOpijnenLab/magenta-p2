#!/usr/bin/perl -w

#Margaret Antonio 15.12.26

#DESCRIPTION: Filter to identify which "essential" regions are actually "cold spots", 
# here defined as regions with 2000 bp before and after the designated region that also
# have low insertion rate (i.e. they are also defined as "essential" with significant 
# p-values)

use strict;
use Getopt::Long;
use warnings;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR

#ASSIGN INPUTS TO VARIABLES
our ($infile,$h,$size);
GetOptions(
'in:s' => \$infile,
'usage' => \$h,
'b'=>\$bracket,
'size'=>\$size,
'step'=>\$step,
);

if (!$size) { $size=500 };   #set the default sliding window size to 500
if (!$step) {$step=10};
if (!$bracket) { $bracket=2000 };
#GET DATA OUT OF FILE AND INTO 2D ARRAY
my @outArray;
open(DATA, '<', $infile) or die "Could not open '$infile' Make sure input .csv files are entered in the command line\n";
my $dummy=<DATA>;
while (my $entry = <DATA>) {
	chomp $entry;
	my @fields=split(",",$entry);
	push (@outArray,\@fields); #should be ref array
}
close DATA;

#Quick and dirty way: check back and forward 2000/step=2000/10=200 windows
#Better way to do it in case file is not from sliding windows: rec. until 2000 bef and aft

my $lastLine=scalar @outArray;
for (my $i=0;$i<$lastLine;$i++){
	my @roi=$outArray[$i];
	my ($backRes)=backRec(@outArray,$i);
	my ($forwardRes)=forwardRes(@outArray,$i);
	
	foreach (@roi) {
		foreach (@$_){
			print $_, "\t";
		}
		print "\n";
	}
	}
		
	

=begin comment

Step1: Do backwards and forwards checks
$roiStart=$start or $keys in %startHash
$roiEnd=$roiStart+$size

For backwards check all keys where $roiStart-2000 < $key < $roiStart
For forwards check all keys where $roiStart+$steo

Step2: Append result of both checks to %outHash value array for keys $start and $end corresponding
	
=end comment
=cut

