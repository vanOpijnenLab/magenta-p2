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
our ($in,$h,$out,$rIndex,$rCond);
GetOptions(
'i:s' => \$in,
'o:s' =>\$out,
'usage' => \$h,
'r:i' => \$rIndex,
'c:i' => \$rCond,
);


if (!$out) {$out="removalFilt.csv"};
if (!$rIndex) {$rIndex=2}; #Removal of index
if (!$rCond) {$rCond=0};    #If entry at index meets this condition

#GET DATA OUT OF FILE AND INTO 2D ARRAY
my @outArray;

sub cleaner{
	my $line=$_[0];
	chomp($line);
	$line =~ s/\x0d{0,1}\x0a{0,1}\Z//s; 
	return $line;
	}
	
my $csv = Text::CSV->new({ binary => 1, auto_diag => 1, eol => "\n"}) or die "Cannot use CSV: " . Text::CSV->error_diag();
open (my $outfile,'>',"$out")or die "Could not open '$out' \n";

open(DATA, '<', "$in") or die "Could not open '$in' \n";
	
my $line=<DATA>;
#print "This is the line: ",$line,"stoooooop";
$line=cleaner($line); #gets rid of carriage returns (^M)
my @header=split(',',$line);
push (@header,'csTest');
my $tick=0;
$csv->print($outfile, \@header); #header

while (my $entry = <DATA>) {
	#$tick+=1;
	#print $tick, "\t";
	$entry=cleaner($entry);
	#chomp($entry);
	#$entry =~ s/\x0d{0,1}\x0a{0,1}\Z//s; 
	my @fields=split(',',$entry);
    if ($fields[$rIndex]!=$rCond){
        #push (@outArray,\@fields); #should be ref array
        $csv->print($outfile, \@fields);
    }
}

close DATA;
close $outfile;