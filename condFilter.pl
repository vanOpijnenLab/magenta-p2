#!/usr/bin/perl -w

#Margaret Antonio Began: 15.12.26 Updated: 16.01.11

#DESCRIPTION: Filter to remove windows (lines) that satisfy the specified condition
# Where i=number, conditions include <i, >i <=i, >=i, ==i, !=i, etc.

#allow multiple ri rc!!!!

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use diagnostics;

#USAGE from /8-essFilters
#perl ../Blueberries/test.pl --in slidingWindowF.csv

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR

#ASSIGN INPUTS TO VARIABLES
our ($in,$h,$out,$kIndex,$kCond,$kVal);
GetOptions(
'i:s' => \$in,
'o:s' =>\$out,
'help' => \$h,
'ki:i' => \$kIndex,  #keep index #any number that is an index
'kc:s' => \$kCond,   #keep operator condition  Can be e,l,g,le,ge
'kv:i' => \$kVal,     #keep value #any number
);
if ($kCond eq "e"){
    $kCond="==";
}
elsif ($kCond eq "l"){
    $kCond="<";
}
elsif ($kCond eq "g"){
    $kCond=">";
}
elsif ($kCond eq "le"){
    $kCond="<=";
}
elsif ($kCond eq "ge"){
    $kCond=">=";
}
elsif ($kCond eq "ne"){
    $kCond="!=";
}


#if ($h){
    # print "Help information coming soon..." and die "Goodbye..."; #must fix this!!!
    # removal filter can be < > <= >= == or !=
    # }

#SET DEFAULTS IF NOTHING SPECIFIED FOR REQUIRED OPTIONS
if (!$out) {$out="filt.csv"};
if (!$kIndex) {$kIndex=2};      #Removal of index
if (!$kCond) {$kCond="==0"};    #If entry at index meets this condition

#GET DATA OUT OF FILE AND INTO 2D ARRAY
my @outArray;

sub cleaner{
	my $line=$_[0];
    #print "This is it ", $_[0];
	chomp($line);
	$line =~ s/\x0d{0,1}\x0a{0,1}\Z//s; 
	return $line;
	}

open (OUT,'>',"$out")or die "Could not open '$out' \n";

open(DATA, '<', "$in") or die "Could not open '$in' \n";
	
my $head=<DATA>;
$head=cleaner($head); #gets rid of carriage returns (^M)
my @header=split(',',$head);
print OUT (join(",",@header),"\n"); #header

while (my $entry = <DATA>) {
    #print $entry;
	$entry=cleaner($entry);
	my @fields=split(',',$entry);
    my $testVal=$fields[$kIndex];
    if (eval($testVal.$kCond.$kVal)){
        #print ($testVal.$kCond.$kVal);
        print OUT (join(",",@fields),"\n");
    }
}
close DATA;
close OUT;

