#!/usr/bin/perl -w

#Margaret Antonio
#Began: 15.12.26 Updated: 16.01.17
#DESCRIPTION: Filter to keep (lines) that satisfy the specified condition
#Improvements to make: allow multiple indices and conditions for filtering

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use diagnostics;

#ASSIGN INPUTS TO VARIABLES
our ($h,$out,$kIndex,$kOper,$kVal);
GetOptions(
'o:s' =>\$out,
'h' => \$h,
'ki:i' => \$kIndex,   #keep index #any number that is an index (i for integer)
'ko:s' => \$kOper,    #keep operator condition. Can be e,l,g,le,ge,ne. (s for string)
'kv:f' => \$kVal,     #keep value #any number (f for float)
);


sub print_usage() {
    print "\n";
    print "condFilter.pl: Conditionally filter out lines in a csv file\n\n";
    print "DESCRIPTION: Takes a single file and keeps only lines that satisfy condition\n";
    print "the difference in mean fitness, the pval for each gene.\n";
    print "Example: Keep all lines in a comparison file where the difference between\n";
    print "\tmean fitnesses is greater than .15\n";
    print "\nUSAGE: perl condFilter.pl -ki <index of col> -ko <e,l,g,le,ge,or ne>\n";
    print "\t-kv <keep value> <options> <file.csv>\n\n";
    print "OPTIONS:\n\n";
    print "-ki\t Index of column to test keep value and keep operator. Indices start at 0.\n";
    print "-ko\t String version of operator to use (, l, g, le, ge, ne).\n";
    print "-kv\t Value to use in conditional test.\n";
    print " -h\t Prints usage and quits\n";
    print " -o\t Output file for filtered data. Default: filt.csv\n";
    print "\n\n";
}
#IF HELP FLAG (-h) WAS SPECIFIED THEN PRINT USAGE AND QUIT
if ($h){
    print_usage();
    exit;
}

#INTERPRET KEEP CONDITION FROM STRING TO OPERATOR. AND CHECK FOR OTHER PARTS OF CONDTIONAL TEST.
#Can't directly use operator because > is interpreted as output into (> file)

if ($kOper eq "e"){
    $kOper="==";
}
elsif ($kOper eq "l"){
    $kOper="<";
}
elsif ($kOper eq "g"){
    $kOper=">";
}
elsif ($kOper eq "le"){
    $kOper="<=";
}
elsif ($kOper eq "ge"){
    $kOper=">=";
}
elsif ($kOper eq "ne"){
    $kOper="!=";
}
else {
    print "\nPlease enter one of the following for the keep operator (-ko): e, l, g, le, ge, ne\n\n";
    print_usage();
    exit;
}
if (!$kVal){
    print "\nPlease enter the value to use in the conditional test for filtering. -kv <value> \n\n";
    print_usage();
    exit;
}
if (!$kIndex){
    print "\nPlease enter the column index to use for filtering. -ki <index> \n\n";
    print_usage();
    exit;
}
if (!$out) {$out="filt.csv"};

sub cleaner{
    my $line=$_[0];
    chomp($line);
    $line =~ s/\x0d{0,1}\x0a{0,1}\Z//s;
    return $line;
}

#READ INPUT FILE, FILTER, AND IMMEDIATELY WRITE TO THE OUTPUT FILE
my @outArray;

my $in=$ARGV[0];
open (OUT,'>',"$out")or (print_usage() and die "Could not open '$out' \n");
open(DATA, '<', "$in") or (print_usage() and die "Could not open '$out' \n");

#REUSE the header (column names) from the input file for the output file
my $head=<DATA>;
$head=cleaner($head);
my @header=split(',',$head);
print OUT (join(",",@header),"\n");

while (my $entry = <DATA>) {
	$entry=cleaner($entry);
	my @fields=split(',',$entry);
    my $testVal=$fields[$kIndex];
    if (eval($testVal.$kOper.$kVal)){
        print OUT (join(",",@fields),"\n");
    }
}
close DATA;
close OUT;

