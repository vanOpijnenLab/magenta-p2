#!/usr/bin/perl -w

#Margaret Antonio 15.03.27

use strict; use warnings; use autodie;
my $files = $#ARGV; # get number of files - 1
while (my $file = shift @ARGV) {
    open my $fh, "<", $file;
    <$fh> unless $files == @ARGV; # discard header unless first file
    print while <$fh>; # output the rest
}



#In terminal, use command $ perl ../script/mergeCSV.pl results/L1_2394eVI_PennG.csv results/L2_2394eVI_PennG.csv results/L3_2394eVI_PennG.csv results/L4_2394eVI_PennG.csv results/L5_2394eVI_PennG.csv results/L6_2394eVI_PennG.csv >compiledPennG.csv

#http://stackoverflow.com/questions/17993919/merge-multiple-csv-files-using-perl


#Then sort the merged file:
#$ sort -t, -k1 -n compiledPennG.csv > sortedCompiledPennG.csv

#This can be used for the sliding window calculation script (window.pl)


#IF MULTIPLE FILES WERE GIVEN--->NEED TO COMPILE THEM INTO ONE AND THEN SORT THE COMPILED FILE

 my $mergefile;
 my $files = $#ARGV; # get number of files - 1
 my $csvMerge = Text::CSV->new();
 my @comRows;
 use constant POS => 0;
 my $handle;
 while ($mergefile = shift @ARGV) {
 open $handle, "<", $mergefile;
 <$handle>; #no header---messes up with sorting--> #unless $files == @ARGV; # discard header unless first file
 while(my $line_ref=$csvMerge->getline($handle)){
 push @comRows, $line_ref;
 }
 print while <$handle>; # output the rest
 }
 
 #NOW SORT THE MERGED FILE,
 @comRows = sort { $a->[POS] <=> $b->[POS] } @comRows;
 for my $line_ref (@comRows) {
 $csvMerge->combine(@$line_ref);
 print $csvMerge->string(), "\n";
 }
 close $handle;