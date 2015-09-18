#!/usr/bin/perl -w

#Margaret Antonio 15.03.27

use strict; use warnings; use autodie;
my $files = $#ARGV; # get number of files - 1
while (my $file = shift @ARGV) {
    open my $fh, "<", $file;
    <$fh> unless $files == @ARGV; # discard header unless first file
    print while <$fh>; # output the rest
}



#In terminal, use command $ perl ../script/mergeCSV.pl results/L1_2394eVI_PennG.csv results/L2_2394eVI_PennG.csv results/L3_2394eVI_PennG.csv results/L4_2394eVI_PennG.csv results/L5_2394eVI_PennG.csv results/L6_2394eVI_PennG.csv >compiled.csv