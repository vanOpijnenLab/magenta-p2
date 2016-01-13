#!/usr/bin/perl -w

#Margaret Antonio 15.12.26

#DESCRIPTION: After s
use Data::Dumper;
use strict;
use Getopt::Long;
use warnings;
use File::Path;
use File::Basename;

#ASSIGN INPUTS TO VARIABLES
our ($indir,$h,$out,$sortby,$round,$labels);
GetOptions(
'indir:s' => \$indir,
'o:s' =>\$out,
's:i' => \$sortby,
'r:i'=> \$round,
'l:s'=> \$labels,
);

if (!$out) {$out="compGenes.csv"}
if (!$round){$round='%.4f'}



sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}

my @files;
if ($indir){
    my $directory="$indir";
    opendir(DIR, $directory) or die "couldn't open $directory: $!\n";
    my @direct= readdir DIR;
    my $tail=".csv";
    foreach (@direct){
        if (index($_, $tail) != -1){
            $_=$indir.$_;
            push (@files,$_);
        }
    }
    closedir DIR;
}
else{
    @files=@ARGV;
}
my $num=(scalar @files);


my @header;
my @sample;
push (@header,"gene");
if ($labels){
    @sample=split(',',$labels);
}
else{
    #get rid of .csv part in file names and use as labels
    foreach (@files){
        my @temp=split('\\.',$_);
        my $colName=$temp[0];
        push (@sample,$colName);
    }
}
push (@header,@sample);
push (@header,"abs(diffMean)");


my %all;
for (my $i=0; $i<$num; $i++){   #Read files from ARGV
    print "File #",$i+1,"\t";
    
    my $file=$files[$i];
    print $file,"\n";
    
    open(DATA, '<', $file) or die "Could not open '$file' Make sure input .csv files are entered in the command line\n";
    my $dummy=<DATA>; #the header in the file
    #my %hash;
    while (my $entry = <DATA>) {
        chomp $entry;
        my @line=split(",",$entry);
        my $gene=$line[0];
        chomp($gene);
        my $mean = sprintf("$round",$line[1]);
        
        ############################## hash assignment
        if(!exists $all{$gene}){
            my @means;
            push (@means,$mean);
            $all{$gene}=\@means;
        }
        else{
            my @means=@{$all{$gene}};
            #for (my $j=scalar @means+1;$j<$i;$j++){
            #    push (@means,"NA");
            #}
            push (@means,$mean);
            my $diff=sprintf("$round",abs($means[0]-$means[1]));
            push (@means,$diff);
            $all{$gene}=\@means;
        }
       
    }
    
    close DATA;
}


#PRINT TO OUT FILE

#Set up the header (column names)


my @unsorted;

foreach my $entry (sort {$a cmp $b}  keys %all) {
    my @means=@{$all{$entry}};
    #print OUT $entry,",";      #if unhashed, must go after OPEN
    #print OUT (join(",",@means),"\n");
    my @temp=($entry);
    push (@temp,@means);
    push (@unsorted,\@temp);
}
my @sorted = sort { $b->[3] <=> $a->[3] } @unsorted;

open OUT, '>',"$out";
print OUT (join(",",@header),"\n");
foreach (@sorted){
    print OUT (join(",",@{$_}),"\n");
}
close OUT;


