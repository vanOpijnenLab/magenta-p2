#!/usr/bin/perl -w

#Margaret Antonio 16.01.13

#USAGE: perl compGenes.pl <options> <aggregateFile1.csv aggregateFile2.csv> OR -indir <indir/>

use Data::Dumper;
use strict;
use Getopt::Long;
use warnings;
use diagnostics;
use File::Path;
use File::Basename;
use Statistics::Distributions;

#ASSIGN INPUTS TO VARIABLES USING FLAGS
our ($indir,$h,$out,$sortkey,$round,$l);
GetOptions(
'd:s' => \$indir,
'h' => \$h,
'o:s' =>\$out,
's:i' => \$sortkey,
'r:i'=> \$round,
'l:s'=> \$l,
);

sub print_usage() {
    print "\n";
    print "compGenes.pl: COMPARE GENES FROM A TN-SEQ EXPERIMENT USING AGGREGATE FILES\n\n";
    print "DESCRIPTION: Takes two aggregate.pl outputs and compares them by calculating\n";
    print "the difference in mean fitness, the pval for each gene.\n";
    print "Example: compare control vs antibiotic, where both from same strain (genome).\n";
    print "For different strains/genomes, use compStrains.pl\n";
    print "\nUSAGE: perl compGenes.pl <options> <aggregateFile1.csv aggregateFile2.csv>\n";
    print "       or -d <directory/> \n\n";
    print "OPTIONS:\n\n";
    print " -h\tPrints usage and quits\n";
    print " -d\tDirectory containing input files. Make sure / is included after name\n";
    print " -o\tOutput file for comparison data. Default: compFile1File2.csv\n";
    print " -s\tSort output by this index of the file (indices begin at 0). Default: by mean\n";
    print " -r\tRound final output numbers to this number of decimals\n";
    print " -l\tLabels for for compared files. Used for column names and default output file name.\n";
    print "   \tTwo strings, comma separated (i.e. -l gluc,dapto). Order should match file order.\n";
    print "\n\n";
}
if ($h){
    print_usage();
    exit;
}

#THE @files ARRAY WILL CONTAIN INPUT FILE NAMES, EXTRACTED FROM A DIRECTORY (-indir) OR ARGV
my @files;
if ($indir){
    my $directory="$indir";
    opendir(DIR, $directory) or die "Couldn't open $directory: $!\n";
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

#GET LABELS: USE (-l) OR USE FILNEAMES AS LABELS FOR COLUMNS IN OUTPUT FILE

my @labels;
if ($l){
    @labels=split(',',$l);
}
else{
    foreach (@files){
       my @temp=split('\\.',$_);
       my $colName=$temp[0];
       push (@labels,$colName);
   }
}

#CHECK IF REQ. VARIABLES WERE DEFINED USING FLAGS. IF NOT THEN USE DEFAULT VALUES

if (!$out) {$out="comp-".$labels[0].$labels[1].".csv"}
if (!$round){$round='%.4f'}

#READ INPUT FILES AND PUT GENE INFO FROM EACH LINE INTO HASH %all WHERE KEY=GENE_ID AND
#VALUE=GENE INFO (MEAN FITNESS, # INSERTIONS, ETC.)

my %all;
my @header;
for (my $i=0; $i<2; $i++){
    print "File #",$i+1,"\t";
    my $file=$files[$i];
    print $file,"\n";
    
    open(DATA, '<', $file) or die "Could not open '$file'\n";
    
    #EXTRACT THE HEADER (COLUMN NAMES) OF THE FILE AND KEEP FOR OUTPUT HEADER
    #APPEND FILE NAME OR GIVEN LABEL (-l) TO COLUMN NAME SO ITS DIFFERENT FROM OTHER INPUT FILE
    
    my $head=<DATA>;
    my @cols=split(',',$head);
    @cols=@cols[0,1,2,3,4,5,6];
    for (my $j=0;$j<scalar @cols;$j++){
        $cols[$j]=$cols[$j].'-'.$labels[$i];
    }
    push (@header,@cols);
    while (my $entry = <DATA>) {
        chomp $entry;
        my @line=split(",",$entry);
        if (!$line[5]){
            $line[5]="NA";
        }
        if (!$line[6]){
            $line[6]=0;
        }
        @line=@line[0,1,2,3,4,5,6];
        my $gene=$line[0];
        chomp($gene);
        
        #PUT GENE AND INFO INTO THE HASH FOR EXISTING KEY (1ST FILE) OR CREATE NEW KEY (2ND FILE)
        
        if(!exists $all{$gene}){
            my @info;
            push (@info,@line);
            $all{$gene}=\@info;
        }
        else{
            my @info=@{$all{$gene}};
            push (@info,@line);
            my $diff=sprintf("$round",($info[1]-$info[8]));
            my $total1=$info[6];
            my $total2=$info[13];
            my $sd1=$info[3];
            my $se1=$info[4];
            my $sd2=$info[10];
            my $se2=$info[11];
            my $df=$total1+$total2-2;
            my $tdist;
            my $pval;
            
            # CHECK TO MAKE SURE ALL VARIABLES IN TDIST,PVAL CALCULATIONS ARE NUMBERS AND NO
            # ZEROS (0) IN THE DENOMINATOR
            
            if ($se1 eq "X" or $se2 eq "X" or $sd1 eq "X" or $sd2 eq "X" or $total1==0 or $total2==0 or $sd1==0 or $sd2==0 or $df<=0){
                ($tdist,$pval)=("NA","NA");
            }
            else{
                $tdist=sqrt((($diff)/(sqrt((($sd1**2)/$total1)+(($sd2**2)/$total2))))**2);
                $pval=Statistics::Distributions::tprob($df,$tdist);
            }
            push (@info,$diff,$df,$tdist,$pval);
            $all{$gene}=\@info;
        }
    }
    close DATA;
}

#READ ALL HASH CONTENTS INTO 2D ARRAY FOR EASY SORTING AND PRINTING TO OUT FILE
my @unsorted;

foreach my $entry (keys %all) {
    my @info=@{$all{$entry}};
    my @temp;
    push (@temp,@info);
    push (@unsorted,\@temp);
}

#SORT GENES BY PVALUE OR FITNESS DEPENDING ON WHICH FLAG WAS USED IN COMMANDLINE
if (!$sortkey){
    $sortkey=14; #default: sort by difference of means
}
my @sorted = sort { $b->[$sortkey] <=> $a->[$sortkey] } @unsorted;

#ADD NEW FIELDS TO HEADER (COLUMN NAMES)
my $field="Mean".$labels[0].'-'.$labels[1];
push (@header,$field,"DOF","TDIST","PVALUE");

#PRINT TO OUT FILE
open OUT, '>',"$out";
print OUT (join(',',@header),"\n");
foreach (@sorted){
    my @woo=@{$_};
    print OUT join(',',@woo),"\n";
    }
close OUT;


