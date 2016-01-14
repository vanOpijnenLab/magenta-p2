#!/usr/bin/perl -w

#Margaret Antonio 16.01.13

#DESCRIPTION: Takes two aggregate.pl outputs and compares them using mean difference, pval for each
#gene. Can compare, for example, control vs antibiotic.

#Improvements to be made: default output should b

#USAGE: perl compGenes.pl <aggregateFile1.csv aggregateFile2.csv> OR -indir <indir/>

use Data::Dumper;
use strict;
use Getopt::Long;
use warnings;
use File::Path;
use File::Basename;
use Statistics::Distributions;

#ASSIGN INPUTS TO VARIABLES
our ($indir,$h,$out,$sortkey,$round,$l);
GetOptions(
'indir:s' => \$indir,
'o:s' =>\$out,
's:i' => \$sortkey,
'r:i'=> \$round,
'l:s'=> \$l,
);

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}
i
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
push (@header,"id");
my @labels;
if ($l){
    @labels=split(',',$l);
}
else{
    #get rid of .csv part in file names and use as labels
    foreach (@files){
       my @temp=split('\\.',$_);
       my $colName=$temp[0];
       push (@labels,$colName);
   }
}

if (!$out) {$out="comp-".$labels[0].$labels[1].".csv"}
if (!$round){$round='%.4f'}



#push (@header,@sample);

my %all;
for (my $i=0; $i<$num; $i++){   #Read files from ARGV
    print "File #",$i+1,"\t";
    
    my $file=$files[$i];
    print $file,"\n";
    
    open(DATA, '<', $file) or die "Could not open '$file' Make sure input .csv files are entered in the command line\n";
    
    #extract the header
    my $head=<DATA>; #the header in the file
    my @cols=split(',',$head);
     @cols=@cols[1,2,3,4,5,6]; #get rid of gene name
    for (my $j=0;$j<scalar @cols;$j++){
        $cols[$j]=$cols[$j].'-'.$labels[$i];
    }
    push (@header,@cols);
    #my %hash;
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
        #my $mean = sprintf("$round",$line[1]);
        shift @line; #pop off the "gene name"
        
        ############################## hash assignment
        if(!exists $all{$gene}){
            my @info;
            push (@info,@line);
            $all{$gene}=\@info;
        }
        else{
            my @info=@{$all{$gene}};
            push (@info,@line);
            my $diff=sprintf("$round",($info[0]-$info[6]));
            my $total1=$info[5];
            my $total2=$info[11];
            my $sd1=$info[2];
            my $se1=$info[3];
            my $sd2=$info[8];
            my $se2=$info[9];
            my $df=$total1+$total2-2;
            my $tdist;
            my $pval;
            #print $gene, "\t",$total1,"\t",$total2,"\n";
            if ($se1 eq "X" or $se2 eq "X" or $sd1 eq "X" or $sd2 eq "X" or $total1==0 or $total2==0 or $sd1==0 or $sd2==0){
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


#PRINT TO OUT FILE

#Set up the header (column names)

my @unsorted;

foreach my $entry (keys %all) {
    my @info=@{$all{$entry}};
    #print $entry,",";      #if unhashed, must go after OPEN
    #print (join(",",@info),"\n");
    my @temp;
    push (@temp,$entry);
    push (@temp,@info);
    push (@unsorted,\@temp);
}

if (!$sortkey){
    $sortkey=13; #for mean difference
}
my @sorted = sort { $b->[$sortkey] <=> $a->[$sortkey] } @unsorted;

#Finish up the header
my $field="Mean".$labels[0].'-'.$labels[1];
push (@header,$field,"DOF","TDIST","PVALUE");

open OUT, '>',"$out";
print OUT (join(',',@header),"\n");
foreach (@sorted){
    my @woo=@{$_};
    #print scalar @woo,"\t";
    print OUT join(',',@woo),"\n";
    }

    #print OUT (join(',',@foo),"\n");

close OUT;


