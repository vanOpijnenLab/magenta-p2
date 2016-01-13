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
our ($indir,$h,$out,$sortby,$round,$l);
GetOptions(
'indir:s' => \$indir,
'o:s' =>\$out,
's:i' => \$sortby,
'r:i'=> \$round,
'l:s'=> \$l,
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
#push (@header,@sample);


my $sortkey;
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
            push (@info,$diff);
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
    $sortkey=(scalar @temp)-1;
    push (@unsorted,\@temp);
}


my @sorted = sort { $b->[$sortkey] <=> $a->[$sortkey] } @unsorted;

#Finish up the header
push (@header,"abs(diffMean)");

open OUT, '>',"$out";
print OUT (join(',',@header),"\n");
foreach (@sorted){
    my @woo=@{$_};
    #print scalar @woo,"\t";
    print OUT join(',',@woo),"\n";
    }

    #print OUT (join(',',@foo),"\n");

close OUT;


