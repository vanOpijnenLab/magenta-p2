#!/usr/bin/perl -w

#Margaret Antonio 16.01.13

#DESCRIPTION: Takes two aggregate.pl outputs and compares them using mean difference, pval for each
#gene. Can compare, for example, control vs antibiotic. DIFFERENT GENOMES (ie. diff. strains).
#Requires CONVERSION FILE

#USAGE: perl compGenes.pl <aggregateFile1.csv aggregateFile2.csv> OR -indir <indir/> -c <convFile>

use Data::Dumper;
use strict;
use Getopt::Long;
use warnings;
use File::Path;
use File::Basename;
use Statistics::Distributions;

#ASSIGN INPUTS TO VARIABLES
our ($indir,$h,$out,$sortkey,$round,$l,$cfile);
GetOptions(
'indir:s' => \$indir,
'o:s' =>\$out,
's:i' => \$sortkey,
'r:i'=> \$round,
'l:s'=> \$l,
'c:s'=>\$cfile,
);

if (!$cfile){ print "Must input conversion file under the flag -c\n" and die;}

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


#SET LABELS

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

#OPEN AGGREGATE GENE FILES INPUTTED INTO ARGV or INDIR and put them in hashes


####### FILE 1
my @header;
my %one; #Hash to store KEY: geneID VALUE: everything in the line
open (F1,'<',$files[0]);

#STORE COLUMN NAMES FOR HEADER AND APPEND LABELS
my $head=<F1>; #the header in the file
my @cols=split(',',$head);
@cols=@cols[0,1,2,3,4,5,6]; #get rid of blank columns
for (my $j=0;$j<scalar @cols;$j++){
    $cols[$j]=$cols[$j].'-'.$labels[0];   #mark each column name with file it comes from
}
push (@header,@cols);

while (my $line=<F1>){
    chomp $line;
    my @info=split(",",$line);
    @info=@info[0,1,2,3,4,5,6];
    if (!$info[5]){
        $info[5]="NA";
    }
    if (!$info[6]){
        $info[6]=0;
    }
    #make hash entry with KEY: geneID and VALUE:array of geneID,meanFit,...,totalInsertion
    $one{$info[0]}=\@info;
}
close F1;

####### FILE 2

my %two; #Hash to store KEY: geneID VALUE: everything in the line
open (F2,'<',$files[1]);

#STORE COLUMN NAMES FOR HEADER AND APPEND LABELS
$head=<F2>; #the header in the file
@cols=split(',',$head);
@cols=@cols[0,1,2,3,4,5,6]; #get rid of blank columns
for (my $j=0;$j<scalar @cols;$j++){
    $cols[$j]=$cols[$j].'-'.$labels[1];   #mark each column name with file it comes from
}
push (@header,@cols);

while (my $line=<F2>){
    chomp $line;
    my @info=split(",",$line);
    @info=@info[0,1,2,3,4,5,6];
    #For field 5 and field 6, can't have empty field so fill in for geneName and totalInsertion
    if (!$info[5]){
        $info[5]="NA";
    }
    if (!$info[6]){
        $info[6]=0;
    }
    #make hash entry with KEY: geneID and VALUE:array of geneID,meanFit,...,totalInsertion
    $two{$info[0]}=\@info;
}
close F2;


#READ CONVERSION FILE (.CSV) INTO ARRAY. Meanwhile if homologs exist then take info from hashes (above)
#conversion file must have strain 1 for file 1 in column 1 (index 0) and strain 2 for file 2 in column 2 (index 1)
#conversion file must be tab delimited with no NA fields

my @all; #store all homologs in this hash
open (CONV,'<',$cfile);
while (my $line=<CONV>){
    chomp $line;
    my @genes=split("\t",$line);   #Array @genes will contain two genes (SP_0000,SPT_0000)
    if (scalar @genes==2 and $genes[0] ne "" and $genes[1] ne ""){
        my @info;
        my @oneArray=@{$one{$genes[0]}};
        my @twoArray=@{$two{$genes[1]}};
        push (@info,@oneArray,@twoArray);
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
        #Do calculations
        if ($se1 eq "X" or $se2 eq "X" or $sd1 eq "X" or $sd2 eq "X" or $total1==0 or $total2==0 or $sd1==0 or $sd2==0){
            ($tdist,$pval)=("NA","NA");
        }
        else{
            $tdist=sqrt((($diff)/(sqrt((($sd1**2)/$total1)+(($sd2**2)/$total2))))**2);
            $pval=Statistics::Distributions::tprob($df,$tdist);
        }
        push (@info,$diff,$df,$tdist,$pval);
        #foreach (@info){
        #    print $_,"\t";
        #}
        #print "\n";
        push (@all,\@info);
    }
}
close CONV;


#PRINT TO OUT FILE

if (!$sortkey){
    $sortkey=14; #for mean difference
}
my @sorted = sort { $b->[$sortkey] <=> $a->[$sortkey] } @all;

#Finish up the header
my $field="MeanDiff(".$labels[0].'-'.$labels[1].")";
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


