#!/usr/bin/perl -w

#Margaret Antonio 16.01.13

#DESCRIPTION: Takes two aggregate.pl outputs and compares them using mean difference, pval for each
#gene. Can compare, for example, 19F in glucose and TIGR4 in glucose.
#DIFFERENT GENOMES (ie. diff. strains).
#Requires CONVERSION FILE

#USAGE: perl compStrains.pl -c <conversion.csv> <options>
    #[<aggregateFile1.csv aggregateFile2.csv> OR -indir <indir/>]

use Data::Dumper;
use strict;
use Getopt::Long;
use warnings;
use File::Path;
use File::Basename;
use Statistics::Distributions;

#ASSIGN INPUTS TO VARIABLES USING FLAGS
our ($indir,$h,$out,$sortkey,$round,$h,$l,$cfile);
GetOptions(
'd:s' => \$indir,
'h' => \$h,
'o:s' =>\$out,
's:i' => \$sortkey,
'r:i'=> \$round,
'l:s'=> \$l,
'cfile:s'=> \$cfile,
);

sub print_usage() {
    print "\n";
    print "compStrains.pl: COMPARE GENES FROM A TWO TN-SEQ EXPERIMENT\n";
    print "\t(OF DIFFERENT ORGANISMS/STRAINS (GENOMES) USING AGGREGATE FILES\n\n";
    print "DESCRIPTION: Takes two aggregate.pl outputs and compares them by calculating\n";
    print "the difference in mean fitness, the pval for each gene.\n";
    print "Example: two strains tested under same condition.\n";
    print "For same strains (genomes), conversion file not needed. Use compGenes.pl\n";
    print "\nUSAGE: perl compStrains.pl -c <conversion.csv> <options> \n";
    print "\t<aggregateFile1.csv aggregateFile2.csv> or -d <directory/> \n\n";
    print "OPTIONS:\n\n";
    print " -h\tPrints usage and quits\n";
    print " -d\tDirectory containing input files. Make sure / is included after name\n";
    print " -o\tOutput file for comparison data. Default: compFile1File2.csv\n";
    print " -s\tSort output by this index of the file (indices begin at 0). Default: by mean\n";
    print " -r\tRound final output numbers to this number of decimals\n";
    print " -l\tLabels for for compared files. Used for column names and default output file name.\n";
    print "   \tTwo strings, comma separated (i.e. -l gluc,dapto). Order should match file order.\n";
    print " -c\tConversion file: two columns with homologs for genome1 and genome 2\n";
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

#OPEN INPUTTED AGGREGATE GENE FILES AND STORE THEIR CONTENTS INTO TWO HASHES
#FILE1 GOES INTO HASH %ONE AND FILE2 GOES INTO HASH %TWO.

#FILE1 OPENING ---> %one WHERE KEY:VALUE IS GENE_ID:(GENE_ID,INSERTIONS,MEAN,ETC.)
my @header;
my %one;

open (F1,'<',$files[0]);

#STORE COLUMN NAMES (FIRST LINE OF FILE1) FOR HEADER AND APPEND LABELS
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
    #Only keep the first 7 columns (Ones about blanks aren't needed for comparisons)
    @info=@info[0,1,2,3,4,5,6];
    #Sometimes genes that don't have a gene name can't be blank, so fill with NA
    if (!$info[5]){
        $info[5]="NA";
    }
    #If there are no insertions in the column "total", then make it =0 rather than blank
    if (!$info[6]){
        $info[6]=0;
    }
    $one{$info[0]}=\@info;
}
close F1;

#FILE2 OPENING ---> %two WHERE KEY:VALUE IS GENE_ID:(GENE_ID,INSERTIONS,MEAN,ETC.)

my %two;
open (F2,'<',$files[1]);

#STORE COLUMN NAMES (FIRST LINE OF FILE2) FOR HEADER AND APPEND LABELS
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
    if (!$info[5]){
        $info[5]="NA";
    }
    if (!$info[6]){
        $info[6]=0;
    }
    $two{$info[0]}=\@info;
}
close F2;


#READ CONVERSION FILE INTO ARRAY.
#Conversion file must have strain 1 for file 1 in column 1 (index 0) and
    #strain 2 for file 2 in column 2 (index 1)
    #conversion file must be tab delimited with no NA fields
#If homologs (exist then take info from hashes (%one and %two) by referring to gene_id in KEY

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
        #TDIST, PVAL calculations with fail if standard dev, error, or counts are not real numbers
        #or if 0 ends up in denominator
        if ($se1 eq "X" or $se2 eq "X" or $sd1 eq "X" or $sd2 eq "X" or $total1==0 or $total2==0 or $sd1==0 or $sd2==0){
            ($tdist,$pval)=("NA","NA");
        }
        else{
            $tdist=sqrt((($diff)/(sqrt((($sd1**2)/$total1)+(($sd2**2)/$total2))))**2);
            $pval=Statistics::Distributions::tprob($df,$tdist);
        }
        push (@info,$diff,$df,$tdist,$pval);
        push (@all,\@info);
    }
}
close CONV;

#SORT THE HOMOLOGS BY THE SORTKEY OR BY DEFAULT DIFFERENCE IN MEAN FITNESSES
if (!$sortkey){
    $sortkey=14; #for mean difference
}
my @sorted = sort { $b->[$sortkey] <=> $a->[$sortkey] } @all;

#FINISH THE HEADER BY ADDING COLUMN NAMES FOR MEAN-DIFF, DOF, TDIST, AND PVALUE
my $field="MeanDiff(".$labels[0].'-'.$labels[1].")";
push (@header,$field,"DOF","TDIST","PVALUE");

#PRINT MATCHED HOMOLOG INFORMATION INTO A SINGLE OUTPUT FILE
open OUT, '>',"$out";
print OUT (join(',',@header),"\n");
foreach (@sorted){
    my @woo=@{$_};
    print OUT join(',',@woo),"\n";
    }

close OUT;


