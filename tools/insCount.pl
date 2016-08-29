#!/usr/bin/perl -w

#Margaret Antonio updated 15.10.18


#perl ../singleFit.pl <all .csv input files L1_2394....>

use strict;
use Getopt::Long;
use warnings;
use Bio::SeqIO;
use Scalar::Util;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
sub print_usage() {
    print "\nDescription:\n";
    print "Takes multiple results files and returns average fitness or insertion count per insertion. Can ouput files as text, csv, wig.\n";
    print "\nCommand line: singleVal.pl <OPTIONS> <REQ OUTPUT TYPE(S)> <INPUT FILE(S)>\n\n";
    print "\nRequired:\n";
    print "In the command line (without a flag), input the name(s) of the file(s) containing fitness values for individual insertion mutants.\n";
    print "--h for help\n\n";
    
}


#ASSIGN INPUTS TO VARIABLES
our ($txt,$txtg,$cutoff,$h,$val, $wig,$ref_genome,$infile, $csv, $step,$indir, $size);
GetOptions(
'wig:s' => \$wig,
'ref:s' => \$ref_genome,
'cutoff:i'=>\$cutoff,
'i:s' => \$infile,
'csv:s'  => \$csv,
'step:i' => \$step,
'size:i' => \$size,
'txtg:s' => \$txtg,
'txt:s' => \$txt,
'd:s' => \$indir,
'h'=>\$h,
'v:s' => \$val, # fit for avg fitness val at site, count for # insertions
);


sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}
# Just to test out the script opening

print "\n";
if ($h){
    print_usage();
    exit;
}

if (!$cutoff){$cutoff=15};
if (!$val){$val="count"};


#CREATE AN ARRAY OF DATA FROM INPUT CSV FILE(S)
print "\nStart input array ",get_time(),"\n";

my $rowCount=-1;
my $last;
my @unsorted;

######################### INPUT FILES
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

########################## READ FILES

my $num=scalar @files;
print "\nNumber of files in csv: ", $num,"\n";

my %select; 
my %input;

for (my $i=0; $i<$num; $i++){   
    my $file=$files[$i];
    print "File #",$i+1,"\t",$file,"\n";
    
    open(DATA,'<', $file) or die "Could not open '$file'";
    
    my $dummy=<DATA>;
    while(my $line=<DATA>){
    chomp($line);
    my @fields=split(",",$line);
    my $pos=int($fields[0]);
    my $c2=int($fields[3]);
    if (exists $input{$pos}){
    	$input{$pos}+=$c2;
    	}
    else{
    	$input{$pos}=$c2;
    	}
    }
    close DATA; 
}

###################### OUTPUT A WIG FILE WTH POS AND COUNT

open (OUT,'>',"output_test.wig");

print OUT "variableStep\tchrom=chrN\n";
for my $pos(sort {$a<=>$b} keys %input){
	print OUT $pos, "\t",$input{$pos},"\n";
}
close OUT;

die;


