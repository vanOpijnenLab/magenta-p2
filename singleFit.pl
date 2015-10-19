#!/usr/bin/perl -w

#Margaret Antonio updated 15.09.15


#perl ../filterMerge.pl <all .csv input files L1_2394....>

use strict;
use Getopt::Long;
use warnings;
use Bio::SeqIO;
use Scalar::Util;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
sub print_usage() {
    print "\nDescription:\n";
    print "Integrates multiple files of transposon insertion data and outputs aggregate fitness within a sliding window (specified by size and step). Can ouput files as text, csv, wig.\n";
    print "\nCommand line: windowFit.pl <OPTIONS> <REQ OUTPUT TYPE(S)> <INPUT FILE(S)>\n\n";
    print "\nRequired:\n";
    print "In the command line (without a flag), input the name(s) of the file(s) containing fitness values for individual insertion mutants.\n";
    print "--h for help\n\n";
    
}


#ASSIGN INPUTS TO VARIABLES
our ($txt,$txtg,$cutoff,$h, $wig,$ref_genome,$infile, $csv, $step, $size);
GetOptions(
'wig:s' => \$wig,
'ref:s' => \$ref_genome,
'cutoff:i'=>\$cutoff,
'in:s' => \$infile,
'csv:s'  => \$csv,
'step:i' => \$step,
'size:i' => \$size,
'txtg:s' => \$txtg,
'txt:s' => \$txt,
'h'=>\$h,
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


#CREATE AN ARRAY OF DATA FROM INPUT CSV FILE(S)
print "\nStart input array ",get_time(),"\n";

my $rowCount=-1;
my $last;
my @unsorted;
my $num=$#ARGV+1;
print "\nNumber of files in csv: ", $num,"\n";

my %select; #this is a hash to store insertion position and fitness value at that position for each file (lane)

for (my $i=0; $i<$num; $i++){   #Read files from ARGV---loops once for each file in ARGV
    my $file=$ARGV[$i];
    open(DATA,'<', $file) or die "Could not open '$file' Make sure input .csv files are entered in the command line\n";
    print "opened!";
    print "\t",$file,"\n";
    my $dummy=<DATA>;
    print "$dummy\t\t";
    while (my $line = <DATA>) {
        		my @fields=split(",",$line);
        my $w = $fields[12];
        if (!$w or $w eq"\n"){
            next;} # For blanks
        else{
            my $c1 = $fields[2];
            my $c2 = $fields[3];
            my $avg = ($c1+$c2)/2;
            if ($avg < $cutoff) {
                next;
                
            } # Skip cutoff genes.
            else {  #This is a good value and should be added to the hash
            	my $pos=$fields[0];
            	$w=sprintf("%.2f",$w);
            	if(!exists $select{$pos}){
            	#print $pos;
            		my @fits;
            		push(@fits,"1");
            		for(my $j=0;$j<$i-1;$j++){
            			push (@fits,"0");
            		}
            		push (@fits,$w);
            		$select{$pos}=\@fits;
            	}
            	else{
            		my @fits=@{$select{$pos}};
					for (my $j=scalar @fits+1;$j<$i;$j++){
						push (@fits,"0");
					}
            		push (@fits,$w);
                    $select{$pos}=\@fits;
            	}
            }
        }
    }
    close DATA;
}

foreach my $val(values %select){
	if (scalar @$val <$num+1){
		for (my $j=scalar @$val;$j<$num+1;$j++){
			push (@$val,"0");
		}
	}
}
		

#IF MAKING A REGULAR TEXT FILE fields: [chrom,start,end,fitness,count]
open TXT,'>', "heatPrep151004.txt";
    #print TXT "Start text file creation time: ",get_time(),"\n";
print TXT "seqnames\tstart\tend\tcontrol\tL1\tL3\tL4\tL5\tL6\n";
foreach my $entry (sort keys %select) {
    print TXT "NC003028\t",$entry,"\t",$entry+1,"\t";
    my @entryFits=@{$select{$entry}};
    foreach (@entryFits){
    	print TXT $_, "\t";
    }
    print TXT "\n";
}
close TXT;


open ALL,'>', "singleFit.txt";
foreach my $entry (sort keys %select) {
    print ALL "NC003028\t",$entry,"\t",$entry+1,"\t";
    my @entryFits=@{$select{$entry}};
    my $sum=0; my $count=0;    
    foreach (@entryFits){
    	if ($_!=0){
    		$count++;
    		$sum+=$_;
    	}
    }
    my $avgFit=sprintf("%.2f",$sum/$count);
    print ALL $avgFit,"\n";
}
close ALL;



