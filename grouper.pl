#!/usr/bin/perl -w

#Margaret Antonio 15.12.26

#DESCRIPTION: After sliding windows have been filtered and ordered, can be grouped
#must decide to group by fitness or insertion representation (one or the other)

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use List::Util qw(sum);
use List::BinarySearch qw( :all );
use List::BinarySearch::XS;
use List::MoreUtils qw(uniq);
use feature qw/say/;
use autodie;

sub print_usage() {
    print "\n";
    print "grouper.pl: Takes slidingWindow.csv file and groups consecutive like windows\n\n";
    print "Example: Group consecutive windows by fitness. If fitness is the same then expand.\n";
    print "For different strains/genomes, use compStrains.pl\n";
    print "\nUSAGE: perl grouper.pl <options> <slidingWindow.csv>\n\n";
    print "OPTIONS:\n\n";
    print " -h\tPrints usage and quits\n";
    print " -o\tOutput file for grouped windows. Default: groupedWindows.csv\n";
    print " -s\tSort output by this index of the file (indices begin at 0).\n";
    print " -sig\t Group consecutive windows of same pvalue\n";
    print " -fit\t Group consecutive windows of same fitness value\n";
    print " -r\tRound final output numbers to this number of decimals\n";
    print "\n\n";
}
if ($h){
    print_usage();
    exit;
}


#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
# 0:start,1:end,2:fitness,3:mutant_count,4:insertions,5:TA_sites,6:ratio,7:p-value

#ASSIGN INPUTS TO VARIABLES
our ($h,$out,$sortby,$round,$sig,$fit);
GetOptions(
'o:s' =>\$out,
'h' => \$h,
's:i' => \$sortby,
'r:i'=> \$round,
'sig'=>\$sig,
'fit'=>\$fit,

);

#Assign defaults
if (!$out) {$out="groupedWindows.csv"}
if (!$round){$round='%.3f'}
#if sortby was not specified then use -fit flag or default by significance
if (!$sortby){
    if ($fit){
        $sortby=7;
    }
    else{
        $sortby=10;
    }
}

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}
sub mean {
    return sum(@_)/@_;
}
sub cleaner{
    my $line=$_[0];
    chomp($line);
    $line =~ s/\x0d{0,1}\x0a{0,1}\Z//s;
    return $line;
}

#READ INPUT FILE INTO A 2D ARRAY @windows

my $in=$ARGV[0]; #get input file from first entry in ARGV
my @windows;
open(DATA, '<', $in) or die "Could not open '$in' \n";
#Store first line for use as column names of output file
my $line=<DATA>;
$line=cleaner($line);
my @header=split(',',$line);
while (my $entry = <DATA>) {
	$entry=cleaner($entry);
	my @fields=split(',',$entry);
	push (@windows,\@fields);     
}
close DATA;

my $count=0; #number of windows added to cumm that will be used for avg calc
my $sumFit=0;
my $sumRatio=0;

	
#SORT ALL REGIONS BY ABS. DIFF. BTWN REGION'S MEAN DIFF AND AVG MEAN FITNESS (index 8).
#If sort by fitness, most interesting regions have high mean differece so sort largest to smallest
if ($fit){
    @windows= sort {$b->[$sortby]<=>$a->[$sortby] || $a->[0]<=>$b->[0] } @expWindows;
}
#If sort by significance, most interesting regions are low pvals, so sort smallest to largest
else{
    @windows= sort {$a->[$sortby]<=>$b->[$sortby] || $a->[0]<=>$b->[0] } @expWindows;
}

open my $gOut, '>', "$out" or die $!;

#PRINT COLUMN NAMES (HEADER) TO THE OUTPUT FILE
my $string = join(",", @header);
print $gOut "$string\n"; 

#Start cummulative array to begin for loop
my $lastLine=scalar @windows;
my @cumu=@{$windows[0]};
my ($cstart,$cend,$cfit,$cmcount,$cins,$cta,$cratio,$cpval,$cdev)=@cumu;

for (my $i=1;$i<$lastLine;$i++){
	my @field=@{$windows[$i]};
	my ($fstart,$fend,$ffit,$fmcount,$fins,$fta,$fratio,$fpval,$fdev)=@field;
    
    #Either this window needs to be added or need to start new cummulative
	
	#Add field window (@field) if overlaps with cumulative window (@cumu)
	if (($cend>=$fstart and $cstart<=$fstart) or ($cend>=$fend and $cstart<=$fend)){ 
		#Keep cstart as it is
		$cend=$fend; #change the end coordinate
		$sumFit+=($ffit*$fins);
		$cins+=$fins;
		$cmcount+=$fmcount;
	}
    #need to output this cumm region with average calcs
	else{
		$cratio=sprintf("$round",($cins/$cta));
		if ($cins !=0){
			$cfit=sprintf('%.2f',($sumFit/$cins)); #not accurate
			}
		else{
			$cfit=0;
			}
		@cumu=($cstart,$cend,$cfit,$cmcount,$cins,$cta,$cratio,$cpval,$cdev);
		print $gOut ($_,",") for @cumu;
		print $gOut ("\n");
		#Set up current entry as new cumulative
		@cumu=@field;
	    ($cstart,$cend,$cfit,$cmcount,$cins,$cta,$cratio,$cpval,$cdev)=@cumu;
	    $sumFit=0;
	    $sumRatio=0;
	}
}
close $gOut;