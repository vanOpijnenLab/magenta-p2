#!/usr/bin/perl -w

#Margaret Antonio 15.12.26

#DESCRIPTION: After sliding windows have been filtered and ordered, can be grouped

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use List::Util qw(sum);

use Bio::SeqIO;
use Data::Random qw(:all);
use List::Util qw(sum);
use List::BinarySearch qw( :all );
use List::BinarySearch::XS;
use List::MoreUtils qw(uniq);
use File::Path;
use File::Basename;
use feature qw/say/;
use autodie;


#USAGE from /8-essFilters
#perl ../Blueberries/grouper.pl -i apples.csv

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR


# 0:start,1:end,2:fitness,3:mutant_count,4:insertions,5:TA_sites,6:ratio,7:p-value
#print "Sort by options:\n0:Start coord.\n1:End coord\n2:Fitness\n3: mutant count\n4: insertions\n5:TA_sites\n6:ratio\n7: p-value\n8: deviation from mean fitness\n";

#ASSIGN INPUTS TO VARIABLES
our ($in,$h,$size, $bracket,$step,$out,$sortby,$round,$sig,$fit);
GetOptions(
'i:s' => \$in,
'o:s' =>\$out,
'h' => \$h,
'b'=>\$bracket,
'size'=>\$size,
'step'=>\$step,
's:i' => \$sortby,
'r:i'=> \$round,
'sig'=>\$sig,
'fit'=>\$fit,

);

if (!$size) {$size=500}   #set the default sliding window size to 500
if (!$step) {$step=10}
if (!$out) {$out="orderedGroup.csv"}
if (!$round){$round='%.3f'}
#if sortby was not specified then a quicker sortby fitness using -fit flag or default by significance
if (!$sortby){
    if ($fit){
        $sortby=8;
    }
    else{
        $sortby=7;
    }
}




sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}

#GET DATA OUT OF FILE AND INTO 2D ARRAY

sub mean {
	return sum(@_)/@_;
}
my @windows;

sub cleaner{
	my $line=$_[0];
	chomp($line);
	$line =~ s/\x0d{0,1}\x0a{0,1}\Z//s; 
	return $line;
	}
	
open(DATA, '<', $in) or die "Could not open '$in' \n";
	
my $line=<DATA>;
$line=cleaner($line); #gets rid of carriage returns (^M)
my @header=split(',',$line);
while (my $entry = <DATA>) {
	$entry=cleaner($entry);
	my @fields=split(',',$entry);
	push (@windows,\@fields);     
}
close DATA;

#my @finalOutput;
#my $lastLine=scalar @outArray;

my $count=0; #number of windows added to cumm that will be used for avg calc
my $sumFit=0;
my $sumRatio=0;

#print "Start grouped txt file creation time: ",get_time(),"\n";


#What's the mean fitness value for all of the windows?
my @allFits = map $_->[2], @windows;
my $meanFit=mean(@allFits);

#Add the absolute deviation from mean fitness to each window array (use this to sort)
my @expWindows=();
foreach (@windows){
	my @entry=@{$_};	
	my $absdev=sprintf('%.1f',abs($entry[2]-$meanFit));
	push @entry,$absdev;
	push @expWindows,\@entry;	
	} 
	
#Now sort @expWindows by the abs. dev. of fitness (index 8). Reuse old @windows variable

#For sorting by fitness, want top of file to show most interesting regions---high dev. from mean
if ($fit){
    @windows= sort {$b->[$sortby]<=>$a->[$sortby] ||
	$a->[0]<=>$b->[0] } @expWindows;
}
#For sorting by significance, most interesting region are low pvals. so sort with smallest values at top
else{
    @windows= sort {$a->[$sortby]<=>$b->[$sortby] ||
    $a->[0]<=>$b->[0] } @expWindows;
}

open my $gOut, '>', "$out" or die $!;

#column fields from sliding window input: 
# 0:start,1:end,2:fitness,3:mutant_count,4:insertions,5:TA_sites,6:ratio,7:p-value

#print the header
push(@header, "abs(diff_mean)");
my $string = join(",", @header);
print $gOut "$string\n"; 

#Seed entry to begin grouping
my $lastLine=scalar @windows;
my @cumu=@{$windows[0]};
my ($cstart,$cend,$cfit,$cmcount,$cins,$cta,$cratio,$cpval,$cdev)=@cumu;

for (my $i=1;$i<$lastLine;$i++){
	my @field=@{$windows[$i]};
	my ($fstart,$fend,$ffit,$fmcount,$fins,$fta,$fratio,$fpval,$fdev)=@field;
    #either this window needs to be added or need to start new cumm
	
	#Add field window (@field) if overlaps with cumulative window (@cumu)
	if (($cend>=$fstart and $cstart<=$fstart) or ($cend>=$fend and $cstart<=$fend)){ 
		#Keep cstart as it is
		$cend=$fend; #change the end coordinate
		$sumFit+=($ffit*$fins);
		$cins+=$fins;
		$cmcount+=$fmcount;
	}
	else{ #need to output this cumm region with average calcs
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