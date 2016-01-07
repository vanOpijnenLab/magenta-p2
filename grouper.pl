#!/usr/bin/perl -w

#Margaret Antonio 15.12.26

#DESCRIPTION: After sliding windows have been filtered and ordered, can be grouped

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;


#USAGE from /8-essFilters
#perl ../Blueberries/grouper.pl --in apples.csv

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR

#ASSIGN INPUTS TO VARIABLES
our ($in,$h,$size, $bracket,$step,$defEss,$out);
GetOptions(
'i:s' => \$in,
'o:s' =>\$out,
'h' => \$h,
'b'=>\$bracket,
'size'=>\$size,
'step'=>\$step,
'ess'=>\$defEss,
);

if (!$size) { $size=500 };   #set the default sliding window size to 500
if (!$step) {$step=10};
if (!$out) {$out="orderedGroup.csv"};


sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}

#GET DATA OUT OF FILE AND INTO 2D ARRAY


my @newWindows;

sub cleaner{
	my $line=$_[0];
	chomp($line);
	$line =~ s/\x0d{0,1}\x0a{0,1}\Z//s; 
	return $line;
	}
	
open(DATA, '<', $in) or die "Could not open '$in' \n";
	
my $line=<DATA>;
#print "This is the line: ",$line,"stoooooop";
$line=cleaner($line); #gets rid of carriage returns (^M)
my @header=split(',',$line);
push (@header,'csTest');
my $tick=0;
while (my $entry = <DATA>) {
	#$tick+=1;
	#print $tick, "\t";
	$entry=cleaner($entry);
	#chomp($entry);
	#$entry =~ s/\x0d{0,1}\x0a{0,1}\Z//s; 
	my @fields=split(',',$entry);
	push (@newWindows,\@fields); #should be ref array    
	#freezes after line 214820
}
close DATA;

#my @finalOutput;
#my $lastLine=scalar @outArray;

my $div=0; #number of windows added to cumm that will be used for avg calc
my $sumFit=0;
my $sumRatio=0;

print "Start grouped txt file creation time: ",get_time(),"\n";

open my $gOut, '>', "$out" or die $!;

my $lastLine=scalar @newWindows;
my @cumm=@{$newWindows[0]};

for (my $i=1;$i<$lastLine;$i++){
	my @field=@{$newWindows[$i]};
	$sumFit=$cumm[3];
	$div=1;
	$sumRatio=$cumm[7];
	
    #either this window needs to be added or need to start new cumm
	if (($cumm[2]>=$field[1] and $cumm[2]<=$field[2]) or ($field[2]>=$cumm[1] and $field[2]<=$cumm[2])){ 
	#Add window if overlapping windows
		$cumm[2]=$field[2]; #change the end coordinate
		$sumFit+=$field[3];
		$div++;
		$sumRatio+=$field[7];
	}
	else{ #need to output this cumm region with average calcs
		$cumm[3]=$sumFit/$div;
		$cumm[7]=$sumRatio/$div;
		shift @cumm;
		print $gOut ($_,",") for @cumm;
		print $gOut ("\n");
		@cumm=@field;
	}
	
}
close $gOut;