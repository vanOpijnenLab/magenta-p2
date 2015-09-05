#!/usr/bin/perl -w

#Margaret Antonio
#VERSION 11: 15.04.06 9:41 pm | This version reads the input file into an array for access when gathering values that fall in the window. The start, end, fitness, and count values for each window are put in a row array and then appendend to a 2D array of all the windows made from the file. That array is then written all at once into the ouput csv file
#Use this commandline prompt for testing:   ../script/windowFit11.pl  --in results/sortedCompiledPennG.csv --out windowFit/12CwithoutZero.csv --step 10 --size 500
#Compiled the .csv files using mergeCSV.pl and then sorted them by GNU sort in terminal from /results/
# $ sort -t, -k1 -n compiledPennG.csv > sortedCompiledPennG.csv

#../script/windowFit11.pl --cutoff 15 --in results/sortedCompiledPennG.csv --out windowFit/11Xwhole.csv --step 10 --size 500

#../script/windowFit11.pl  --wig gview/12G.wig --ref=NC_003028b2.gbk --cutoff 15 --in results/sortedCompiledPennG.csv --out windowFit/12G-withWig.csv --step 10 --size 500

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use Text::CSV_XS;
use Bio::SeqIO;

#AVAILABLE OPTIONS. WILL PRINT UPON ERRO
sub print_usage() {
    print "\nRequired:\n";
    print "--size \t The size of the sliding window \n";
    print "--step \t The window spacing \n";
    print "--in \t The name of the file containing fitness values for individual insertion mutants.\n";
    print "--out \t Name of a file to enter the .csv output for sliding windows.\n";
    print "--cutoff \tCutoff: Don't include fitness scores with average counts (c1+c2)/2 < x (default: 0)\n\n";
    print "\nOptional:\n";
    print "--wig\tCreate a wiggle file for viewing in a genome browser. Provide a filename. Also provide genome under --ref\n";
    print "--ref\tThe name of the reference genome file, in GenBank format. Needed for wig file creation\n";
    
}

#ASSIGN INPUTS TO VARIABLES
our ($cutoff,$wig,$ref_genome,$infile, $outfile, $step, $size);
GetOptions(
'wig:s' => \$wig,
'ref:s' => \$ref_genome,
'cutoff:i'=>\$cutoff,
'in:s' => \$infile,
'out:s'  => \$outfile,
'step:i' => \$step,
'size:i'   => \$size,

);

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}
# Just to test out the script opening
print print_usage(),"\n";
print "\n";
print "Input file: ", $infile,"\n";
print "Output file: ", $outfile,"\n";
print "Window size: $size\n";
print "Step value: $step\n";

#CHECKING PARAMETERS: Check to make sure required option inputs are there and if not then assign default
if (!$outfile) {
    print "\nThere must be an input out output file assigned\n";
    print_usage();
    die; }
if (!$size) { $size=500 };   #set the default sliding window size to 500
if (!$step) { $step = 10 };   #set the default step size to 10

###############################################################################################

print "\n---------This is the sliding window fitness calculation part--------\n\n";

print "Start input array ",get_time(),"\n";
#Create an array of data from input csv file
my @data;
my $csvIN=Text::CSV->new;   #({sep_char => ','})
open my $fh, '<', $infile or die "Could not open";
$csvIN->getline($fh);  #get rid of header
while (my $line = $csvIN->getline($fh)){      #While there is data in the line
    my $w = $line->[12];
    #if ($w and $w eq 'nW') {next;}
    if (!$w){next;} # For blanks
    else{
    my $c1 = $line->[2];
    my $c2 = $line->[3];
    my $avg = ($c1+$c2)/2; # Later: change which function to use? C1? AVG(C1+C2)?
    if ($avg < $cutoff) { next; } # Skip cutoff genes.
    else {
        push(@data,$line);
        }
        
    }
}
close $fh;
print "Finished input array ",get_time(),"\n";

sub OneWindow{
    
    my $Wstart=shift @_;
    my $Wend=shift@_;
    my $WwindowAvg=0;
    my $Wcount=0;
    my $Wsum=0;
    
    foreach my $row(@data){
        my @fields=@$row;   #for retrieving info from the input 2D array @data with $row (which references @fields--12 columns or input csv)
        my @window=();         #an empty array to contain start, end, fitness, and count for each window calculation
        # Check if positions fall in window parameters:
        
        if (($fields[12] eq 0) or (!$fields[12])){next;}
        else{
            if (($fields[0] >= $Wstart) and ($fields[0]<=$Wend) and ($fields[12])){
                $Wsum+=$fields[12];
                $Wcount+=1;
            }
            if ($fields[0]>$Wend){         #if finished with that window, then:
                if ($Wcount!=0){
                    $WwindowAvg=$Wsum/$Wcount;
                    my @Wwindow=($Wstart,$Wend,$WwindowAvg,$Wcount);
                    print @Wwindow,"\n";
                    return (\@Wwindow);
                }
                else{return -1}
            
                
            }
        }
        }
    }

print "Start calculation time: ",get_time(),"\n";
my $start=1;
my $end=$size;
my $windowNum=0;
my @allWindows; #will be a 2D array containing all window info to be written into output file

while ($end<300000){  #100,000bp-->9,950 windows--> only 8500 windows in csv because 0
    my($window)=OneWindow($start,$end);
    if ($window!=-1){
        push (@allWindows,$window);
    }
    $start=$start+$step;
    $end=$end+$step;
    print "$end  ";
}
print "End calculation time: ",get_time(),"\n\n";

print "Start wig file creation: ",get_time(),"\n";
#make wig file
if ($wig){
my $in = Bio::SeqIO->new(-file=>$ref_genome);
my $refseq = $in->next_seq;
my $refname = $refseq->id;
    open WIG, ">$wig";
    print WIG "track type=wiggle_0 name=$wig\n";
    print WIG "variableStep chrom=$refname\n";
    foreach my $wigLine(@allWindows){
        my @wigFields=@$wigLine;
        print  WIG $wigFields[0]," ",$wigFields[2],"\n";
    }
    close WIG;
}
print "End wig file creation: ",get_time(),"\n\n";
#Convert wig to big wig with terminal command: wigToBigWig gview/12G.wig tigr4_chrom.txt >BigWig/12G.bw
#Make output CSV file
print "Start csv ouput file creation: ",get_time(),"\n";
my $csv = Text::CSV->new({ binary => 1, auto_diag => 1, eol => "\n"}) or die "Cannot use CSV: " . Text::CSV->error_diag();  # open in append mode
open my $file, ">", "$outfile" or die "Failed to open file";
$csv->print($file, [ "start", "end","fitness","mutant_count" ]); #header
for my $winLine(@allWindows){
    $csv->print($file,$winLine);
}
close $file or die "$outfile: $!";
print "End csv ouput file creation: ",get_time(),"\n\n";
