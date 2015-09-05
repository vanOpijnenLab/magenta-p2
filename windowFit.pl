#!/usr/bin/perl -w

#Margaret Antonio
#VERSION 5: This version of the windowFit script can output to a .csv file with the step option but lacks wig file creation for BigWig and edited out unnecessary Perl modules

#Use this commandline prompt for testing:
#../script/windowFit5.pl  --in results/sortedCompiledPennG.csv --out windowFit/5B-time100000.csv --step 10 --size 500

#Compiled the .csv files using mergeCSV.pl and then sorted them by GNU sort in terminal from /results/
# $ sort -t, -k1 -n compiledPennG.csv > sortedCompiledPennG.csv

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use DBI;
use DBD::CSV;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
sub print_usage() {
    
    print "\nRequired:\n";
    print "--size \t The size of the sliding window \n";
    print "--step \t The window spacing \n";
    print "--in \t The name of the file containing fitness values for individual insertion mutants.\n";
    print "--out \t Name of a file to enter the .csv output for sliding windows.\n\n";
}

#ASSIGN INPUTS TO VARIABLES

our ($infile, $outfile, $step, $size);

GetOptions(
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
print get_time(),"\n";
print "\n";
print "This is the input file: ", $infile,"\n";
print "THis is the output file: ", $outfile,"\n";
print "This is the window size chosen: $size\n";
print "This is the step value: $step\n";

#CHECKING PARAMETERS: Check to make sure required option inputs are there and if not then assign default
if (!$outfile) {
    #if there is no ref_genome, t1, t2, or out file in designated in the commandline, then: print all the user options and end the script
    print "\nThere must be an input out output file assigned\n";
    print_usage();
    die; }
if (!$size) { $size=500 };   #set the default sliding window size to 500
if (!$step) { $step = 10 };   #set the default step size to 10

###############################################################################################

print "\n---------This is the sliding window fitness calculation part--------\n\n";

my $csv = Text::CSV->new({ binary => 1, auto_diag => 1, eol => "\n"}) or die "Cannot use CSV: " . Text::CSV->error_diag();  # open in append mode
open my $TheOutFile, ">>", "$outfile" or die "Failed to open file: $!";
$csv->print($TheOutFile, [ "start", "end","fitness","mutant_count" ]); #header

sub OneWindow{
    my $Wstart=shift @_;
    my $Wend=shift@_;
    my $WwindowAvg=0;
    my $Wcount=0;
    my $Wsum=0;
    
    my $csv=Text::CSV->new({sep_char => ','});
    my $file=$infile or die "Need to get CSV file on the command line\n";
    open(my $data, '<', $file) or die "Could not open";
    $csv->getline ($data);   #Skips the header line in the .csv file
    my @data;
    #Going through each data-containing line in the .csv file
    while (my $line = <$data>){            #While there is data in the line
        chomp $line;                       #Chomp the line and parse it
        my @fields = split(/,/, $line);
        push @data, \@fields;
        
        # Check if positions fall in window parameters:
        if ($fields[0] eq "position"){next;}
        else{
            if (($fields[0] >= $Wstart) and ($fields[0]<=$Wend) and ($fields[12]) and ($fields[12]!=0)){
                $Wsum+=$fields[12];
                $Wcount+=1;
            }
            
            if ($fields[0]>$Wend){         #if finished with that window, then:
                if ($Wcount!=0){
                    $WwindowAvg=$Wsum/$Wcount;
                }
                else{
                    $WwindowAvg=0;
                }
                return ($WwindowAvg,$Wcount);
            }
        }
    }
    close $file;
}
print "Starting time: ",get_time(),"\n";
my $start=0;
my $end=$size;
my $windowNum=0;
print "This is how many windows have been completed: ";

while ($end<100000){
    my($windowAvg,$count)=OneWindow($start,$end);
    $csv->print($TheOutFile, [ "$start","$end","$windowAvg","$count"]);
    $start=$start+$step;
    $end=$end+$step;
    $windowNum+=1;
    print "$windowNum  ";
    
}
close $TheOutFile;
print "End time:",get_time();

#NOTES ###############################################################################################

#IMPROVEMENTS:
#Hashify columns for step option https://www.biostars.org/p/61526/
#Going to have to make new .wig file for BigWig for IGV
#Compile option for multiple file inputs
# Eventually should make sorting part of the script http://stackoverflow.com/questions/4221424/perl-sort-csv-on-a-certain-column


###############################################################################################

