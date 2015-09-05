
#!/usr/bin/perl -w

#Margaret Antonio
#VERSION 13: This version reads the input file into an array for access when gathering values that fall in the window. Thought it would be faster than chomping each line but may not be.

#Use this commandline prompt for testing:
#../script/13window.pl  --in results/L3_2394eVI_PennG.csv --out windowFit/Awesome10.csv --step 10 --size 500

#Compiled the .csv files using mergeCSV.pl and then sorted them by GNU sort in terminal from /results/
# $ sort -t, -k1 -n compiledPennG.csv > sortedCompiledPennG.csv

use strict;
use Getopt::Long;
use warnings;
#use Text::CSV;
#use Text::CSV::Hashify;
use Text::CSV_XS;
#use Scalar::Util;


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
print "Starting time 1: ",get_time(),"\n";

# Just to test out the script opening
print print_usage(),"\n";
print get_time(),"\n";
print "\n";
print "This is the input file: ", $infile,"\n";
print "This is the output file: ", $outfile,"\n";
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



my @data;
my $csvIN=Text::CSV->new;   #({sep_char => ','})
open my $fh, '<', $infile or die "Could not open";
while (my $line = $csvIN->getline($fh)){      #While there is data in the line
    
    push(@data,$line)
}
close $fh;

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
        if ($fields[0] eq "position"){next;}
        if (($fields[12] eq 0) or (!$fields[12])){next;}
        else{
            if (($fields[0] >= $Wstart) and ($fields[0]<=$Wend) and ($fields[12])){
                $Wsum+=$fields[12];
                $Wcount+=1;
            }
            if ($fields[0]>$Wend){         #if finished with that window, then:
                if ($Wcount!=0){
                    $WwindowAvg=$Wsum/$Wcount;
                }
                else{
                    next;
                    #$WwindowAvg=0;
                }
                $window[0]=$Wstart;
                $window[1]=$Wend;
                $window[2]=$WwindowAvg;
                $window[3]=$Wcount;
                
                return (\@window);
            }
        }
    }
}





print "Starting time: ",get_time(),"\n";
my $start=0;
my $end=$size;
my $windowNum=0;

#print "This is how many windows have been completed: ";

my @allWindows; #will be a 2D array containing all window info to be written into output file

while ($end<6000){     #($end<3000000){
    my($window)=OneWindow($start,$end);
    push (@allWindows,$window);
    $start=$start+$step;
    $end=$end+$step;
    $windowNum+=1;
    print "$windowNum  ";
}


use Tie::Array::CSV;
#my $aoa = csv (in => @allWindows);
#Make and write output csv file all at once
print "Writing the new file\n";
#my $csv = Text::CSV->new({ binary => 1, auto_diag => 1, eol => "\n"}) or die "Cannot use CSV: " . Text::CSV->error_diag();  # open in append mode
#open my $TheOutFile, ">>", "$outfile" or die "Failed to open file: $!";
#$csv->print($TheOutFile, [ "start", "end","fitness","mutant_count" ]); #header
#csv (in => $aoa, out => "$outfile", sep_char=> ";");
#$csv->print ($TheOutFile, $_) for @rows;
tie @allWindows, 'Tie::Array::CSV',$outfile;


print get_time();



#NOTES #############################################################################

#IMPROVEMENTS:
#Going to have to make new .wig file for BigWig for IGV
#Compile option for multiple file inputs
# Eventually should make sorting part of the script http://stackoverflow.com/questions/4221424/perl-sort-csv-on-a-certain-column


####################################################################################

