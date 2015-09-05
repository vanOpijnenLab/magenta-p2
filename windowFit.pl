#!/usr/bin/perl -w

#Margaret Antonio 15.04.11

#VERSION 14: Added two new options onto version 12: 1)--txt s ouputs a textfile of the data, could be bed file and 2)--txtg groups consecutive windows of same fitness into one text file (or bed) and ouputs [start,end,fitness,count] (of course, end is adjusted to last window added and count is cummulative across all windows that fit into the cumulative. 3)made 4 different output options (optional) wig,csv,txt,or txtg

#Compiled the .csv files using mergeCSV.pl and then sorted them by GNU sort in terminal from /results/
# $ sort -t, -k1 -n compiledPennG.csv > sortedCompiledPennG.csv

#../script/windowFit14.pl --ref=NC_003028b2.gbk --cutoff 15 --in results/sortedCompiledPennG.csv --csv windowFit/14B.csv --step 10 --size 500 --txtg ../viewer/14Bgrouped.txt --txt ../viewer/14B.txt --wig 14B.wig

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use Text::CSV_XS;
use Bio::SeqIO;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
sub print_usage() {
    print "\nRequired:\n";
    print "--size \t The size of the sliding window \n";
    print "--step \t The window spacing \n";
    print "--in \t The name of the file containing fitness values for individual insertion mutants.\n";
    print "--csv \t Name of a file to enter the .csv output for sliding windows.\n";
    print "--cutoff \tCutoff: Don't include fitness scores with average counts (c1+c2)/2 < x (default: 0)\n";
    
    print "\nMust choose at least one type of output:\n";
    print "--wig\tCreate a wiggle file for viewing in a genome browser. Provide a filename. Also provide genome under --ref\n";
    print "--txt\t Output all data [start,end,W,count] into a text of bed file.\n";
    print "--txtg\t If consecutive windows have the same value, then group them into one window. Ouput into txt file or bed file.\n";
    print "--ref\tThe name of the reference genome file, in GenBank format. Needed for wig and txt file creation\n";
}

#ASSIGN INPUTS TO VARIABLES
our ($txt,$txtg,$cutoff,$wig,$ref_genome,$infile, $csv, $step, $size);
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

);

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}
# Just to test out the script opening
print print_usage(),"\n";
print "\n";
print "Input file: ", $infile,"\n";
if ($csv){print "Output file: ", $csv,"\n";}
if ($txt){print "Text file for window data: $txt\n";}
if ($txtg){print "Text file for grouped windows: $txtg\n";}


#CHECKING PARAMETERS: Check to make sure required option inputs are there and if not then assign default
if (!$size) { $size=500 };   #set the default sliding window size to 500
if (!$step) { $step = 10 };   #set the default step size to 10
if (!$cutoff) {$cutoff=0};
print "Window size: $size\n";
print "Step value: $step\n";
print "Cutoff: $cutoff\n";

if ((!$csv) and (!$txt) and (!$txtg) and (!$wig)) {
    print "\nThere must be some kind of output file assigned (--csv, --txt, or --txtg)\n";
    print_usage();
    die;
}


###############################################################################################

print "\n---------This is the sliding window fitness calculation part--------\n\n";

print "Start input array ",get_time(),"\n";
#CREATE AN ARRAY OF DATA FROM INPUT CSV FILE

my @data;
my $csvIN=Text::CSV->new;   #({sep_char => ','})
open my $fh, '<', $infile or die "Could not open";
$csvIN->getline($fh);
my $rowCount=-1;
my $last;
while (my $line = $csvIN->getline($fh)){      #While there is data in the line
    my $w = $line->[12];
    #if ($w and $w eq 'nW') {next;}
    if (!$w){next;} # For blanks
    else{
        my $c1 = $line->[2];
        my $c2 = $line->[3];
        my $avg = ($c1+$c2)/2;
        if ($avg < $cutoff) { next; } # Skip cutoff genes.
        else {
            my @select=($line->[0],$line->[12]);
            my $select=\@select;
            push(@data,$select);
            $rowCount+=1;
            $last=$select->[0];
        }
    }
}
close $fh;
print "Finished input array ",get_time(),"\n\n";

my $index=-1;
my $marker=0;

#SUBROUTINE FOR EACH WINDOW CALCULATION
sub OneWindow{
    my $Wstart=shift @_;
    my $Wend=shift@_;
    my $WwindowAvg=0;
    my $Wcount=0;
    my $Wsum=0;
    my $i;
    for ($i=$marker;$i<$rowCount;$i++){
        my @fields=@{$data[$i]};
        if ($fields[0]<($Wstart+$step)){
            $index++;
        }
        if ($fields[0]<$Wstart){  #if deleted, error shows up
            next;
        }
        if ($fields[0]<=$Wend){
            if ($fields[0]<=($Wstart+$step)) {$marker=$index;}
            $Wsum+=$fields[1];
            $Wcount+=1;
        }
        else{
            if ($fields[0]>$Wend){         #if finished with that window, then:
                if ($Wcount!=0){
                    $WwindowAvg=sprintf("%.2f",$Wsum/$Wcount);
                    my @Wwindow=($Wstart,$Wend,$WwindowAvg,$Wcount);
                    #print @Wwindow;
                    return (\@Wwindow);
                }
                else{return -1}  #Becuse count=0 (i.e. there were no insertion mutants in that window)
                
            }
        }
    }
}

print "Start calculation: ",get_time(),"\n";
my $start=1;
my $end=$size+1;
my $windowNum=0;
my @allWindows=(); #will be a 2D array containing all window info to be written into output file

#WHILE LOOP TO CALL THE ONE WINDOW SUBROUTINE FOR CALCULATIONS===INCREMENTS START AND END VALUES OF THE WINDOW
while ($end<=$last-$size){  #100,000bp-->9,950 windows--> only 8500 windows in csv because 0
    my($window)=OneWindow($start,$end);
    if ($window!=-1){
        push (@allWindows,$window);
    }
    $start=$start+$step;
    $end=$end+$step;
    #print "$end  ";
}
print "End calculation: ",get_time(),"\n\n";

#MAKE OUTPUT CSV FILE WITH WINDOW CALCULATIONS
if ($csv){
    print "Start csv ouput file creation: ",get_time(),"\n";
    my $csvBIG = Text::CSV->new({ binary => 1, auto_diag => 1, eol => "\n"}) or die "Cannot use CSV: " . Text::CSV->error_diag();  # open in append mode
    open my $file, ">", "$csv" or die "Failed to open file";
    $csvBIG->print($file, [ "start", "end","fitness","mutant_count" ]); #header
    foreach my $winLine(@allWindows){
        $csvBIG->print($file,$winLine);
    }
    close $file;
    print "End csv ouput file creation: ",get_time(),"\n\n";
}

#MAKE WIG FILE---->later make BW--->IGV
if ($wig){
    print "Start wig file creation: ",get_time(),"\n";
    my $in = Bio::SeqIO->new(-file=>$ref_genome);
    my $refseq = $in->next_seq;
    my $refname = $refseq->id;
    open WIG, ">$wig";
    print WIG "track type=wiggle_0 name=$wig\n";
    print WIG "variableStep chrom=$refname\n";
    foreach my $wigLine(@allWindows){
        my @wigFields=@$wigLine;
        my $position=$wigFields[0];
        #while ($position<=$wigFields[1]){
        print WIG $position," ",$wigFields[2],"\n";
        #$position=$position+1;
        #}
        #print  WIG $wigFields[0]," ",$wigFields[2],"\n";
    }
    close WIG;
    print "End wig file creation: ",get_time(),"\n\n";
    print "If this wig file needs to be converted to a Big Wig, then use USCS program wigToBigWig in terminal: \n \t./wigToBigWig gview/12G.wig organism.txt BigWig/output.bw \n\n";
}

#IF GOING TO MAKE A TEXT FILE FOR BED CONVERSION TO BIGBED, NEED CHROM # IN COLUMN 0
my @cummulative;
if ($txtg or $txt){
    for my $line(@allWindows){
        unshift($line, "NC_003028");
    }
}
#IF MAKING A REGULAR TEXT FILE fields: [chrom,start,end,fitness,count]
if ($txt){
    print "Start text file creation time: ",get_time(),"\n";
    open my $TXT, '>', "$txt" or die $!;
    print $TXT (join("\t",@$_),"\n") for @allWindows;
    close $TXT;
    print "End text file creation: ",get_time(),"\n\n";
}


#IF MAKING A TEXT FILE OF GROUPED CONSECUTIVE WINDOWS WITH SAME FITNESS
if ($txtg){
    print "Start grouped txt file creation time: ",get_time(),"\n";
    open my $TXTg, '>', "$txtg" or die $!;
    for my $line(@allWindows){
        my @field=@$line;
        if (!@cummulative){
            @cummulative=@field;
            if ($cummulative[4]>1000){$cummulative[4]=1000}
        }
        else{
            if ($cummulative[3]==$field[3]){
                $cummulative[2]=$field[2];
                if ($cummulative[4]<=1000-$field[4]){
                    $cummulative[4]+=$field[4];
                }
            }
            else{
                print $TXTg ($_,"\t") for @cummulative;
                print $TXTg ("\n");
                @cummulative=@field;
            }
        }
    }
    close $TXTg;
    print "End grouped text file creation: ",get_time(),"\n\n";
    if ($txt or $txtg){
        print "\nTo make a BigBed file from this text file, rename file to .bed and use USCS program bedToBigBed in terminal \n\t\n";
    }
}





