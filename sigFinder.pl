#!/usr/bin/perl -w

#Margaret Antonio 16.01.08

#DESCRIPTION: Exapnding windows to identify max region of significance (essential or cold spot)

#cpanm Getopt::Long Data::Random List::Util File::Path File::Basename
#cpanm List::BinarySearch List::BinarySearch::XS List::MoreUtils


use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
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

use Text::CSV;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
sub print_usage() {
    print "\nRequired:\n";
    print "In the command line (without a flag), input the name(s) of the file(s) containing fitness values for individual insertion mutants.\n";
    print "\n slidingWindow.pl <OPTIONS> <OUTPUT> <--essentials AND --genome> <INPUT CSV FILES>\n\n";
    
    print "\nOPTIONS:\n";
    print "--size \t The size of the sliding window(default=500) \n";
    print "--step \t The window spacing (default=10) \n";
    #print "--csv \t Name of a file to enter the .csv output for sliding windows.\n";
    print "--cutoff \tCutoff: Don't include fitness scores with average counts (c1+c2)/2 < x (default: 0)\n";
    print "--essentials \t Calculate genome region essentiality based on transposon insertion representation\n";
    print "--outdir\tSpecify name of new directory for all output files\n";
    print "--log\t Send all output to a log file instead of the terminal\n";
    print "--indir\t Directory containing all input files (results files from calc fitness script\n";
    print "--usage\t Print usage\n";
    print "--tan\t Expected max number of TA sites in each window---used for creating null distribution library (default: 200)";
    print "--d\t Size of null distributions (default: 1000)";
    
    
    print "\nREQUIRED: Must choose at least one type of output:\n";
    print "--wig\tCreate a wiggle file for viewing in a genome browser. Provide a filename. Also provide genome under --ref\n";
    print "--txt\t Output all data [start,end,W,count] into a text of bed file.\n";
    print "--txtg\t If consecutive windows have the same value, then group them into one window. Ouput into txt file or bed file.\n";
    print "--ref\tThe name of the reference genome file, in GenBank format. Needed for wig and txt file creation\n";
    
}


#ASSIGN INPUTS TO VARIABLES
our ($round,$random,$txt,$txtg,$cutoff,$wig,$infile, $csv, $step, $h, $outdir,$size,
	$fasta, $log, $ref_genome,$tan,$indir,$d,$minSize,$maxSize,$cond);
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
'random:s' =>\$random,
'round:i' =>\$round,
'fasta:s' => \$fasta,
'outdir:s' =>\$outdir,
'log' => \$log,
'usage' => \$h,
'tan'=>\$tan,
'indir:s'=>\$indir,
'd:i' =>\$N,
'min:i' =>\$minSize,
'max:i' =>\$maxSize,
'inc:i'=>\$inc,
);

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}
# Just to test out the script opening
if ($h){
    print print_usage(),"\n";
    print "\n";
    if ($csv){print "CSV output file: ", $csv,"\n";}
    if ($txt){print "Text file for window data: $txt\n";}
    if ($txtg){print "Text file for grouped windows: $txtg\n";}
}
if (!$round){$round='%.3f';}
if (!$outdir){
    $outdir="F-151228";
}
mkpath($outdir);

if ($log){
    print "\nSending all output to log file\n";
    # redirect STDOUT to log.txt
    open (STDOUT, ">>$outdir/log.txt");
}


#CHECKING PARAMETERS: Check for option inputs, if none then assign default
if (!$size) { $size=500 };   #set the default sliding window size to 500
if (!$step) { $step=10 };   #set the default step size to 10
if (!$cutoff) {$cutoff=15};
if (!$tan){$tan=200};
if (!$N){$N=10000};
if (!$minSize){$minSize=500};
if (!$maxSize){$maxSize=10000};   #if all works efficiently, default max should be ~5000 bp
if (!$sigCut){$sigCut=.05};
if (!$inc){$inc=2000};

print "Window size: $size\n";
print "Step value: $step\n";
print "Cutoff: $cutoff\n";







######################## PART 1: Read input files into manageable array #####################

#Create an array of data from input csv files(s)

my $rowCount=-1;
my $last=0;
my @unsorted;
my @insertPos; #array to hold all positions of insertions. Use this later to match up with TA sites

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
my $num=(scalar @files);

my $datestring = localtime();
print "Local date and time $datestring\n";

print "\n---------Importing files--------\n";
print "\tStart input array ",get_time(),"\n";
print "\tNumber of csv files: ", $num,"\n";


#Go through each file from the commandline (ARGV array) and read each line as an array
#into select array if values satisfy the cutoff
my $index=0;

for (my $i=0; $i<$num; $i++){   #Read files from ARGV
    print "File #",$i+1,"\t";
    
    my $file=$files[$i];
    print $file,"\n";
    
    open(DATA, '<', $file) or die "Could not open '$file' \n";
    my $dummy=<DATA>;
    while (my $entry = <DATA>) {
        chomp $entry;
        my @line=split(",",$entry);
        my $w = $line[12];
        if (!$w){next;} # For blanks
        else{
            my $c1 = $line[2];
            my $c2 = $line[3];
            my $avg = ($c1+$c2)/2;
            if ($avg > $cutoff) {
                my @select=($index,$line[0],$line[12]);
                my $select=\@select;
                push(@unsorted,$select);
                push(@insertPos,$line[0]);   #keep track of actual insertion site position
                $last=$select[0];
                $rowCount++;
                $index++;
            }
        }
    }
    close DATA;
}

my @sorted = sort { $a->[0] <=> $b->[0] } @unsorted;
@insertPos = sort { $a <=> $b } @insertPos;
@insertPos= uniq @insertPos;
print "\n\tFinished input array ",get_time(),"\n";





########################## PART 2: Whether insertions happened at TA sites #####################

#Counting the number of TA sites in the genome and whether an insertion occurred there or not

print "\n---------Assessing essentiality of genome region in each window--------\n\n";

my @sites;

#First read fasta file into a string
my $seqio = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');
my $prev;
my $total=0;
while(my $seq = $seqio->next_seq) {
    $fasta = $seq->seq;
}
#Get number of "TA" sites, regardless of insertion---this is just the fasta file
my $x="TA";
my @c = $fasta =~ /$x/g;
my $countTA = @c;

#At what positions in the genome do TA sites occur?
my $pos=0;
my $countInsert=0;
my @selector;   #use this array to hold all 1 and 0 of whether insertion happened or not.

#Go through all TA sites identified above and see if an insertion occurs there.
#Push results onto two arrays a 2d array with the position and a 1D array with just the 0 or 1

my @unmatched;  #hold all unmatched ta sites
my @allTAsites; #2d array to hold all occurences of TA sites in genome
my $unmatchedCount=0;
my $offset=0;
my $genPos = index($fasta,'TA',$offset);

while (($genPos != -1) and ($pos!=scalar @insertPos)) { #as long as the TA site is foun
    my $res=0;
    if ($genPos>$insertPos[$pos]){
        push @unmatched,$insertPos[$pos];
        $unmatchedCount++;
        $pos++;
    }
    if ($genPos==$insertPos[$pos]){
        $res=1;
        $countInsert++;
        $pos++;
    }
    my @sites=($genPos,$res);
    push @selector,$res;
    #push the 0 or 1 onto the array @selector---use this to draw random sets for the null distr
    push (@allTAsites,\@sites);
    $offset = $genPos + 1;
    $genPos = index($fasta, 'TA', $offset);
    $countTA++;
}
my $FILE1 = "$outdir/allTAsites.txt";
open (ALL_TA, ">", $FILE1);
foreach my $sit(@allTAsites){
    foreach (@$sit){
        print ALL_TA $_, "\t";
    }
    printf ALL_TA "\n";
}
close ALL_TA;
my $FILE2 = "$outdir/unmatched.txt";
unless(open UNM, ">", $FILE2){
    die "\nUnable to create $FILE2:\n$!";
}
foreach (@unmatched){
    print UNM $_, "\n";
}
close UNM;

print "\n\nTotal of unmatched insertions: $unmatchedCount\n";
print "See unmatched.txt for genome indices of unmatched sites\n";
print "See allTAsites.txt for details on all TA sites and insertions\n\n";

#print "\nTotal: $countInsert insertions in $countTA TA sites.\n";







################# PART 3: Create Null Distributions  ##############################

#Now, have an array for each TA site and if an insertion occurred there.
    #So per site @sites=(position, 0 or 1 for insertion).
#Next step, create null distribution of 10,000 random sets with same number
    #of TA sites as the window and then calculate p-value

#SUBROUTINE FOR MAKING THE NULL DISTRIBUTION SPECIFIC TO THE WINDOW

#DISTRIBUTION'S STANDARD DEVIATION

sub mean {
    return sum(@_)/@_;
}
sub stdev{
    my @data = @{shift @_};
    my $average=shift @_;
    my $sqtotal = 0;
    foreach(@data) {
        $sqtotal += ($average-$_) ** 2;
    }
    my $std = ($sqtotal / ($N-1)) ** 0.5;
    return $std;
}

#MAKE LIBRARY OF NULL DISTRIBUTIONS:
print "Making a library of null distributions.\n For information on distribution library, see nullDist.txt\n";

my @distLib;

my $FILE3 = "$outdir/nullDist.txt";
unless(open DIST, ">", $FILE3){
    die "\nUnable to create $FILE3:\n$!";
}
printf DIST "Sites\tSize(n)\tMin\tMax\tMean\tstdev\tvar\tMinPval\n";

#Loop once for each distribution in the library
#(i.e. distribution of 10,000 sites each with 35 TA sites, then a distribution of 10,000 sites each with 36 TA sites, etc)

sub pvalue{
    
    #takes in window count average (countAvg) and number of TAsites and makes a null distribution to calculate the pvalue, which it returns
    my $mean=shift@_;
    my $TAsites=shift@_;
    my $N=10000;
    my @specDist=@{$distLib[$TAsites-1]};
    my $rank= binsearch_pos { $a cmp $b } $mean,@specDist;
    my $i=$rank;
    while ($i<scalar(@specDist)-1 and $specDist[$i+1]==$specDist[$rank]){
        $i++;
    }
    $rank=$i;
    my $pval=$rank/$N; #calculate pval as rank/N
    return $pval;
    
}

for (my $sitez=1; $sitez<=$tan;$sitez++){
    #print "In the first for loop to make a null distribution\n";
    my @unsorted;
    my $count=0;
    my $sum=0;
    
    for (my $i=1; $i<=$N; $i++){
        my @random_set = rand_set( set => \@selector, size => $sitez);
        my $setMean=mean(@random_set);
        push (@unsorted, $setMean);
        #print "$i:\t", "$setMean\n";
        $count++;
        $sum+=$setMean;
    }
    my @nullDist= sort { $a <=> $b } @unsorted;
    my $nullMean=sprintf("$round",($sum/$count));
    my $standev =sprintf("$round",stdev(\@nullDist, $nullMean));
    my $variance=sprintf("$round",($standev*$standev));
    my $min=sprintf("$round",$nullDist[0]);
    my $max=sprintf("$round",$nullDist[scalar @nullDist-1]);
    my $setScalar=scalar @nullDist;
    push (@distLib,\@nullDist);
    my $minp=pvalue(0,$sitez);
    printf DIST "$sitez\t$N\t$min\t$max\t$nullMean\t$standev\t$variance\t$minp\n";
}
close DIST;







#########################  PART 4:  Functions  #########################



#SigCheck()
#extendFront()
#extendBack()

sub pvalue{
    
    #takes in window count average (countAvg) and number of TAsites and makes a null distribution to calculate the pvalue, which it returns
    my $mean=shift@_;
    my $TAsites=shift@_;
    my $N=10000;
    my @specDist=@{$distLib[$TAsites-1]};
    my $rank= binsearch_pos { $a cmp $b } $mean,@specDist;
    my $i=$rank;
    while ($i<scalar(@specDist)-1 and $specDist[$i+1]==$specDist[$rank]){
        $i++;
    }
    $rank=$i;
    my $pval=sprintf("$round",$rank/$N); #calculate pval as rank/N
    return $pval;
    
}

sub extendFront{
	my ($start,$end,$pval)=shift @_;
	my $nextEnd=$end+$inc;
	my $nextpval=sigCheck($start,$nextEnd);
	#Stop condition
	if (($nextEnd-$start>$maxSize) or ($nextpval>$cond)) {
		return ($start,$end,$pval);
		}
	#Recursion Condition: when ($nextpval<=$cond)
	else{     
		extendFront($start,$newEnd);
		}
	}
		
sub extendBack{
	my ($start,$end,$pval)=shift @_;
	my $nextStart=$start-$inc;
	my $nextpval=sigCheck($start,$nextStart);
	#Stop condition
	if (($end-$nextStart>$maxSize) or ($nextpval>$cond)) {
		return ($start,$end,$pval);
		} #Return without extending back
	#Recursion Condition: when ($nextpval<=$cond)
	else{ 
		$indexMarker+=
		extendBack($start,$nextStart);
		}
	}

my $index=-1;
my $fm=0;
my $lm=0;

#Fields[0]=index   Fields[1]=Position
sub sigCheck{ #WILL RETURN A WINDOW
my $count=0;    #keeps track of the number of insertions in a window
my $start=shift @_;
my $end=shift@_;
my $avg=0;
my $sum=0;
my $lastPos=0;
my $i;

#$i=$firstMarker
	my ($start,$end,$indexMarker)=shift @_;
	for (my $i=$indexMarker;$i<$rowCount;$i++){
        my @fields=@{$sorted[$i]};
        my $pos=$fields[1]; 
        my $index=$fields[0];
        my $fit=$fields[2];
        #Shouldn't need this
        if ($pos<$Wstart){  #if deleted, error shows up
            next;
        }  
        elsif ($pos<=$end) and ($pos>=$start){
            $count++; #For significance
            $fitSum+=$fit; #For importance
            next;
			}
        else{   
            #if finished with that window, then:
       
            if ($count!=0){
                my $avgFit=sprintf("%.2f",$fitSum/$count);
                }
            my $len=
            my $seq = substr($fasta,$start-1,$end-$start);  #start-1 becase $start and $end are positions in genome starting at 1,2,3.... substr(string,start, length) needs indexes
    		my $ta="TA";
   			my @c = $seq =~ /$ta/g;
   			my $TAsites = scalar @c;
   			my $ratio=$count/$TAsites;
   			my $pval=pvalue($ratio,$TAsites);
            my @window=($start,$end,$fm,$lm,$count,$TAsites,$ratio,$pval);
    		return \@window;
    		}
        }
	}

		
	
	
	





#########################  PART 5: Begin window ################################

#Declare these vars: $minSize, $maxSize, $cond

print "Start calculation: ",get_time(),"\n";

my @allWindows=(); #will be a 2D array containing all window info to be written into output file


#WHILE LOOP TO CALL THE ONE WINDOW SUBROUTINE FOR CALCULATIONS===INCREMENTS START AND END VALUES OF THE WINDOW
#fm=0; lm=0;
my ($start,$end,$pval)=(1,$minSize+1,0);
my $genLen=length($fasta);
while ($end <= $genLen){  #100,000bp-->9,950 windows--> only 8500 windows in csv because 0
    my $sCreturn=sigCheck($start,$end,$fm);
    #if min size window is significant, expand to see if surrounding regions can be incl.
    if ($pval<$cond){  
    	($start,$end,$pval)=checkFront($start,$end+$inc,$pval);
    	($start,$end,$pval)=checkBack($start-$inc,$end,$pval);
    	#resulting $start, $end, and $pval should be of max size sig. region
    	my @sigRegion=($start,$end,$pval);
    	push (@allSigWindows,\@sigRegion);
    	}
    $start=$end+$minSize;  #Not sure if this increment makes sense
    $end=$start+$minSize;
}

print "End calculation: ",get_time(),"\n";


#SUBROUTINE FOR PRELIMINARY SIGNIFICANCE CHECK
#Is there underrepresentation of insertions at minimum window size. If yes, then expand window.

my $index=-1;
my $firstMark=0;
my $lastMark=0;
my $totalInsert=0;
my $totalWindows=0;

sub sigCheck{
    my $start=shift @_;
    my $end=shift@_;
    my $Wavg=0;
    my $Wcount=0;
    my $insertion=0;
    my $Wsum=0;
    my $lastPos=0;
    my $i;
    
    for ($i=$marker;$i<$rowCount;$i++){
        my @fields=@{$sorted[$i]};
        if ($fields[0]<$Wstart){  #if deleted, error shows up
            next;
        }
        if ($fields[0]<=$Wend){
            if ($fields[0]<($Wstart+$step)){
                $marker++;
            }
            $Wsum+=$fields[1];
            $Wcount++;
            if ($fields[0]!=$lastPos){
                $insertion+=1;
                $lastPos=$fields[0];
            }
        }
        
        else{   #if ($fields[0]>$Wend){         #if finished with that window, then:
            if ($Wcount!=0){
                $Wavg=sprintf("%.2f",$Wsum/$Wcount);
                my @window=($Wstart,$Wend,$Wavg,$Wcount,$insertion);
                $totalWindows++;
                $totalInsert+=$insertion;
                #print @Wwindow, "\n";
                return (\@window);
            }
            
            else{ #Even if there were no insertions, still want window in file for consistent start/end
                my @window=($Wstart,$Wend,0,0,0);
                return (\@window);
            }  	#Because count=0 (i.e. there were no insertion mutants in that window)
        }
    }
}









#########################  PART 6: Expand Window  ##################################

#Quick and dirty way: check back and forward 2000/step=2000/10=200 windows
#Better way to do it in case file is not from sliding windows: rec. until 2000 bef and aft

sub isEssentialBack{
	my ($currentIndex,$status,$stop)=@_;
	if ($currentIndex==0){
		return $status;
		}
	my @currentRegion=@{$outArray[$currentIndex]};
	my $regionStart=$currentRegion[0];
	my $regionEss=$currentRegion[7]; #p-value for essentiality
	if ($regionStart<$stop){
		return $status; #Problem: for first window where roiStart=2, will return status=1
	}
	elsif ($regionEss>$defEss){return 2}
	else{isEssentialBack($currentIndex-1,$status,$stop)}
}

sub isEssentialForward{
	my ($currentIndex,$status,$stop)=@_;
	#my @outArray=@{$wholeArray};
	my $last=scalar @outArray-1;
	if ($currentIndex==$last){return $status}
	my @currentRegion=@{$outArray[$currentIndex]};
	my $regionEnd=$currentRegion[1];
	my $regionEss=$currentRegion[7]; #p-value for essentiality
	if ($regionEnd>$stop){
		return $status; #Problem: for first window where roiStart=2, will return status=1
	}
	if ($regionEss>$defEss){   #Does this ever execute? Was error that wasn't caught
		return 2;
	}
	else{
		isEssentialForward($currentIndex+1,$status,$stop);
	}
}

my @finalOutput;
my $lastLine=scalar @outArray;

for (my $roiIndex=0;$roiIndex<$lastLine;$roiIndex++){

	my @roi=@{$outArray[$roiIndex]};
	my $backResult=0; #Default no test, then backResult=0
	my $forwardResult=0;
	my $isCold=0; #Region has not been tested / will not be tested
	#only check windows that are "essential" i.e. satisfy $defEss
	if ($roi[7]<$defEss){
		$isCold=1; #Region was tested for cold spot and by default is a cold spot.
		my $backStop=$roi[0]-2000;
		my $forwardStop=$roi[1]+2000;
    	$backResult=isEssentialBack($roiIndex,1,$backStop);
    	$forwardResult=isEssentialForward($roiIndex,1,$forwardStop);
    	if ($backResult!=1 or $forwardResult!=1){ 
			$isCold=2; #Region is not a cold spot given it was tested
		}
	}
			
	push (@roi,$isCold);
	push (@finalOutput,\@roi);
	#print $roiIndex,"\t";
}


my $csv = Text::CSV->new({ binary => 1, auto_diag => 1, eol => "\n"}) or die "Cannot use CSV: " . Text::CSV->error_diag(); 
open (my $OUT, ">$out"); 

$csv->print($OUT, \@header); #header

foreach (@finalOutput){
	$csv->print ($OUT, $_);
	
}
close $OUT;


