#!/usr/bin/perl -w

#Margaret Antonio 15.12.10

#15.12.10: Things To Do for ESSENTIALS:
    #1) If region is essential, then fitness of single insertion doesn't matter. Instead look at 2000 bp before and after to check if it is a cold spot. If many insertions before and after essential region, then it is a TRUE essential region. DONE
    #2) Should output everything including regions that don't have "enough" TA sites DONE

#SLIDING WINDOW: a culmination of fitness calculation, essentiality determination,
    #and TA site and insertion statistics for automatically generated sliding widnow throughout genome

#perl ../Blueberries/slidingWindow.pl --size 500 --fasta tigr4_genome.fasta --indir 2394_PennG/ --log --ref NC_003028b2.gbk

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
    print "--csv \t Name of a file to enter the .csv output for sliding windows.\n";
    print "--cutoff \tCutoff: Don't include fitness scores with average counts (c1+c2)/2 < x (default: 0)\n";
    print "--essentials \t Calculate genome region essentiality based on transposon insertion representation\n";
    print "--outdir\tSpecify name of new directory for all output files\n";
    print "--log\t Send all output to a log file instead of the terminal\n";
    print "--indir\t Directory containing all input files (results files from calc fitness script\n";
    print "--usage\t Print usage\n";
    print "--tan\t Max number of TA sites in each window---used for creating null distribution library (default: 100)";
    
    
    print "\nREQUIRED: Must choose at least one type of output:\n";
    print "--wig\tCreate a wiggle file for viewing in a genome browser. Provide a filename. Also provide genome under --ref\n";
    print "--txt\t Output all data [start,end,W,count] into a text of bed file.\n";
    print "--txtg\t If consecutive windows have the same value, then group them into one window. Ouput into txt file or bed file.\n";
    print "--ref\tThe name of the reference genome file, in GenBank format. Needed for wig and txt file creation\n";
    print "--sort\tIndex of column to sort by\n";

}


#ASSIGN INPUTS TO VARIABLES
our ($round,$random,$txt,$txtg,$cutoff,$wig,$infile, $csv, $step, $h, $outdir,$size,$fasta, $log, $ref_genome,$tan,$indir,$inc,$sortby);
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
'inc:i'=>\$inc,
'sort:i'=>\$sortby,
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
if (!$outdir){$outdir="F-151228";}
if (!$inc){$inc=20;}
mkpath($outdir);

if ($log){
	print "\nSending all output to log file\n";
	# redirect STDOUT to log.txt
    open (STDOUT, ">>$outdir/log.txt");
}
if (!$sortby){
    $sortby=8;
}


#CHECKING PARAMETERS: Check to make sure required option inputs are there and if not then assign default
if (!$size) { $size=500 };   #set the default sliding window size to 500
if (!$step) { $step=10 };   #set the default step size to 10
if (!$cutoff) {$cutoff=15};
if (!$tan){$tan=100};

print "Window size: $size\n";
print "Step value: $step\n";
print "Cutoff: $cutoff\n";

#if ((!$csv) and (!$txt) and (!$txtg) and (!$wig)) {
#print "\nThere must be some kind of output file assigned (--csv, --txt, or --txtg)\n";
#print_usage();
#die;
#}

#CREATE AN ARRAY OF DATA FROM INPUT CSV FILE(S)

my $rowCount=-1;
my $last=0;
my @unsorted;
my @insertPos; #array to hold all positions of insertions. Going to use this later to match up with TA sites

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


#Go through each file from the commandline (ARGV array) and read each line as an array into select array if values satisfy the cutoff
for (my $i=0; $i<$num; $i++){   #Read files from ARGV
	print "File #",$i+1,"\t";
    
    my $file=$files[$i];
    print $file,"\n";
    
    open(DATA, '<', $file) or die "Could not open '$file' Make sure input .csv files are entered in the command line\n";
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
                my @select=($line[0],$line[12]);
                my $select=\@select;
                push(@unsorted,$select);
                push(@insertPos,$line[0]);   #keep track of actual insertion site position
                $last=$select[0];
                $rowCount++;  
            }
        }
    }
    close DATA;
}

my @sorted = sort { $a->[0] <=> $b->[0] } @unsorted;

@insertPos = sort { $a <=> $b } @insertPos;
@insertPos= uniq @insertPos;

#for (my $i=0;$i<30;$i++){
#    foreach my $element ( @{ $sorted[$i] }){
#        print $element,"\t";
#    }
#    print "\n";
#}

print "\n\tFinished input array ",get_time(),"\n";

###############################################################################################

print "\n---------Creating sliding windows across the genome--------\n\n";

my $index=-1;
my $marker=0;
my $totalInsert=0;
my $totalWindows=0;


#SUBROUTINE FOR EACH WINDOW CALCULATION
sub OneWindow{
    my $Wstart=shift @_;
    my $Wend=shift@_;
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




print "Start calculation: ",get_time(),"\n";
my $start=1;
my $end=$size+1;
my $windowNum=0;
my @allWindows=(); #will be a 2D array containing all window info to be written into output file

#WHILE LOOP TO CALL THE ONE WINDOW SUBROUTINE FOR CALCULATIONS===INCREMENTS START AND END VALUES OF THE WINDOW
while ($end<=$last-$size){  #100,000bp-->9,950 windows--> only 8500 windows in csv because 0
    my($window)=OneWindow($start,$end);
    #if ($window!=-1){
        push (@allWindows,$window);
    #}
    $start=$start+$step;
    $end=$end+$step;
}
print "End calculation: ",get_time(),"\n";

my $avgInsert=$totalInsert/$totalWindows;
print "Average number of insertions for $size base pair windows: $avgInsert\n";

#ESSENTIALS: Counting the number of TA sites in the genome and whether an insertion occurred there or not

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
        	#push the 0 or 1 onto the array @selector---going to use this to draw random sets for the null distribution
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

#--------------------------------------------------------------------------------------------------
    
    #Now, have an array for each TA site and if an insertion occurred there. So per site @sites=(position, 0 or 1 for insertion).
    #Next step, create null distribution of 10,000 random sets with same number of TA sites as the window and then calculate p-value
  
#SUBROUTINE FOR MAKING THE NULL DISTRIBUTION SPECIFIC TO THE WINDOW

#DISTRIBUTION'S STANDARD DEVIATION
my $N=10000;

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

#SUBROUTINE TO CALCULATE THE P-VALUE OF THE WINDOW INSERTIONS AGAINST THE NULL DISTRIBUTION    
#---------------------------------------------------------------------------------------------------
    print "\n In case you were wondering....the size of genome is: ", length($fasta), " bp\n";
    
    #Now we have an array called @allTAsites which contains every TAsite position with a 0 next to it for "no insertion".
    #Now just need to replace 0 with 1 if there IS and insertion at that site

    
    my @newWindows=();
    my $printNum=0;

print "Start p-value calculation for individual windows: ",get_time(),"\n\n";

#my $allWindows=\@allWindows;
for (my $i=0;$i<scalar @allWindows;$i++){
    my @win=@{$allWindows[$i]};
    my $starter=$win[0];
    my $ender=$win[1];
    #print "num $printNum -->\tStart pos: $starter\tEnd pos: $ender\n";
    #How many TA sites are there from $genome[$start] to $genome[$end]?

    my $seq = substr($fasta,$starter-1,$size);  #start-1 becase $start and $end are positions in genome starting at 1,2,3.... substr(string,start, length) needs indexes
    my $ta="TA";
    my @c = $seq =~ /$ta/g;
    my $TAsites = scalar @c;
    push(@win,$TAsites);
    my $countAvg=sprintf("$round",$win[4]/$TAsites);
    push (@win,$countAvg);
    my $pval=pvalue($countAvg,$TAsites);
    push (@win,$pval);
    push (@newWindows,\@win);
    $printNum++;
    }

   
 print "End p-value calculation: ",get_time(),"\n\n";

#print "This is the TA count: $count\nTotal genome size is: $total\n\n";


### For all windows, add a field that has the difference between that window's mean fitness and the
#average mean fitness for all of the windows

my @allFits = map $_->[2], @newWindows;
my $meanFit=mean(@allFits);
print "Average fitness for all windows: ",$meanFit,"\n";
my @expWindows=();
foreach (@newWindows){
    my @entry=@{$_};
    my $mean=$entry[2];
    my $absdev=sprintf("$round",abs($mean-$meanFit));
    push @entry,$absdev;        #now becomes index [8], column 7 after pvalue
    push @expWindows,\@entry;
}
@newWindows= sort {$b->[$sortby]<=>$a->[$sortby]} @expWindows;

#-------------------------------------------Essentials OUTPUTS------------------------------------------------------------

#MAKE OUTPUT CSV FILE WITH ESSENTIAL WINDOW CALCULATIONS

    print "Start csv ouput file creation: ",get_time(),"\n";
    my $csvBIG = Text::CSV->new({ binary => 1, auto_diag => 1, eol => "\n"}) or die "Cannot use CSV: " . Text::CSV->error_diag();  
    # open in append mode
    
	open (my $FH8, ">$outdir/slidingWindows.csv");
	
    $csvBIG->print($FH8, [ "start", "end","fitness","mutant_count","insertions","TA_sites","ratio","p-value","fitness-meanFit"]); #header
    foreach my $winLine(@newWindows){
        $csvBIG->print($FH8,$winLine);
    }
    close $FH8;
    print "End csv ouput file creation: ",get_time(),"\n\n";

#printwig();
my $in = Bio::SeqIO->new(-file=>$ref_genome);
   my $refseq = $in->next_seq;
    my $refname = $refseq->id;

#MAKE essentials WIG FILE---->later make BW--->IGV   #####
 sub printwig{
    print "Start wig file creation: ",get_time(),"\n";
    my $in = Bio::SeqIO->new(-file=>$ref_genome);
    my $refseq = $in->next_seq;
    my $refname = $refseq->id;
    my $ewig="essentialWindows.wig";
    my $FILE5 ="$outdir/$ewig";
    open (eWIG, ">$FILE5");
    printf eWIG "track type=wiggle_0 name=$ewig\n";
    printf eWIG "variableStep chrom=$refname\n";
    foreach my $wigLine(@newWindows){
        my @wigFields=@$wigLine;
        my $position=$wigFields[0];
        while ($position<=$wigFields[1]){
        printf eWIG $position, " ",$wigFields[7],"\n";    #7 for pvalue, but 2 for fitness!!!!!!
        $position=$position+1;
        }
    }
    close eWIG;
    print "End wig file creation: ",get_time(),"\n\n";
    print "If this wig file needs to be converted to a Big Wig, then use USCS program wigToBigWig in terminal: \n \t./wigToBigWig gview/12G.wig organism.txt BigWig/output.bw \n\n";
}


#IF GOING TO MAKE A TEXT FILE FOR BED CONVERSION TO BIGBED, NEED CHROM # IN COLUMN 0
my @ecummulative;
my @temp=split("$ref_genome",".");
my $refName=$temp[0];
if ($txtg or $txt){
    for my $line(@newWindows){
        unshift($line, "$refName");
    }
}
#IF MAKING A REGULAR TEXT FILE fields: [chrom,start,end,fitness,count]
#Don't think this is necessary....what would we do with a regular text file?:
if ($txt){
    print "Start text file creation time: ",get_time(),"\n";
    open my $TXT, '>', "$outdir/$txt" or die "\nUnable to create $txt:\n$!";
    print $TXT (join("\t",@$_),"\n") for @newWindows;
    close $TXT;
    print "End text file creation: ",get_time(),"\n\n";
}

#IF MAKING A TEXT FILE OF GROUPED CONSECUTIVE WINDOWS WITH SAME FITNESS


#Edit: what is the 1000?....probably should be a variable

    my $etxtg="groupedWindows.txt";
    print "Start grouped txt file creation time: ",get_time(),"\n";
    open my $eTXTg, '>', "$outdir/$etxtg" or die $!;
    for my $line(@newWindows){
        my @field=@$line;
        if (!@ecummulative){
            @ecummulative=@field;
            if ($ecummulative[4]>1000){$ecummulative[4]=1000}
        }
        else{
            if ($ecummulative[3]==$field[3]){
                $ecummulative[2]=$field[2];
                if ($ecummulative[4]<=1000-$field[4]){
                    $ecummulative[4]+=$field[4];
                }
            }
            else{
                print $eTXTg ($_,"\t") for @ecummulative;
                print $eTXTg ("\n");
                @ecummulative=@field;
            }
        }
    }
    close $eTXTg;
    print "End grouped text file creation: ",get_time(),"\n\n";
    
    print "\nTo make a BigBed file from this text file, rename file to .bed and use USCS program bedToBigBed in terminal \n\t\n";



#----------------OUTPUT REGULAR SLIDING INFORMATION WITHOUT ESSENTIALS CALCULATIONS (P-VALUES)----------------------


#MAKE OUTPUT CSV FILE WITH FITNESS WINDOW CALCULATIONS
	my $fcsv="fitWindows.csv";
    print "Start csv ouput file creation: ",get_time(),"\n";
    my $fcsvBIG = Text::CSV->new({ binary => 1, auto_diag => 1, eol => "\n"}) or die "Cannot use CSV: " . Text::CSV->error_diag();  # open in append mode
    open my $file, ">", "$outdir/$fcsv" or die "Failed to open file";
    $fcsvBIG->print($file, [ "start", "end","fitness","mutant_count" ]); #header
    foreach my $winLine(@allWindows){
        $fcsvBIG->print($file,$winLine);
    }
    close $file;
    print "End csv ouput file creation: ",get_time(),"\n\n";

#MAKE WIG FILE---->later make BW--->IGV
    print "Start wig file creation: ",get_time(),"\n";
    my $fin = Bio::SeqIO->new(-file=>$ref_genome);
    $refseq = $fin->next_seq;
    $refname = $refseq->id;
    my $fwig="fitWindows.wig";
    open fWIG, ">$outdir/$fwig";
    print fWIG "track type=wiggle_0 name=$fwig\n";
    print fWIG "variableStep chrom=$refname\n";
    foreach my $wigLine(@allWindows){
        my @wigFields=@$wigLine;
        my $position=$wigFields[0];
        #while ($position<=$wigFields[1]){
        print fWIG $position," ",$wigFields[2],"\n";
        #$position=$position+1;
        #}
        #print  WIG $wigFields[0]," ",$wigFields[2],"\n";
    }
    close fWIG;
    print "End wig file creation: ",get_time(),"\n\n";
    print "If this wig file needs to be converted to a Big Wig, then use USCS program wigToBigWig in terminal: \n \t./wigToBigWig gview/12G.wig organism.txt BigWig/output.bw \n\n";


#GOING TO MAKE A TEXT FILE FOR BED CONVERSION TO BIGBED, NEED CHROM # IN COLUMN 0
my @cummulative;
if ($txtg or $txt){
    for my $line(@allWindows){
        unshift($line, $refName);
    }
}
#MAKING A REGULAR TEXT FILE fields: [chrom,start,end,fitness,count]

    print "Start text file creation time: ",get_time(),"\n";
    open my $TXT, '>', "$outdir/fitWindows.txt" or die $!;
    print $TXT ("start","end","W","mutant_count","insertion_count\n");
    print $TXT (join("\t",@$_),"\n") for @allWindows;
    close $TXT;
    print "End text file creation: ",get_time(),"\n\n";


#MAKING A TEXT FILE OF GROUPED CONSECUTIVE WINDOWS WITH SAME FITNESS

    print "Start grouped txt file creation time: ",get_time(),"\n";
    open my $TXTg, '>', "$outdir/groupedWindows.txt" or die $!;
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

    #print "\nTo make a BigBed file from this text file, rename file to .bed and use USCS program bedToBigBed in terminal \n\t\n";
    







