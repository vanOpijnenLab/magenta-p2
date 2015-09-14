#!/usr/bin/perl -w

#Margaret Antonio 15.04.15

#essentials.pl: An attempt at calculating window regions with underrepresented number of reads
#essentials update: now adding p-value to each window. For each window of x possible TA sites, generate 10,000 sets of x TA sites and then get ratio for (insertions at those TA sites) /(TA sites).
#essentials update: makes one null distribution "library" of random 10,000 sites (instead of remaking it every time) and uses it for all statistical testing. Much faster.

#../Tn_SeqAnalysisScripts/essentials.pl --excel essentialTest/kill.xls --ref=NC_003028b2.gbk --essential tigr4_genome.fasta --csv essentialTest/1Bessentials.csv results/L1_2394eVI_PennG.csv results/L3_2394eVI_PennG.csv results/L4_2394eVI_PennG.csv results/L5_2394eVI_PennG.csv results/L6_2394eVI_PennG.csv

#../Tn_SeqAnalysisScripts/essentials.pl --ref=NC_003028b2.gbk  --excel essentialTest/kill.xls  --csv essentialTest/1Bessential.csv results/L1_2394eVI_PennG.csv results/L3_2394eVI_PennG.csv

#results/L4_2394eVI_PennG.csv results/L5_2394eVI_PennG.csv results/L6_2394eVI_PennG.csv


use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use Text::CSV_XS;
use Bio::SeqIO;
use Data::Random qw(:all);
use List::Util qw(sum);
use Spreadsheet::WriteExcel;
use List::BinarySearch qw( :all );
use List::BinarySearch::XS;
use List::MoreUtils qw(uniq);

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
sub print_usage() {
    print "\nRequired:\n";
    print "In the command line (without a flag), input the name(s) of the file(s) containing fitness values for individual insertion mutants.\n";
    
    print "\nOptional:\n";
    print "--size \t The size of the sliding window(default=500) \n";
    print "--step \t The window spacing (default=10) \n";
    print "--csv \t Name of a file to enter the .csv output for sliding windows.\n";
    print "--cutoff \tCutoff: Don't include fitness scores with average counts (c1+c2)/2 < x (default: 0)\n";
    
    print "\nMust choose at least one type of output:\n";
    print "--wig\tCreate a wiggle file for viewing in a genome browser. Provide a filename. Also provide genome under --ref\n";
    print "--txt\t Output all data [start,end,W,count] into a text of bed file.\n";
    print "--txtg\t If consecutive windows have the same value, then group them into one window. Ouput into txt file or bed file.\n";
    print "--ref\tThe name of the reference genome file, in GenBank format. Needed for wig and txt file creation\n";
    print "--excel \t Name of an excel file for output data. Currently only takes .xls extension (Default: outputs to the terminal) \n";
}


#ASSIGN INPUTS TO VARIABLES
our ($round,$random,$txt,$txtg,$cutoff,$wig,$ref_genome,$infile, $eWig, $csv, $step, $size,$genome, $excel);
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
'essential:s'=> \$genome,
'random:s' =>\$random,
'round:i' =>\$round,
'excel:s'  => \$excel,
'essentialWig:s'=>\$eWig,

);

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}
# Just to test out the script opening
print print_usage(),"\n";
print "\n";
if ($csv){print "CSV output file: ", $csv,"\n";}
if ($txt){print "Text file for window data: $txt\n";}
if ($txtg){print "Text file for grouped windows: $txtg\n";}
if (!$round){$round='%.3f';}
if (!$cutoff){$cutoff=15;}

#CHECKING PARAMETERS: Check to make sure required option inputs are there and if not then assign default
if (!$size) { $size=500 };   #set the default sliding window size to 500
if (!$step) { $step=10 };   #set the default step size to 10
if (!$cutoff) {$cutoff=0};
print "Window size: $size\n";
print "Step value: $step\n";
print "Cutoff: $cutoff\n";

#if ((!$csv) and (!$txt) and (!$txtg) and (!$wig)) {
#print "\nThere must be some kind of output file assigned (--csv, --txt, or --txtg)\n";
#print_usage();
#die;
#}

#CREATE AN ARRAY OF DATA FROM INPUT CSV FILE(S)
print "\nStart input array ",get_time(),"\n";

my $rowCount=-1;
my $last=0;
my @unsorted;
my @insertPos; #array to hold all positions of insertions. Going to use this later to match up with TA sites
my $num=$#ARGV+1;
print "\nNumber of files in csv: ", $num,"\n";


#Go through each file from the commandline (ARGV array) and read each line as an array into select array if values satisfy the cutoff

for (my $i=0; $i<$num; $i++){   #Read files from ARGV
    my $csvtemp=Text::CSV->new;
    my $file=$ARGV[$i];
    open(my $data, '<', $file) or die "Could not open '$file' Make sure input .csv files are entered in the command line\n";
    $csvtemp->getline($data);
    print"FileNum $i ------------------------\n";
    while (my $line = $csvtemp->getline($data)) {
        chomp $line;
        my $w = $line->[12];
        if (!$w){next;} # For blanks
        else{
            my $c1 = $line->[2];
            my $c2 = $line->[3];
            my $avg = ($c1+$c2)/2;
            if ($avg > $cutoff) {
                my @select=($line->[0],$line->[12]);
                my $select=\@select;
                push(@unsorted,$select);
                push(@insertPos,$line->[0]);   #keep track of actual insertion site position
                $last=$select->[0];
                $rowCount++;
                
            }
        }
    }
    close $data;
}

my @sorted = sort { $a->[0] <=> $b->[0] } @unsorted;

@insertPos = sort { $a <=> $b } @insertPos;
@insertPos= uniq @insertPos;

for (my $i=0;$i<30;$i++){
    foreach my $element ( @{ $sorted[$i] }){
        print $element,"\t";
    }
    print "\n";
}

print "Finished input array ",get_time(),"\n";

###############################################################################################

print "\n---------This is the sliding window fitness calculation part--------\n\n";

my $index=-1;
my $marker=0;


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
                #print @Wwindow, "\n";
                return (\@window);
            }
            else{
                return -1;
            }  #Because count=0 (i.e. there were no insertion mutants in that window)
                
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
}
print "End calculation: ",get_time(),"\n\n";

#-------ESSENTIALS:Counting the number of TA sites in the genome and whether an insertion occurred there or not

#if ($genome){

    my @sites;

    #First read fasta file into a string
    my $seqio = Bio::SeqIO->new(-file => "tigr4_genome.fasta", '-format' => 'Fasta');
    my $prev;
    my $total=0;
    while(my $seq = $seqio->next_seq) {
        $genome = $seq->seq;
    }
    #Get number of "TA" sites, regardless of insertion---this is just the fasta file
    my $x="TA";
    my @c = $genome =~ /$x/g;
    my $countTA = @c;
    
    #At what positions in the genome do TA sites occur?


    my $pos=0;
    my $countInsert=0;
    my @selector;   #use this array to hold all 1 and 0 of whether insertion happened or not.
        #for my $some(@insertPos){print $some," ";}

    
    #Go through all TA sites identified above and see if an insertion occurs there.
    #Push results onto two arrays a 2d array with the position and a 1D array with just the 0 or 1
   
    my @unmatched; #hold all unmatched ta sites
    my @allTAsites; #2d array to hold all occurences of TA sites in genome
    my $unmatchedCount=0;
    my $offset=0;
    my $result = index($genome, 'TA',$offset);

    while (($result != -1) and ($pos!=scalar @insertPos)) {
        my $yes=0;
        if ($result>$insertPos[$pos]){
            push @unmatched,$insertPos[$pos];
            $unmatchedCount++;
            $pos++;
        }
        if ($result==$insertPos[$pos]){
            $yes=1;
            $countInsert++;
            $pos++;
        }

        #print "$result \t $insert[$posit] \t $yes\n";
        #print $insert[$index],"     ";
        my @sites=($result,$yes);
        push @selector,$yes;    #push the 0 or 1 onto the array @selector---going to use this to draw random sets for the null distribution
        push (@allTAsites,\@sites);
        #print $result, " ", $yes,"\n";
        $offset = $result + 1;
        $result = index($genome, 'TA', $offset);
        $countTA++;
        
    }
    #  foreach my $pos(@unmatched){
    #    print $pos, "\n";}
    
    print "Total of unmatched insertions: $unmatchedCount\n";
    #PRINT 2D ARRAY OF TA SITES AND WHETHER THEY HAVE INSERTIONS ONTO EXCEL SHEET
    my $workbook;
    my $wholeSheet;
    $workbook = Spreadsheet::WriteExcel->new($excel);
    
    $wholeSheet = $workbook->add_worksheet();
    $wholeSheet->write('A1', "TA_SITE");
    $wholeSheet->write('B1', "INSERTION?");
    
    $wholeSheet->write_col(1,0, \@allTAsites); # Write a column of data
    
    #PRINT 2D ARRAY TO WIG FILE SO WE CAN SEE IT IN IGV
    
        print "Start wig file creation: ",get_time(),"\n";
        my $in = Bio::SeqIO->new(-file=>$ref_genome);
        my $refseq = $in->next_seq;
        my $refname = $refseq->id;
        open WIG, ">$eWig";
        print WIG "track type=wiggle_0 name=$eWig\n";
        print WIG "variableStep chrom=$refname\n";
        foreach my $entry(@unmatched){
            #my @wigFields=@$wigLine;
            print WIG $entry," ",2,"\n";
        }
        close WIG;
        print "End wig file creation: ",get_time(),"\n\n";
        print "If this wig file needs to be converted to a Big Wig, then use USCS program wigToBigWig in terminal: \n \t./wigToBigWig gview/12G.wig organism.txt BigWig/output.bw \n\n";
 

#print "\nTotal: $countInsert insertions in $countTA TA sites.\n";

#------------------------------------------------------------------------------------------------------------------------------------------------------
    
    #Now, have an array for each TA site and if an insertion occurred there. So per site @sites(position, 0 or 1 for insertion).
    #Next step, create null distribution of 10,000 random sets with same number of TA sites as the window and then calculate p-value
  
    #SUBROUTINE FOR MAKING THE NULL DISTRIBUTION SPECIFIC TO THE WINDOW
    sub mean {
        return sum(@_)/@_;
    }


#MAKE LIBRARY OF NULL DISTRIBUTIONS:
print "About to make a library of null distributions\n";

my @distLib;

my $N=10000;


for (my $sitez=1; $sitez<50;$sitez++){
    
    #print "In the first for loop to make a null distribution\n";
    my @unsorted;
    my $count=0;
    my $sum=0;
    #my $average=0;
    
    for (my $i=0; $i<50; $i++){
        #print "In the second for loop to make a null distribution\n";
        my @random_set = rand_set( set => \@selector, size => $sitez);
        my $setAvg=sprintf("$round",mean(@random_set));
        push (@unsorted, $setAvg);
        #print "$i:\t", "$setAvg\n";
        $count++; $sum+=$setAvg;
    }
    my @nullDist= sort { $a <=> $b } @unsorted;
    
    my $min=$nullDist[0];
    my $max=$nullDist[scalar @nullDist-1];
    my $nullAvg=$sum/$count;
    #print "NULL DIST Min: $min | Max: $max | Average: $nullAvg\n\n";
    #print "Done with null distrubtion #$sitez\n";
    push (@distLib,\@nullDist);
    
}
    
    
    #SUBROUTINE TO CALCULATE THE P-VALUE OF THE WINDOW INSERTIONS AGAINST THE NULL DISTRIBUTION

    sub pvalue{
        
        #takes in window count average (countAvg) and number of TAsites and makes a null distribution to calculate the pvalue, which it returns
        my $countAvg=shift@_;
        my $TAsites=shift@_;
        my $N=10000;
        #my @nullDist=random($TAsites,$N);
        #my compare=compare(\@nullDist,$countAvg);
        #my $rank=ranker(\&compare(\@nullDist,$countAvg), $countAvg,\@nullDist,\@nullDist,$countAvg)+1;      #rank window average among null averages
        my $rank= binsearch_pos { $a cmp $b } $countAvg, @{$distLib[$TAsites+1]};
        #print "rank=$rank\n";
        my $pval=$rank/$N; #calculate pval as rank/N
        # print "$pval\t";
        return $pval;
        
    }
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    
    print "Size of genome is: ", length($genome), "\n";
    
    #Now we have an array called @allTAsites which contains every TAsite position with a 0 next to it for "no insertion".
    #Now just need to replace 0 with 1 if there IS and insertion at that site

    #FOR TESTING: to print subset of the genome
        #my $sub=substr($genome,3780,10);
        #print "\nFrom positions 3780 to 3790: \n  $sub  ";
    
    #Just counting number of TA sites in the window
    
    my @newWindows=();
#my $countAvg;
    my $printNum=0;
#print $printNum,"\t";
print "Creating null distribution and calculating pvalue for windows:\n";
#print "This is the array of all Windows\n";

print "Done testing if array of allWindows prints\n";


print "Start p-value calculation: ",get_time(),"\n\n";

#my $allWindows=\@allWindows;
for (my $i=0;$i<5000;$i++){
    my @win=@{$allWindows[$i]};
    my $starter=$win[0];
    my $ender=$win[1];
    #print "num $printNum -->\tStart pos: $starter\tEnd pos: $ender\n";
    #How many TA sites are there from $genome[$start] to $genome[$end]?


    my $seq = substr($genome,$starter-1,500);  #start-1 becase $start and $end are positions in genome starting at 1,2,3.... substr(string,start, length) needs indexes
    my $ta="TA";
    my @c = $seq =~ /$ta/g;
    my $TAsites = scalar @c;
    push(@win,$TAsites);
    my $countAvg=$win[3]/$TAsites;
    #print "\tCountAVG=$countAvg\n";
    push (@win,$countAvg);
    my $pval=pvalue($countAvg,$TAsites);
    push (@win,$pval);
    push (@newWindows,\@win);
    
    $printNum++;
    }
@allWindows=@newWindows;

   
 print "End p-value calculation: ",get_time(),"\n\n";

#print "This is the TA count: $count\nTotal genome size is: $total\n\n";

#---------------------------------------------------------OUTPUTS-------------------------------------------------------------------


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

#MAKE OUTPUT CSV FILE WITH WINDOW CALCULATIONS
elsif ($csv){
    print "Start csv ouput file creation: ",get_time(),"\n";
    my $csvBIG = Text::CSV->new({ binary => 1, auto_diag => 1, eol => "\n"}) or die "Cannot use CSV: " . Text::CSV->error_diag();  # open in append mode
    open my $file, ">", "$csv" or die "Failed to open file";
    $csvBIG->print($file, [ "start", "end","fitness","mutant_count","insertions","TA_sites","ratio","p-value"]); #header
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
        print WIG $position," ",$wigFields[7],"\n";    #7 for pvalue, but 2 for fitness!!!!!!
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

