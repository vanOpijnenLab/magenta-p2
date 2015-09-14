#!/usr/bin/perl -w

#Margaret Antonio 15.09.12

#../Tn_SeqAnalysisScripts/windowFitStats.pl --normal tigr4_normal.txt --dist_set 7 --dist_size 1000 --excel testerExcel2.xls results/L1_2394eVI_PennG.csv results/L3_2394eVI_PennG.csv results/L4_2394eVI_PennG.csv results/L5_2394eVI_PennG.csv results/L6_2394eVI_PennG.csv --version

use strict;
use Getopt::Long;
use warnings;
use Text::CSV;
use Text::CSV_XS;
use Bio::SeqIO;
use Data::Random qw(:all);
use List::Util qw(sum);
use Spreadsheet::WriteExcel;
#use Excel::Writer::XLSX;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
sub print_usage() {
    print "\nRequired:\n";
    print "Without a flag\t In the command line , input the name(s) of the file(s) containing fitness values for individual insertion mutants.\n";
    print "--normal \t Text file of normalization gene SP numbers (transposon genes)\n";
    print "--dist_set \t Give the number of insertions that should be in each random set of insertions\n";
    print "--dist_size \t Give the number of randomized sets that should be generated\n";
    print "\nOptional:\n";
    print "--cutoff \tCutoff: Don't include fitness scores with average counts (c1+c2)/2 < x (default: 0)\n";
    print "--excel \t Name of an excel file for output data. Currently only takes .xls extension (Default: outputs to the terminal) \n";
    print "--round \t Round all data files to x decimal places in the form of %.3f for 3 decimal places\n";
    print "--random \t Enter either 'transposon' or 'whole' depending on which file should be used for the random distribution"
    
}

#ASSIGN INPUTS TO VARIABLES
our ($normal, $excel,$set,$cutoff,$N,$round, $random, $version);
GetOptions(
'normal:s' => \$normal,
'excel:s'  => \$excel,
'dist_set:i' => \$set,
'dist_size:i' => \$N,
'cutoff:i' =>\$cutoff,
'round:i' =>\$round,
'random:s' =>\$random,
'version' =>\$version,
);

if ($version){print "\nThe transposonStat.pl script generates statistical information about transposon or whole file insertions outputted to the terminal or an excel file.\n";}
sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}
# Just to test out the script opening
print print_usage(),"\n";
print "\n";

#CHECKING PARAMETERS: Check to make sure required option inputs are there and if not then assign default

if (!$cutoff){$cutoff=0;}
if (!$round){$round='%.3f';}
print "\nFor a random distribution of x sets of y insertions in transposon genes: x=$N and y=$set"
    or die "\n A random distribution of transposon insertions requires the set size '$set' and the number of sets in the distribution '$N'\n";
print "\nExcel output file: ", $excel,"\n" or die "\n Need to specify excel ouput file";

open IN, $normal or die "No normalization file given\n";
my $transposon_genes = '';
while (<IN>) { chomp $_; $transposon_genes .= "$_ "; }
close IN;
print "Normalize genes loaded\n";   #: $transposon_genes\n";

#CREATE AN ARRAY OF DATA FROM INPUT CSV FILE(S)
my $rowCount=-1;
my $last;
my @unsorted;
my $num=$#ARGV+1;
print "\nNumber of files in csv: ", $num,"\n";

my $regInsert=0;

my $csvtemp=Text::CSV->new;
my @allFitness;

my @lonely;

for (my $i=0; $i<$num; $i++){   #Read files from ARGV
    my $csvtemp=Text::CSV->new;
    my $file=$ARGV[$i];
    open(my $data, '<', $file) or die "Could not open '$file' Make sure input .csv files are entered in the command line\n";
    $csvtemp->getline($data);
    while (my $line = $csvtemp->getline($data)) {
        chomp $line;
        my $w = $line->[12];
        if (!$w){next;} # For blanks
        else{
            my $c1 = $line->[2];
            my $c2 = $line->[3];
            my $avg = ($c1+$c2)/2;
            if ($avg < $cutoff) { next; } # Skip cutoff genes.
            else {
                my @select=($line->[0],$line->[9],$line->[11],$line->[12]);
                push (@lonely, $line->[12]);
                my $fitness=$line->[12];
                push(@allFitness,$fitness);
                my $select=\@select;
                push(@unsorted,$select);
                $rowCount+=1;
                $last=$select->[0];
                $regInsert++;
            }
        }
    }
    close $data;
    
}

my @sorted = sort { $a->[0] <=> $b->[0] } @unsorted;
@lonely = sort { $a <=> $b } @lonely;   # numerical sort
sub uniq {
my %seen;
grep !$seen{$_}++, @_;
}
@lonely = uniq(@lonely);

#foreach my $ele(@lonely){
    #print $ele, " ";
#}
print "\nTotal number of transposon insertions in the genome: $regInsert\n";

my $workbook;
my $wholeSheet;
my $TnSheet;

if ($excel){
    $workbook = Spreadsheet::WriteExcel->new($excel);
    
    $wholeSheet = $workbook->add_worksheet();
    $wholeSheet->write('A1', "Position");
    $wholeSheet->write('B1', "Gene");
    $wholeSheet->write('C1',"W");
    $wholeSheet->write('D1',"nW");
    
    $wholeSheet->write('F1',"SortednW");
    
    $wholeSheet->write_col(1,0, \@sorted); # Write a column of data
    

    $wholeSheet->write_col(4,0, \@lonely); # Write a column of data
    
    
    $TnSheet = $workbook->add_worksheet();
    $TnSheet->write('A1', "Gene");
    $TnSheet->write('B1', "Old_Fitness");
    $TnSheet->write('C1',"Weighted_Fitness");
}

print "GENE:\t\t\tOld\t\tWeighted\n";
my @transposon;
my $counter=0;
my $current="";
my $i=2;
for my $lines(@sorted){
    if ($lines->[1] ne '' and $transposon_genes =~ $lines->[1] and $lines->[2]) {
        if ($lines->[1] ne $current){
            #print "\n";
            #$current=$lines->[1];
            $TnSheet->write($i,0, $lines->[1]);
        }
        $current=$lines->[1];
        my $singleTn=(sprintf("$round",$lines->[3]));
        if ($excel){
            $TnSheet->write($i,1, sprintf("$round",$lines->[2]));
            $TnSheet->write($i,2, sprintf("$round",$lines->[3]));
            $i++;
        }
        #print $lines->[1], ":\t\t", sprintf("$round",$lines->[2]), "\t\t", sprintf("$round",$lines->[3]), "\n";
        push(@transposon,$singleTn);
        $counter++;
    }
    
}
my $allTnSheet = $workbook->add_worksheet();
$allTnSheet->write_col('A2', \@transposon); # Write a column of data

#Distribution of random sets from insertions in transposon genes; need set size and number of sets for distribution

sub mean {
    return sum(@_)/@_;
}
my $RandomTnSheet = $workbook->add_worksheet();
#print "\n The average fitnesses of $N Random sets of $set insertions:\n";
my @distribution;
for (my $i=0; $i<$N; $i++){
    my @random_set = rand_set( set => \@allFitness, size => $set );
    my $average=sprintf("$round",mean(@random_set));
    push (@distribution, $average);
    #print "$average\n";
}
@distribution = sort { $a <=> $b } @distribution;   # numerical sort
$RandomTnSheet->write('A1', "$N random sets of $set insertions that occured in transposon genes"); # Write a column of data
$RandomTnSheet->write_col('A2', \@distribution); # Write a column of data


#print "\n Total number of insertions in transposons: $counter\n";

my $seqio = Bio::SeqIO->new(-file => "fake.txt", '-format' => 'Fasta');
my $prev;
my $total=0;
my $genome;
while(my $seq = $seqio->next_seq) {
    $genome = $seq->seq;
}
#my $x="TA";
#my @c = $genome =~ /$x/g;
#my $countTA = @c;
my $countTA=0;

my $x = 'TA';
while ($genome =~ /TA/g) { $countTA++ }
print "There are $countTA negative numbers in the string";

       
