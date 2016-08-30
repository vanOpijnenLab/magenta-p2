#!/usr/bin/perl -w

#Margaret Antonio 16.08.29

use strict;
use Getopt::Long;
use Set::Scalar;
use Text::CSV;
use Bio::SeqIO;
use Data::Random qw(:all);
use List::Util qw(sum);
use List::BinarySearch qw( :all );
use List::BinarySearch::XS;
use List::MoreUtils qw(uniq);
use File::Path;
use File::Basename;
use autodie;
use List::MoreUtils qw( minmax );
no warnings;
#use warnings;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
sub print_usage() {

    print "\n###############################################################\n";
    print "\nperl dataOverview.pl -i inputs/ -f genome.fasta -r genome.gbk\n";
        
    print "\nREQUIRED:\n";
    print "-i\tDirectory containing all input files (results files from \n\tcalc fitness script)\n";
    print "\t  OR\n";
    print "\tIn the command line (without a flag), input the name(s) of \n\tthe files containing fitness values for individual \n\tinsertion mutants\n";
    print "-f\tFasta file for genome\n";
    print "-r\tThe name of the reference genome file, in GenBank format\n";

    print "\nOPTIONAL:\n";
    print "-l\tSend all output to a log file instead of the terminal\n";
    print "-h\tPrint usage\n";
    print "-c\tCutoff average(c1+c2)>c. Default: 15\n";
    print "-o\tFilename for output\n";
    print "\n~~~~Always check that file paths are correctly specified~~~~\n";
    print "\n###############################################################\n";

}

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
    }
sub mean {
    return sum(@_)/@_;
}
#ASSIGN INPUTS TO VARIABLES
our ($cutoff,$fastaFile, $outfile,$help,$ref,$indir,$weight_ceiling,$log);
GetOptions(
'r:s' => \$ref,
'f:s' => \$fastaFile,
'i:s'=>\$indir,
'c:i'=>\$cutoff,
'o:s' => \$outfile,
'h'=> \$help,
'w:i' => \$weight_ceiling,
'l' => \$log,
);

# Set defaults
if (!$weight_ceiling){$weight_ceiling=999999;}
if (!$cutoff){$cutoff=15;}

# If help option is specified or required files are not specified:

if ($help) {
    print_usage();
	print "\n";
	exit;
}
if (!$indir and (scalar @ARGV==0)){
	print "\nERROR: Please correctly specify input files or directory\n";
    print_usage();
	print "\n";
	exit;
}
if (!$fastaFile or !$ref){
	print "\nERROR: Please correctly specify reference genome fasta and genbank files\n";
	print "Most genomes (in fasta and gbk format) are available at NCBI\n";
    print_usage();
	print "\n";
	exit;
}

# Redirect STDOUT to log.txt. Anything printed to the terminal will go into the log file
if ($log){
	print "\nSending all output to log file\n";
    open (STDOUT, ">>$outfile");
}

#Not sure if I'll need this but sometimes funky data inputs have hidden characters
sub cleaner{
    my $line=$_[0];
    chomp($line);
    $line =~ s/\x0d{0,1}\x0a{0,1}\Z//s;
    return $line;
}

#Get the input files out of the input directory, or take off of command line
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

print "Gathering data overview for Tn-Seq experiment\n\n";
print "Begin time: ",get_time(),"\n\n";

#CREATE AN ARRAY OF DATA FROM INPUT CSV FILE(S). 
#These are the "results" files from calc_fitness.pl. Insertion location, fitness, etc.
#Go through each file from the commandline (ARGV array) and read each line as an array
#into select array if values satisfy the cutoff


#Store ALL insertion locations in this array. Later, get unique insertions
my @insertions_all;
#Store all genes with valid insertions here
my @genes_insertions;
#all lines that satisfied cutoff
my @unsorted;
#array to hold all positions of insertions. Going to use this later to match up with TA sites
my @insertPos;

#Markers
my $rows=-1;
my $last=0;

print "Library description\n\n";
my @header=("library","file_path","ins","ins.f","genes.ins");
print join ("\t",@header),"\n";

for (my $i=0; $i<$num; $i++){
	#Temp arrays for library
	my(@insertions_all_lib,@genes_insertions_lib,@insertPos_lib);
    my $file=$files[$i];
    open(DATA, '<', $file) or die "Could not open '$file' Make sure input .csv files are entered in the command line\n";
    my $dummy=<DATA>; #read and store column names in dummy variable
    while (my $entry = <DATA>) {
    	chomp $entry;
		my @line=split(",",$entry);
        my $locus = $line[9]; #gene id (SP_0000)
        my $w = $line[12]; #nW
        if (!$w){ $w=0 }   # For blanks
        my $c1 = $line[2];
        my $c2 = $line[3];
        my $coord= $line[0];
        push (@insertions_all_lib,$coord);
         #Average counts must be greater than cutoff (minimum allowed)
        my $avg = ($c1+$c2)/2;
        if ($avg > $cutoff) {
        	my @select=($coord,$w,$avg,$locus);
            my $select=\@select;
            push(@unsorted,$select);
            push(@insertPos_lib,$line[0]);   #keep track of actual insertion site position
            push (@genes_insertions_lib,$locus);
            $last=$select[0];
            $rows++;
        }
        if ($avg >= $weight_ceiling) { $avg = $weight_ceiling } # Maximum weight
    }
    close DATA;
    push (@insertions_all,@insertions_all_lib);
    @genes_insertions_lib= uniq @genes_insertions_lib;
    push (@genes_insertions,@genes_insertions_lib);
    push (@insertPos,@insertPos_lib);
    my @stat=($i+1,$file,scalar @insertions_all_lib,scalar @insertPos_lib,scalar @genes_insertions_lib);
    print join("\t",@stat),"\n";
}

@insertPos = sort { $a <=> $b } @insertPos;
@insertPos= uniq @insertPos;
@genes_insertions= uniq @genes_insertions;
@insertions_all=uniq @insertions_all;
my $totalAll=scalar @insertions_all;
my $total=scalar @insertPos;
my $temp="1-".$num;
my @all_stat=($temp,"NA",$totalAll,$total,scalar @genes_insertions);
print join("\t",@all_stat),"\n";

#Genome description: #TA sites, distance between TA sites, #TA sites in ORFS
print "\n-------------------------\n";
print "\nGenome description\n\n";
print "File for genome: ", $fastaFile,"\n";

my @sites;
#First read fasta file into a string
my $seqio = Bio::SeqIO->new(-file => $fastaFile, '-format' => 'Fasta');
my $fasta;
while(my $seq = $seqio->next_seq) {
	$fasta = $seq->seq;
}
#Just in case $fasta file is in lowercase, change it to uppercase
$fasta=uc $fasta;

#Get genomic coordinate for TA sites:
my $x="TA";
my $offset=0;
my @indices;
my $result=index($fasta,$x,$offset);
while ($result !=-1){
	push (@indices,$result);
	$offset=$result+1;
	$result=index($fasta,$x,$offset);
}
my $countTA=scalar @indices;

#Get longest stretch with no TA sites
my @tempta=@indices;
my $prev=shift @tempta;
my $current=shift @tempta;
my $lg_dist_ta=$current-$prev;
foreach my $site(@tempta){
	$prev=$current;
	$current=$site;
	my $d=$current-$prev;
	if ($d>$lg_dist_ta){
		$lg_dist_ta=$d;
	}
}

#Get longest stretch of with no insertions
my @tempins=@insertPos;
$prev=shift @tempins;
$current=shift @tempins;
my $lg_dist_ins=$current-$prev;
foreach my $site(@tempins){
	$prev=$current;
	$current=$site;
	my $d=$current-$prev;
	if ($d>$lg_dist_ins){
		$lg_dist_ins=$d;
	}
}


my $genSize=length $fasta;
print "$genSize\tGenome size\n";
print "$countTA\tTotal number of TA sites\n\n";

my $sat=sprintf("%.2f", ($total/$countTA)*100);
my $satAll=sprintf("%.2f", ($totalAll/$countTA)*100);
my $inscov=sprintf("%.2f", ($total/$genSize)*100);
my $tacov=sprintf("%.2f", ($countTA/$genSize)*100);

#Get GC content of genome

my $sequence = ' ';
my $Ccount = 0;
my $Gcount = 0;
my $identifier = ' ';

my @nucleotides = split('', $fasta);

foreach my $nuc (@nucleotides) {
	if ($nuc eq 'G') {$Gcount++} 
	elsif ($nuc eq 'C') {$Ccount++}
}
my $sequencelength=length $fasta;

my $GCcontent = sprintf("%.2f",((($Gcount + $Ccount) / $sequencelength) * 100));
my $ATcontent =100-$GCcontent;

print "$GCcontent%\tGC content of this genome\n";
print "$ATcontent%\tAT content of this genome\n";

print "$satAll%\tSaturation of TA sites before cutoff filter (allInsertions/TAsites)\n";
print "$sat%\tSaturation of TA sites after cutoff filter (validInsertions/TAsites)\n";
print "$inscov%\tGenome coverage by insertions (validInsertions/genomeSize)\n";
print "$tacov%\tGenome coverage by TA sites (TAsites/genomeSize)\n";
print "$lg_dist_ta\tLargest distance between TA sites\n";
print "$lg_dist_ins\tLargest distance between insertions\n";
print "\n\nOpen Reading Frames\n\n";

#Store everything to be printed in array
my @table;

#Find open reading frames from fasta file
local $_  = $fasta;
my @orfSize;
my @allc; #numbers of TAs in the ORFS here.
my $blank=0; #ORFS that don't have any TA sites.
my $orfCount=0; #keep track of the number of ORFs found.
my $minSize=0; 
#Read somewhere that 99 is a good min but there is an annotated 86 bp gene for 19F
while ( /ATG/g ) {
   my $start = pos() - 3;
   if ( /T(?:AA|AG|GA)/g ) {
     my $stop = pos;
     my $size=$stop - $start;
     if ($size>=$minSize){
		 push (@orfSize,$size);
		 my $seq=substr ($_, $start, $stop - $start); 
		 my @ctemp = $seq =~ /$x/g;
		 my $countTA = @ctemp;
		 if ($countTA==0){$blank++}
		 push (@allc,$countTA);  
		 $orfCount++;  
	   }
	}
}

print "\nORFs based on Fasta sequence and start (ATG) and end (TAA,TAG,TGA) codons\n";
push (@table,["Set minimum size for an ORF",$minSize]);
print "$orfCount\tTotal number of ORFs found\n";
my ($minORF, $maxORF) = minmax @orfSize;
print "$minORF\tSmallest ORF\n";
print "$maxORF\tLargest ORF\n";
my ($mintaORF, $maxtaORF) = minmax @allc;
print "$mintaORF\tFewest # TA sites in an ORF\n";
print "$maxtaORF\tGreatest # TA sites in an ORF\n";
print "$blank\tNumber of ORFs that don't have any TA sites\n";


print "\nGenes using the genbank annotation file\n\n";
###Get genbank file. Find all start and stop for genes
#See how many insertions fall into genes vs intergenic regions
#Get array of coordinates for all insertions then remove insertion if it is
#within a gene region
my $gb = Bio::SeqIO->new(-file=>$ref);
my $refseq = $gb->next_seq;

#store number of insertions in a gene here
my @geneIns;
my @allLengths;
my $blankGene=0; #Number of genes that don't have any insertions in them
my @genomeSeq=split('',$fasta);


#keep a copy to remove insertions that are in genes
my @insertPosCopy=@insertPos;

my @features = $refseq->get_SeqFeatures(); # just top level
foreach my $feature ( @features ) {
	if ($feature->primary_tag eq "gene"){
		my $start=$feature->start;
		my $end=$feature->end;
		my $length=$end-$start;
		push (@allLengths,$length);
		#turn this into a for loop
		my $i=0;
		my $ins=0;
		my $current=$insertPos[$i];;
		while ($current<=$end && $i!=scalar @insertPos){
			if ($current>=$start){
				splice(@insertPosCopy, $i, 1);
				$ins++;			
			}
			$current=$insertPos[$i++];
		}
		if ($ins==0){$blankGene++}
		push (@geneIns,$ins);
	}
}
my $avgLength=sprintf("%.2f",mean(@allLengths));
my ($minLength, $maxLength) = minmax @allLengths;
my $avgInsGene=sprintf("%.2f",mean(@geneIns));
my ($minInsGene, $maxInsGene) = minmax @geneIns;
my $nonGeneIns=scalar @insertPosCopy;
my $totalIns=scalar @insertPos;
my $percNon=sprintf("%.2f",($nonGeneIns/$totalIns)*100);
print "Length of a gene\n";
print "$avgLength\tAverage","\n$minLength\tMininum","\n$maxLength\tMaximum\n";
print "\nFor insertions in a gene:\n";
print "$avgInsGene\tAverage","\n$minInsGene\tMininum","\n$maxInsGene\tMaximum\n";
print "Number of genes that do not have any insertions: ",$blankGene,"\n";
print "\n$nonGeneIns\tInsertions that are not in genes","\n$percNon% of all insertions\n";
#How many insertions are in genes and how many are in non-gene regions?


print "End time: ",get_time();

print "\n\n-----THE END-------\n\n";
