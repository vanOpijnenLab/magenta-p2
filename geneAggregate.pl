#!/usr/bin/perl -w

#This script (aggregate_mla.pl) is a modified version of aggregate.pl. It includes the following changes:
#1 Differentiated between boolean and argument flags. So now random argument isn't needed for -w
#2 Since we always run the following parameters for flags, made these the defaults
    # -w 1 -x 10 -l 50 -b 0
#3 Added directory input so individual files don't have to be specified in the command line


use strict;
use Getopt::Std;
use Bio::SeqIO;
use Set::Scalar;

use Getopt::Long;
use Set::Scalar;
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
use Data::Dumper;


#perl ../Blueberries/aggregate_MLA.pl -d data/ -o testAgg_160509.csv -w -l 50 -x 10 -r ../0-genome/NC_003028b2.gbk -f ../0-genome/NC_003028b2.fasta -n ../9-daptomycin/nullDist/nullDist_3000t4dapto.csv

#perl ../Blueberries/aggregate_MLA.pl -d ../9-daptomycin/data/T4Dapto/results/ -o testAgg_160508.csv -w -l 50 -x 10 -r ../0-genome/NC_003028b2.gbk -f ../0-genome/NC_003028b2.fasta -n ../9-daptomycin/nullDist/nullDist_3000t4dapto.csv
sub print_usage() {
    print "\n";
    print "Usage: ./aggregate.pl -o outfile -d indirectory -r genbankFile -f fastaFile\n\n";
    print "Description: Finds average fitness of insertion mutations within annotated genes\n\n";
    print "Option List:\n\n";
    print " -o\tOutput file for aggregated data. (Required)\n";
    print " -d\tDirectory containing input files. Make sure / is included after name\n";
    print " -r\tCheck for missing genes in the data set - provide a reference genome in\n";
    print " \tgenbank format. Missing genes will be sent to stdout.\n";
    print " -m\tPlace a mark in an extra column for this set of genes. Provide a file\n";
    print " \twith a list of genes seperated by newlines.\n";
    print " -x\tCutoff: Don't include fitness scores with average counts (c1+c2)/2 < x (default: 0)\n";
    print " -b\tBlanks: Exclude -b % of blank fitness scores (scores where c2 = 0) (default: 0 = 0%)\n";
    print " -w\tUse weighted algorithm to calculate averages, variance, sd, se\n";
    print " -l\tWeight ceiling: maximum value to use as a weight (default: 999,999)\n";
    print "-g\tGene list with coordinates so essentiality (pvalue) can be calculated\n";
    print "-f\tFasta file\n";
    print "-r\tRound to this number of decimals\n";
    print "-n\tNull distribution";
    print "\n";
}

#Global options
my %options;
getopts('n:u:r:f:g:s:m:x:wl:b:o:d:a:n:',\%options);
my $out = $options{o};
my $ref = $options{r};
my $marked = $options{m};
my $cutoff = $options{x} || 10;
my $blank_pc = $options{b} || 0;
my $weight_ceiling = $options{l} || 50;
my $indir=$options{d};
my $coords=$options{g};
my $fastaFile=$options{f};
my $round=$options{s} || '%.3f';
my $maxLength=$options{a};
my $library=$options{n};


print "\nOutfile: $options{o}\n";
# other things found on the command line
if (!$options{o}) { &print_usage; exit; }

#Read in the fasta file so it's stored as a string
my $seqio = Bio::SeqIO->new(-file => $fastaFile, '-format' => 'Fasta');
my $fasta;
while(my $seq = $seqio->next_seq) {
	$fasta = $seq->seq;
}

#ADDED BY MLA allows input to be directory---good for inputting L1-L6
my @files;

if ($indir){
    my $directory=$indir;
    print "Input directory: ", $directory,"\n";

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


#If calculating essentiality (pvalue) then read list of gene coordinates into a hash

if ($coords){
    my %geneCoords;
    open (COORDS,'<',$coords);
    my $header=<COORDS>;
    while (my $line = <COORDS>) {
        chomp $line;
        my @line=split(",",$line);
        my $gene_id=$line[5];
        my $start=$line[1];
        my $end=$line[2];
        my @value=($start,$end);
        $geneCoords{$gene_id}=\@value;
    }

}


# Returns mean, variance, sd, se
sub average {
    my $scoreref = shift @_;
    my @scores = @{$scoreref};
    my ($variance,$sd,$se);

    my $sum=0;
    my $num=0;

    # Get the average.
    foreach my $w (@scores) {
        $sum += $w;
        $num++;
    }
    my $average= $sum/$num;
    my $xminusxbars = 0;

    # Get the variance.
    foreach my $w (@scores) {
        $xminusxbars += ($w-$average)**2;
    }
    if ($num>1){
   		$variance = (1/($num-1)) * $xminusxbars;
   		$sd = sqrt($variance);
   		$se = $sd / sqrt($num);
   	 }
   	 else{
   	 	$variance=0;
   	 	$sd=0;
   	 	$se=0;
   	 	}

    return ($average, $variance, $sd, $se);

}

# Takes two parameters, both hashrefs to lists.
# 1) hashref to list of scores
# 2) hashref to list of weights, equal in length to the scores.
sub weighted_average {

    my $scoreref = shift @_;
    my $weightref = shift @_;
    my @scores = @{$scoreref};
    my @weights = @{$weightref};

    my $sum=0;
    my ($weighted_average, $weighted_variance)=(0,0);
    my $v2;

    # Go through once to get total, calculate V2.
    for (my $i=0; $i<@weights; $i++) {
       $sum += $weights[$i];
       $v2 += $weights[$i]**2;
    }
    if ($sum == 0) { return 0; } # No scores given?

    my $scor = join (' ', @scores);
    my $wght = join (' ', @weights);
    #print "Scores: $scor\n";
    #print "Weights: $wght\n";

    # Now calculated weighted average.
    my ($top, $bottom) = (0,0);
    for (my $i=0; $i<@weights; $i++) {
        $top += $weights[$i] * $scores[$i];
        $bottom += $weights[$i];
    }
    $weighted_average = $top/$bottom;
    #print "WA: $weighted_average\n";

    ($top, $bottom) = (0,0);
    # Now calculate weighted sample variance.
    for (my $i=0; $i<@weights; $i++) {
       $top += ( $weights[$i] * ($scores[$i] - $weighted_average)**2);
       $bottom += $weights[$i];
    }
    $weighted_variance = $top/$bottom;
    #print "WV: $weighted_variance\n";

    my $weighted_stdev = sqrt($weighted_variance);
    my $weighted_stder = $weighted_stdev / sqrt(@scores);  # / length scores.

    #print "$weighted_average, $weighted_variance, $weighted_stdev\n";
    return ($weighted_average, $weighted_variance, $weighted_stdev, $weighted_stder);
}

my $N=10000;

sub mean {
    if (scalar @_ ==0){
        return 0;
    }
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


my @distLib;

my $FILE3 = "nullDist.txt";
unless(open DIST, ">", $FILE3){
    die "\nUnable to create $FILE3:\n$!";
}
printf DIST "Sites\tSize(n)\tMin\tMax\tMean\tstdev\tvar\tMinPval\n";

#Loop once for each distribution in the library
#(i.e. distribution of 10,000 sites each with 35 TA sites, then a distribution of 10,000 sites each with 36 TA sites, etc)

my %distLibrary;
my @selector;
my @nullDist;
my @unsorted;
my $sum;
my $count;

if ($library){
	#If a library was already made  using the MakeNullDist tool then read in result
	open LIB,'<',$library;
	while (<LIB>){
		chomp $_;
		my @dist=split (",",$_);
		my $key=$dist[0];
		shift @dist;
		$distLibrary{$key}=\@dist;
		#print $key,"\t",scalar @dist,"\n";
	}
	close LIB;
	print "size\tN\tmin\tmax\tnullMean\tstdev\tvar\tminp\n";
	#Print stats on each array in the distribution library
	foreach (sort keys %distLibrary){
	    my $sites=$_;
		@unsorted=@{$distLibrary{$_}};
		@nullDist= sort { $a <=> $b } @unsorted;	
		my $size=scalar @nullDist;
		#my ($nullMean, $var, $stdev, $sterr)=average(@nullDist);
		my $counter=scalar @nullDist;
		my $nullMean=sprintf("$round",mean(@nullDist));
		my $stdev =sprintf("$round",stdev(\@nullDist, $nullMean));
		my $var=sprintf("$round",($stdev*$stdev));
		my $min=sprintf("$round",$nullDist[0]);
		my $max=sprintf("$round",$nullDist[scalar @nullDist-1]);
		my $minp=pvalue(0,$size);
		print "$sites\t$N\t$min\t$max\t$nullMean\t$stdev\t$var\t$minp\n";

}
}


else{
	for (my $sitez=1; $sitez<=$maxLength;$sitez++){
		#print "In the first for loop to make a null distribution\n";
	    @unsorted=();
		$count=0;
		$sum=0;
	
		for (my $i=1; $i<=$N; $i++){
				my @random_set = rand_set( set => \@selector, size => $sitez);
				my $setMean=mean(@random_set);
				push (@unsorted, $setMean);
			}
		@nullDist= sort { $a <=> $b } @unsorted;	
		$distLibrary{$sitez}=\@nullDist;
		my ($nullMean, $nullVar, $nullsd, $nullse)=average(@nullDist);
		my $counter=scalar @nullDist;
		$nullMean=sprintf("$round",($sum/$count));
		my $stdev =sprintf("$round",stdev(\@nullDist, $nullMean));
		my $variance=sprintf("$round",($stdev*$stdev));
		my $min=sprintf("$round",$nullDist[0]);
		my $max=sprintf("$round",$nullDist[scalar @nullDist-1]);
		my $minp=pvalue(0,$sitez);
		printf DIST "$sitez\t$N\t$min\t$max\t$nullMean\t$stdev\t$variance\t$minp\n";
	}
	close DIST;

}

sub pvalue{
    
    #takes in window count average (countAvg) and number of TAsites and makes a null distribution to calculate the pvalue, which it returns
    my $mean=shift@_;
    my $TAsites=shift@_;
    my @specDist=$distLibrary{$TAsites};
    my $rank= binsearch_pos { $a cmp $b } $mean,@specDist;
    my $i=$rank;
    while ($i<scalar(@specDist)-1 and $specDist[$i+1]==$specDist[$rank]){
        $i++;
    }
    $rank=$i;
    my $pval=$rank/$N; #calculate pval as rank/N
    return $pval;
    
}
############################################### End of subroutines

my @marked;
my $marked_set = Set::Scalar->new;
if ($marked) {
   open MARKED, $marked;
   while (<MARKED>) {
      chomp $_;
      $marked_set->insert($_);
   }
   close MARKED;
   #print "Set members: ", $marked_set->members, "\n";
}
my @insertPos; #keep track of ALL insertion positions. Need this for essential distributions.

# Create Gene Summary File.
# Contains Gene => [w1,w2,w3,w4,etc.]
my %data;
foreach my $filename (@files) {
   print $filename,"\n";
   open (IN,'<', $filename);
   my $headerLine=<IN>; #read the header (column names) line and store in dummy variable
   while (my $line = <IN>) {
      chomp $line;
      my @lines = split(/,/,$line);
      my $locus = $lines[9]; #gene id (SP_0000)
      my $w = $lines[12];    #nW
      my $pos = $lines[0];
      push (@insertPos,$pos);
      if (!$w) { $w = 0; }  # For blanks
      my $c1 = $lines[2];  # count_1
      my $c2 = $lines[3];  # count_2
      my $avg = ($c1+$c2)/2; #Later: change which function to use? C1? AVG(C1+C2)?
       
      if ($avg < $cutoff) { next; } # Skip cutoff genes.
      if ($avg >= $weight_ceiling) { $avg = $weight_ceiling; } # Maximum weight.
      
      my @empty;
	  
	  #IF a this is the first insertion for this locus then start empty arrays for values
      if (!$data{$locus}) {
        $data{$locus}{w} = [@empty];
        $data{$locus}{s} = [@empty];
        $data{$locus}{p} = [@empty];
      }
      #IF we've seen this an insertion in this gene already, then just add to the value arrays
      $data{$locus}{w} = [@{$data{$locus}{w}}, $w];  # List of Fitness scores.
      $data{$locus}{s} = [@{$data{$locus}{s}}, $avg]; # List of counts used to generate those fitness scores.
      $data{$locus}{p} = [@{$data{$locus}{p}}, $pos]; #List of insertion positions in the gene
          #later can get UNIQUE insertion positions and total of that mutant
   }
   close IN;
}

#print Dumper(\%data);
    @insertPos = sort { $a <=> $b } @insertPos;
    @insertPos = uniq @insertPos;

my $x="TA";

my @complete; #2d array to hold everything that will be printed

#Get fasta file as string so we can look up the number of TA sites given genomic coordinates
my $in = Bio::SeqIO->new(-file=>$ref); #FindMissing is the genbank file.
my $refseq = $in->next_seq;
my @features = $refseq->get_SeqFeatures;
my $mainheader="locus,mean,var,sd,se,gene,Total,Blank,Not Blank,Blank Removed,M\n";

my %hashtest;

my @header=("locus","start","stop","length","strand","gene", "countTA","insertions","essPval","avgFit","variance","stdev","stderr","blank_ws","num","removed","marked");


foreach my $feature (@features) {
	#print $feature,"\n";
    # Check if it's a gene.
    if ($feature->primary_tag eq 'gene') {

		my $start=$feature->start;
		my $stop=$feature->end;
		my $strand=$feature->strand;
		my $length=$stop-$start;
		my $seq=substr ($fasta, $start, $stop - $start); 
		my @ctemp = $seq =~ /$x/g;
		my $countTA = @ctemp;
        
		my @locus = $feature->get_tagset_values('locus_tag');
		my $locus = $locus[0] || '';
		#print "$locus\t$countTA\n";

		my @gene = $feature->get_tagset_values('gene');
		my $gene = $gene[0] || '';

		my $sum=0;
		my $num=0;
		my $avgsum = 0;
		my $blank_ws = 0;
		
		# Count blanks
        my @summary; #Array to hold all info for a gene that will be printed to out file

		#Check if we have tn-seq fitness information about that gene
		if ($data{$locus}) {
			my $i=0;
			#Count blank fitness scores.
			foreach my $w (@{$data{$locus}{w}}) {
				if (!$w) {
					$blank_ws++;
			  	}
			  	else {
				 	$sum += $w;
				 	$num++;
			  	}
			  	$i++;
			}
			my $count = $num + $blank_ws;

			#Remove blanks from scores if we need to.
		   	my $removed=0;
		   	my $to_remove = int($blank_pc * $count);
		   	if ($blank_ws > 0) {
			  	for (my $i=0; $i < @{$data{$locus}{w}}; $i++) {
					if ($removed == $to_remove) {last}
					my $w = ${$data{$locus}{w}}[$i];
					if (!$w) {
						$removed++;
						splice( @{$data{$locus}{w}}, $i, 1);
						splice( @{$data{$locus}{s}}, $i, 1);
						$i-=1;
					}
			  	}
			}
			
			#First do insertion representation assessment
			my $avgInsert=sprintf("$round",$count/$countTA);
			my $pval=pvalue($avgInsert,$countTA);
				
				
			#Get statistics on gene: average, standard deviation, variance, standard error, etc
			
			#Can't do stats if no observations (insertions)
			if ($num == 0 ) {
				   @summary=($locus,$start,$stop,$length,$strand,$gene,$countTA,$count,$pval,"NA","NA","NA","NA",$blank_ws,$num,$removed);
			}
			#If there are observations (insertions), then do stats
			else {
				
			   my ($average, $variance, $stdev, $stderr);
			   # If weighted was not specified then get regular average
			   if (!$options{w}) {
				   ($average, $variance, $stdev, $stderr) = &average($data{$locus}{w});
			   }
			   else {
				   ($average, $variance, $stdev, $stderr)= &weighted_average($data{$locus}{w},$data{$locus}{s});
			   }
			   
			   
			   
			   #Store all of these variables into a summary for the gene
			   @summary=($locus,$start,$stop,$length,$strand,$gene,$countTA,$count,$pval,$average,$variance,$stdev,$stderr, $blank_ws,$num,$removed); 
			}
		} #End of routine for when we have Tn-Seq data (insertions) for gene
        
        #TEST DELETE
        #In the case that there are no insertions for the gene:
		else {
		   @summary=($locus,$start,$stop,$length,$strand,$gene,$countTA,"NA","NA","NA","NA","NA","NA","NA","NA");
		}

		# If this gene is in the optional list of genes to be marked, then mark it
		if ($marked && $marked_set->contains($locus)) {
			   push(@summary,"M");
		}
		
		#Now that we have all available data on the gene, store it in the 2d array

		push (@complete,\@summary);
		#Try storing the data in a hash instead of a 2d array
		$hashtest{$locus}=\@summary;

	} #end of routine for features that are genes
}


#my @header=("locus","start","stop","length","strand","gene","countTA","avgFit","variance","stdev","stderr","count","blank_ws","num","removed","marked");

open OUT, ">",$out;
print OUT join(",",@header),"\n";
foreach my $key (sort keys %hashtest){
	#print $key,"\n";
	my @vals=@{$hashtest{$key}};
	print OUT join(",",@vals);
	print OUT "\n";
	}
close OUT;
