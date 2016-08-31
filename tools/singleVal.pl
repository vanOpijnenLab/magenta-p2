#!/usr/bin/perl -w

#Margaret Antonio updated 15.10.18

use strict;
use Getopt::Long;
use warnings;
use Bio::SeqIO;
use Scalar::Util;

#AVAILABLE OPTIONS. WILL PRINT UPON ERROR
sub print_usage() {
	print "\n####################################################################\n";
    print "singleVal: creates a wig file with position and insertion count OR fitness\n";
    print "\nDESCRIPTION: ";
    print "Integrates multiple files of transposon insertion data and outputs\n";
    print "aggregate fitness within a sliding window (specified by size and step).\n";
    
    print "\nUSAGE:\n";
    print "perl windowFit.pl <OPTIONS> <REQ OUTPUT TYPE(S)> <INPUT FILE(S)>\n\n";
    
    print "\nREQUIRED:\n";
    print " -d\tDirectory containing all input files (files from\n";
    print "   \taggregate script)\n";
    print "   \tOR\n";
    print "   \tIn the command line (without a flag), input the name(s) of\n";
    print "   \tfiles containing gene fitness values (output of calcFit). \n\n";
    
    print "OPTIONAL:\n";
    print " -h\tPrints usage and exits program\n";
    print " -o\tOutput file for comparison data. Default: singleVal.wig\n";
    print " -v\tString value for output: 'fit' for fitness OR 'count' for count\n";
    print "   \tDefault: fit for fitness\n";
    print " -n\tName of the reference genome, to be included in the wig header\n";
    print "   \tDefault: genome\n";

    print " \n~~~~Always check that file paths are correctly specified~~~~\n";
    print "\n##################################################################\n";
    
}


#ASSIGN INPUTS TO VARIABLES
our ($cutoff,$help,$ref_genome,$indir,$out,$name,$val);
GetOptions(
'cutoff:i'=>\$cutoff,
'd:s' => \$indir,
'h'=>\$help,
'v:s' =>\$val,
'o:s' =>\$out,
'n:s' =>\$name,
);

sub get_time() {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return "$hour:$min:$sec";
}

if ($help){
    print_usage();
    exit;
}
if (!$indir and (scalar @ARGV==0)){
	print "\nERROR: Please correctly specify input files or directory\n";
    print_usage();
	print "\n";
	exit;
}

#SET DEFAULTS
if (!$cutoff){$cutoff=15}
if (!$name){$name="genome"}

if ($val eq "count"){

		my $rowCount=-1;
		my $last;
		my @unsorted;

		# EXTRACT INPUT FILES
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

		# READ INPUT FILES

		my $num=scalar @files;
		print "\nNumber of files in csv: ", $num,"\n";

		my %select; 
		my %input;

		for (my $i=0; $i<$num; $i++){   
			my $file=$files[$i];
			print "File #",$i+1,"\t",$file,"\n";
	
			open(DATA,'<', $file) or die "Could not open '$file'";
	
			my $dummy=<DATA>;
			while(my $line=<DATA>){
			chomp($line);
			my @fields=split(",",$line);
			my $pos=int($fields[0]);
			my $c2=int($fields[3]);
			if (exists $input{$pos}){
				$input{$pos}+=$c2;
				}
			else{
				$input{$pos}=$c2;
				}
			}
			close DATA; 
		}

		# OUTPUT A WIG FILE WTH POS AND COUNT
        if (!$out){$out="singleCount.wig"}
        
		open (OUT,'>',$out);

		print OUT "variableStep\tchrom=chrN\n";
		for my $pos(sort {$a<=>$b} keys %input){
			print OUT $pos, "\t",$input{$pos},"\n";
		}
		close OUT;
}
##########################################################################

else{

		#CREATE AN ARRAY OF DATA FROM INPUT CSV FILE(S)
		print "\nStart input array ",get_time(),"\n";

		my $rowCount=-1;
		my $last;
		my @unsorted;

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

		my $num=scalar @files;
		print "\nNumber of files in csv: ", $num,"\n";

		my %select; #this is a hash to store insertion position and fitness value at that position for each file (library)

		for (my $i=0; $i<$num; $i++){   #Read files from ARGV---loops once for each file in ARGV
			my $file=$files[$i];
			open(DATA,'<', $file) or (print "Could not open '$file'\n" and print_usage() and exit);
			print "File #",$i+1,"\t",$file,"\n";
			my $dummy=<DATA>;
			while (my $line = <DATA>) {
				chomp $line;
				my @fields=split(",",$line);
				my $w = $fields[12];
				if (!$w or $w eq"\n"){
					next;
				} 
				# For blanks
				else{
					my $c1 = $fields[2];
					my $c2 = $fields[3];
					my $avg = ($c1+$c2)/2;
					if ($avg < $cutoff) {
						next;
				
					} # Skip cutoff genes.
					else {  #This is a good value and should be added to the hash
						my $pos=int($fields[0]);
						$w=sprintf("%.2f",$w);
						if(!exists $select{$pos}){
						#print $pos;
							my @fits;
							push(@fits,"1");
							for(my $j=0;$j<$i-1;$j++){
								push (@fits,"0");
							}
							push (@fits,$w);
							$select{$pos}=\@fits;
						}
						else{
							my @fits=@{$select{$pos}};
							for (my $j=scalar @fits+1;$j<$i;$j++){
								push (@fits,"0");
							}
							push (@fits,$w);
							$select{$pos}=\@fits;
						}
					}
				}
			}
			close DATA;
		}
		##### FINISHED READING INPUT FILES

		foreach my $val(values %select){
			if (scalar @$val <$num+1){
				for (my $j=scalar @$val;$j<$num+1;$j++){
					push (@$val,"0");
				}
			}
		}

		#### DEFAULT OUTPUT	
		if (!$out){
			$out="singleFit.wig";
			}
	
		#### OUTPUT SELECTED VALUE
		open OUT2,'>', $out;
	    print OUT2 "variableStep\tchrom=chrN\n";

		for my $entry (sort {$a<=>$b} keys %select) {
			print OUT2 $entry,"\t",$entry+1,"\t";
			my @entryFits=@{$select{$entry}};
			my $sum=0; my $count=0;    
			foreach (@entryFits){
				if ($_!=0){
					$count++;
					$sum+=$_;
				}
			}
			my $avgFit=sprintf("%.2f",$sum/$count);
			print OUT2 $avgFit,"\n";
		}
		close OUT2;
	

}

