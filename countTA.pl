#!/usr/bin/perl -w

#Margaret Antonio 15.04.13

#my $genome=shift@_; #Fasta file
#$lgenome=length(genome);
#$count=0;



use Bio::SeqIO;

my $count=0;

my $seqio = Bio::SeqIO->new(-file => "tigr4_genome.fasta", '-format' => 'Fasta');
my $prev;
my $total=0;
while(my $seq = $seqio->next_seq) {
    my $string = $seq->seq;
    # do stuff with $string
   for my $c (split //, $string){
       if (($c eq'A') and ($prev eq'T')){
            $count++;
        }
       $prev=$c;
       $total++;
    }
}

print "This is the TA count: $count\nTotal genome size is: $total\n\n";

my $genseq = Bio::SeqIO->new(-file => "NC_003028.gbk", '-format' => 'Fasta');