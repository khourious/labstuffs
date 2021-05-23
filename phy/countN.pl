#!/usr/bin/perl -w 
use strict;
use Bio::SeqIO;

my $seqio = Bio::SeqIO->new(-file => "$ARGV[0]", -format => 'fasta');
while(my $seq = $seqio->next_seq) {
  my $s = $seq->seq;
  $s =~ s/n/N/g;
  my $id  = $seq->display_id; # sample id
  my $len = $seq->length; # length
  my %count;
  $count{$_}++ foreach split //, $s; # count bases
  $count{"N"} = 0 if !exists($count{"N"});
  my $ratio = ($count{"N"} / $len * 100); # percent of N
  my @output = ($id, $len, $count{"N"}, $ratio);
  print join("\t", @output), "\n";
}
