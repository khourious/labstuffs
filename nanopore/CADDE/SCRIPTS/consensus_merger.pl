#!/usr/bin/perl

use warnings;
use strict;

# Silly consensus merger #
# usage: perl consensus_merger.pl ref.fa cons.fa
# ref.fa -> reference file in fasta format
# cons.fa -> fasta file containing consensus sequences to be merged

my $line;
my $head;
my $refseq;
my $k;
my $i;
my $ref = $ARGV[0];
my $cons = $ARGV[1];
my %data;


open IN,"<$ref";
while($line = <IN>){
	chomp($line);
	if($line =~  />/){
		$head = $line;
	}else{
		$refseq .= uc($line);
	}
}
close IN;


open IN,"<$cons";
while($line = <IN>){
	chomp($line);
	if($line =~ />/){
		$k++;
	}else{
		$data{$k} .= uc($line);
	}
}
close IN;

my @sitesref = split('',$refseq);
my @sites1 = split('',$data{"1"});
my @sites2 = split('',$data{"2"});


print "$head\n";
for($i = 0; $i <= length($refseq) - 1; $i++){
		if($sites1[$i] ne "N"){
			print "$sites1[$i]";
		}elsif($sites2[$i] ne "N"){
			print "$sites2[$i]";
		}else{
			print "N";
		}
}
print "\n";























