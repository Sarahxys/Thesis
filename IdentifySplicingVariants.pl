#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
#no autovivification;

################################################################################
# goal: identify total exons, max exons per transcript, min exons per transcript, location of variants (ends or middle)
# method:
# - firstly, run Blast with setof assembled transcripts flagged with " no error" or "contain scafold" against itself.
# - then use the blast result and coordinate information from previous genomic blast to achieve goal
# critiria
# - trnascript that are just repeat region of each other does not count as splcing variants
# - transcripts that identified as splicing variants should have genomic coordinate close to each other
# - transcripts that identifiedas splcing variants have to have two exons mapped to each other or have one large mapped exons that is greater than 200bp (changable variable)
#######################################

#open the input file
my ($inputFile)=@ARGV;
my %transcript_hash;

open(INPUT,"<", "$inputFile") or die "could not open the input file";

while ( my $line = <INPUT>) {
	chomp($line);
	@line = split (/\s+/,$line); #split the line with spaces (\s = space; and \s+ = more than one spece);

	if (exists $transcript_hash{$line[0]}){ 
		$transcript_hash{matched transcript id}{fragment number}[array number]
	}
	else{

	}