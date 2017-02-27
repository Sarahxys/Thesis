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
my @line;

open(INPUT,"<", "$inputFile") or die "could not open the input file";

sub addtight{
	my $fragtid = scalar keys %{$transcript_hash{$line[1]} +1;
    $transcript_hash{$line[1]}{$fragtid}[1]= $line[6]+0; 
    $transcript_hash{$line[1]}{$fragtid}[2]= $line[7]+0; 
    $transcript_hash{$line[1]}{$fragtid}[3] = $line[8]+0;
    $transcript_hash{$line[1]}{$fragtid}[4]= $line[9]+0;
    $transcript_hash{$line[1]}{$fragtid}[5]= $line[10]+0;
    $transcript_hash{$line[1]}{$fragtid}[6]= $line[11];
    $transcript_hash{$line[1]}{$fragtid}[8]= $line[3];
}

sub IdentifySplicingVariant {


	#delete the matched transcript name with only one fragment
	foreach my $sseqID (sort {$a cmp $b} keys %transcript_hash)
		if (scalar keys %{$transcript_hash{$sseqID} < 2){
			delete $transcript_hash{$sseqID};
		}
		else{
			
		}
	}

	if (scalar keys %transcript_hash){
		foreach my $sseqID (sort {$a cmp $b} keys %transcript_hash)
			foreach my $fragment (sort {$a <=> $b} keys %{$transcript_hash{$sseqID}}){
	        	print "$oldfragid\t$fragment\t";
		        foreach my $i (0..$#{$transcript_hash{$oldfragid}{$fragment}}){
		          print "$transcript_hash{$sseqID}{$fragment}[$i]\t";
		        }
		         print "\n";
		    }
		}
		return %transcript_hash;
	}
	else{
		return %transcript_hash;
	}

}




#---------------------------------------------------main body--------------------------------------------------------
while ( my $line = <INPUT>) {
	chomp($line);
	@line = split (/\s+/,$line); #split the line with spaces (\s = space; and \s+ = more than one spece);

	if ($line[0] == $oldfragid){ 
		addfragment;
	}
	else{
		if ($oldfragid){
			IdentifySplicingVariant; 
			%transcript_hash=();
		}
		addfragment;
		$oldfragid = $line[0];
	}

}
close INPUT;