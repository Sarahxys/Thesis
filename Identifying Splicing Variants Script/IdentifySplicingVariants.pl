#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
#no autovivification;

################################################################################
# goal: identify total exons, max exons per transcript, min exons per transcript, location of variants (ends or middle)
# method:
# - firstly, run Blast with setof assembled transcripts flagged with " no error" or "contain scaffold" against itself.
# - then use the blast result and coordinate information from previous genome blast to achieve goal
# criteria
# - transcript that are just repeat region of each other does not count as splicing variants
# - transcripts that identified as splicing variants should have genome coordinate close to each other
# - transcripts that identified as splicing variants have to have two exons mapped to each other or have one large mapped exons that is greater than 200bp (changable variable)
#######################################

#open the input file
my ($inputFile1,$inputFile2)=@ARGV; #inputfile 1 is the detailfile(ex, 3909_Detail_OLL5); inputfile2 is the blastn output generated from blastn of extractedtranscriptome to itself(ex, blastout_3909);
my %transcript_hash;
my %transcriptsDetail_hash;
my %transcriptsChrID_hash;
my %isASVariants;
my %ASVariants;
my %GenomeCoordinate;
my @line;
my $oldfragid;
my $numberoftranscript = 0;
my $ToleratedProximaty = 10000; 

open(INPUT1,"<", "$inputFile1") or die "could not open the input file";
open(INPUT2,"<", "$inputFile2") or die "could not open the input file";

sub addfragment{
	my $fragtid = (scalar keys %{$transcript_hash{$line[1]}}) + 1;
    ############
	#blast format 6
	#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
	############

    $transcript_hash{$line[1]}{"F$fragtid"}[0]= $line[6]+0; #qstart
    $transcript_hash{$line[1]}{"F$fragtid"}[1]= $line[7]+0; #qend 
    $transcript_hash{$line[1]}{"F$fragtid"}[2] = $line[8]+0; #sstart
    $transcript_hash{$line[1]}{"F$fragtid"}[3]= $line[9]+0; #send

    if ($transcript_hash{$line[1]}{"F$fragtid"}[2] < $transcript_hash{$line[1]}{"F$fragtid"}[3]){ #determine the orientation of the fragment: if forward = 1; if reverse = -1;
        $transcript_hash{$line[1]}{"F$fragtid"}[4] = 1;
    }
    else {
        $transcript_hash{$line[1]}{"F$fragtid"}[4] = -1;
    }
    return %transcript_hash;

}


sub IsAlternativeSplcingVariant {
	# the sub will evalute is the aligned/matched transcript is alternative splicing variant for the query transcript
	my $sseqID=shift;
	#re-store the fragments in ascending order based on qstart position
	my @qstart_array;
	my @sortedqstart;
	foreach my $fragmentid (sort keys %{$transcript_hash{$sseqID}}){
		push (@qstart_array, $transcript_hash{$sseqID}{$fragmentid}[0]);
	}

	foreach my $fragmentid (sort keys %{$transcript_hash{$sseqID}}){        
	    my $counter = 1;
	    my $compareqstart = $transcript_hash{$sseqID}{$fragmentid}[0];
	    foreach my $i (@qstart_array){        
		     if ($compareqstart > $i){
			     $counter ++;
		     }
	    }
	    $transcript_hash{$sseqID}{$counter} = delete $transcript_hash{$sseqID}{$fragmentid};
	}

	#assuming transcripts that are alternative splicing variants of each other should be matching to each other with asending order.
	my $oldqend=0;
	my $oldsend =0;
	foreach my $fragmentid (sort {$a <=>$b} keys %{$transcript_hash{$sseqID}}){ 
		my $direction = $transcript_hash{$sseqID}{$fragmentid}[4];

		if ($fragmentid > 1){ #
			my $qstart = $transcript_hash{$sseqID}{$fragmentid}[0];
			my $sstart = $transcript_hash{$sseqID}{$fragmentid}[2];
			if (($qstart < $oldqend && $oldqend - $qstart > 10) ||
				($direction == 1 && $oldsend - $sstart >10)){
				delete $transcript_hash{$sseqID};
				return %transcript_hash;
			}
			elsif ($direction == -1){
				delete $transcript_hash{$sseqID};
				return %transcript_hash;
			}
			else{
				$oldqend = $transcript_hash{$sseqID}{$fragmentid}[1];
				$oldsend = $transcript_hash{$sseqID}{$fragmentid}[3];
			}


		}
		else{ #it is the first fragemnt, have nothing to compare to 
			
			if ($direction == -1){
				delete $transcript_hash{$sseqID};
				return %transcript_hash;
			}

			$oldqend = $transcript_hash{$sseqID}{$fragmentid}[1];
			$oldsend = $transcript_hash{$sseqID}{$fragmentid}[3];
		}
	}

	#check if the two transcripts are close to each other on the genomic coordinate
	if ($transcript_hash{$sseqID}){
		my $lastFragId = (scalar keys %{$transcript_hash{$sseqID}}) +1;
		my $qstart = $transcript_hash{$sseqID}{"1"}[0];
		my $qend = $transcript_hash{$sseqID}{$lastFragId}[1];
		my $sstart = $transcript_hash{$sseqID}{"1"}[2];
		my $send = $transcript_hash{$sseqID}{$lastFragId}[3];

		#calculate genomic location for the hit
		my @sClosestCoordinate;
		$sClosestCoordinate[0] = 10e10;
		my @sGenomeCoordinate;
		foreach my $fragmentid(sort {$a <=> $b} keys %{$GenomeCoordinate{$sseqID}}){
			if (abs($GenomeCoordinate{$sseqID}{$fragmentid}[0] - $qstart) < $sClosestCoordinate[0]){
				$sClosestCoordinate[0] = $qstart - $GenomeCoordinate{$sseqID}{$fragmentid}[0];
				$sClosestCoordinate[1] = $GenomeCoordinate{$sseqID}{$fragmentid}[2]; #get the genomic start coordinate
				$sClosestCoordinate[2] = $GenomeCoordinate{$sseqID}{$fragmentid}[3];#get the genomic end coordiante; may not be neccessary
			}

		}
		
		if ($sClosestCoordinate[0] == 0){
			$sGenomeCoordinate[0]=$sClosestCoordinate[1];
			$sGenomeCoordinate[1]=$sClosestCoordinate[2];
		}
		elsif($sClosestCoordinate[0] >0){
			$sGenomeCoordinate[0]=$sClosestCoordinate[1] + abs($sClosestCoordinate[0]);
			$sGenomeCoordinate[1]=$sClosestCoordinate[2] + abs($sClosestCoordinate[0]);
		}
		elsif($sClosestCoordinate[0] <0){
			$sGenomeCoordinate[0]=$sClosestCoordinate[1] - abs($sClosestCoordinate[0]);
			$sGenomeCoordinate[1]=$sClosestCoordinate[2] - abs($sClosestCoordinate[0]);
		}

		#calculate the genome location for the query
		my @qClosestCoordinate;
		$qClosestCoordinate[0] = 10e10;
		my @qGenomeCoordinate;
		foreach my $fragmentid(sort {$a <=> $b} keys %{$GenomeCoordinate{$oldfragid}}){
			if (abs($GenomeCoordinate{$oldfragid}{$fragmentid}[0] - $qstart) < $qClosestCoordinate[0]){
				$qClosestCoordinate[0] = abs($GenomeCoordinate{$oldfragid}{$fragmentid}[0] - $qstart);
				$qClosestCoordinate[1] = $GenomeCoordinate{$oldfragid}{$fragmentid}[2]; #get the genomic start coordinate
				$qClosestCoordinate[2] = $GenomeCoordinate{$oldfragid}{$fragmentid}[3];#get the genomic end coordiante; may not be neccessary
			}

		}
		if ($qClosestCoordinate[0] == 0){
			$qGenomeCoordinate[0]=$qClosestCoordinate[1];
			$qGenomeCoordinate[1]=$qClosestCoordinate[2];
		}
		elsif($qClosestCoordinate[0] >0){
			$qGenomeCoordinate[0]=$qClosestCoordinate[1] + abs($qClosestCoordinate[0]);
			$qGenomeCoordinate[1]=$qClosestCoordinate[2] + abs($qClosestCoordinate[0]);
		}
		elsif($qClosestCoordinate[0] <0){
			$qGenomeCoordinate[0]=$qClosestCoordinate[1] - abs($qClosestCoordinate[0]);
			$qGenomeCoordinate[1]=$qClosestCoordinate[2] - abs($qClosestCoordinate[0]);
		}

		#compare genomic coordinate to determine if they are alternative splicing variants
		if(($qGenomeCoordinate[0] <= $sGenomeCoordinate[0] && $sGenomeCoordinate[0] <= $qGenomeCoordinate[1])||
			($qGenomeCoordinate[0] <= $sGenomeCoordinate[1] && $sGenomeCoordinate[1] <= $qGenomeCoordinate[1])||
			($sGenomeCoordinate[0] <= $qGenomeCoordinate[0] && $qGenomeCoordinate[0] <= $sGenomeCoordinate[1])||
			($sGenomeCoordinate[0] <= $qGenomeCoordinate[1] && $qGenomeCoordinate[1] <= $sGenomeCoordinate[1])) { #case 1: their genomic coordinate are overlapping or are exact the same;
			return %transcript_hash;
		}
		elsif((abs($sGenomeCoordinate[0] - $qGenomeCoordinate[0]) <= $ToleratedProximaty)||
			(abs($sGenomeCoordinate[1] - $qGenomeCoordinate[1]) <= $ToleratedProximaty)){
				return %transcript_hash;
		}
		else {
			delete $transcript_hash{$sseqID};
			return %transcript_hash;
		}
	}
	return %transcript_hash;
}

sub IdentifySplicingVariant{
	
	#eliminate aligned transcripts that have only one fragment or doesn't have the same ChrID as the query transcript 
	foreach my $sseqID (sort keys %transcript_hash){
		#delete the matched transcript name with only one fragment
		if (scalar keys %{$transcript_hash{$sseqID}} < 2 || $sseqID eq $oldfragid){
			delete $transcript_hash{$sseqID};
		}
		
		if (scalar keys %transcript_hash <=1){%transcript_hash =(); return %transcript_hash;}

		#check if they are on the same chr or scaffold, if not, delete them
		my $delete = 1; #set delete to 1(true) if no match found; once a match of chrID found, set it to 0(false) and sseqID will not be eliminated.
		
		foreach my $sseqIDchrID (@{$transcriptsChrID_hash{$sseqID}}){
			foreach my $QueryChrID(@{$transcriptsChrID_hash{$oldfragid}}){
				if ($QueryChrID eq $sseqIDchrID){
					$delete = 0;
				}
			}
		}

		if ($delete){delete $transcript_hash{$sseqID};}
		if (scalar keys %transcript_hash <=1){%transcript_hash =(); return %transcript_hash;}
	}


	if (scalar keys %transcript_hash >=1){
		foreach my $sseqID (sort {$a cmp $b} keys %transcript_hash){
			IsAlternativeSplcingVariant($sseqID);
			if ($transcript_hash{$sseqID}){
				$numberoftranscript ++;
				$isASVariants{$oldfragid}{$sseqID} = delete $transcript_hash{$sseqID}
			}
		}

		# #print out the detail of transcripts that passed the checks so that I can decide how to precess for the next step
		# foreach my $sseqID (sort {$a cmp $b} keys %{$isASVariants{$oldfragid}}){
		# 	foreach my $fragment (sort {$a <=> $b} keys %{$isASVariants{$oldfragid}{$sseqID}}){
	 #        	print "$oldfragid\t$sseqID\t$fragment\t";
		#         foreach my $i (@{$isASVariants{$oldfragid}{$sseqID}{$fragment}}){
		#           print "$i\t";
		#         }
		#          print "\n";
		#     }
		# }
		return %transcript_hash;
	}
	else{
		%transcript_hash =();
		return %transcript_hash;
	}

}




#---------------------------------------------------main body--------------------------------------------------------
my $count=0;
while (my $line = <INPUT1>){
	chomp($line);
	@line = split (/\s+/,$line); #split the line with spaces (\s = space; and \s+ = more than one spece);

	#format of transcript in Detail file
	#0=qseqID; 1= fragment number; 2=chrID; 3= qstart; 4=qend; 5=sstart; 6=send; 7=e-value; 8=bit score; 9=direction(1 or -1); 10=length
	push (@{$transcriptsChrID_hash{$line[0]}}, $line[2]);

	my $fragNum = scalar keys %{$GenomeCoordinate{$line[0]}}; 
	$GenomeCoordinate{$line[0]}{$fragNum}[0]= $line[3];
	$GenomeCoordinate{$line[0]}{$fragNum}[1]= $line[4];
	$GenomeCoordinate{$line[0]}{$fragNum}[2]= $line[5];
	$GenomeCoordinate{$line[0]}{$fragNum}[3]= $line[6];
	$count++;
}
close INPUT1;
@line=();
print "the number of push:$count\n";

open(INPUT2,"<", "$inputFile2") or die "could not open the input file";


while ( my $line = <INPUT2>) {
	chomp($line);
	@line = split (/\s+/,$line); #split the line with spaces (\s = space; and \s+ = more than one spece);

	if ($oldfragid && $line[0] eq $oldfragid){ 
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
close INPUT2;

my $noMatch = 1;
my @seqID_array;
foreach my $sseqID(sort keys %isASVariants){
	push (@seqID_array, $sseqID);
	foreach my $hitID(sort keys %{$isASVariants{$sseqID}}){
		push (@seqID_array, $hitID);
	}
	
	my $chrID = $transcriptsChrID_hash{$sseqID};
	if (exists $ASVariants{$chrID}){
		foreach my $groupNum(sort {$a cmp $b} keys %{$ASVariants{$chrID}}){
			if ( !array_diff(@seqID_array, @{$ASVariants{$chrID}{$groupNum}})){
	            #move the item that is in @seqID_array to 
	            my @diff = array_diff(@seqID_array, @{$ASVariants{$chrID}{$groupNum}});
	            push (@{$ASVariants{$chrID}{$groupNum}}, @diff);
	            $noMatch = 0;
	    	}
	    }

	    if ($noMatch){
	    	my $newgroupID = (scalar keys %{$ASVariants{$chrID}}) +1;
	    	@{$ASVariants{$chrID}{$newgroupID}} = @seqID_array;
	    }
	}
	else{
		@{$ASVariants{$chrID}{"1"}} = @seqID_array;
	}


}


#calculate the total number of ASV group
my $totalASVgroup =0;
foreach my $chrID(sort keys %ASVariants){
	$totalASVgroup = $totalASVgroup + (scalar keys %{$ASVariants{$chrID}});
	foreach my $groupID(sort {$a <=> $b} keys %{$ASVariants{$chrID}}){
		print "$chrID\t$groupID\t";
		foreach my $i (@{$ASVariants{$chrID}{$groupID}}){
			print "$i\t";
		}
		print "\n";
	}
}


print "the total number of alternative splicing variants group are: $totalASVgroup\n";
