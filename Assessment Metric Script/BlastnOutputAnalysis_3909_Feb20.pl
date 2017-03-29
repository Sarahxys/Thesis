#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
#no autovivification;

#####
#This program will analysis the blastn output of trinity assembly. -> this version will store the tiess and output the number of sequences with tiess ->
#-> so that I can know which method to precess
#This script will read into the result file of blastn(in fasta format) and for each transcriptome sequence, aligned fragments will be filter based on e-value.
#After filtering, fragments will be re-order based on query coordinates. 
#The script will then identify the possible error in each transcritome sequence.
#VersionJan11_printlist: this version will be motified to select the optimal path when 2 chrid both have the max occurence for tiess. 
# this version will also motified to have transcript_hash contain fragments/tiess for one query id only. The resulting transcript-id will be stored at error_hash and print out to printlist(output #1); 
#infomation including coordinates, chrid, and etc will be print out to printinfor(output #2)  
############
#blast format 6
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
############
#progress
#The current script can select and store fragment in a hash(transcript_hash) based on evalue;
#done: foreach loop to compare the new fragment to all the fragment store in the hash
############
#Decisions 
#hash structure: qseqid, fragment number, ties number, array for stats storage(last level)
#last level array order: 0=sseqid; 1=qsatart; 2=qend; 3=sstard; 4= send; 5=evalue; 6=bitsore; 7=forward(1)/reverse(-1); 8=length;
#3 main type of errors: are they the same chromosom; order of the chromosome; the orientation of the fragments(forward or reverse)
#use e-value as the index; always keep the fragment with the lowest e-value(closer to 0)
#the toleranted overlapping region is 5bp
#I will modify transcript_hash as a lexical(global/my) variable in subroutine instead of passing it to subroutine as a reference;(it will be interesting to change it to the reference/dereference method and do a time-analysis)
#there will be two different subroutines for sorting tiess for tropicalis and laevis since they have difference genomic characteristics
###########
#improvement
#store the one that is ties and make decision later: 
#split the blastnoutput file into multiple smaller file and do parallel running
#blast the gene id of single reads into xenbase and see where is that located

#open the input file
my ($inputFile)=@ARGV;
my @line;
my %transcript_hash;
my %error_hash;
my $numfragment = 0;
my $DoneLoading = 0;
my $OverlapLengthLimit = 5;
my %tiesseq_hash;
my $tiesOccur = 0;
my $check_bothdir = 0;
my $oldfragid = 0;

my $printlist = "3909_Printlist_OLL$OverlapLengthLimit";
my $printdetail = "3909_PrintDetail_OLL$OverlapLengthLimit";
my $printother = "3909_Printother_OLL$OverlapLengthLimit";

GetOptions(
	'oll=i' => \$OverlapLengthLimit,
	"out1=s" => \$printlist,
	"out2=s" => \$printdetail
	);

$printlist = "3909_Printlist_OLL$OverlapLengthLimit";
$printdetail = "3909_PrintDetail_OLL$OverlapLengthLimit";
$printother = "3909_Printother_OLL$OverlapLengthLimit";



open(INPUT,"<", "$inputFile") or die "could not open the input file";
open(OUTPUT2, ">","$printdetail") or die "could not open the output printdetail";
open(OUTPUT1, ">","$printlist") or die "could not open the output printlist";
open(OUTPUT3, ">","$printother") or die "could not open the output printother";


#--------------------------------------------subroutine-------------------------------------------------------------
sub addfragment {
     #transcript_hash stored the data/info of each fragment as an array in the order of: 0=sseqid; 1=qsatart; 2=qend; 3=sstard; 4= send; 5=evalue; 6=bitsore; 7=forward(1)/reverse(-1); 8=length;
    $transcript_hash{$line[0]}{"fragment$numfragment"}{"T0"}[0]= $line[1]; 
    $transcript_hash{$line[0]}{"fragment$numfragment"}{"T0"}[1]= $line[6]+0; 
    $transcript_hash{$line[0]}{"fragment$numfragment"}{"T0"}[2]= $line[7]+0; 
    $transcript_hash{$line[0]}{"fragment$numfragment"}{"T0"}[3] = $line[8]+0;
    $transcript_hash{$line[0]}{"fragment$numfragment"}{"T0"}[4]= $line[9]+0;
    $transcript_hash{$line[0]}{"fragment$numfragment"}{"T0"}[5]= $line[10]+0;
    $transcript_hash{$line[0]}{"fragment$numfragment"}{"T0"}[6]= $line[11];
    $transcript_hash{$line[0]}{"fragment$numfragment"}{"T0"}[8]= $line[3];
   
    if ($transcript_hash{$line[0]}{"fragment$numfragment"}{"T0"}[3] < $transcript_hash{$line[0]}{"fragment$numfragment"}{"T0"}[4]){ #determine the orientation of the fragment: if forward = 1; if reverse = -1;
        $transcript_hash{$line[0]}{"fragment$numfragment"}{"T0"}[7] = 1;
    }
    else {
        $transcript_hash{$line[0]}{"fragment$numfragment"}{"T0"}[7] = -1;
    }
    return %transcript_hash;
}


sub addties{
	my ($fragtid)=shift;
	my $tiesnumber = scalar keys %{$transcript_hash{$line[0]}{$fragtid}};
	$transcript_hash{$line[0]}{$fragtid}{"T$tiesnumber"}[0]= $line[1]; 
    $transcript_hash{$line[0]}{$fragtid}{"T$tiesnumber"}[1]= $line[6]+0; 
    $transcript_hash{$line[0]}{$fragtid}{"T$tiesnumber"}[2]= $line[7]+0; 
    $transcript_hash{$line[0]}{$fragtid}{"T$tiesnumber"}[3] = $line[8]+0;
    $transcript_hash{$line[0]}{$fragtid}{"T$tiesnumber"}[4]= $line[9]+0;
    $transcript_hash{$line[0]}{$fragtid}{"T$tiesnumber"}[5]= $line[10]+0;
    $transcript_hash{$line[0]}{$fragtid}{"T$tiesnumber"}[6]= $line[11];
    $transcript_hash{$line[0]}{$fragtid}{"T$tiesnumber"}[8]= $line[3];
    
    if ($transcript_hash{$line[0]}{$fragtid}{"T$tiesnumber"}[3] < $transcript_hash{$line[0]}{$fragtid}{"T$tiesnumber"}[4]){ #determine the orientation of the fragment: if forward = 1; if reverse = -1;
        $transcript_hash{$line[0]}{$fragtid}{"T$tiesnumber"}[7] = 1;
    }
    else {
        $transcript_hash{$line[0]}{$fragtid}{"T$tiesnumber"}[7] = -1;
    }
    return %transcript_hash;
}

sub printdetail{
	my $qseqID = shift;
	#print the detail to output2
	if ($error_hash{$qseqID} eq "no error"||$error_hash{$qseqID} eq "Contain scaffold"){
		if ($tiesseq_hash{$qseqID}){#print detail of transcripts with tiess
			foreach my $fragment (sort keys %{$transcript_hash{$qseqID}}){
		        foreach my $ties(sort keys %{$transcript_hash{$qseqID}{$fragment}}){
		        	print OUTPUT2 "$qseqID\t$fragment\t";
			        foreach my $i (0..$#{$transcript_hash{$qseqID}{$fragment}{$ties}}){
			          print OUTPUT2 "$transcript_hash{$qseqID}{$fragment}{$ties}[$i]\t";
			        }
			         print OUTPUT2 "\n";
			    }
		    }
		}
		else{#print detail of transcripts without tiess
			foreach my $fragment (sort {$a <=> $b} keys %{$transcript_hash{$qseqID}}){
	        	print OUTPUT2 "$qseqID\t$fragment\t";
		        foreach my $i (0..$#{$transcript_hash{$qseqID}{$fragment}}){
		          print OUTPUT2 "$transcript_hash{$qseqID}{$fragment}[$i]\t";
		        }
		         print OUTPUT2 "\n";
		    }
		}
	}
	return;
}

sub checkerror_tropicalis{
	my ($qseqID) = shift;
	#re-store the fragments in ascending order based on qstart position
	my @qstart_array;
	foreach my $fragmentid (sort keys %{$transcript_hash{$qseqID}}){
		push (@qstart_array, $transcript_hash{$qseqID}{$fragmentid}{"T0"}[1]);
	}

	foreach my $fragmentid (sort keys %{$transcript_hash{$qseqID}}){        
	    my $counter = 1;
	    my $compareqstart = $transcript_hash{$qseqID}{$fragmentid}{"T0"}[1];
	    foreach my $i (@qstart_array){        
		     if ($compareqstart > $i){
			     $counter ++;
		     }
	    }
	    $transcript_hash{$qseqID}{$counter} = $transcript_hash{$qseqID}{$fragmentid}{"T0"};
	    delete $transcript_hash{$qseqID}{$fragmentid};
	}

	#identify the errors      
	$error_hash{$qseqID} = "no error";

	#check if the fragments are on the same chromosome
	my $firstfrag;
	$firstfrag = $transcript_hash{$qseqID}{"1"}[0];
	my $firstfragname = $firstfrag; 
	$firstfragname =~ s/[^a-zA-Z]//g;#extrat chrID name - whetther it is a scaffold or chromosome
	my $firstfragnum = $firstfrag; 
	$firstfragnum =~ s/\D//g;#extrat the chromosome number
	
	foreach my $fragmentid (sort {$a<=>$b} keys %{$transcript_hash{$qseqID}}){
		my $comparefrag;
		my $comparefragname;
		my $comparefragnum;
		$comparefrag = $transcript_hash{$qseqID}{$fragmentid}[0];
		$comparefragname = $comparefrag; 
		$comparefragname =~ s/[^a-zA-Z]//g;#extrat chrID name - whetther it is a scaffold or chromosome
	    $comparefragnum = $comparefrag; 
	    $comparefragnum =~ s/\D//g;#extrat the chromosome number
		if ($firstfrag ne $comparefrag){
		    if ($firstfragnum > 10 && $comparefragnum <=10){ #case 1: scaffold >10 compare with chr/scaffold 1-10; no error but document it as containing scaffold; meanwhile, prevent it proceed with orientation check.  
		     	 $error_hash{$qseqID} = "Contain scaffold";
		     	 $firstfragnum = $comparefragnum; #if the first scaffold > 10 and now encounter a scaffold/chr 1-10, change $firstfrag to scaffold/chr1-10. 
		     	 $firstfragname =$comparefragname;#because scaffold>10 will be ignored in comparison and will cause the un-documentation of error case 
		     	 $firstfrag = $comparefrag;
		    }
		    elsif(($firstfragnum <= 10 && $comparefragnum <=10)&&($firstfragnum == $comparefragnum)&&($firstfragname ne $comparefragname)){#case 2: scaffold 1-10 compare with chr/scaffold 1-10 and have the same number; no error but document it as containing scaffold; meanwhile, prevent it proceed with orientation check.
		    	 $error_hash{$qseqID} = "Contain scaffold";
		    }
		    elsif(($firstfragnum <= 10 && $comparefragnum <=10)&&($firstfragnum != $comparefragnum)){#case 3: scaffold/chr1-10 compare with caffold/chr1-10 and it doesnt match; it is an error
			     $error_hash{$qseqID} = "Problem: different chromosome ID";					    
			     last;
			}

		 }
	}

	#check the direction of the fragments; why not checking the order first? because it is meaningless to compare the coordinates of send and sstart if there are fragment in reverse and forward direction
	if($error_hash{$qseqID} eq "no error"){		
		foreach my $fragmentid (sort keys %{$transcript_hash{$qseqID}}){        
			foreach my $comparefragmentid (sort keys %{$transcript_hash{$qseqID}}){        
				 if ($transcript_hash{$qseqID}{$fragmentid}[7] != $transcript_hash{$qseqID}{$comparefragmentid}[7]){
					 $error_hash{$qseqID} = "Problem: different orientation";						
					 last;
				 }
			}
		}
	}

	#ccheck if fragments are assemblied in order;*****
	if ($error_hash{$qseqID} eq "no error"){		
	   if ((scalar keys %{$transcript_hash{$qseqID}}) ==1){
	   }
	   else {
		   foreach my $counter (sort {$a <=> $b} keys %{$transcript_hash{$qseqID}}){         
				 my $counter1 = $counter - 1;
				 if ($counter > 1){
					 if ((($transcript_hash{$qseqID}{$counter}[7] == 1) && ($transcript_hash{$qseqID}{$counter}[3] < $transcript_hash{$qseqID}{$counter1}[3]))||
					      (($transcript_hash{$qseqID}{$counter}[7] == -1) && ($transcript_hash{$qseqID}{$counter}[3] > $transcript_hash{$qseqID}{$counter1}[3]))){
						 	$error_hash{$qseqID} = "Problem: not in order";				
						 	last;
					 }
				 }
		    }
		}
	} 
}

sub sortties_tropicalis{
	my ($qseqID)=shift;
	my %tempties_hash;
	my %count_hash;
	my $ChrID;
	my $max=0;
	my $maxChrID;
	my %keepChrID; #the list of chrID to be keep when maxChrID does not exist in every fragment
	$error_hash{$qseqID} = "no error";

	#store the name of chromosome/scafold in the tempties_hash;
	foreach my $fragID(sort {$a cmp $b} keys %{$transcript_hash{$qseqID}}){
		foreach my $tiesID(sort {$a cmp $b} keys %{$transcript_hash{$qseqID}{$fragID}}){
			$ChrID = $transcript_hash{$qseqID}{$fragID}{$tiesID}[0];
			push (@{$tempties_hash{$ChrID}{$fragID}}, $tiesID); 
			
			#determine which ChrID occurs the most of the fragment(note: this is counting how many fragID have ChrID)
			if (scalar keys %{$tempties_hash{$ChrID}} > $max){
				$max=scalar keys %{$tempties_hash{$ChrID}}; #the occurance of $maxChrID
				$maxChrID = $ChrID;
			}
		}


	}

		####################test if it is possible that two chr occur the most equally############
		my @test;
		foreach my $ChrID(sort keys %tempties_hash){
		 	if (scalar keys %{$tempties_hash{$ChrID}} == $max){ #the occurance of $ChrID equal to $max, meaning it is the one of the maxChr that would generate optimal path
		 		push (@test,$ChrID);
		 	}
		}

		if (scalar @test > 1){
			selectOptimal($qseqID, $max, \@test); #send the qseqID and the list of maxChrID to the subroutine selectOptiaml
			return; #exist the subroutine sortties_tropicalis
		}
		#########################test only##################################

	#determine the best route and keep tiess that fit the best route	
	if (scalar keys %{$transcript_hash{$qseqID}} == $max){#there is a ChrID that occurs in all fragment
		foreach my $ChrID(sort keys %tempties_hash){ #delete all the tiess in the transcription_hash{$qseqID} that doesnt have ChrID
			if ($ChrID ne $maxChrID){
				foreach my $fragID(sort keys %{$tempties_hash{$ChrID}}){	
					foreach my $tiesID(@{$tempties_hash{$ChrID}{$fragID}}){
						delete $transcript_hash{$qseqID}{$fragID}{$tiesID};
					}
				}
			}
		}
	}
	else{#ChrID does not exist in every fragment, need to sort throught the fragment **********very badly structured!!!!! change it if have time
		foreach my $f1(sort keys %{$transcript_hash{$qseqID}}){
			if ($tempties_hash{$maxChrID}{$f1}){#if this fragment contain maxChrID, keep tiess with ChrID only and delete the rest.
					$keepChrID{$f1} = $maxChrID; 					
			}			
			else{
				my $tempChrID;
				my $priority = 4;
				foreach my $t1(sort keys %{$transcript_hash{$qseqID}{$f1}}){#select the one that will lead to less issue
 					my $tempChrID = $transcript_hash{$qseqID}{$f1}{$t1}[0];
 					my $tempChrIDnum = $tempChrID;
 					$tempChrIDnum =~ s/\D//g;
 					my $maxChrIDnum= $maxChrID; 
 					$maxChrIDnum =~ s/\D//g;
					if ($priority > 1 && $tempChrID eq $maxChrID){
						$keepChrID{$f1} = $tempChrID;
						$priority = 1;
					}
					elsif($priority>2 && $maxChrIDnum ~~ [1..10] && $tempChrIDnum==$maxChrIDnum){
						$keepChrID{$f1} = $tempChrID;
						$priority = 2;
					}
					elsif($priority>3 && $tempChrIDnum>10){
						$keepChrID{$f1} = $tempChrID;
						$priority = 3;
					}
					else{
						if (exists $keepChrID{$f1}){
						}
						else{
							$keepChrID{$f1} = $tempChrID;
							$priority =4;
						}
					}
				}	
			}
		}

		#delete tiess that don't have $keepChrID
		foreach my $f1(sort keys %{$transcript_hash{$qseqID}}){
			foreach my $t1(sort keys %{$transcript_hash{$qseqID}{$f1}}){
				if ($transcript_hash{$qseqID}{$f1}{$t1}[0] ne $keepChrID{$f1}){
					delete $transcript_hash{$qseqID}{$f1}{$t1};
				}
			}

		}
	}
	
	############after selecting which ties to keep, check if it has chromosome issue###########
	#check if the fragments are on the same chromosome
	my $firstfrag;
	$firstfrag = $maxChrID;
	my $firstfragname = $firstfrag =~ s/[^a-zA-Z]//g;#extrat chrID name - whetther it is a scaffold or chromosome
	my $firstfragnum = $firstfrag =~ s/\D//g;#extrat the chromosome number

	if ($firstfragname =~ "scaffold"){
		$error_hash{$qseqID} = "Contain scaffold";
	}
	
	foreach my $fragmentid (sort keys %keepChrID){
		my $comparefrag = $keepChrID{$fragmentid};
		my $comparefragname = $comparefrag; 
		$comparefragname =~ s/[^a-zA-Z]//g;#extrat chrID name - whetther it is a scaffold or chromosome
	    my $comparefragnum = $comparefrag; 
	    $comparefragnum =~ s/\D//g;#extrat the chromosome number
		if ($firstfrag ne $comparefrag){
		    if ($firstfragnum > 10 && $comparefragnum <=10){ #case 1: scaffold >10 compare with chr/scaffold 1-10; no error but document it as containing scaffold; meanwhile, prevent it proceed with orientation check.  
		     	 $error_hash{$qseqID} = "Contain scaffold";
		     	 $firstfragnum = $comparefragnum; #if the first scaffold > 10 and now encounter a scaffold/chr 1-10, change $firstfrag to scaffold/chr1-10. 
		     	 $firstfragname =$comparefragname;#because scaffold>10 will be ignored in comparison and will cause the un-documentation of error case 
		     	 $firstfrag = $comparefrag;
		    }
		    elsif(($firstfragnum <= 10 && $comparefragnum <=10)&&($firstfragnum == $comparefragnum)&&($firstfragname ne $comparefragname)){#case 2: scaffold 1-10 compare with chr/scaffold 1-10 and have the same number; no error but document it as containing scaffold; meanwhile, prevent it proceed with orientation check.
		    	 $error_hash{$qseqID} = "Contain scaffold";
		    }
		    elsif(($firstfragnum <= 10 && $comparefragnum <=10)&&($firstfragnum != $comparefragnum)){#case 3: scaffold/chr1-10 compare with caffold/chr1-10 and it doesnt match; it is an error
			     $error_hash{$qseqID} = "Problem: different chromosome ID";
			     last;
			}

		}
	}
			
	#########if it doesnt have chromosome issue, check if they have orientation issue#########
	#if there is orientation issue, keep the one has the best choice and document it as orientation issue########
	my $orientationID;
	if($error_hash{$qseqID} eq "no error"){	
		%tempties_hash = ();#empty the tempties_hash so that i can use it here without declare another temp hash
		
		foreach my $fragID(sort keys %{$transcript_hash{$qseqID}}){
			foreach my $tiesID(sort keys %{$transcript_hash{$qseqID}{$fragID}}){
				$orientationID = $transcript_hash{$qseqID}{$fragID}{$tiesID}[7];
				push (@{$tempties_hash{$orientationID}{$fragID}}, $tiesID); 			
			}
		}

 		$orientationID = 0; 
		if (scalar keys %{$transcript_hash{$qseqID}} == scalar keys %{$tempties_hash{"-1"}}){
			$orientationID = -1;
			##just a check, how many ties case with tiess occurs in both direction for all fragments
			if (scalar keys %{$tempties_hash{"-1"}} == scalar keys %{$tempties_hash{"1"}}){
				$check_bothdir ++;
			}

		}
		elsif (scalar keys %{$transcript_hash{$qseqID}} == scalar keys %{$tempties_hash{"1"}}) {#if there is one direction that exist in all fragment 
			$orientationID = 1;
		}
		else{
			$error_hash{$qseqID} = "Problem: different orientation";
		}

		#if $orientationID is not 0, deleted the one that doesnt have the right orientateion id $orientationID
		if ($orientationID){
			foreach my $f1(sort keys %{$transcript_hash{$qseqID}}){
				foreach my $t1(sort keys %{$transcript_hash{$qseqID}{$f1}}){
					if ($transcript_hash{$qseqID}{$f1}{$t1}[7] ne $orientationID){
						delete $transcript_hash{$qseqID}{$f1}{$t1};
					}
				}

			}
		}
	}

	#############if it doesnt have orientation issue, check if they have order issue#########
	if ($error_hash{$qseqID} eq "no error"){	
		##########if it doestnt have chromosome issue, re-order them based on qstart and qend#################
		foreach my $fragID (sort keys %{$transcript_hash{$qseqID}}){        
		    my $counter = 1;
		    my $firstkey = "empty";
		    $firstkey = ((keys %{$transcript_hash{$qseqID}{$fragID}})[0]); #try to get one existing ties(any key) in the fragtment using: ((keys %h)[0])
		    foreach my $comparefragmentid (sort keys %{$transcript_hash{$qseqID}}){ 
			    my $comparetiesID; 
			    $comparetiesID = ((keys %{$transcript_hash{$qseqID}{$comparefragmentid}})[0]);
			     if ($transcript_hash{$qseqID}{$fragID}{$firstkey}[1]>$transcript_hash{$qseqID}{$comparefragmentid}{$comparetiesID}[1]){
				     $counter ++;
			     }
		    }
		    
		    $transcript_hash{$qseqID}{$counter} = delete $transcript_hash{$qseqID}{$fragID};
		}

		########compare the sstart and send to access if they are in order############
		#though: they are a set of range(ex, forward direction). Take the first range and sort th ranges with range 1 and set 2
		my %range;
		my @sortedRange;
		my $start = 0;
		my $end = 0;
		foreach my $fragID (sort {$a <=> $b} keys %{$transcript_hash{$qseqID}}){
			foreach my $tiesID(sort {$a cmp $b} keys %{$transcript_hash{$qseqID}{$fragID}}){
				$start = $transcript_hash{$qseqID}{$fragID}{$tiesID}[3];
				$end = $transcript_hash{$qseqID}{$fragID}{$tiesID}[4];
				$range{$fragID}{$start}=$end; #store s-start ($start) and s-end ($end) in %range
			}
		}

		my $compareCoor;
		my $found = 0;
  
		foreach my $fragmentid(sort {$a <=> $b} keys %range){
			$found = 0; #re-set $found to be 0 when start a new fragment check
			if ($orientationID == 1){
				foreach my $i (sort {$a <=> $b} keys %{$range{$fragmentid}}){ #$i is s-start
					if ($fragmentid == 1){
						$compareCoor = $range{$fragmentid}{$i}; #if it is the first fragment, assign $compareCoor to be s-end of first range of fragment 1 
						$found = 1;
						last; #if it is the first fragment, exist the inner loop to skip to fragment 2;
					}
					if ($compareCoor < $i){
						$compareCoor = $range{$fragmentid}{$i};
						$found = 1;
					}
				}
				if ($found == 0){$error_hash{$qseqID}="Problem: not in order";}#if no s-end in the list is greater than $selected, it means they are not in order.}
			}
			elsif($orientationID == -1){
				foreach my $i (sort {$b <=> $a} keys %{$range{$fragmentid}}){ #$i is s-start; reverse sorting for reverse direction
					if ($fragmentid == 1){
						$compareCoor = $i; #if it is the first fragment, assign $compareCoor to be s-end of first range of fragment 1 
						$found = 1;
						last; #if it is the first fragment, exist the inner loop to skip to fragment 2;
					}
					if ($compareCoor < $range{$fragmentid}{$i}){
						$compareCoor = $i;
						$found = 1;
					}
				}
				if ($found == 0){$error_hash{$qseqID}="Problem: not in order";}#if no s-end in the list is greater than $selected, it means they are not in order.}
			}
		}
	} 	
	
}


sub selectOptimal{
	my $qseqID=shift;
	my $max = shift;
	my $maxChrIDlist = shift;
	my @maxChrIDlist = @$maxChrIDlist;
	my %temptranscript_hash;
	my %tempties_hash;
	my %keepChrID;
	my %patherror_hash;
	my $pathnumber;

	#determine the best route and keep tiess that fit the best route	
	foreach my $maxChrID (@maxChrIDlist){
		$pathnumber++;
		$patherror_hash{$pathnumber} = "no error";
		if (scalar keys %{$transcript_hash{$qseqID}} == $max){#there is a ChrID that occurs in all fragment
			foreach my $fragID (sort {$a cmp $b} keys %{$transcript_hash{$qseqID}}){
				foreach my $tiesID(sort {$a cmp $b} keys %{$transcript_hash{$qseqID}{$fragID}}){
					if ($transcript_hash{$qseqID}{$fragID}{$tiesID}[0]eq $maxChrID){
						foreach my $value (@{$transcript_hash{$qseqID}{$fragID}{$tiesID}}){
							push (@{$temptranscript_hash{$pathnumber}{$qseqID}{$fragID}{$tiesID}},$value); 
						}
					}
				}
			}
		}
		else{#ChrID does not exist in every fragment, need to sort throught the fragment **********very badly structured!!!!! change it if have time
			foreach my $f1(sort keys %{$transcript_hash{$qseqID}}){
				my $tempChrID;
				my $priority = 4;
				foreach my $t1(sort keys %{$transcript_hash{$qseqID}{$f1}}){#select the one that will lead to less issue
 					my $tempChrID = $transcript_hash{$qseqID}{$f1}{$t1}[0];
 					my $tempChrIDnum = $tempChrID;
 					$tempChrIDnum =~ s/\D//g;
 					my $maxChrIDnum= $maxChrID; 
 					$maxChrIDnum =~ s/\D//g;
					if ($priority > 1 && $tempChrID eq $maxChrID){
						$keepChrID{$f1} = $tempChrID;
						$priority = 1;
					}
					elsif($priority>2 && $maxChrIDnum ~~ [1..10] && $tempChrIDnum==$maxChrIDnum){
						$keepChrID{$f1} = $tempChrID;
						$priority = 2;
					}
					elsif($priority>3 && $tempChrIDnum>10){
						$keepChrID{$f1} = $tempChrID;
						$priority = 3;
					}
					else{
						if (exists $keepChrID{$f1}){
						}
						else{
							$keepChrID{$f1} = $tempChrID;
							$priority =4;
						}
					}
				}	
				
			}

			#delete tiess that don't have $keepChrID
			foreach my $f1(sort keys %{$transcript_hash{$qseqID}}){
				foreach my $t1(sort keys %{$transcript_hash{$qseqID}{$f1}}){
					if ($transcript_hash{$qseqID}{$f1}{$t1}[0] eq $keepChrID{$f1}){
						foreach my $value (@{$transcript_hash{$qseqID}{$f1}{$t1}}){
							push (@{$temptranscript_hash{$pathnumber}{$qseqID}{$f1}{$t1}},$value); 
						}						
					}
				}

			}
		}

		############after selecting which ties to keep, check if it has chromosome issue###########
		#check if the fragments are on the same chromosome
		my $firstfrag;
		$firstfrag = $maxChrID;
		my $firstfragname = $firstfrag =~ s/[^a-zA-Z]//g;#extrat chrID name - whetther it is a scaffold or chromosome
		my $firstfragnum = $firstfrag =~ s/\D//g;#extrat the chromosome number
		
			if ($firstfragname =~ "scaffold"){
				$patherror_hash{$pathnumber} = "Contain scaffold";
			}

		foreach my $fragmentid (sort keys %keepChrID){
			my $comparefrag = $keepChrID{$fragmentid};
			my $comparefragname = $comparefrag; 
			$comparefragname =~ s/[^a-zA-Z]//g;#extrat chrID name - whetther it is a scaffold or chromosome
		    my $comparefragnum = $comparefrag; 
		    $comparefragnum =~ s/\D//g;#extrat the chromosome number
			if ($firstfrag ne $comparefrag){
			    if ($firstfragnum > 10 && $comparefragnum <=10){ #case 1: scaffold >10 compare with chr/scaffold 1-10; no error but document it as containing scaffold; meanwhile, prevent it proceed with orientation check.  
			     	 $patherror_hash{$pathnumber} = "Contain scaffold";
			     	 $firstfragnum = $comparefragnum; #if the first scaffold > 10 and now encounter a scaffold/chr 1-10, change $firstfrag to scaffold/chr1-10. 
			     	 $firstfragname =$comparefragname;#because scaffold>10 will be ignored in comparison and will cause the un-documentation of error case 
			     	 $firstfrag = $comparefrag;
			    }
			    elsif(($firstfragnum <= 10 && $comparefragnum <=10)&&($firstfragnum == $comparefragnum)&&($firstfragname ne $comparefragname)){#case 2: scaffold 1-10 compare with chr/scaffold 1-10 and have the same number; no error but document it as containing scaffold; meanwhile, prevent it proceed with orientation check.
			    	 $patherror_hash{$pathnumber} = "Contain scaffold";

			    }
			    elsif(($firstfragnum <= 10 && $comparefragnum <=10)&&($firstfragnum != $comparefragnum)){#case 3: scaffold/chr1-10 compare with caffold/chr1-10 and it doesnt match; it is an error
				     $patherror_hash{$pathnumber} = "Problem: different chromosome ID";
				     delete $temptranscript_hash{$pathnumber};
				     last;
				}

			}
		}
		

		#########if it doesnt have chromosome issue, check if they have orientation issue#########
		#if there is orientation issue, keep the one has the best choice and document it as orientation issue########
		my $orientationID;
		if($patherror_hash{$pathnumber} eq "no error"){	
			%tempties_hash = ();#empty the tempties_hash so that i can use it here without declare another temp hash
			
			foreach my $fragID(sort keys %{$temptranscript_hash{$pathnumber}{$qseqID}}){
				foreach my $tiesID(sort keys %{$temptranscript_hash{$pathnumber}{$qseqID}{$fragID}}){
					$orientationID = $temptranscript_hash{$pathnumber}{$qseqID}{$fragID}{$tiesID}[7];
					push (@{$tempties_hash{$orientationID}{$fragID}}, $tiesID); 			
				}
			}

	 		$orientationID = 0; 
			if (scalar keys %{$temptranscript_hash{$pathnumber}{$qseqID}} == scalar keys %{$tempties_hash{"-1"}}){
				$orientationID = -1;
				##just a check, how many ties case with tiess occurs in both direction for all fragments
				if (scalar keys %{$tempties_hash{"-1"}} == scalar keys %{$tempties_hash{"1"}}){
					$check_bothdir ++;
				}

			}
			elsif (scalar keys %{$temptranscript_hash{$pathnumber}{$qseqID}} == scalar keys %{$tempties_hash{"1"}}) {#if there is one direction that exist in all fragment 
				$orientationID = 1;
			}
			else{
				$patherror_hash{$pathnumber} = "Problem: different orientation";
				delete $temptranscript_hash{$pathnumber};
			}

			#if $orientationID is not 0, deleted the one that doesnt have the right orientateion id $orientationID
			if ($orientationID){
				foreach my $f1(sort keys %{$temptranscript_hash{$pathnumber}{$qseqID}}){
					foreach my $t1(sort keys %{$temptranscript_hash{$pathnumber}{$qseqID}{$f1}}){
						if ($temptranscript_hash{$pathnumber}{$qseqID}{$f1}{$t1}[7] ne $orientationID){
							delete $temptranscript_hash{$pathnumber}{$qseqID}{$f1}{$t1};
						}
					}

				}
			}
		}


		#############if it doesnt have orientation issue, check if they have order issue#########
		if ($patherror_hash{$pathnumber} eq "no error"){	
			##########if it doestnt have chromosome issue, re-order them based on qstart and qend#################
			foreach my $fragID (sort keys %{$temptranscript_hash{$pathnumber}{$qseqID}}){        
			    my $counter = 1;
			    my $firstkey = "empty";
			    $firstkey = ((keys %{$temptranscript_hash{$pathnumber}{$qseqID}{$fragID}})[0]); #try to get one existing ties(any key) in the fragtment using: ((keys %h)[0])
			    foreach my $comparefragmentid (sort keys %{$temptranscript_hash{$pathnumber}{$qseqID}}){ 
				    my $comparetiesID; 
				    $comparetiesID = ((keys %{$temptranscript_hash{$pathnumber}{$qseqID}{$comparefragmentid}})[0]);
				     if ($temptranscript_hash{$pathnumber}{$qseqID}{$fragID}{$firstkey}[1]>$temptranscript_hash{$pathnumber}{$qseqID}{$comparefragmentid}{$comparetiesID}[1]){
					     $counter ++;
				     }
			    }
			    
			    $temptranscript_hash{$pathnumber}{$qseqID}{$counter} = delete $temptranscript_hash{$pathnumber}{$qseqID}{$fragID};
			}

			########compare the sstart and send to access if they are in order############
			#though: they are a set of range(ex, forward direction). Take the first range and sort th ranges with range 1 and set 2
			my %range;
			my @sortedRange;
			my $start = 0;
			my $end = 0;
			foreach my $fragID (sort {$a <=> $b} keys %{$temptranscript_hash{$pathnumber}{$qseqID}}){
				foreach my $tiesID(sort {$a cmp $b} keys %{$temptranscript_hash{$pathnumber}{$qseqID}{$fragID}}){
					$start = $temptranscript_hash{$pathnumber}{$qseqID}{$fragID}{$tiesID}[3];
					$end = $temptranscript_hash{$pathnumber}{$qseqID}{$fragID}{$tiesID}[4];
					$range{$fragID}{$start}=$end; #store s-start ($start) and s-end ($end) in %range
				}
			}

			my $compareCoor;
			my $found = 0;
	  
			foreach my $fragmentid(sort {$a cmp $b} keys %range){
				$found = 0; #re-set $found to be 0 when start a new fragment check
				if ($orientationID == 1){
					foreach my $i (sort {$a <=> $b} keys %{$range{$fragmentid}}){ #$i is s-start
						if ($fragmentid == 1){
							$compareCoor = $range{$fragmentid}{$i}; #if it is the first fragment, assign $compareCoor to be s-end of first range of fragment 1 
							$found = 1;
							last; #if it is the first fragment, exist the inner loop to skip to fragment 2;
						}
						if ($compareCoor < $i){
							$compareCoor = $range{$fragmentid}{$i};
							$found = 1;
						}
					}
					if ($found == 0){
						$patherror_hash{$pathnumber}="Problem: not in order"; 
						delete $temptranscript_hash{$pathnumber};
					}#if no s-end in the list is greater than $selected, it means they are not in order
				}
				elsif($orientationID == -1){
					foreach my $i (sort {$b <=> $a} keys %{$range{$fragmentid}}){ #$i is s-start; reverse sorting for reverse direction
						if ($fragmentid == 1){
							$compareCoor = $i; #if it is the first fragment, assign $compareCoor to be s-end of first range of fragment 1 
							$found = 1;
							last; #if it is the first fragment, exist the inner loop to skip to fragment 2;
						}
						if ($compareCoor < $range{$fragmentid}{$i}){
							$compareCoor = $i;
							$found = 1;
						}
					}
					if ($found == 0){
						$patherror_hash{$pathnumber}="Problem: not in order"; 
						delete $temptranscript_hash{$pathnumber};
					}#if no s-end in the list is greater than $selected, it means they are not in order.}
				}
			}
		}#end of order check

	}#end of foreach loop for each MaxChr ID

	#now, we need to select the optimal path based on the error message stored in %patherror_hash. we will go throught the hash. if a path is marked as no error, it will be selected.
	#if no path was marked as no error, the first one that is marked as "contain scafold" will be selected. 
	my $priority = 4;
	my @pathinfo;
	foreach my $path (sort {$a <=> $b} keys%patherror_hash){
		if ($priority >1 && $patherror_hash{$path} eq "no error"){
			$priority =1;
			$pathinfo[0]=1;
			$pathinfo[1]="no error";
			$pathinfo[2]=$path;

		}
		elsif($priority > 2 && $patherror_hash{$path} eq "Contain scaffold"){
			$priority = 2;
			$pathinfo[0]=1;
			$pathinfo[1]="Contain scaffold";
			$pathinfo[2]=$path;

		}
		elsif($priority >3){
			$pathinfo[0]=0;
			$error_hash{$qseqID} = $patherror_hash{$path};
		}
	}

	if ($pathinfo[0] > 0){
		$error_hash{$qseqID} = $pathinfo[1];
		$transcript_hash{$qseqID} = ();
		#transfer the information of $temptranscript_hash{$pathinfo[2]} to $transcript_hash
		foreach my $fragID (sort {$a cmp $b} keys %{$temptranscript_hash{$pathinfo[2]}{$qseqID}}){
			foreach my $tiesID(sort {$a cmp $b} keys %{$temptranscript_hash{$pathinfo[2]}{$qseqID}{$fragID}}){
				foreach my $value (@{$temptranscript_hash{$pathinfo[2]}{$qseqID}{$fragID}{$tiesID}}){
					push (@{$transcript_hash{$qseqID}{$fragID}{$tiesID}},$value); 
				}
			}
		}
		return %transcript_hash;
	}
}

#------------------------------------------main body-------------------------------------------------------------------------------------------------------------------------------------------------------------------
while ( my $line = <INPUT>) {
	chomp($line);
	@line = split (/\s+/,$line); #split the line with spaces (\s = space; and \s+ = more than one spece);
	
	#compare new line (%temptranscript_hash) and old line (%transcript_hash)
	if (exists $transcript_hash{$line[0]}){ #if the transcript with seqid already exist, assign it to a temperary hash for later comparison
					
		#declare variables: store qend, qstart, and e-value of temp fragment in variable for easier visual; 
		my $tempqstart = $line[6];
		my $tempqend = $line[7];
		my $tempEvalue= $line[10];
		my %OverlappedFrag;
		
		foreach my $fragt (sort {$a cmp $b} keys %{$transcript_hash{$line[0]}}){ #compare the new fragment(temp fragment) to all the stored fragment (with same qseqid) in t_hash;
		    #asign value of qend, qstart and evalue to scalar variables for easier visual
		    my $tqstart = $transcript_hash{$line[0]}{$fragt}{"T0"}[1];
		    my $tqend = $transcript_hash{$line[0]}{$fragt}{"T0"}[2];
		    my $tEvalue = $transcript_hash{$line[0]}{$fragt}{"T0"}[5];	

		   if (($tqstart <= $tempqstart)&& ($tqend >= $tempqend)){ #fragment comparison case#1: temp fragment is a repeated region of existing fragment in hash or vice versa-> compare e-value;
		    	if ($tEvalue < $tempEvalue){ #the e-value of temp is larger of if have same e-value, temp is shorter -> dont store it;                     
					@line=(); #emptying the @line array 
		    	}
		    	elsif($tEvalue > $tempEvalue){ #e-value of temp is smaller or  -> store it; 
		    		delete $transcript_hash{$line[0]}{$fragt};
		    	}
		    	else{ #$tEvalue = $tempEvalue
					if($line[3] == $transcript_hash{$line[0]}{$fragt}{"T0"}[8]){#if have same e-value and same length
						#$tiesfragtid = $fragt;
						addties($fragt);
						$tiesseq_hash{$line[0]} = "1"; #the ties array store the qseqid of the gene that have tiess
						$tiesOccur ++;
						@line=();
					}
					else{
							@line=(); #emptying the @line array
					}
		    	}
		    	last; #exiting the foreach loop
		    }
		   elsif ((($tqstart > $tempqstart)&&($tqstart >= $tempqend)) ||
		         (($tqend <= $tempqstart)&&($tqend < $tempqend))) { #fragment comparison case #2: no overlap at all, do nothing to @line;
		    }
		   elsif((($tqstart == $tempqstart)&&($tqend < $tempqend))||
		   		(($tqstart > $tempqstart)&&($tqend == $tempqend))||
		   		(($tqstart > $tempqstart)&&($tqend < $tempqend))){ #fragment comparison case #4: temp cover the region of the existing fragment, store the id of the fragment in case it is long fragment that overlap with 3+ existing fragment -> compare e-value
		   		if ($tEvalue < $tempEvalue) {
             			@line=(); #emptying the @line array 
		    	    	last; #exiting the foreach loop
             	}
             	elsif ($tEvalue >= $tempEvalue){
             			$OverlappedFrag{$fragt}[0]=1;
             	}
		   }
		   elsif ((($tqstart < $tempqstart)&&($tqend < $tempqend)) || 
		         (($tqstart > $tempqstart)&&($tqend > $tempqend))){ #fragment comparison case #3: not a repeated fragment but have some overlapping part -> compare e-value
             		my $OverlapLength = 0;
             		if ($tqstart > $tempqstart){
             			$OverlapLength = $tempqend - $tqstart;
             		}
             		elsif($tqstart < $tempqstart){
             			$OverlapLength = $tqend - $tempqstart;
             		}

             		if ($OverlapLength > $OverlapLengthLimit){ #the overlapping part > 5
             			if ($tEvalue < $tempEvalue){ #if e-value of temp is larger, disgar temp
             				@line=(); #emptying @line array
		    	    		last; #exiting the foreach loop
             			}
             			else{ #if e-value of temp is smaller: keep it but set the value of $OverlappedFrag{$fragt} to 1 so that the existing fragment can be deleted when store tempfrag
             				$OverlappedFrag{$fragt}[0]=1;
             			}
             		}
             		else { #overlapping part is lower than the tolerated number (=5 at Oct 12, 2016). the temp fragment will be added to the hash without replacing the existing fragment   
             			$OverlappedFrag{$fragt}[0]=0;   			
             		}
		   }		   		    
		}#end of foreach loop 		
		
       	#make decision for overlapping fragment
       	if ((scalar keys %OverlappedFrag > 0)&& (@line)){
       		foreach my $OlFrag(sort keys %OverlappedFrag){
       			if ($OverlappedFrag{$OlFrag}[0] == 1){ #the condition is true when the value stored is = 1 -> delete the existing fragment; 
       				delete $transcript_hash{$line[0]}{$OlFrag}; 
       			}
       		}
        }
        
        if (@line){ #add the temp fragment to t_hash if it passed all the checks
        	$numfragment++;
        	addfragment(@line);
        }

        #if it is the last line of the file, do the check
	    if (eof){
	    	if ((exists $tiesseq_hash{$oldfragid}) && ($tiesseq_hash{$oldfragid} eq "1")){ #if it is a transcript with ties, sort it in subroutine and error will be identified in subroutine too 
				sortties_tropicalis ($oldfragid);
				printdetail($oldfragid);
				%transcript_hash=();
			}
			else{	
				checkerror_tropicalis($oldfragid);
				printdetail($oldfragid);
				%transcript_hash=();		
			}
	    }


    }
	else { #case 0: if the transcript with seqid does not exist, it means it is the first fragment(with the lowest e-value) of the transcript ->  store it.
            if ($oldfragid){
            	if ((exists $tiesseq_hash{$oldfragid}) && ($tiesseq_hash{$oldfragid} eq "1")){ #if it is a transcript with ties, sort it in subroutine and error will be identified in subroutine too 
					sortties_tropicalis ($oldfragid);

					##########TEST ONLY, PRINTING TRANSCRIPT FLAGGED AS ERROR TO SEE IF THEY ARE REALLY ERRORS OR NOT##########
					if ($error_hash{$oldfragid} eq "no error"||$error_hash{$oldfragid} eq "Contain scaffold"){ 
					}
					else{
						print "$oldfragid: $error_hash{$oldfragid}\n";
						print "This transcript contains tiess \n";
						foreach my $fragment (sort keys %{$transcript_hash{$oldfragid}}){
					        foreach my $ties(sort keys %{$transcript_hash{$oldfragid}{$fragment}}){
					        	print  "$oldfragid\t$fragment\t";
						        foreach my $i (0..$#{$transcript_hash{$oldfragid}{$fragment}{$ties}}){
						          print "$transcript_hash{$oldfragid}{$fragment}{$ties}[$i]\t";
						        }
						         print "\n";
						    }
					    } 
					}
					##############################PRINT END##################
					
					printdetail($oldfragid);
					%transcript_hash=();
				}
				else{	
					checkerror_tropicalis($oldfragid);

					##########TEST ONLY, PRINTING TRANSCRIPT FLAGGED AS ERROR TO SEE IF THEY ARE REALLY ERRORS OR NOT##########
					if ($error_hash{$oldfragid} eq "no error"||$error_hash{$oldfragid} eq "Contain scaffold"){ 
					}
					else{
						print "$oldfragid: $error_hash{$oldfragid}\n";
						if ($tiesseq_hash{$oldfragid}){#print detail of transcripts with tiess
							#print OUTPUT2 "$qseqID\t$error_hash{$qseqID}\n";
							foreach my $fragment (sort keys %{$transcript_hash{$oldfragid}}){
						        foreach my $ties(sort keys %{$transcript_hash{$oldfragid}{$oldfragid}}){
						        	print  "$oldfragid\t$fragment\t";
							        foreach my $i (0..$#{$transcript_hash{$oldfragid}{$fragment}{$ties}}){
							          print  "$transcript_hash{$oldfragid}{$fragment}{$ties}[$i]\t";
							        }
							         print  "\n";
							    }
						    }
						}
						else{#print detail of transcripts without tiess
							foreach my $fragment (sort {$a <=> $b} keys %{$transcript_hash{$oldfragid}}){
					        	print "$oldfragid\t$fragment\t";
						        foreach my $i (0..$#{$transcript_hash{$oldfragid}{$fragment}}){
						          print "$transcript_hash{$oldfragid}{$fragment}[$i]\t";
						        }
						         print "\n";
						    }
						} 
					}
					##############################PRINT END##################

					printdetail($oldfragid);
					%transcript_hash=();
				}
            }
            $numfragment = 1;
            addfragment (@line); 
            $oldfragid = $line[0];
	}

}
close INPUT;
$DoneLoading = 1;

######################################################Printing result summaries#################################################################################

	my $ChromosomeError = 0;
	my $OrientationError = 0;
	my $OrderError = 0;
	my $tiesCounter = 0;
	my $noerror = 0;

print OUTPUT3 "list of transcript and its error flag\n";
foreach my $qseqid (sort {$a cmp $b } keys %error_hash){
	if ($error_hash{$qseqid} eq "no error"||$error_hash{$qseqid} eq "Contain scaffold"){
		$noerror ++;
		print OUTPUT1 "$qseqid\n"; 
		print OUTPUT3 "$qseqid\t$error_hash{$qseqid}\n";
	}
	elsif($error_hash{$qseqid} eq "Problem: different chromosome ID"){
		$ChromosomeError++;
		print OUTPUT3 "$qseqid\t$error_hash{$qseqid}\n";
	}
	elsif($error_hash{$qseqid} eq "Problem: different orientation"){
		$OrientationError++;
		print OUTPUT3 "$qseqid\t$error_hash{$qseqid}\n";
	}
	elsif($error_hash{$qseqid} eq "Problem: not in order"){
		$OrderError++;
		print OUTPUT3 "$qseqid\t$error_hash{$qseqid}\n";
	}
}

print OUTPUT3 "list of ties\n";
foreach my $qseqid (sort {$a cmp $b } keys %tiesseq_hash){

	print OUTPUT3 "$qseqid\n";
}



print OUTPUT1 "Summary\n";
print OUTPUT1 "Total number of genes:", scalar keys %error_hash, " \n";
print OUTPUT1 "Number of error free assembled transcripts: $noerror (", ($noerror/scalar keys %error_hash)*100, "%)\n";
print OUTPUT1 "There are $ChromosomeError genes that were assembled with different Chomosomes (", ($ChromosomeError/scalar keys %error_hash)*100, "%)\n";
print OUTPUT1 "There are $OrientationError genes that have orientation problems (", ($OrientationError/scalar keys %error_hash)*100, "%) \n";
print OUTPUT1 "There are $OrderError genes that have order problem (", ($OrderError/scalar keys %error_hash)*100, "%)\n";


close OUTPUT1;
close OUTPUT2;
close OUTPUT3;
