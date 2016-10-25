Oct 15, 2016
```
#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

#####
#This program will analysis the blastn output of trinity assembly. 
#This script will read into the result file of blastn(in fasta format) and for each transcriptome sequence, aligned fragments will be filter based on e-value.
#After filtering, fragments will be re-order based on query coordinates. 
#The script will then identify the possible error in each transcritome sequence. 
############
#blast format 6
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
############
#progress
#The current script can select and store fragment in a hash(transcript_hash) based on evalue;
#done: foreach loop to compare the new fragment to all the fragment store in the hash
############
#decisions 
#hash storage order: 0=sseqid; 1=qsatart; 2=qend; 3=sstard; 4= send; 5=evalue; 6=bitsore; 7=forward(1)/reverse(-1); 8=length;
#3 main type of errors: are they the same chromosom; order of the chromosome; the orientation of the fragments(forward or reverse)
#use e-value as the index; always keep the fragment with the lowest e-value
#the toleranted overlapping region is 5bp
###########

#open the input file
# my $inputFile = "/Users/lenduensis/Desktop/XueSS/Testfile/testfile5m";
my $inputFile = "C:/Users/Sarah/Desktop/PerlScript/testfile1m.txt";
#my ($inputFile)=@ARGV;
my @line;
my %temptranscript_hash;
my %transcript_hash;
my %error_hash;
my $numfragment = 0;
my $DoneLoading = 0;
my $OverlapLength = 5;

GetOptions('OverlapLength=i' => $OverlapLength);
 
open(INPUT,"<", "$inputFile") or die "could not open the input file";

while ( my $line = <INPUT>) {
	chomp($line);
	@line = split (/\s+/,$line); #split the line with spaces (\s = space; and \s+ = more than one spece);	
		
	#compare new line (%temptranscript_hash) and old line (%transcript_hash)
	if (exists $transcript_hash{$line[0]}){ #if the transcript with seqid already exist, assign it to a temperary hash for later comparison
		%temptranscript_hash =(); #reset/empty %temptranscript 
		$temptranscript_hash{$line[0]}[0]= $line[1]; #seqid
		$temptranscript_hash{$line[0]}[1]= $line[6]+0; #qstart
		$temptranscript_hash{$line[0]}[2]= $line[7]+0; #qend
		$temptranscript_hash{$line[0]}[3] = $line[8]+0; #sstart
		$temptranscript_hash{$line[0]}[4]= $line[9]+0; #send
		$temptranscript_hash{$line[0]}[5]= $line[10]+0; #evalue
		$temptranscript_hash{$line[0]}[8]= $line[3]+0; #length
				
		#declare variables: store qend, qstart, and e-value of temp fragment in variable for easier visual; 
		my $tempqstart = $temptranscript_hash{$line[0]}[1];
		my $tempqend = $temptranscript_hash{$line[0]}[2];
		my $tempEvalue= $temptranscript_hash{$line[0]}[5];
		my %OverlappedFrag;
		
						
		foreach my $fragt (sort keys %{$transcript_hash{$line[0]}}){ #compare the new fragment(temp fragment) to all the stored fragment (with same qseqid) in t_hash;
		    #asign value of qend, qstart and evalue to scalar variables for easier visual
		    my $tqstart = $transcript_hash{$line[0]}{$fragt}[1];
		    my $tqend = $transcript_hash
		    {$line[0]}{$fragt}[2];
		    my $tEvalue= $transcript_hash{$line[0]}{$fragt}[5];	
		
		   if (($tqstart <= $tempqstart)&& ($tqend >= $tempqend)){ #fragment comparison case#1: temp fragment is a repeated region of existing fragment in hash or vice versa-> compare e-value;
		    	if ($tEvalue < $tempEvalue){ #the e-value of temp is larger of if have same e-value, temp is shorter -> dont store it;                     
				%temptranscript_hash =(); #emptying temp_hash 
		    	}
		    	elsif($tEvalue > $tempEvalue){ #e-value of temp is larger or if have same e-value, temp is longer -> store it; 
		    		delete $transcript_hash{$line[0]}{$fragt};
		    		addfragment(@line);
		    	}
		    	else{ #$tEvalue = $tempEvalue
				if($temptranscript_hash{$line[0]}[8] <= $transcript_hash{$line[0]}{$fragt}[8]){
					%temptranscript_hash =(); #emptying temp_hash;
				}
		    	}
		    	last; #exiting the foreach loop
		    }
		   elsif ((($tqstart > $tempqstart)&&($tqstart >= $tempqend)) ||
		         (($tqend <= $tempqstart)&&($tqend < $tempqend))) { #fragment comparison case #2: no overlap at all, do nothing to temp_hash;
		    }
		   elsif(($tqstart >= $tempqstart)&&($tqend <= $tempqend)){ #fragment comparison case #4: temp cover the region of the existing fragment, store the id of the fragment in case it is long fragment that overlap with 3+ existing fragment -> compare e-value
		   		if ($tEvalue < $tempEvalue) {
             			%temptranscript_hash =(); #emptying temp_hash 
		    	    	last; #exiting the foreach loop
             		}
             		elsif ($tEvalue > $tempEvalue){
             			$OverlappedFrag{$fragt}[0]=1;
             	}
		   }
		   elsif ((($tqstart < $tempqstart)&&($tqend < $tempqend)) || 
		         (($tqstart > $tempqstart)&&($tqend > $tempqend))){ #fragment comparison case #3: not a repeated fragment but have some overlapping part -> compare e-value
             		if ((($tqstart < $tempqstart)&&(($tqend - $tempqstart)>$OverlapLength))
             		   ||(($tqstart > $tempqstart)&&(($tempqend - $tqstart) > $OverlapLength))){ #the overlapping part > 5
             			if ($tEvalue < $tempEvalue){ #if e-value of temp is larger, disgar temp
             				%temptranscript_hash =(); #emptying temp_hash 
		    	    		last; #exiting the foreach loop
             			}
             			else{ #if e-value of temp is smaller: keep it but set the value of $OverlappedFrag{$fragt} to 1 so that the existing fragment can be deleted when store tempfrag
             				$OverlappedFrag{$fragt}[0]=1;
             			}
             		}
             		else { #overlapping part is lower than the toleranted number (=5 at Oct 12, 2016). the temp fragment will be added to the hash without replacing the existing fragment   
             			$OverlappedFrag{$fragt}[0]=0;   			
             		}
		   }		   		    
		}#end of foreach loop 		
		
       	#make decision for overlapping fragment
       	if (keys %OverlappedFrag > 0){
       		foreach my $OlFrag(sort keys %OverlappedFrag){
       			if ($OverlappedFrag{$OlFrag}[0]){ #the condition is true when the value stored is = 1 -> delete the existing fragment; 
       				delete ($transcript_hash{$line[0]}{$OlFrag}); 
       			}
       		}
        }
        
        if (keys %temptranscript_hash){ #add the temp fragment to t_hash if it passed all the checks
        	$numfragment++;
        	addfragment(@line);
        }
    }
	else { #case 0: if the transcript with seqid does not exist, it means it is the first fragment(with the lowest e-value) of the transcript ->  store it.
            $numfragment = 1;
            addfragment (@line);                                  
	}
}
close INPUT;
$DoneLoading = 1;

######################################################Finished storing the fragment; Analysis start#################################################################################

my $ChromosomeError = 0;
my $OrientationError = 0;
my $OrderError = 0;

if ($DoneLoading){
    #re-store the fragments in ascending order based on qstart position
    foreach my $qseqid (sort keys %transcript_hash){ #print the content stored in the transcript_hash to the screen; can modify to output the result to a output file
		foreach my $fragmentid (sort keys %{$transcript_hash{$qseqid}}){        
		    my $counter = 1;
		    foreach my $comparefragmentid (sort keys %{$transcript_hash{$qseqid}}){        
		     if ($transcript_hash{$qseqid}{$fragmentid}[1]>$transcript_hash{$qseqid}{$comparefragmentid}[1]){
			     $counter ++;
		     }
		    }
		    $transcript_hash{$qseqid}{"Frag$counter"}=delete $transcript_hash{$qseqid}{$fragmentid};
		}
    }

    foreach my $qseq (sort keys %transcript_hash){ #print the content stored in the transcript_hash to the screen; can modify to output the result to a output file
	    foreach my $fragment (sort keys %{$transcript_hash{$qseq}}){
	        print "$qseq        $fragment     ";
	        foreach my $i (0..$#{$transcript_hash{$qseq}{$fragment}}) {
	          print $transcript_hash{$qseq}{$fragment}[$i], "\t";
	        }
	         print "\n";
	    }
	     print "\n";
	}
  
 #identify the errors
    my $firstfrag;
    foreach my $qseqid (sort keys %transcript_hash){ #print the content stored in the transcript_hash to the screen; can modify to output the result to a output file
		$error_hash{$qseqid} = "no error";
		#check if the fragments are on the same chromosome
		$firstfrag = $transcript_hash{$qseqid}{"Frag1"}[0];
		foreach my $fragmentid (sort keys %{$transcript_hash{$qseqid}}){
			if ($firstfrag ne $transcript_hash{$qseqid}{$fragmentid}[0]){
			     $error_hash{$qseqid} = "Different Chromosome";
			     $ChromosomeError ++; 
			     last;
			 }
		}

		#check the direction of the fragments; why not checking the order first? because it is meaningless to compare the coordinates of send and sstart if there are fragment in reverse and forward direction
			if($error_hash{$qseqid} eq "no error"){		
				foreach my $fragmentid (sort keys %{$transcript_hash{$qseqid}}){        
					foreach my $comparefragmentid (sort keys %{$transcript_hash{$qseqid}}){        
					 if ($transcript_hash{$qseqid}{$fragmentid}[7] != $transcript_hash{$qseqid}{$comparefragmentid}[7]){
						 $error_hash{$qseqid} = "Orientation Problem";
						 $OrientationError ++;
						 last;
					 }
					}
				}
			}

		#check the order of the chromosome if they are in the same orientation;*****
		if ($error_hash{$qseqid} eq "no error"){		
		   if (scalar keys %{$transcript_hash{$qseqid}}==1){

		   }
		   else {
			   foreach my $counter (1..scalar keys %{$transcript_hash{$qseqid}}) {         
				 my $counter2 = $counter+1;
				 if (!$transcript_hash{$qseqid}{"Fragment$counter2"}[3]){
				 	last;
				 }
				 if ($transcript_hash{$qseqid}{"Fragment$counter"}[3] > $transcript_hash{$qseqid}{"Fragment$counter2"}[3]){
					 $error_hash{$qseqid} = "Order Problem";
					 $OrderError ++;
					 last;
				 }
			    }
			}
		} 
	}
	
	#count the number of each errors;
	print "There are ", scalar keys %transcript_hash, " total number of gene\n";
	print "There are $ChromosomeError genes that assembled with different Chomosome\n";
	print "There are $OrientationError genes that have orientation problems\n";
	print "There are $OrderError genes that have order problem\n";
}

#########################################################################Subroutine#################################################################################################
sub addfragment {
     #transcript_hash stored the data/info of each fragment as an array in the order of: 0=sseqid; 1=qsatart; 2=qend; 3=sstard; 4= send; 5=evalue; 6=bitsore; 7=forward(1)/reverse(-1); 8=length;
     $transcript_hash{$line[0]}{"fragment$numfragment"}[0]= $line[1]; 
    $transcript_hash{$line[0]}{"fragment$numfragment"}[1]= $line[6]+0; 
    $transcript_hash{$line[0]}{"fragment$numfragment"}[2]= $line[7]+0; 
    $transcript_hash{$line[0]}{"fragment$numfragment"}[3] = $line[8]+0;
    $transcript_hash{$line[0]}{"fragment$numfragment"}[4]= $line[9]+0;
    $transcript_hash{$line[0]}{"fragment$numfragment"}[5]= $line[10]+0;
     $transcript_hash{$line[0]}{"fragment$numfragment"}[6]= $line[11];
     $transcript_hash{$line[0]}{"fragment$numfragment"}[8]= $line[3];
    
    if ($transcript_hash{$line[0]}{"fragment$numfragment"}[3] < $transcript_hash{$line[0]}{"fragment$numfragment"}[4]){ #determine the orientation of the fragment: if forward = 1; if reverse = -1;
        $transcript_hash{$line[0]}{"fragment$numfragment"}[7] = 1;
    }
    else {
        $transcript_hash{$line[0]}{"fragment$numfragment"}[7] = -1;
    }
    return %transcript_hash;
}

```
