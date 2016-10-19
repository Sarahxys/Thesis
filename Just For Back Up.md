Oct 19, 2016
```
#!/usr/bin/perl

use strict;
use warnings;

#####
#This program will identify if the transcript is assembled in the right order.
############
#blast format 6
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
############
#hash storage order
#0=sseqid; 1=qsatart; 2=qend; 3=sstard; 4= send; 5=evalue; 6=bitsore; 7=forward(1)/reverse(-1); 8=length;
############
#progress
#The current script can select and store fragment in a hash(transcript_hash) based on evalue;
#done: foreach loop to compare the new fragment to all the fragment store in the hash
############
#General Guide 
#use e-value as the index; always keep the fragment with the lowest e-value
###########
#improvement needed
#done: create a larger database to allow better testing
#done: add a selection for slight overlap(tolerant 5bp)
#done: create a error hash and store the error in it
#3 main type of errors: are they the same chromosom; order of the chromosome; the orientation of the fragments(forward or reverse)
   

#open the input file
#my $infile = "/Users/lenduensis/Desktop/XueSS/Testfile/testfile1m";
my $inputFile = "C:/Users/Sarah/Desktop/Perl Practice/testfile1m.txt"; 
open(INPUT,"<", "$inputFile") or die "could not open the input file";


#store the fragment information in a hash: seqid, sseqid, qstart, qend, sstart, send, evalue, bitscore
my @line;
my %temptranscript_hash;
my %transcript_hash;
my %error_hash;
my $numfragment = 0;
my $DoneLoading = 0;

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
		$temptranscript_hash{$line[0]}[5]= $line[10]; #evalue
		$temptranscript_hash{$line[0]}[8]= $line[3]+0; #length
				
		#declare variables: store qend, qstart, and e-value of temp fragment in variable for easier visual; 
		my $tempqstart = $temptranscript_hash{$line[0]}[1];
		my $tempqend = $temptranscript_hash{$line[0]}[2];
		my $tempEvalue= $temptranscript_hash{$line[0]}[5];
		my $OverlapLength = 5;
		my %OverlappedFrag;
		
						
		foreach my $fragt (sort keys %{$transcript_hash{$line[0]}}){ #compare the new fragment to all the sotred fragment (with same qseqid) in hash;
		    #asign value of qend, qstart and evalue to scalar variables for easier visual
		    my $tqstart = $transcript_hash{$line[0]}{$fragt}[1];
		    my $tqend = $transcript_hash{$line[0]}{$fragt}[2];
		    my $tEvalue= $transcript_hash{$line[0]}{$fragt}[5];	
		
		   if (($tqstart <= $tempqstart)&& ($tqend >= $tempqend)){ #fragment comparison case#1: temp fragment is a repeated region of existing fragment in hash or vice versa-> compare e-value;
		    	if (($tEvalue < $tempEvalue)||
		    		(($tEvalue = $tempEvalue)&&($temptranscript_hash{$line[0]}[8]>$transcript_hash{$line[0]}{$fragt}[8]))){ #the e-value of temp is larger of if have same e-value, temp is shorter -> dont store it;                     
				%temptranscript_hash =(); #emptying temp_hash 
		    	  
		    	}
		    	elsif(($tEvalue > $tempEvalue)|| #e-value of temp is larger or if have same e-value, temp is longer -> store it;
		    		 (($tEvalue = $tempEvalue)&&($temptranscript_hash{$line[0]}[8]>$transcript_hash{$line[0]}{$fragt}[8]))){
		    		delete $transcript_hash{$line[0]}{$fragt};
		    		addfragment(@line);
		    		
		    	}
		    	last; #exiting the foreach loop
		    }
		   elsif ((($tqstart > $tempqstart)&&($tqstart >= $tempqend)) ||
		         (($tqend <= $tempqstart)&&($tqend < $tempqend))) { #fragment comparison case #2: no overlap at all, add the fragment to hash; 
		           $numfragment ++;
		           addfragment (@line);
		           last;#exit the foreach loop		
		    }
		   elsif ((($tqstart <= $tempqstart)&&($tqend <= $tempqend)) || #fragment comparison case #3: not a repeated fragment but have some overlapping part -> compare e-value
		         (($tqstart >= $tempqstart)&&($tqend >= $tempqend))){
             		if ((($tqend - $tempqstart)>$OverlapLength)||(($$tempqend - $tqstart) > $OverlapLength)){ #the overlapping part > 5
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
		   elsif(($tqstart >= $tempqstart)&&($tqend <= $tempqend)){ #fragment comparison case #4: temp cover the region of the existing fragment, store the id of the fragment in case it is long fragment that overlap with 3+ existing fragment -> compare e-value
		   	if ($tEvalue < $tempEvalue) {
             			%temptranscript_hash =(); #emptying temp_hash 
		    	    	last; #exiting the foreach loop
             		}
             		elsif ($tEvalue > $tempEvalue){
             			$OverlappedFrag{$fragt}[0]=1;
             		}
		   }		    
		}#end of foreach loop 

       	#make decision for overlapping fragment
       	if (keys %OverlappedFrag){
       		foreach my $OlFrag(sort keys %OverlappedFrag){
       			if ($OverlappedFrag{$OlFrag}[0]){ #the condition is true when the value stored is = 1 -> delete the existing fragment; 
       				delete ($transcript_hash{$line[0]}{$OlFrag}); 
       			}
       		}
           	addfragment(@line);
        }
        }
	else { #case 0: if the transcript with seqid does not exist, it means it is the first fragment(with the lowest e-value) of the transcript ->  store it. 
            $numfragment = 1;
            addfragment (@line);
            print "the qend for the first fragment stored in the hash with qseqid ", $line[0], "is: ", $transcript_hash{$line[0]}{"fragment$numfragment"}[4], "\n"; #check if hash is funtioning
	}	    
}
	
close INPUT;
$DoneLoading = 1;

foreach my $qseqid (sort keys %transcript_hash){ #print the content stored in the transcript_hash to the screen; can modify to output the result to a output file
    foreach my $fragmentid (sort keys %{$transcript_hash{$qseqid}}){
          print "$qseqid        $fragmentid     ";
        
        foreach my $i (0..$#{$transcript_hash{$qseqid}{$fragmentid}}) {
          print $transcript_hash{$qseqid}{$fragmentid}[$i], "\t";
          }
        
         print "\n";
    }
}

######################################################Finished storing the fragment; Analysis start#################################################################################
#order of the chromosome; the orientation of the gene(exons, individual exon) (forward or reverse)

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
		    $transcript_hash{$qseqid}{"Fragment$counter"}=delete $transcript_hash{$qseqid}{$fragmentid};
		}
    }

    
    #identify the errors
    foreach my $qseqid (sort keys %transcript_hash){ #print the content stored in the transcript_hash to the screen; can modify to output the result to a output file
		$error_hash{$qseqid} = "no error";
		
		#check if the fragments are on the same chromosome, directions, 
		foreach my $fragmentid (sort keys %{$transcript_hash{$qseqid}}){        
		    foreach my $comparefragmentid (sort keys %{$transcript_hash{$qseqid}}){			     
			     if ($transcript_hash{$qseqid}{$fragmentid}[0] ne $transcript_hash{$qseqid}{$comparefragmentid}[0]){
				     $error_hash{$qseqid} = "Different Chromosome";
				     $ChromosomeError ++; 
				     last;
			     }			     
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

		#check the order of the chromosome if they are in the same orientation;
		if ($error_hash{$qseqid} eq "no error"){		
		    foreach my $i (1..scalar keys %{$transcript_hash{$qseqid}}) {         
			 my $ii = $i+1;
			 if ($transcript_hash{$qseqid}{"Fragment$i"}[3] > $transcript_hash{$qseqid}{"Fragment$ii"}[3]){
				 $error_hash{$qseqid} = "Order Problem";
				 $OrderError ++;
				 last;
			 }
		    }
		}
	}
	
	#count the number of each errors;
	print "There are $ChromosomeError genes that assembled with different Chomosome/n";
	print "There are $OrientationError genes that have orientation problems";
	print "There are $OrderError genes that have order problem";

}



#########################################################################Subroutine#################################################################################################
sub addfragment {
     #transcript_hash stored the data/info of each fragment as an array in the order of: 0=sseqid; 1=qsatart; 2=qend; 3=sstard; 4= send; 5=evalue; 6=bitsore; 7=forward(1)/reverse(-1); 8=length;
     $transcript_hash{$line[0]}{"fragment$numfragment"}[0]= $line[1]; 
    $transcript_hash{$line[0]}{"fragment$numfragment"}[1]= $line[6]+0; 
    $transcript_hash{$line[0]}{"fragment$numfragment"}[2]= $line[7]+0; 
    $transcript_hash{$line[0]}{"fragment$numfragment"}[3] = $line[8]+0;
    $transcript_hash{$line[0]}{"fragment$numfragment"}[4]= $line[9]+0;
    $transcript_hash{$line[0]}{"fragment$numfragment"}[5]= $line[10];
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
